import os
from osgeo import gdal
from Py6S import *
import numpy as np
import subprocess
import tarfile
import Py6S
import json
import xml.dom.minidom


def untar(fname, dirs):
    try:
        t = tarfile.open(fname)
    except Exception as e:
        print("文件%s打开失败" % fname)
    t.extractall(path=dirs)


def rpc_ortho(srcfile, dem, dstfile, scale=None, **kwargs):
    isNorth = 1 if os.path.basename(srcfile).split('_')[3][0] == 'N' else 0
    zone = str(int(float(os.path.basename(srcfile).split('_')[2][1:]) / 6) + 31)
    zone = int('326' + zone) if isNorth else int('327' + zone)
    rpc_ortho_cmd = "gdalwarp -t_srs EPSG:{zone} -rpc -to RPC_DEM={RPC_DEM} -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -r cubic -dstnodata {nodata} -of GTiff -overwrite -q {srcfile} {dstfile}".format(zone=zone, RPC_DEM=dem, nodata=0.00, srcfile=srcfile, dstfile=dstfile)
    subprocess.run([x for x in rpc_ortho_cmd.split(" ") if x != ""])


def pansharpening(spectral, pancromatic, output):
    pansharp_cmd = "gdal_pansharpen.py {} {} {} -r cubic -of GTiff -nodata 0.00 -threads ALL_CPUS -co NUM_THREADS=ALL_CPUS -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -co BIGTIFF=IF_SAFER -q".format(pancromatic, spectral, output)
    subprocess.run([x for x in pansharp_cmd.split(" ") if x != ""])


def compress(input, output):
    compress_cmd = "gdal_translate -of GTiff -co NUM_THREADS=ALL_CPUS -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -co BIGTIFF=IF_SAFER -q {} {}".format(input, output)
    subprocess.run([x for x in compress_cmd.split(" ") if x != ""])
    os.remove(input)


def mosaic(mosaic_file, files_txt, vrt, resample="cubic"):
    buildvrt_cmd = "gdalbuildvrt -resolution highest -r {resample} -q -input_file_list {files_txt} {vrt}".format(resample=resample, files_txt=files_txt, vrt=vrt)
    translate_cmd = "gdal_translate -of GTiff -co NUM_THREADS=ALL_CPUS -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 -co BIGTIFF=IF_SAFER -q {vrt} {mosaic_file}".format(vrt=vrt, mosaic_file=mosaic_file)
    subprocess.run([x for x in buildvrt_cmd.split(" ") if x != ""])
    subprocess.run([x for x in translate_cmd.split(" ") if x != ""])


def save_txt(txt_path, contents):
    fh = open(txt_path, 'a')  # 'a'为新增，'w'为覆盖
    fh.write('%s\n' % contents)
    fh.close()


def build_pyramid(file):
    build_pyramid_cmd = "gdaladdo {file} -r nearest -ro --config COMPRESS_OVERVIEW DEFLATE --config BIGTIFF_OVERVIEW IF_SAFER".format(
        file=file)
    subprocess.run([x for x in build_pyramid_cmd.split(" ") if x != ""])


class correction():
    def __init__(self, input_raster, metadata=None, scale=10):
        super(correction, self).__init__()
        self.IDataSet = gdal.Open(input_raster)
        self.cols = self.IDataSet.RasterXSize
        self.rows = self.IDataSet.RasterYSize
        self.bands = self.IDataSet.RasterCount

        # 读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        self.config = json.load(open(os.path.join(os.path.dirname(__file__), "data", "RadiometricCorrectionParameter.json")))
        # 获取json文件索引信息
        filename_split = os.path.basename(input_raster).split("_")
        self.SatelliteID = filename_split[0]
        self.SensorID = filename_split[1]
        self.Year = filename_split[4][:4]
        self.ImageType = "MUX" if "MUX" in filename_split[-1] else "MSS"
        self.ImageType = "PAN" if "PAN" in filename_split[-1] else self.ImageType
        self.metedata = metadata
        self.scale = scale
        self.rcfile = input_raster.replace(".tiff", "_rc_temp.tiff")
        self.sixsfile = input_raster.replace(".tiff", "_6s_temp.tiff")

    def radiometric(self):
        # 设置输出波段
        Driver = self.IDataSet.GetDriver()
        trans = self.IDataSet.GetGeoTransform()
        proj = self.IDataSet.GetProjection()
        outDataset = Driver.Create(self.rcfile, self.cols, self.rows, self.bands, gdal.GDT_UInt16)
        outDataset.SetGeoTransform(trans)
        outDataset.SetProjection(proj)
        # 分别读取8个波段
        for m in range(1, self.bands+1):
            ImgBand = self.IDataSet.GetRasterBand(m)
            outband = outDataset.GetRasterBand(m)
            outband.SetNoDataValue(0)
            # 获取对应波段的增益gain和偏移bias
            Gain, Bias = self.RadiometricCalibration(m)
            try:
                Image = ImgBand.ReadAsArray(0, 0, self.cols, self.rows)
                Image = np.where(Image != 0, (Image * Gain + Bias) * self.scale, 0)
                outband.WriteArray(Image, 0, 0)
                outband.FlushCache()
            except:
                pass
        return self.rcfile

    def atmospheric(self):
        # 设置输出波段
        Driver = self.IDataSet.GetDriver()
        trans = self.IDataSet.GetGeoTransform()
        proj = self.IDataSet.GetProjection()
        outDataset = Driver.Create(self.sixsfile, self.cols, self.rows, self.bands, gdal.GDT_UInt16)
        outDataset.SetGeoTransform(trans)
        outDataset.SetProjection(proj)
        # 获取大气校正系数
        AtcCof_list = self.get_ac_cof()
        # 分别读取8个波段
        for m in range(1, self.bands+1):
            ImgBand = self.IDataSet.GetRasterBand(m)
            outband = outDataset.GetRasterBand(m)
            outband.SetNoDataValue(0)
            AtcCofa, AtcCofb, AtcCofc = AtcCof_list[m - 1]
            try:
                Image = ImgBand.ReadAsArray(0, 0, self.cols, self.rows)
                Image = np.where(Image != 0, AtcCofa * (Image / self.scale) - AtcCofb, 0)
                Image = np.where(Image != 0, (Image / (1 + AtcCofc * Image)) * 1000, 0)
                outband.WriteArray(Image, 0, 0)
                outband.FlushCache()
            except:
                pass
        return self.sixsfile

    def MeanDEM(self, pointUL, pointDR):
        '''
        计算影像所在区域的平均高程.
        '''
        try:
            DEM = os.path.join(os.path.dirname(__file__), 'data', 'GMTED2km.tif')
            DEMIDataSet = gdal.Open(DEM)
        except Exception as e:
            pass

        DEMBand = DEMIDataSet.GetRasterBand(1)
        geotransform = DEMIDataSet.GetGeoTransform()
        # DEM分辨率
        pixelWidth = geotransform[1]
        pixelHight = geotransform[5]

        # DEM起始点：左上角，X：经度，Y：纬度
        originX = geotransform[0]
        originY = geotransform[3]

        # 研究区左上角在DEM矩阵中的位置
        yoffset1 = int((originY - pointUL['lat']) / pixelWidth)
        xoffset1 = int((pointUL['lon'] - originX) / (-pixelHight))

        # 研究区右下角在DEM矩阵中的位置
        yoffset2 = int((originY - pointDR['lat']) / pixelWidth)
        xoffset2 = int((pointDR['lon'] - originX) / (-pixelHight))

        # 研究区矩阵行列数
        xx = xoffset2 - xoffset1
        yy = yoffset2 - yoffset1
        # 读取研究区内的数据，并计算高程
        DEMRasterData = DEMBand.ReadAsArray(xoffset1, yoffset1, xx, yy)

        MeanAltitude = np.mean(DEMRasterData)
        return MeanAltitude

    def RadiometricCalibration(self, BandId):
        # 如果当年定标系数未发布，使用最近一年的定标系数
        while self.Year not in self.config["Parameter"][self.SatelliteID][self.SensorID]:
            self.Year = str(int(self.Year) - 1)
        if self.SensorID[0:3] == "WFV":
            Gain = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year]["gain"][BandId - 1]
            Bias = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year]["offset"][BandId - 1]
        else:
            Gain = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year][self.ImageType]["gain"][BandId - 1]
            Bias = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year][self.ImageType]["offset"][BandId - 1]

        return Gain, Bias

    # 6s大气校正
    def get_ac_cof(self):
        # 读取头文件
        dom = xml.dom.minidom.parse(self.metedata)

        # 6S模型
        s = SixS()

        # 传感器类型 自定义
        s.geometry = Geometry.User()
        s.geometry.solar_z = 90 - float(dom.getElementsByTagName('SolarZenith')[0].firstChild.data)
        s.geometry.solar_a = float(dom.getElementsByTagName('SolarAzimuth')[0].firstChild.data)
        s.geometry.view_z = 0
        s.geometry.view_a = 0
        # 日期
        DateTimeparm = dom.getElementsByTagName('StartTime')[0].firstChild.data
        DateTime = DateTimeparm.split(' ')
        Date = DateTime[0].split('-')
        s.geometry.month = int(Date[1])
        s.geometry.day = int(Date[2])

        # 中心经纬度
        TopLeftLat = float(dom.getElementsByTagName('TopLeftLatitude')[0].firstChild.data)
        TopLeftLon = float(dom.getElementsByTagName('TopLeftLongitude')[0].firstChild.data)
        TopRightLat = float(dom.getElementsByTagName('TopRightLatitude')[0].firstChild.data)
        TopRightLon = float(dom.getElementsByTagName('TopRightLongitude')[0].firstChild.data)
        BottomRightLat = float(dom.getElementsByTagName('BottomRightLatitude')[0].firstChild.data)
        BottomRightLon = float(dom.getElementsByTagName('BottomRightLongitude')[0].firstChild.data)
        BottomLeftLat = float(dom.getElementsByTagName('BottomLeftLatitude')[0].firstChild.data)
        BottomLeftLon = float(dom.getElementsByTagName('BottomLeftLongitude')[0].firstChild.data)

        ImageCenterLat = (TopLeftLat + TopRightLat + BottomRightLat + BottomLeftLat) / 4

        # 大气模式类型
        if ImageCenterLat > -15 and ImageCenterLat < 15:
            s.atmos_profile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.Tropical)

        if ImageCenterLat > 15 and ImageCenterLat < 45:
            if s.geometry.month > 4 and s.geometry.month < 9:
                s.atmos_profile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeSummer)
            else:
                s.atmos_profile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.MidlatitudeWinter)

        if ImageCenterLat > 45 and ImageCenterLat < 60:
            if s.geometry.month > 4 and s.geometry.month < 9:
                s.atmos_profile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.SubarcticSummer)
            else:
                s.atmos_profile = Py6S.AtmosProfile.PredefinedType(Py6S.AtmosProfile.SubarcticWinter)

        # 气溶胶类型大陆
        s.aero_profile = Py6S.AtmosProfile.PredefinedType(AeroProfile.Continental)

        # 下垫面类型
        s.ground_reflectance = Py6S.GroundReflectance.HomogeneousLambertian(0.36)

        # 550nm气溶胶光学厚度,对应能见度为40km
        s.aot550 = 0.14497

        # 通过研究去区的范围去求DEM高度。
        pointUL = dict()
        pointDR = dict()
        pointUL["lat"] = max(TopLeftLat, TopRightLat, BottomRightLat, BottomLeftLat)
        pointUL["lon"] = min(TopLeftLon, TopRightLon, BottomRightLon, BottomLeftLon)
        pointDR["lat"] = min(TopLeftLat, TopRightLat, BottomRightLat, BottomLeftLat)
        pointDR["lon"] = max(TopLeftLon, TopRightLon, BottomRightLon, BottomLeftLon)
        meanDEM = (self.MeanDEM(pointUL, pointDR)) * 0.001

        # 研究区海拔、卫星传感器轨道高度
        s.altitudes = Altitudes()
        s.altitudes.set_target_custom_altitude(meanDEM)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)
        # 获取s.wavelength参数
        sensor = self.SatelliteID + "-" + self.SensorID
        support_sensors = ["GF1-PMS1", "GF1-PMS2", "GF1-WFV1", "GF1-WFV2", "GF1-WFV3", "GF1-WFV4", "GF1B-PMS", "GF1C-PMS", "GF1D-PMS", "GF2-PMS1", "GF2-PMS2", "GF6-WFV"]

        if sensor == support_sensors[-1]:
            print("Warning: A minimum of 64G RAM is required for GF6-WFV atmospheric correction")
            wavelengths = [(0.45, 0.52), (0.52, 0.59), (0.63, 0.69), (0.77, 0.89), (0.69, 0.73), (0.73, 0.77), (0.40, 0.45), (0.59, 0.63)]
            s_wl_list = self.get_py6s_wavelength(wavelengths)
        elif sensor in support_sensors[:-1]:
            wavelengths = [(0.45, 0.52), (0.52, 0.59), (0.63, 0.69), (0.77, 0.89)]
            s_wl_list = self.get_py6s_wavelength(wavelengths)
        elif sensor not in support_sensors and self.bands == 4:
            # Using bilitin spectral reflectance function of Sentinel-2
            s_wl_list = [Wavelength(PredefinedWavelengths.S2A_MSI_02),
                         Wavelength(PredefinedWavelengths.S2A_MSI_03),
                         Wavelength(PredefinedWavelengths.S2A_MSI_04),
                         Wavelength(PredefinedWavelengths.S2A_MSI_08)]
        else:
            print(f"ERROR: Atmospheric correction for {sensor} is not supported yet")

        AtcCof_list = list()
        for s_wl in s_wl_list:
            s.wavelength = s_wl
            # 运行6s大气模型
            s.run()
            xa = s.outputs.coef_xa
            xb = s.outputs.coef_xb
            xc = s.outputs.coef_xc
            AtcCof_list.append((xa, xb, xc))
        return AtcCof_list

    def get_py6s_wavelength(self, wavelengths):
        s_wl_list = list()
        for i, wl in enumerate(wavelengths):
            SRFband = self.config["Parameter"][self.SatelliteID][self.SensorID]["SRF"][f"{str(i+1)}"]
            s_wl = Wavelength(wl[0], wl[1], SRFband)
            s_wl_list.append(s_wl)
        return s_wl_list