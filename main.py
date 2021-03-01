import os
import sys
import osgeo
from osgeo import gdal
from Py6S import *
import numpy as np
import subprocess
import tarfile
import Py6S
import json
import xml.dom.minidom
import rasterio as rio

def untar(fname, dirs):
    try:
        t = tarfile.open(fname)
    except Exception as e:
        print("文件%s打开失败" % fname)
    t.extractall(path=dirs)


def rpc_ortho(srcfile, dem, dstfile, pixel_size=None, utm=True, **kwargs):
    pixel_size = f"-tr {pixel_size} {pixel_size} -tap" if pixel_size else ""
    if utm:
        isNorth = 1 if os.path.basename(srcfile).split('_')[3][0] == 'N' else 0
        zone = str(int(float(os.path.basename(srcfile).split('_')[2][1:]) / 6) + 31)  # int(longitude / 6) + 31
        zone = int('326' + zone) if isNorth else int('327' + zone)
        target_coords = f"EPSG:{zone}"
    else:
        target_coords = "EPSG:4326" # WGS84
    rpc_ortho_cmd = "gdalwarp -t_srs {target_coords} {pixel_size} -rpc -to RPC_DEM={RPC_DEM} -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -r cubic -dstnodata {nodata} -of GTiff -overwrite -q {srcfile} {dstfile}".format(target_coords=target_coords, pixel_size=pixel_size, RPC_DEM=dem, nodata=0.00, srcfile=srcfile, dstfile=dstfile)
    subprocess.run([x for x in rpc_ortho_cmd.split(" ") if x != ""])


def reproject(input_raster, output_raster, target_coords="EPSG:4326"):
    reproj_cmd = f"gdalwarp -t_srs {target_coords} -r cubic -of GTiff -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -overwrite -q {input_raster} {output_raster}"
    subprocess.run([x for x in reproj_cmd.split(" ") if x != ""])


def pansharpening(spectral, pancromatic, output):
    if sys.platform == "win32":
        # Put gdal_pansharpen.py alongside this main.py
        pansharpen_py = os.path.join(os.path.dirname(__file__), 'gdal_pansharpen.py')
        version = osgeo.gdal.__version__
        contents = f"__requires__ = 'GDAL=={version}'\n__import__('pkg_resources').run_script('GDAL=={version}', 'gdal_pansharpen.py')"
        save_txt(pansharpen_py, contents, "w")
        sys.path.append(os.getcwd())
        from gdal_pansharpen import gdal_pansharpen
        gdal_pansharpen(["", pancromatic, spectral, output, "-r", "cubic", "-of", "GTiff", "-threads", "ALL_CPUS", "-co", "NUM_THREADS=ALL_CPUS", "-co", "COMPRESS=LZW", "-co", "PREDICTOR=2", "-co", "BIGTIFF=IF_SAFER", "-q"])

    elif "linux" in sys.platform:
        pansharp_cmd = "gdal_pansharpen.py {} {} {} -r cubic -of GTiff -threads ALL_CPUS -co NUM_THREADS=ALL_CPUS -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER".format(pancromatic, spectral, output)
        subprocess.run([x for x in pansharp_cmd.split(" ") if x != ""])


def clip_raster_by_mask_layer(input_raster, mask, output_raster, pixel_size=None, target_coords=None):
    pixel_size = f"-tr {pixel_size} {pixel_size} -tap" if pixel_size else ""
    target_coords = f"-t_srs {target_coords}" if target_coords else ""
    mask_name = os.path.splitext(os.path.basename(mask))[0]
    clip_cmd = f"gdalwarp -of GTiff {target_coords} {pixel_size} -cutline {mask} -cl {mask_name} -crop_to_cutline -r cubic -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -overwrite {input_raster} {output_raster}"
    subprocess.run([x for x in clip_cmd.split(" ") if x != ""])


def compress(input, output):
    compress_cmd = "gdal_translate -of GTiff -co NUM_THREADS=ALL_CPUS -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -q {} {}".format(input, output)
    subprocess.run([x for x in compress_cmd.split(" ") if x != ""])
    os.remove(input)


def mosaic(input_files, mosaic_file, resample="cubic"):
    paths = ""
    for file in input_files:
        paths = paths + f"{file}\n"
    vrt = input_files[0].replace(".tif", ".vrt")
    input_file_list = input_files[0].replace(".tif", ".txt")
    txt = open(input_file_list, "w")
    txt.write(paths)
    txt.close()
    buildvrt_cmd = "gdalbuildvrt -resolution highest -r {resample} -q -input_file_list {input_file_list} {vrt}".format(resample=resample, input_file_list=input_file_list, vrt=vrt)
    translate_cmd = "gdal_translate -of GTiff -co NUM_THREADS=ALL_CPUS -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -q {vrt} {mosaic_file}".format(vrt=vrt, mosaic_file=mosaic_file)
    subprocess.run([x for x in  buildvrt_cmd.split(" ") if x != ""])
    subprocess.run([x for x in translate_cmd.split(" ") if x != ""])
    os.remove(input_file_list)
    os.remove(vrt)


def save_txt(txt_path, contents, mode="a"):
    fh = open(txt_path, mode)  # 'a'为新增，'w'为覆盖
    fh.write('%s\n' % contents)
    fh.close()


def build_pyramid(file):
    build_pyramid_cmd = "gdaladdo {file} -r nearest -ro --config COMPRESS_OVERVIEW LZW --config BIGTIFF_OVERVIEW IF_SAFER".format(file=file)
    subprocess.run([x for x in build_pyramid_cmd.split(" ") if x != ""])


def stack(rasters, out_raster):
    # Read in metadata
    profile = rio.open(rasters[0], 'r').profile
    count = 0
    for raster in rasters:
        count += rio.open(raster, 'r').count
    profile.update(
        count=count,
    )
    out_raster = rio.open(out_raster, 'w', **profile)
    band_index = 0
    for raster in rasters:
        for band in rio.open(raster, "r").read():
            band_index+=1
            out_raster.write(band, band_index)
    out_raster.close()

class correction():
    def __init__(self, input_raster, metadata=None):
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
        self.metadata = metadata
        self.rcfile = input_raster.replace(".tiff", "_rc_temp.tiff")
        self.sixsfile = input_raster.replace(".tiff", "_6s_temp.tiff")

    def radiometric(self, TOA=False):
        # 设置输出波段
        Driver = self.IDataSet.GetDriver()
        trans = self.IDataSet.GetGeoTransform()
        proj = self.IDataSet.GetProjection()
        outDataset = Driver.Create(self.rcfile, self.cols, self.rows, self.bands, gdal.GDT_Float32)
        outDataset.SetGeoTransform(trans)
        outDataset.SetProjection(proj)
        # 分别读取8个波段
        for m in range(1, self.bands+1):
            ImgBand = self.IDataSet.GetRasterBand(m)
            outband = outDataset.GetRasterBand(m)
            outband.SetNoDataValue(0)
            # 获取对应波段的太阳辐照度ESUN，增益gain，偏移bias
            params = self.RadiometricCalibration(m, TOA)
            try:
                Image = ImgBand.ReadAsArray(0, 0, self.cols, self.rows)
                Image = np.where(Image != 0, (Image * params[-2] + params[-1]), 0)
                if TOA:
                    # Equation: TOA = (π∗Lλ∗d2)/(ESUNλ∗cosθs), d=1
                    # https://semiautomaticclassificationmanual-v4.readthedocs.io/en/latest/remote_sensing.html#top-of-atmosphere-toa-reflectance
                    cosine_solar_zenith = np.cos(np.deg2rad(float(xml.dom.minidom.parse(self.metadata).getElementsByTagName('SolarZenith')[0].firstChild.data)))
                    Image = (np.pi * Image) / (params[0] * cosine_solar_zenith)
                outband.WriteArray(Image, 0, 0)
                outband.FlushCache()
            except Exception as e:
                pass
        return self.rcfile

    def atmospheric(self):
        # 设置输出波段
        Driver = self.IDataSet.GetDriver()
        trans = self.IDataSet.GetGeoTransform()
        proj = self.IDataSet.GetProjection()
        outDataset = Driver.Create(self.sixsfile, self.cols, self.rows, self.bands, gdal.GDT_Float32)
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
                Image = np.where(Image != 0, AtcCofa * Image - AtcCofb, 0)
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

    def RadiometricCalibration(self, BandId, TOA):
        # 如果当年定标系数未发布或太懒没写，使用最近一年的定标系数
        dct = self.config["Parameter"][self.SatelliteID][self.SensorID]
        self.Year = str(min([int(k) for k in dct if k.isdigit()], key=lambda x: abs(x - int(self.Year))))

        params = list()
        if self.SensorID[0:3] == "WFV":
            Gain = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year]["gain"][BandId - 1]
            Bias = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year]["offset"][BandId - 1]
            if TOA:
                ESUN = self.config["Parameter"][self.SatelliteID][self.SensorID]["2020"]["ESUN"][BandId - 1]
                params.append(ESUN)
        else:
            Gain = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year][self.ImageType]["gain"][BandId - 1]
            Bias = self.config["Parameter"][self.SatelliteID][self.SensorID][self.Year][self.ImageType]["offset"][BandId - 1]
            if TOA:
                ESUN = self.config["Parameter"][self.SatelliteID][self.SensorID]["2020"][self.ImageType]["ESUN"][BandId - 1]
                params.append(ESUN)
        params.append(Gain)
        params.append(Bias)

        return params

    # 6s大气校正
    def get_ac_cof(self):
        # 读取头文件
        dom = xml.dom.minidom.parse(self.metadata)

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
            #print("Warning: A minimum of 64G RAM is required for GF6-WFV atmospheric correction")
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