import os
import re
import argparse
import traceback
import uuid
import sys
import numpy as np
import subprocess
import tarfile
import yaml
import shutil
import xml.dom.minidom
import rasterio
import rasterio.control
import rasterio.crs
import rasterio.sample
import rasterio.vrt
import rasterio._features
from itertools import product
from tqdm import tqdm
from os.path import join
from glob import glob
from loguru import logger


def untar(filePath, extractFolderPath):
    logger.info("解压中...")
    try:
        # Determine the mode based on file extension
        if filePath.endswith('.tar.gz'):
            mode = 'r:gz'
        elif filePath.endswith('.tar.xz'):
            mode = 'r:xz'
        else:
            mode = 'r'
        
        # Extract the tar file directly to extractFolderPath
        with tarfile.open(filePath, mode) as tar:
            tar.extractall(path=extractFolderPath)
        
        # Check for extra nested folders and adjust if necessary
        extracted_folders = [f for f in os.listdir(extractFolderPath) if os.path.isdir(join(extractFolderPath, f))]
        for folder in extracted_folders:
            nested_path = join(extractFolderPath, folder)
            if os.path.exists(nested_path):
                for item in os.listdir(nested_path):
                    shutil.move(join(nested_path, item), extractFolderPath)
                os.rmdir(nested_path)
    except Exception:
        logger.error(f"文件{filePath}打开失败: \n", traceback.format_exc())


def rpc_ortho(srcfile, dem, dstfile, pixel_size=None, utm=True):
    logger.info("RPC正射校正中...")
    pixel_size = f"-tr {pixel_size} {pixel_size} -tap" if pixel_size else ""
    if utm:
        isNorth = 1 if os.path.basename(srcfile).split('_')[3][0] != 'S' else 0
        for string in os.path.basename(srcfile).split('_'):
            if re.match(r'^E\d', string) or re.match(r'^W\d', string):
                zone = str(int(float(string[1:]) / 6) + 31)  # int(longitude / 6) + 31
        zone = int('326' + zone) if isNorth == 1 else int('327' + zone)
        target_coords = f"EPSG:{zone}"
    else:
        target_coords = "EPSG:4326"  # WGS84
    gdalwarp = join(otbPath, r"bin\gdalwarp.exe")
    rpc_ortho_cmd = f"{gdalwarp} -t_srs {target_coords} {pixel_size} -rpc -to RPC_DEM={dem} -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -r cubic -dstnodata 0.00 -of GTiff -overwrite {srcfile} {dstfile}"
    subprocess.run([x for x in rpc_ortho_cmd.split(" ") if x != ""])


def reproject(input_raster, output_raster, target_coords="EPSG:4326"):
    gdalwarp = join(otbPath, r"bin\gdalwarp.exe")
    reproj_cmd = f"{gdalwarp} -t_srs {target_coords} -r cubic -of GTiff -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -overwrite {input_raster} {output_raster}"
    subprocess.run([x for x in reproj_cmd.split(" ") if x != ""])


def pansharpening_otb(spectral, pancromatic, output):
    logger.info("融合中...")
    os.environ['OTB_appPath'] = join(otbPath, r"lib\otb\applications")
    otb_bat = join(otbPath, r"bin\otbcli_BundleToPerfectSensor.bat")
    output_temp = join(os.path.dirname(output), f"temp_{str(uuid.uuid4())[:5]}.tif")
    dtype = 'uint16' if args.level == "DN" else 'float'
    cmd = f"{otb_bat} -inp {pancromatic} -inxs {spectral} -out {output_temp} {dtype}"
    subprocess.run(join(otbPath, "otbenv.bat"))
    subprocess.run(cmd)
    compress(output_temp, output)


def clip_raster_by_mask_layer(input_raster, mask, output_raster, pixel_size=None, target_coords=None):
    pixel_size = f"-tr {pixel_size} {pixel_size} -tap" if pixel_size else ""
    target_coords = f"-t_srs {target_coords}" if target_coords else ""
    mask_name = os.path.splitext(os.path.basename(mask))[0]
    gdalwarp = join(otbPath, r"bin\gdalwarp.exe")
    clip_cmd = f"{gdalwarp} -of GTiff {target_coords} {pixel_size} -cutline {mask} -cl {mask_name} -crop_to_cutline -r cubic -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER -multi -wo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS -overwrite {input_raster} {output_raster}"
    subprocess.run([x for x in clip_cmd.split(" ") if x != ""])


def compress(input, output):
    logger.info("压缩中...")
    gdal_translate = join(otbPath, r"bin\gdal_translate.exe")
    compress_cmd = f"{gdal_translate} -of GTiff -co NUM_THREADS=ALL_CPUS -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER {input} {output}"
    subprocess.run([x for x in compress_cmd.split(" ") if x != ""])
    os.remove(input)


def mosaic(input_files, mosaic_file, resample="cubic"):
    logger.info("镶嵌中...")
    paths = ""
    for file in input_files:
        paths = paths + f"{file}\n"
    vrt = input_files[0].replace(".tif", ".vrt")
    input_file_list = input_files[0].replace(".tif", ".txt")
    txt = open(input_file_list, "w")
    txt.write(paths)
    txt.close()
    gdalbuildvrt = join(otbPath, r"bin\gdalbuildvrt.exe")
    gdal_translate = join(otbPath, r"bin\gdal_translate.exe")
    buildvrt_cmd = f"{gdalbuildvrt} -resolution highest -r {resample} -input_file_list {input_file_list} {vrt}"
    translate_cmd = f"{gdal_translate} -of GTiff -co NUM_THREADS=ALL_CPUS -co COMPRESS=LZW -co PREDICTOR=2 -co BIGTIFF=IF_SAFER {vrt} {mosaic_file}"
    subprocess.run([x for x in  buildvrt_cmd.split(" ") if x != ""])
    subprocess.run([x for x in translate_cmd.split(" ") if x != ""])
    os.remove(input_file_list)
    os.remove(vrt)


def save_txt(txt_path, contents, mode="a"):
    fh = open(txt_path, mode)  # 'a'为新增，'w'为覆盖
    fh.write('%s\n' % contents)
    fh.close()


def build_pyramid(file):
    logger.info("建金字塔中...")
    gdaladdo = join(otbPath, r"bin\gdaladdo.exe")
    build_pyramid_cmd = f"{gdaladdo} {file} -r nearest -ro --config COMPRESS_OVERVIEW LZW --config BIGTIFF_OVERVIEW IF_SAFER"
    subprocess.run([x for x in build_pyramid_cmd.split(" ") if x != ""])
    return


def stack(rasters, out_raster):
    logger.info("波段堆叠中...")
    count = sum([rasterio.open(raster).count for raster in rasters])

    with rasterio.open(rasters[0], 'r') as src:
        profile = src.profile
        profile.update(count=count, BIGTIFF='IF_SAFER')

        with rasterio.open(out_raster, 'w', **profile) as dst:
            band_index = 0
            for raster in rasters:
                with rasterio.open(raster, "r") as src_raster:
                    for idx in range(1, src_raster.count + 1):
                        band = src_raster.read(indexes=idx)
                        band_index += 1
                        dst.write(band, band_index)


class correction():
    def __init__(self, input_raster, metadata):
        self.input_raster = input_raster
        self.config = yaml.safe_load((open(join(appPath, "data", "AtmosphericCorrectionParameter.yaml"))))
        filename_split = os.path.basename(input_raster).split("_")
        self.SatelliteID = filename_split[0]
        self.SensorID = filename_split[1]
        self.Year = filename_split[4][:4]
        self.ImageType = "MUX" if "MUX" in filename_split[-1] else "MSS"
        self.ImageType = "PAN" if "PAN" in filename_split[-1] else self.ImageType
        self.metadata = metadata

        self.supported_sensor = ['GF1-PMS1/2', 
                                 'GF1-WFV1/2/3/4',
                                 'GF1B/C/D-PMS',
                                 'GF2-PMS1/2',
                                 'GF6-PMS',
                                 'GF6-WFV',
                                 'GF7-DLC',
                                 'GF7-BWD',
                                 'ZY303-TMS']

        self.is_valid = True
        if 'ESUN' not in str(self.config.get(self.SatelliteID, {}).get(self.SensorID, {})):
            logger.warning(f"缺少{self.SatelliteID}_{self.SensorID}的太阳辐照度（ESUN）参数，大气校正目前仅支持如下传感器：{self.supported_sensor}")
            logger.warning("退出大气校正！")
            self.is_valid = False  # Set is_valid to False if condition is met
        
    def apply_corrections(self):
        with rasterio.open(self.input_raster) as src:
            profile = src.profile
            profile.update(dtype=np.float32, blockxsize=100, blockysize=100, BIGTIFF='IF_SAFER')
            rasterTOA = self.input_raster.replace(".tif", f"_TOA.tif")
            with rasterio.open(rasterTOA, 'w', **profile) as dst:
                windows = [window for _, window in dst.block_windows()]
                for window in tqdm(windows, desc="大气表观反射率计算中..."):
                    data = src.read(window=window).astype(np.float32)  # Ensure float32
                    data = self.radiometric_block(data, self.metadata)
                    dst.write(data.astype(np.float32), window=window)
        
        dark_object_dns = self.calculate_dark_object(rasterTOA)
        with rasterio.open(rasterTOA) as src:
            rasterSR = self.input_raster.replace(".tif", f"_SR.tif")
            with rasterio.open(rasterSR, 'w', **profile) as dst:
                windows = [window for _, window in dst.block_windows()]
                for window in tqdm(windows, desc="地表反射率计算中..."):
                    data = src.read(window=window).astype(np.float32)  # Ensure float32
                    data = self.DOS_correction(data, dark_object_dns)
                    dst.write(data.astype(np.float32), window=window)
        os.remove(self.input_raster)
        os.remove(rasterTOA)
        os.rename(rasterSR, self.input_raster)

    def calculate_dark_object(self, inImage, patch_size=50, percentile=1):
        logger.info("计算暗物体反射率...")
        with rasterio.open(inImage, 'r') as src:
            nodata = src.profile.get("nodata", None)
            num_bands = src.count
            min_values_per_band = [[] for _ in range(num_bands)]
            # Generate windows based on patch_size if the image is not inherently blocked
            nrows, ncols = src.shape
            offsets = product(range(0, nrows, patch_size), range(0, ncols, patch_size))
            windows = [((row_start, min(row_start + patch_size, nrows)), 
                        (col_start, min(col_start + patch_size, ncols))) for row_start, col_start in offsets]
            
            for window in windows:
                X = src.read(window=window)
                if nodata is not None:
                    mask = (X != nodata) & (X != 0)
                else:
                    mask = (X != 0)
                for band in range(num_bands):
                    X_masked = X[band][mask[band]]
                    if X_masked.size > 0:
                        min_value = np.min(X_masked)
                        min_values_per_band[band].append(min_value)

        # Compute the dark object DN for the entire image for each band using the min values
        dark_object_dns_per_band = [np.percentile(min_values, percentile).astype(np.float32) for min_values in min_values_per_band]
        
        return dark_object_dns_per_band
    
    def radiometric_block(self, block, metadata):
        """
        Apply radiometric correction on the input block.
        """
        block = block.astype(np.float32)
        bands, _, _ = block.shape
        cosine_solar_zenith = np.cos(np.deg2rad(float(xml.dom.minidom.parse(metadata).getElementsByTagName('SolarZenith')[0].firstChild.data)))
        for m in range(bands):
            params = self.RadiometricCalibration(m+1)
            # TOA calculation
            block[m] = (block[m] * params[-2] + params[-1])
            block[m] = (np.pi * block[m]) / (params[0] * cosine_solar_zenith)
        return block

    def DOS_correction(self, data, dark_object_dns):
        """
        Apply Dark Object Subtraction (DOS) correction on the input data.
        """
        data = data.astype(np.float32)  # Convert to float32 for calculations
        for band in range(data.shape[0]):
            data[band] -= dark_object_dns[band]
            data[band] = np.clip(data[band], 0, None)
        return data

    def RadiometricCalibration(self, bandID):
        # 读取增益、偏移、太阳辐照度参数
        if self.SensorID.startswith('WFV'):
            esun = self.config[self.SatelliteID][self.SensorID]['ESUN'][bandID - 1]
            dct = self.config[self.SatelliteID][self.SensorID]['gain']
            self.Year = str(min([int(k) for k in dct if k.isdigit()], key=lambda x: abs(x - int(self.Year))))
            gain = self.config[self.SatelliteID][self.SensorID]['gain'][self.Year][bandID - 1]
            offset = self.config[self.SatelliteID][self.SensorID]['offset'][self.Year][bandID - 1]
        else:
            esun = self.config[self.SatelliteID][self.SensorID][self.ImageType]['ESUN'][bandID - 1]
            dct = self.config[self.SatelliteID][self.SensorID][self.ImageType]['gain']
            self.Year = str(min([int(k) for k in dct if k.isdigit()], key=lambda x: abs(x - int(self.Year))))
            gain = self.config[self.SatelliteID][self.SensorID][self.ImageType]['gain'][self.Year][bandID - 1]
            offset = self.config[self.SatelliteID][self.SensorID][self.ImageType]['offset'][self.Year][bandID - 1]
        return [esun, gain, offset]
    

def main(outputDirectory, dem, level, pansharpen, pyramid, files):
    for i, file in enumerate(files):
        # logger.warning(f"处理第{i+1}景数据中，共{len(files)}景数据")
        logger.info(f"开始处理{file}")
        assert " " not in file, "请确保输入数据的路径不包含空格！"
        assert ".tar" in file, "仅支持原始tar压缩包，请不要自行解压！"

        # Extract to outputDirectory
        dataPath = join(outputDirectory, os.path.basename(file)[:-7])
        try:
            untar(file, dataPath)
        except:
            continue

        # Get spectral rasters and pancromatic rasters
        M, P = get_rasters_name(dataPath)

        # RPC Orthorectification for spectral rasters
        m_out_list = list()
        for m in M:
            logger.info("处理多光谱数据中...")
            try:
                m_in = join(dataPath, m)
                m_out = join(outputDirectory, m)
                rpc_ortho(m_in, dem, m_out)
                m_out_list.append(m_out)
            except:
                continue

        def ato_cor(m_out, metadata):
            correction_obj = correction(m_out, metadata)
            if correction_obj.is_valid:
                correction_obj.apply_corrections()
            return m_out
        
        if "AHSI" in M[0]:
            # Stack GF5(GF5B, ZY1E)-VN, SW
            # sort in [VN, SW] order
            m_out_list.sort(reverse=True)
            m_out = re.sub("_SW.geotiff|_VN.geotiff|_SW.tif|_VN.tif", ".tif", m_out_list[0])
            stack(m_out_list, m_out)
            metadata = re.sub("_SW.geotiff|_VN.geotiff", ".xml", m_in)
            m_out = ato_cor(m_out, metadata) if level == "Surface_Reflectance" else m_out
        elif "GF5B_VIMI" in M[0]:
            # Stack GF5B-B1-6, B7-12
            # sort in [B1-6, B7-12] order
            m_out_list.sort(reverse=False)
            m_out = re.sub("_B1_B6.tif[f]", ".tiff", m_out_list[0])
            stack(m_out_list, m_out)
            metadata = re.sub(".tiff", ".xml", m_in)
            m_out = ato_cor(m_out, metadata) if level == "Surface_Reflectance" else m_out
        elif "GF6_WFV" in M[0]:
            m_out = re.sub("-[0-9].tiff", ".tiff", join(outputDirectory, M[0]))
            mosaic(m_out_list, m_out)
            [os.remove(m) for m in m_out_list]
            metadata = re.sub("-[0-9].tiff", ".xml", m_in)
            m_out = ato_cor(m_out, metadata) if level == "Surface_Reflectance" else m_out
        else:
            for m_out in m_out_list:
                metadata = re.sub(".tiff", ".xml", join(dataPath, os.path.basename(m_out)))
                m_out = ato_cor(m_out, metadata) if level == "Surface_Reflectance" else m_out

        # 先大气校正再融合的顺序通常更为推荐，因为这样可以确保影像在融合前已经反映了地物的真实信息，并且可以保持影像的光谱一致性。
        fusion_list = []
        if pansharpen and len(P) > 0:
            for p in P:
                logger.info("处理全色数据中...")
                try:
                    p_in = join(dataPath, p)
                    p_out = join(outputDirectory, p)
                    rpc_ortho(p_in, dem, p_out)
                    p_out = ato_cor(p_out, metadata) if level == "Surface_Reflectance" else p_out
                    fusion = p_out.replace("PAN", "FUS")
                    # Trying to find corresponding spectral dataset
                    m = glob(p_out.rsplit("PAN", 1)[0] + "M*" + p_out.rsplit("PAN", 1)[1])[0]
                    pansharpening_otb(m, p_out, fusion)
                    fusion_list.append(fusion)
                    # os.remove(m)
                    # os.remove(p_out)
                except Exception:
                    logger.error("Error occurred: \n", traceback.format_exc())
                    pass

        # Build pyramid, we only build for the final raster
        if pyramid:
            if "GF6_WFV" in M[0] or "GF5B_VIMI" in M[0] or "AHSI" in M[0]:
                build_pyramid(m_out)
            elif pansharpen and len(fusion_list) > 0:
                [build_pyramid(fusion) for fusion in fusion_list]
            else:
                [build_pyramid(m_out) for m_out in m_out_list]


def get_rasters_name(dataPath):
    all_file = next(os.walk(dataPath))[2]
    # M is list of spectral rasters' filenames, P is list of pancromatic rasters' filenames
    M = list(filter(lambda x: re.match('.*MUX.*.tiff|.*MSS.*.tiff|.*WFV.*.tiff|GF4_PMS.*.tiff|GF5_AHSI.*geotiff|GF5B_AHSI.*_VN.tif.*|GF5B_AHSI.*_SW.tif.*|GF5B_VIMI.*[0-9].tif.*|ZY1[E-F]_AHSI.*_VN.tif.*|ZY1[E-F]_AHSI.*_SW.tif.*|ZY1F_IRS_NSR.*.tif.*|CB04A_WPM.*MSS.tiff|HJ2[A-B]_CCD.*.tiff', x) != None, all_file))
    P = list(filter(lambda x: re.match('.*PAN.*.tiff', x) != None, all_file))
    assert len(M) > 0, "未找到数据！仅支持如下卫星型号：[CB04A_WPM, GF1_PMS, GF1_WFV, GF1B_PMS, GF1C_PMS, GF1D_PMS, GF2_PMS, GF4_PMI, GF5_AHSI, GF5B_AHSI, GF5B_VIMI, GF6_PMS, GF6_WFV, GF7_BWD, GF7_DLC, HJ2A_CCD, HJ2B_CCD, ZY1E_VNIC, ZY1F_AHSI, ZY303_TMS]"
    return M, P


if __name__ == "__main__":
    # 获取编译后exe或py文本的根目录
    if getattr(sys, 'frozen', False):
        appPath = os.path.dirname(sys.executable)
    else:
        appPath = os.path.dirname(__file__)

    parser = argparse.ArgumentParser()
    parser.add_argument('--outputDirectory', dest='outputDirectory',
                        help='Folder to store preprocessed files',
                        type=str, default=None)
    parser.add_argument('--dem', dest='dem',
                        help='path for dem',
                        type=str, default=join(appPath, 'data', 'GMTED2km.tif'))
    parser.add_argument('--otbPath', dest='otbPath',
                        help='path for Orfeo Toolbox',
                        type=str, default=glob(join(appPath, 'OTB-*-Win64'))[0])
    parser.add_argument('--level', dest='level',
                        help='Processing level, options: DN, Surface_Reflectance',
                        type=str, default='DN')                    
    parser.add_argument('--pansharpen', dest='pansharpen',
                        help='Whether to perform pansharpening',
                        action='store_true', default=False)
    parser.add_argument('--pyramid', dest='pyramid',
                        help='Whether to build pyramid',
                        action='store_true', default=False)
    parser.add_argument('--cache', dest='cache',
                        help='Set GDAL raster block cache size, may speed up processing with higher percentage, default is 5% of usable physical RAM',
                        type=str, default='5%')
    parser.add_argument('--files', dest='files',
                        help='List of files to process',
                        nargs='+', default=None)
    args = parser.parse_args()

    if args.pansharpen:
        # 确保路径为全英文
        for f in args.files:
            assert f.isascii(), "镶嵌前，请确保输入数据的路径不包含中文！"

    otbPath = args.otbPath
    os.environ['GDAL_DATA'] = join(otbPath, r"share\data")
    os.environ['PROJ_LIB'] = join(otbPath, r"share\proj")
    os.environ['GDAL_DRIVER_PATH'] = 'disable'
    os.environ["GDAL_CACHEMAX"] = args.cache
    main(args.outputDirectory, args.dem, args.level, args.pansharpen, args.pyramid, args.files)
