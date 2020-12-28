# -*- coding:utf-8 -*-
import os
from os.path import join
import re
import glob
import argparse
import random
import PySimpleGUI as sg
from main import untar, rpc_ortho, build_pyramid, pansharpening, correction, compress, save_txt, mosaic


def main(InputFolder, OutputFolder, DEM, _6S, pansharpen, pyramid, pid, cache, InputFile_list):

    if InputFolder:
        files = glob.glob(join(InputFolder, "GF*"))
    elif InputFile_list:
        files = InputFile_list
    else:
        raise ValueError("Please specify --InputFolder or --InputFile_list")

    os.environ["GDAL_CACHEMAX"] = cache
    # Batch preprocess begins
    for i, file in enumerate(files):
        print("进程%d: 处理第%d景数据中，共%d景数据" % (pid, i+1, len(files)))
        try:
            # Extract archive
            if ".tar.gz" in file:
                # Extract to current directory
                datapath = file.replace(".tar.gz", "")
                try:
                    untar(file, datapath)
                except:
                    continue
            else:
                datapath = file

            # Get spectral rasters and pancromatic rasters
            M, P = get_rasters_name(datapath)

            # RPC Orthorectification for spectral rasters
            m_out_list = list()
            for m in M:
                try:
                    m_in = join(datapath, m)
                    m_out = join(OutputFolder, m)
                    rpc_ortho(m_in, DEM, m_out)
                    temp = correction(m_out).radiometric()
                    # Compress and delete temp file
                    compress(temp, m_out)
                    m_out_list.append(m_out)
                except:
                    continue

            # Mosaic GF6_WFV-1,2,3
            if "GF6_WFV" in M[0]:
                m = join(OutputFolder, M[0])
                mosaic_file = re.sub("-[0-9].tiff", "-mosaic.tiff", m)
                vrt = re.sub("-[0-9].tiff", ".vrt", m)
                txt = re.sub("-[0-9].tiff", ".txt", m)
                [save_txt(txt, m_out) for m_out in m_out_list]
                mosaic(mosaic_file, txt, vrt)
                [os.remove(m_out) for m_out in m_out_list]

            # 6S atmospheric correction
            if _6S:
                if "GF6_WFV" in M[0]:
                    metadata = re.sub("-[0-9].tiff", ".xml", join(datapath, M[0]))
                    temp = correction(mosaic_file, metadata).atmospheric()
                    compress(temp, mosaic_file)
                else:
                    for m_out in m_out_list:
                        metadata = re.sub(".tiff", ".xml", join(datapath, os.path.basename(m_out)))
                        temp = correction(m_out, metadata).atmospheric()
                        compress(temp, m_out)

            # Pansharpening
            if pansharpen and len(P)>0:
                fusion_list = list()
                for p in P:
                    p_in = join(datapath, p)
                    p_out = join(OutputFolder, p)
                    rpc_ortho(p_in, DEM, p_out)
                    temp = correction(p_out).radiometric()
                    compress(temp, p_out)
                    fusion = p_out.replace("PAN", "FUS")
                    fusion_list.append(fusion)
                    # Trying to find corresponding spectral dataset
                    m = glob.glob(p_out.rsplit("PAN", 1)[0] + "M*" + p_out.rsplit("PAN", 1)[1])[0]
                    pansharpening(m, p_out, fusion)
                    os.remove(m)
                    os.remove(p_out)

            # Build pyramid, we only build for the final raster
            if pyramid:
                if "GF6_WFV" in M[0]:
                    build_pyramid(mosaic_file)
                elif pansharpen and len(P) > 0:
                    [build_pyramid(fusion) for fusion in fusion_list]
                else:
                    [build_pyramid(m_out) for m_out in m_out_list]
        except:
            continue


def get_rasters_name(datapath):
    all_file = next(os.walk(datapath))[2]
    # M is list of spectral rasters' filenames, P is list of pancromatic rasters' filenames
    M = list(filter(lambda x: re.match('.*MUX.*.tiff|.*MSS.*.tiff|.*WFV.*.tiff|GF4_PMS.*.tiff', x) != None, all_file))
    P = list(filter(lambda x: re.match('.*PAN.*.tiff', x) != None, all_file))
    return M, P


def GUI():
    # GUI界面，随机主题
    theme = random.choice(sg.theme_list())
    sg.theme(theme)  # A touch of color

    # All the stuff inside your window.
    layout = [[sg.Text('Input Folder', font=("Helvetica", 25),
                       pad=((10, 50), (50, 50))), sg.InputText()],
              [sg.Text('Output Folder', font=("Helvetica", 25),
                       pad=((10, 50), (50, 50))), sg.InputText()],
              [sg.Text('DEM(optional)', font=("Helvetica", 25),
                       pad=((10, 50), (50, 50))), sg.InputText()],
              [sg.Checkbox('6S', default=False, font=("Helvetica", 25),
                           pad=((10, 25), (0, 35))),
               sg.Checkbox('Pansharpening', default=False, font=("Helvetica", 25),
                           pad=((10, 25), (0, 35))),
               sg.Checkbox('Build Pyramid', default=False, font=("Helvetica", 25),
                           pad=((10, 25), (0, 35)))],
              [sg.Button('OK', font=("Helvetica", 25),
                         pad=((275, 50), (0, 0))),
               sg.Button('CANCEL', font=("Helvetica", 25),
                         pad=((30, 50), (0, 0)))]]

    # Create the Window
    window = sg.Window('GaoFen-Batch', layout, default_element_size=(30, 1), size=(1150, 800))

    # Get the "values" of the inputs
    event, values = window.read(close=True)
    if event == sg.WIN_CLOSED or event == 'CANCEL':  # if user closes window or clicks 取消
        os._exit(1)

    return [values[i] for i in values]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--InputFolder', dest='InputFolder',
                        help='Folder which containing raw data, only one argument should be specified between InputFile_list and InputFolder',
                        default=None)
    parser.add_argument('--OutputFolder', dest='OutputFolder',
                        help='Folder to store preprocessed files',
                        type=str, default=None)
    parser.add_argument('--DEM', dest='DEM',
                        help='path for DEM',
                        type=str, default=os.path.join(os.path.dirname(__file__), 'data', 'GMTED2km.tif'))
    parser.add_argument('--_6S', dest='_6S',
                        help='Whether to perform 6S atmospheric correction',
                        action='store_true', default=False)
    parser.add_argument('--pansharpen', dest='pansharpen',
                        help='Whether to perform pansharpening',
                        action='store_true', default=False)
    parser.add_argument('--pyramid', dest='pyramid',
                        help='Whether to build pyramid',
                        action='store_true', default=False)
    parser.add_argument('--pid', dest='pid',
                        help='process id, for parallel mode',
                        type=int, default=1)
    parser.add_argument('--cache', dest='cache',
                        help='Set GDAL raster block cache size, may speed up processing with higher percentage, default is 5% of usable physical RAM',
                        type=str, default='5%')
    parser.add_argument('--InputFile_list', dest='InputFile_list',
                        help='List of files to process, reserved for parallel_batch.py',
                        nargs='+', default=None)
    args = parser.parse_args()

    # If there is no input from console, open graphic user interface
    if args.InputFolder or args.InputFile_list and args.OutputFolder:
        main(args.InputFolder, args.OutputFolder, args.DEM, args._6S, args.pansharpen, args.pyramid, args.pid, args.cache, args.InputFile_list)
    else:
        InputFolder, OutputFolder, DEM, _6S, pansharpen, pyramid = GUI()
        DEM = DEM if DEM else os.path.join(os.path.dirname(__file__), 'data', 'GMTED2km.tif')
        main(InputFolder, OutputFolder, DEM, _6S, pansharpen, pyramid, pid=1, cache="5%", InputFile_list=None)


