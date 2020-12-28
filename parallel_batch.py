import os
import glob
import math
import time
import argparse
import subprocess
import multiprocessing
import numpy
import random
import PySimpleGUI as sg


def main(InputFolder, OutputFolder, DEM, _6S, pansharpen, pyramid, n_parallel, cache):
    os.environ["GDAL_CACHEMAX"] = cache
    batch_file = os.path.join(os.path.dirname(__file__), 'batch.py')
    files_groups = split_datasets(InputFolder, n_parallel)

    start = time.time()
    processes = []
    _6S = "--_6S" if _6S else ""
    pansharpen = "--pansharpen" if pansharpen else ""
    pyramid = "--pyramid" if pyramid else ""
    for i, files in enumerate(files_groups):
        pid = i + 1
        batch_cmd = f"python {batch_file} --InputFile_list {files} --OutputFolder {OutputFolder} --DEM {DEM} {_6S} {pansharpen} {pyramid} --pid {pid}"
        p = multiprocessing.Process(target=_subprocess, args=(batch_cmd,))
        p.start()
        time.sleep(5)
        processes.append(p)

    for process in processes:
        process.join()

    elapsed_h = (time.time() - start) / 60
    print("批处理总用时： " + format(elapsed_h, ".1f") + " min")


def split_datasets(InputFolder, n_parallel):
    # Split datasets to several groups
    files = glob.glob(os.path.join(InputFolder, "GF*"))
    files_groups = [x for x in range(n_parallel)]
    for i in range(n_parallel):
        files_groups[i] = " ".join(files[math.floor(i / n_parallel * len(files)):math.floor((i + 1) / n_parallel * len(files))])
    return files_groups


def _subprocess(cmd):
    cmd = [x for x in cmd.split(" ") if x != ""]
    subprocess.run(cmd)


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
              [sg.Text('Parallel Number', font=("Helvetica", 25),
                       pad=((10, 50), (50, 50))),
               sg.Slider(range=(1, 8), default_value=2, size=(20, 15),
                         orientation='horizontal', font=('Helvetica', 25),
                         pad=((30, 50), (0, 50)))],
              [sg.Checkbox('6S', default=False, font=("Helvetica", 25),
                           pad=((10, 25), (0, 35))),
               sg.Checkbox('Pansharpening', default=True, font=("Helvetica", 25),
                           pad=((10, 25), (0, 35))),
               sg.Checkbox('Build Pyramid', default=True, font=("Helvetica", 25),
                           pad=((10, 25), (0, 35)))],
              [sg.Button('OK', font=("Helvetica", 25),
                         pad=((275, 50), (0, 0))),
               sg.Button('CANCEL', font=("Helvetica", 25),
                         pad=((30, 50), (0, 0)))]]

    # Create the Window
    window = sg.Window('GaoFen Preprocess', layout, default_element_size=(30, 1), size=(1000, 900))

    # Get the "values" of the inputs
    event, values = window.read(close=True)
    if event == sg.WIN_CLOSED or event == 'CANCEL':  # if user closes window or clicks CACNCEL
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
    parser.add_argument('--n_parallel', dest='n_parallel',
                        help='Number of processes, e.g., 3 means run 3 files in parallel',
                        type=int, default=3)
    parser.add_argument('--cache', dest='cache',
                        help='Set GDAL raster block cache size, may speed up processing with higher percentage, default is 5% of usable physical RAM',
                        type=str, default='5%')
    args = parser.parse_args()

    # If there is no input from console, open graphic user interface
    if args.InputFolder and args.OutputFolder:
        main(args.InputFolder, args.OutputFolder, args.DEM, args._6S, args.pansharpen, args.pyramid, args.n_parallel, args.cache)
    else:
        InputFolder, OutputFolder, DEM, n_parallel, _6S, pansharpen, pyramid = GUI()
        DEM = DEM if DEM else os.path.join(os.path.dirname(__file__), 'data', 'GMTED2km.tif')
        main(InputFolder, OutputFolder, DEM, _6S, pansharpen, pyramid, int(n_parallel), "5%")
