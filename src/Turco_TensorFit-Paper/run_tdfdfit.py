import pandas as pd
import numpy as np
from utils import load_simulated_data
import subprocess
import multiprocessing
from multiprocessing.pool import ThreadPool
from subprocess import DEVNULL, call
import glob
import os
import time
import math
import argparse

def remove_files(folder="tdfdfit_tmp", alles=False):
    keep = [x for x in glob.glob(f"{folder}/*.bin.mat")]
    keep2 = [x for x in glob.glob(f"{folder}/*.bin.fit")]
    remove = [x for x in glob.glob(f"{folder}/test_*")]
    for name in remove:
        if alles or (name not in keep and name not in keep2):
            os.remove(name)

def run(cmd):
    return cmd, call(cmd, stdout=DEVNULL, stderr=DEVNULL)

def run_parallel_tdfdfit(prior_knowledge, fit_only=None, folder="resuts"):
    names = [os.path.basename(x) for x in glob.glob(f"{folder}/test_*.bin")]
    if fit_only is not None:
        names = [f"test_{fit_only}.bin"]
    #TODO: replace for 
    tdfdfit = "tdfdfit2017.exe"
    folder = f"{folder}"+os.sep
    cmd_list = [[tdfdfit, folder, filename, prior_knowledge] for filename in names]
    for cmd, rc in ThreadPool(multiprocessing.cpu_count()-1).imap_unordered(run, cmd_list):
        if rc != 0:
            print("Failed!")

def main():
    parser = argparse.ArgumentParser(description='Argument parser')
    parser.add_argument('-i', '--input', metavar='<input filename>', default="Dataset"+os.sep+"Simulated_dataset.csv", help='filename (.csv) with the spectra dataset. Default: Dataset'+os.sep+'Simulated_dataset.csv')
    parser.add_argument('--size', type=int, help="Amount of Spectra to fit, integer")

    args = parser.parse_args()

    prior_knowledge = "Models"+os.sep+"semiLaser135BW1200_3T"+os.sep+"MasterTDFDFitData"+os.sep+"semiLaser135BW1200_3T_master.par"

    dataset_size = math.inf
    if args.size:
        dataset_size = args.size

    filename_in = args.input
    parameters, Y = load_simulated_data(filename_in, parameters_nro=13, dataset_size=dataset_size)
    run_tdfdfit(Y, prior_knowledge)


def run_tdfdfit(Y, prior_knowledge, clear_all=True):
    dataset_size = len(Y)
    folder = "tdfdfit_tmp"
    Voxels = [i for i in range(0, len(Y))]

    if not os.path.exists(folder):
        os.makedirs(folder)
    remove_files(folder, alles=True)
    for i in range(0, len(Y)):
        Y[i] = Y[i] - np.mean(Y[i][0:200])
        Y_out = np.concatenate(([Y[i].real], [Y[i].imag]), axis=0).T.reshape(2048)
        Y_out.astype('float32').tofile(f"{folder}/test_{Voxels[i]}.bin")
    start = time.time()
    run_parallel_tdfdfit(prior_knowledge, fit_only=None, folder=folder)
    end = time.time()
    remove_files(folder, alles=clear_all)

    if clear_all:
        msg = "|    time to fit (TDFDFit) {} spectrum: {} seconds.    |".format(dataset_size, round(end-start, 3))
        print("\n\n")
        print('    '+'-'*len(msg))
        print('    '+msg)
        print('    '+'-'*len(msg))
    return end-start

if __name__=="__main__":
    main()