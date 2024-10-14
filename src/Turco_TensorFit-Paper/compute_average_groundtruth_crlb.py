import math
import torch
import matplotlib.pyplot as plt
import numpy as np
import json
from BaseTorchLayer import *
from ReadBaselines import *
from utils import load_simulated_data, conver_to_TD, convert_to_FD, compute_crlb_optimized, correct_offset
import copy
from tqdm import tqdm
from torch.autograd import grad
import time

torch.backends.cudnn.benchmark=False

def process_spectra(Y, net_params):
    y = np.array(Y)

    if net_params["remove_offset"]:
        for i in range(0, y.shape[0]):
            y[i] = y[i] - y[i][net_params["noise_limit_min"]:net_params["noise_limit_max"]].mean()

    y_trun = np.zeros_like(y)
    if net_params["use_td"]:
        for i in range(0, len(y)):
            y[i] = conver_to_TD(y[i])
            y_trun[i] = y[i]
    elif net_params["truncate"]:
        for i in range(0, len(y)):
            y_trun[i] = conver_to_TD(y[i])
            y_trun[i] = [y_trun[i][j] if j<net_params["truncate_index"] else 0 for j in range(0, len(y_trun[i]))]
            y_trun[i] = convert_to_FD(y_trun[i])

    y = y.T
    y_trun = y_trun.T
    x = torch.linspace(0, net_params["input_size"]-1, net_params["input_size"]).to(net_params["device"])
    x = x.resize(x.shape[0],1).to(net_params["device"])
    x = x.repeat([1, net_params["number_of_spectra"]])

    y_new_trun = torch.from_numpy(y_trun)
    y_new_trun = y_new_trun.to(net_params["device"])
    y_new = torch.from_numpy(y)
    y_new = y_new.to(net_params["device"])
    return x, y_new, y_new_trun
#-----------------------
##Training

def format_dict(parameters, dictio, net_params):
    names = {
        "Choline(Cho)":"num_pat_Choline(Cho)_area",
        "Glutamate(Glu)":"num_pat_Glutamate(Glu)_area",
        "Lactate(Lac)":"num_pat_Lactate(Lac)_area",
        "Glutamine_noNH2(Gln)":"num_pat_Glutamine_noNH2(Gln)_area",
        "creatine_3T_TE135_BW1200_47ppm":"num_pat_creatine_3T_TE135_BW1200_47ppm_area",
        "Aspartate(Asp)":"num_pat_Aspartate(Asp)_area",
        "NAcetylAspartate(NAA)":"num_pat_NAcetylAspartate(NAA)_area",
        "myo_inositol_3T_TE135_BW1200_47ppm":"num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area",
        "NAcetylAspartate(NAA)-2_6ppmmultiplet":"num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area"
    }
    for meta in names.keys():
        dictio[meta] = torch.tensor(np.array(parameters[names[meta]].values)).to(net_params["device"]).float()[:,None]

    dictio["offset"] = torch.tensor(np.array(parameters["overall_offset"].values)).to(net_params["device"]).float()[:,None]
    dictio["width"] = torch.tensor(np.array(parameters["overall_lorentz_width"].values)).to(net_params["device"]).float()[:,None]
    dictio["phase"] = torch.tensor(np.array(parameters["overall_zero_phase"].values)).to(net_params["device"]).float()[:,None]
    return dictio


def main():
    parser = argparse.ArgumentParser(description='Argument parser')
    parser.add_argument('-i', '--input', metavar='<input filename>', default="Dataset\\Simulated_dataset.csv", help='filename (.csv) with the spectra dataset. Default: Dataset\\Simulated_dataset.csv')
    parser.add_argument('--size', type=int, help="Amount of Spectra to fit")
    
    args = parser.parse_args()

    input_file = args.input

    model_name = "semiLaser135BW1200_3T"

    #Todo improve this to directly read the full (or relative) path to the model.
    prior = tdfdfit_parser(model_name = model_name).read_file()
    dataset_size = math.inf
    if args.size:
        dataset_size = int(args.size)
    run(input_file, dataset_size, prior)

def run(input_file, dataset_size, prior):
    net_params = {}
    net_params["input_size"] = 1024
    net_params["epochs_nro"] = 10000
    net_params["remove_offset"] = True
    net_params["truncate"] = False
    net_params["truncate_index"] = 1024
    net_params["device"] = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    net_params["use_td"] = False
    net_params["noise_limit_min"] = 0
    net_params["noise_limit_max"] = 200

    parameters, Y = load_simulated_data(input_file, parameters_nro=13, dataset_size=dataset_size)
    
    total_size = len(Y)
    corr = correct_offset(Y, prior, avoid=True)
    net_params["offset_corr"] = torch.from_numpy(corr)
    net_params["number_of_spectra"] = total_size
    x, y, y_trun = process_spectra(Y, net_params)

    x_new = x.t().contiguous()
    x_newnew = x_new[0:1,:].contiguous()

    model = Spectra_model(prior, 
                  net_params["number_of_spectra"], 
                  net_params["device"], 
                  x_new, 
                  net_params["input_size"], 
                  initial_offset=net_params["offset_corr"][0:net_params["number_of_spectra"]], 
                  use_td=net_params["use_td"], 
                  truncate=net_params["truncate_index"]).to(net_params["device"])

    dictio = model.get_tensor_params()
    dictio = format_dict(parameters, dictio, net_params)
    model.set_params(dictio)
    crlb, names_list = compute_crlb_optimized(model, x_newnew, y)

    names_dict = {"common.Dense_offset" : "frequency shift",
                  "common.Dense_width" : "lorentizian width",
                  "common.Dense_phase" : "phase",
                  "metabolites.0.Dense_ampl" : "Choline(Cho)",
                  "metabolites.1.Dense_ampl" : "Glutamate(Glu)",
                  "metabolites.2.Dense_ampl" : "Lactate(Lac)",
                  "metabolites.3.Dense_ampl" : "Glutamine_noNH2(Gln)",
                  "metabolites.4.Dense_ampl" : "creatine_3T_TE135_BW1200_47ppm",
                  "metabolites.5.Dense_ampl" : "Aspartate(Asp)",
                  "metabolites.6.Dense_ampl" : "NAcetylAspartate(NAA)",
                  "metabolites.7.Dense_ampl" : "myo_inositol_3T_TE135_BW1200_47ppm",
                  "metabolites.8.Dense_ampl" : "NAcetylAspartate(NAA)-2_6ppmmultiplet"}

    means = np.mean(crlb, axis=0)
    if __name__=="__main__":
        print("Average Percentual CRLB:")
        [print(names_dict[names_list[i]], ":", round(means[i], 3), "%") for i in range(3, len(names_list))]
    return means

if __name__ == "__main__":
    main()