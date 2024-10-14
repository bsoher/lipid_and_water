import time
import math
import torch
import os
import numpy as np
from BaseTorchLayer import *
from ReadBaselines import *
from utils import load_simulated_data, conver_to_TD, convert_to_FD
import copy
import argparse

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

def fit_specta(x, y, y_trun, prior, net_params, epsi):
    max_patience = net_params["max_patience"]
    list_loss = []
    list_time = []
    start_total = time.time()
    delta = net_params["delta"]
    x_new = x.t().contiguous()
    x_newnew = x_new[0:1,:].contiguous()
    for i in range(0,net_params["nro_of_runs"]):
        start = time.time()
        best_model = torch.ones(1, device=net_params["device"])

        y_full = y[:,i*net_params["number_of_spectra"]:(i+1)*net_params["number_of_spectra"]].t().to(torch.complex64).contiguous()
        y_trunc = y_trun[:, i*net_params["number_of_spectra"]:(i+1)*net_params["number_of_spectra"]].t().to(torch.complex64).contiguous()

        patience = max_patience
        model = Spectra_model(prior, 
                              net_params["number_of_spectra"], 
                              net_params["device"], 
                              x_new, 
                              net_params["input_size"], 
                              initial_offset=net_params["offset_corr"][i*net_params["number_of_spectra"]:(i+1)*net_params["number_of_spectra"]], 
                              use_td=net_params["use_td"],
                              truncate=net_params["truncate_index"]).to(net_params["device"])
        
        criterion = torch.nn.MSELoss(reduction="mean")
        optimizer = torch.optim.Rprop(model.parameters(), lr=net_params["lr"])

        original_max_limit = net_params["limit_max"]
        loss = None
        trunc_patience = max_patience
        best_loss = 99999999.9*torch.ones(1).to(net_params["device"])

        if net_params["use_td"]==False and net_params["truncate"]==True:
            y_out = y_trunc
        else:
            y_out = y_full

        for t in range(net_params["epochs_nro"]):
            def closure_test():
                nonlocal best_loss, patience, best_model, y_out, x_newnew, delta, optimizer, criterion, model, trunc_patience, loss, y_full
                optimizer.zero_grad()
                output = model(x_newnew, trunc=net_params["truncate"])
                if net_params["truncate"] and net_params["use_td"]:
                    net_params["limit_max"] = net_params["truncate_index"]
                else:
                    net_params["limit_max"] = original_max_limit
                if net_params["remove_offset"]:
                    output = output - torch.mean(output[:,net_params["noise_limit_min"]:net_params["noise_limit_max"]], dim=1)[:,None]
                loss = criterion(torch.view_as_real(output[:, net_params["limit_min"]:net_params["limit_max"]]), torch.view_as_real(y_out[:, net_params["limit_min"]:net_params["limit_max"]]))
                if 100*(best_loss - loss)/best_loss > delta:
                    patience = max_patience
                    trunc_patience = int(max_patience)
                    best_loss = loss
                    best_model = copy.deepcopy(model.state_dict())
                else:
                    trunc_patience = trunc_patience - 1
                    patience = patience - 1
                if trunc_patience==0 and net_params["truncate"]:
                    net_params["truncate"] = False
                    y_out = y_full
                    best_loss = 99999999.9*torch.ones(1).to(net_params["device"])
                    patience = max_patience
                if patience==0:
                    print("Final epoch: ",t)
                    raise KeyboardInterrupt
                if t%10==0:
                    print("Epoch: ", t, "loss: ", loss.detach().cpu().numpy())
                loss.backward()
                return loss
            try:
                optimizer.step(closure_test)
            except KeyboardInterrupt:
                break

        stop = time.time()
        list_loss.append(best_loss.cpu().detach().numpy())
        list_time.append(stop-start)

    final_total = time.time()
    print("Total Fitting time: ", final_total-start_total)
    model.load_state_dict(best_model)
    return model, final_total-start_total, best_loss.cpu().detach().numpy(), x_newnew


def correct_offset(Y, prior, avoid=False):
    model_spectr = prior.show_model(ret=True).real
    def get_real(x):
        return x.real

    start = time.time()

    spectra = model_spectr[470:950]
    window = 20
    spectr = np.vectorize(get_real)(Y)
    spectr = spectr[:,470-window:950+window]

    corr = []
    for i in range(0, len(spectr)):
        conv = np.correlate(spectr[i], spectra, mode="valid")
        corr.append(np.argmax(np.abs(conv)))

    if avoid:
        return np.zeros(np.array(corr).shape, dtype=np.float32)

    return np.array(corr)-window


def format_dict(parameters, dictio, scale, net_params):
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
        dictio[meta] = torch.tensor(np.array(parameters[names[meta]].values/scale)).to(net_params["device"]).float()[:,None]

    dictio["offset"] = torch.tensor(np.array(parameters["overall_offset"].values)).to(net_params["device"]).float()[:,None]
    dictio["width"] = torch.tensor(np.array(parameters["overall_lorentz_width"].values)).to(net_params["device"]).float()[:,None]
    dictio["phase"] = torch.tensor(np.array(parameters["overall_zero_phase"].values)).to(net_params["device"]).float()[:,None]
    return dictio


def main():
    epsi = False
    parser = argparse.ArgumentParser(description='Argument parser')
    parser.add_argument('--cpu', action='store_true', help='If used, the fitting will be performed in time domain, default is frequency domain.')
    parser.add_argument('-i', '--input', metavar='<input filename>', default="Dataset"+os.sep+"Simulated_dataset.csv", help='filename (.csv) with the spectra dataset. Default: Dataset"+os.sep+"Simulated_dataset.csv')
    parser.add_argument('--size', type=int, help="Amount of Spectra to fit, integer")

    args = parser.parse_args()
    input_file = args.input
    model_name = "semiLaser135BW1200_3T"
    
    prior = tdfdfit_parser(model_name = model_name).read_file()

    dataset_size = math.inf
    if args.size:
        dataset_size = args.size

    net_params = {}
    net_params["input_size"] = 1024
    net_params["epochs_nro"] = 10000
    net_params["remove_offset"] = True
    net_params["truncate"] = True
    net_params["truncate_index"] = 256
    net_params["use_td"] = False

    net_params["device"] = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if args.cpu:
        net_params["device"] = torch.device("cpu")

    net_params["lr"] = 2

    #Values for stop criteria
    net_params["max_patience"] = 20
    net_params["delta"] = 0.01

    #Range to consider in loss computations.
    net_params["limit_max"] = 1024
    net_params["limit_min"] = 0

    #Parameters used to substract offset.
    net_params["noise_limit_min"] = 0
    net_params["noise_limit_max"] = 200

    start = time.time()
    parameters, Y = load_simulated_data(input_file, parameters_nro=13, dataset_size=dataset_size)

    dataset_size = len(Y)
 
    corr = correct_offset(Y, prior, avoid=False)
    net_params["offset_corr"] = torch.from_numpy(corr)
    net_params["number_of_spectra"] = min(dataset_size, dataset_size)
    net_params["nro_of_runs"] = math.ceil(dataset_size/net_params["number_of_spectra"])
    x, y, y_trun = process_spectra(Y, net_params)

    start = time.time()

    model, time_aux, loss_aux, x_newnew = fit_specta(x, y, y_trun, prior, net_params, epsi)

    time_tensorfit = time.time()-start

    msg = "|    time to fit (TensorFit) {} spectrum: {} seconds.    |".format(dataset_size, round(time_tensorfit, 3))
    print("\n\n")
    print('    '+'-'*len(msg))
    print('    '+msg)
    print('    '+'-'*len(msg))

    print("\nRunning TDFDFit...\n")
    from run_tdfdfit import run_tdfdfit
    time_tdfdfit = run_tdfdfit(Y, "Models"+os.sep+"semiLaser135BW1200_3T"+os.sep+"MasterTDFDFitData"+os.sep+"semiLaser135BW1200_3T_master.par", True)

    print("Total Speed-up:", int(time_tdfdfit/time_tensorfit),"x")


if __name__ == "__main__":
    main()