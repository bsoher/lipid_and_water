import time
import math
import torch
import matplotlib.pyplot as plt
import numpy as np
import json
from BaseTorchLayer import *
from ReadBaselines import *
from utils import load_simulated_data, conver_to_TD, convert_to_FD, clean_outliers, compute_quest_error, read_tdfdfit_fit
import copy
from tqdm import tqdm
import pandas as pd
import glob
from torch.autograd import grad

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

def fit_specta(x, y, y_trun, prior, net_params):
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
        print("net param", net_params["truncate"])
        
        model = Spectra_model(prior, 
                              net_params["number_of_spectra"], 
                              net_params["device"], 
                              x_new, 
                              net_params["input_size"], 
                              initial_offset=net_params["offset_corr"][i*net_params["number_of_spectra"]:(i+1)*net_params["number_of_spectra"]], 
                              use_td=net_params["use_td"],
                              truncate=net_params["truncate_index"],
                              add_gaussian = net_params["use_gaussian"]).to(net_params["device"])

        criterion = torch.nn.MSELoss(reduction="mean")
        optimizer = torch.optim.Rprop(model.parameters(), lr=net_params["lr"])
        #optimizer = torch.optim.Adam(model.parameters(), lr=net_params["lr"])
        #import torch_optimizer as optim_2
        #optimizer = optim_2.Adahessian(model.parameters(), lr=net_params["lr"])
        #optimizer = torch.optim.SGD(model.parameters(), lr=10)
        #optimizer = torch.optim.LBFGS(model.parameters(), lr=0.1, max_iter=20, line_search_fn="strong_wolfe")
        if net_params["use_decay"]:
            schedule_lr = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.8, patience=5)

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
                    print("Truncated converged. Fitting full spectra")
                    net_params["truncate"] = False
                    y_out = y_full
                    best_loss = 99999999.9*torch.ones(1).to(net_params["device"])
                    patience = max_patience
                if patience==0:
                    print("Finish at epoch: ",t)
                    raise KeyboardInterrupt
                if t%10==0:
                    print("epoch: ",t, ".Loss: ", np.round(loss.detach().cpu().numpy(),3))
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


def save_model(model, filename, use_gauss):
    dictio = {
        'offset':'overall_offset',
        'phase':'overall_zero_phase',
        'width':'overall_lorentz_width',
        'Choline(Cho)':'num_pat_Choline(Cho)_area',
        'Glutamate(Glu)':'num_pat_Glutamate(Glu)_area',
        'Lactate(Lac)':'num_pat_Lactate(Lac)_area',
        'Glutamine_noNH2(Gln)':'num_pat_Glutamine_noNH2(Gln)_area',
        'creatine_3T_TE135_BW1200_47ppm':'num_pat_creatine_3T_TE135_BW1200_47ppm_area',
        'Aspartate(Asp)':'num_pat_Aspartate(Asp)_area',
        'NAcetylAspartate(NAA)':'num_pat_NAcetylAspartate(NAA)_area',
        'myo_inositol_3T_TE135_BW1200_47ppm':'num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area',
        'NAcetylAspartate(NAA)-2_6ppmmultiplet':'num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area'
    }
    if use_gauss:
        dictio['gausw']='overall_gaussw'

    fit_params = model.get_params()
    
    df = pd.DataFrame()
    for param in dictio.items():
        fitted_params = fit_params['{}'.format(param[0])]
        df[param[1]] = fitted_params
    return df

def fit_statistics(parameters, model, use_gaussian):
    dictio = {
        'offset':'overall_offset',
        'phase':'overall_zero_phase',
        'width':'overall_lorentz_width',
        'Choline(Cho)':'num_pat_Choline(Cho)_area',
        'Glutamate(Glu)':'num_pat_Glutamate(Glu)_area',
        'Lactate(Lac)':'num_pat_Lactate(Lac)_area',
        'Glutamine_noNH2(Gln)':'num_pat_Glutamine_noNH2(Gln)_area',
        'creatine_3T_TE135_BW1200_47ppm':'num_pat_creatine_3T_TE135_BW1200_47ppm_area',
        'Aspartate(Asp)':'num_pat_Aspartate(Asp)_area',
        'NAcetylAspartate(NAA)':'num_pat_NAcetylAspartate(NAA)_area',
        'myo_inositol_3T_TE135_BW1200_47ppm':'num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area',
        'NAcetylAspartate(NAA)-2_6ppmmultiplet':'num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area'
    }
    if use_gaussian:
        dictio.udpate({'gausw': 'overall_gaussw'})

    fit_params = model.get_params()
    error = []
    names = []
    scatter_size = 4
    i = 0
    ret_list = []
    from sklearn.metrics import r2_score
    for param in dictio.items():
        target_param = parameters.loc[:, param[1]].values
        fitted_param = fit_params['{}'.format(param[0])]
        names.append(param[0])
        ret_list.append(100*np.abs(target_param-fitted_param)/np.abs(target_param))
        error.append(100*np.mean(np.abs(target_param-fitted_param)/np.abs(target_param)))

        i = i+1
    for i in range(0, len(error)):
        print(names[i], ":", str(error[i]).replace(".",","))

    return error, names, ret_list


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
    parser = argparse.ArgumentParser(description='Argument parser')
    parser.add_argument('--include_gauss', action='store_true', help='Determine if independiently of the model, it will use also gaussian decay')
    parser.add_argument('-t', '--time', action='store_true', help='If used, the fitting will be performed in time domain, default is frequency domain.')
    parser.add_argument('-o', '--output', metavar='<output filename>', help='Filename with path were the results csv will be saved.')
    parser.add_argument('-i', '--input', metavar='<input filename>', help='filename (.csv) with the spectra dataset. Default: Dataset"+os.sep+"Simulated_dataset.csv', default="Dataset"+os.sep+"Simulated_dataset.csv")
    parser.add_argument('--size', metavar='<number of spectra>', help='Number of spectra to fit. Default, the size of the input .csv file')
    
    args = parser.parse_args()

    use_gauss = False
    if args.include_gauss:
        use_gauss = True

    if not args.output:
        output_file = None

    input_file = args.input
    model_name = "semiLaser135BW1200_3T"

    if args.output is not None:
        folder_path = os.path.dirname(args.output)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        filename, file_extension = os.path.splitext(args.output)
        if file_extension.lower() != '.csv':
            args.output = args.output + ".csv"

    prior = tdfdfit_parser(model_name = model_name).read_file()

    input_size = 1024

    dataset_size = math.inf
    if args.size:
        dataset_size = int(args.size)

    net_params = {}
    net_params["input_size"] = 1024
    net_params["epochs_nro"] = 10000
    net_params["remove_offset"] = True
    net_params["use_gaussian"] = use_gauss
    net_params["truncate_index"] = 1024
    net_params["truncate"] = True
    if net_params["truncate"]:
        net_params["truncate_index"] = 256
    net_params["max_patience"] = 20
    net_params["device"] = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    net_params["lr"] = 2
    net_params["delta"] = 0.001
    net_params["use_decay"] = False
    net_params["use_td"] = args.time
    if net_params["use_td"]:
        net_params["limit_max"] = 1024
        net_params["limit_min"] = 0
        net_params["noise_limit_min"] = 800
        net_params["noise_limit_max"] = 1024
    else:
        net_params["limit_max"] = 1024
        net_params["limit_min"] = 0
        net_params["noise_limit_min"] = 0
        net_params["noise_limit_max"] = 200

    parameters, Y = load_simulated_data(input_file, parameters_nro=13, dataset_size=dataset_size)

    total_size = len(Y)
    corr = correct_offset(Y, prior, avoid=False)
    net_params["offset_corr"] = torch.from_numpy(corr)
    net_params["number_of_spectra"] = min(total_size, total_size)
    net_params["nro_of_runs"] = math.ceil(total_size/net_params["number_of_spectra"])
    x, y, y_trun = process_spectra(Y, net_params)

    model, time_aux, loss_aux, x_newnew = fit_specta(x, y, y_trun, prior, net_params)

    print("\nRunning TDFDFit...\n")
    from run_tdfdfit import run_tdfdfit
    run_tdfdfit(Y, "Models"+os.sep+"semiLaser135BW1200_3T"+os.sep+"MasterTDFDFitData"+os.sep+"semiLaser135BW1200_3T_master.par", False)

    from compute_average_groundtruth_crlb import run as run_crlb
    print("\nComputing CRLB...\n")
    crlb = run_crlb(input_file, dataset_size, prior)

    #compute_quest_error(parameters, filename, tensorfit)
    df_tensorfit = save_model(model, None, False)

    names = [x for x in glob.glob(f"tdfdfit_tmp"+os.sep+"*.bin.mat")]
    df_tdfdfit = read_tdfdfit_fit(names)
    df_tensorfit, df_tdfdfit, parameters = clean_outliers(df_tensorfit, df_tdfdfit, parameters)

    error_net = []
    error_tdfdfit = []

    for param in parameters.columns:
        #There is no gaussian on this dataset
        if param!='overall_gaussw':
            target_param = parameters.loc[:, param].values
            #When TDFDFit fits, the end-result have a factor 1024 extra in the areas parameters. thus, the next /1024
            scale = 1
            if "_area" in param:
                scale = 1024

            fitted_param = df_tdfdfit.loc[:, param]/scale
            NetFit_param = df_tensorfit.loc[:, param]

            error_net_aux = np.mean(np.abs((target_param - NetFit_param)/(target_param)))
            error_tdfdfit_aux = np.mean(np.abs((target_param - fitted_param)/target_param))

            error_tdfdfit.append(error_tdfdfit_aux*100)
            error_net.append(error_net_aux*100)

    i = 3
    for param in parameters.columns[4:]:
        print(param, ":\n\tAveraged Percentual Error TensorFit: ", round(error_net[i], 3), "%\n\tPercentual Error TDFDFit: ", round(error_tdfdfit[i]), "%\n\tCRLB for groundtruth: ", round(crlb[i], 3), "%")
        i = i + 1

if __name__ == "__main__":
    main()