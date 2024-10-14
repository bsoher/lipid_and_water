import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import math
import torch
from tqdm import tqdm
import glob


def flip_spectra(spectr):
	return np.flip(np.conj(spectr))


def convert_to_FD(spectra):
	return np.fft.fft(spectra)


def conver_to_TD(spectra):
	return np.fft.ifft(spectra)


def load_simulated_data(filename, parameters_nro=13, conjugate=False, dataset_size=math.inf):
	df = pd.read_csv(filename, engine="pyarrow")
	if 'Unnamed: 0' in df.columns:
		df.drop(['Unnamed: 0'], axis=1, inplace=True)
	if 'voxel_ID' in df.columns:
		df.drop(['voxel_ID'], axis=1, inplace=True)
	predictions = []
	y = []
	for i in range(0, min(df.shape[0], dataset_size)):
		y_real = df.loc[i].values[parameters_nro:1024+parameters_nro]
		y_imag = df.loc[i].values[parameters_nro+1024:2048+parameters_nro]
		if conjugate:
			y.append(y_real - 1j*y_imag)
		else:
			y.append(y_real + 1j*y_imag)
	parameters = df.iloc[0:min(dataset_size, df.shape[0]), 0:parameters_nro]
	return parameters, y


def compute_crlb_optimized(model, x_newnew, data):
    Y = model(x_newnew, trunc=False)
    data = data.t().to(torch.complex64)
    states = model.state_dict()

    names_list = np.array([values for values in states])
    states = np.array([states[values].cpu().detach().numpy().squeeze() for values in states])

    # Calculate gradients
    gradie = gradients(model, Y)
    gradie_cpu = gradie.detach()

    crlb_list = []
    for i in tqdm(range(Y.shape[0])):
        sig2 = torch.var(data[i][0:150]).detach()
        gra = gradie[:, i, :].squeeze().t()
        grad2 = torch.matmul(gra.conj().t(), gra)
        inv = torch.linalg.inv(grad2.real/sig2)
        crlb = torch.sqrt(torch.abs(torch.diag(inv)))
        crlb_list.append(100 * np.abs(crlb.cpu().detach().numpy().real / states[:,i].squeeze()))

    return crlb_list, names_list

def gradients(model, Y):
    batch_size = Y.shape[1]
    full_vectors = Y.sum(dim=0)
    gradients = []
    for j in tqdm(range(batch_size)):
        full_vector = full_vectors[j]
        full_vector_complex = 1j*full_vectors[j]
        gra_re = torch.autograd.grad(full_vector, model.parameters(), retain_graph=True, create_graph=True)
        gra_im = torch.autograd.grad(full_vector_complex, model.parameters(), retain_graph=True, create_graph=True)
        gradient_re = torch.stack(gra_re, dim=0).detach().squeeze()[:, :, None]
        gradient_im = torch.stack(gra_im, dim=0).detach().squeeze()[:, :, None]
        gradient = gradient_re+1j*gradient_im
        gradients.append(gradient)

    return torch.cat(gradients, dim=2)


def correct_offset(Y, prior, avoid=False):
    model_spectr = prior.show_model(ret=True).real
    def get_real(x):
        return x.real

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


def read_tdfdfit_fit(filenames):
	df_results = pd.DataFrame(columns=["id", "overall_offset","overall_gaussw","overall_zero_phase","overall_lorentz_width","num_pat_Choline(Cho)_area","num_pat_Glutamate(Glu)_area","num_pat_Lactate(Lac)_area","num_pat_Glutamine_noNH2(Gln)_area","num_pat_creatine_3T_TE135_BW1200_47ppm_area","num_pat_Aspartate(Asp)_area","num_pat_NAcetylAspartate(NAA)_area","num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area","num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area"])
	for filename in filenames:
		matrix = []
		file = open(filename, "r")
		for lines in file.readlines():
			matrix.append(lines.split())
		df_aux = pd.DataFrame(data = [[int(os.path.basename(filename).replace(".bin.mat","").replace("test_","")), float(matrix[6][1]), 0.0, float(matrix[6][3]), float(matrix[6][4]),  float(matrix[0][0]), float(matrix[1][0]), float(matrix[2][0]), float(matrix[3][0]), float(matrix[4][0]), float(matrix[5][0]), float(matrix[6][0]), float(matrix[7][0]), float(matrix[8][0])]])
		df_aux.columns = df_results.columns
		df_results = pd.concat([df_results, df_aux]).reset_index(drop=True)
	df_results = df_results.sort_values(by=["id"]).reset_index(drop=True)
	return df_results


def clean_outliers(df_tensorfit, df_tdfdfit, params, ignore=False):
	if not ignore:
		#df_tdfdfit = df_tdfdfit[(df_tdfdfit.overall_offset>254) & (df_tdfdfit.overall_offset<310)]
		df_tdfdfit = df_tdfdfit[((df_tdfdfit.overall_offset-params.overall_offset)<=10) & ((df_tdfdfit.overall_offset-params.overall_offset) >= -10)]
		df_tensorfit = df_tensorfit[df_tensorfit.index.isin(df_tdfdfit.index)]
		df_tensorfit = df_tensorfit[((df_tensorfit.overall_offset-params.overall_offset)<=10) & ((df_tensorfit.overall_offset-params.overall_offset) >= -10)]
		params = params[params.index.isin(df_tdfdfit.index)]
		params = params[params.index.isin(df_tensorfit.index)]
		df_tdfdfit = df_tdfdfit[df_tdfdfit.index.isin(df_tensorfit.index)]

	return df_tensorfit.reset_index(drop=True), df_tdfdfit.reset_index(drop=True), params.reset_index(drop=True)


def compute_quest_error(parameters, filename, netfit):
	#Order:
	#NACETYLASPARTATE(NAA)#
	#ASPARTATE(ASP)#
	#CHOLINE(CHO)#
	#CREATINE_3T_TE135_BW1200_47PPM#
	#GLUTAMATE(GLU)#
	#GLUTAMINE_NONH2(GLN)#
	#LACTATE(LAC)#
	#MYO_INOSITOL_3T_TE135_BW1200_47PPM#
	#NACETYLASPARTATE(NAA)-2_6PPMMULTIPLET
	#TODO: Change order according to .results file
	dict_names = ['num_pat_NAcetylAspartate(NAA)_area',
				  'num_pat_Aspartate(Asp)_area', 
				  'num_pat_Choline(Cho)_area', 
				  'num_pat_creatine_3T_TE135_BW1200_47ppm_area', 
				  'num_pat_Glutamate(Glu)_area', 
				  'num_pat_Glutamine_noNH2(Gln)_area', 
				  'num_pat_Lactate(Lac)_area', 
				  'num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area',
				  'num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area'
				  ]
	try:
		lista = np.loadtxt(open(filename, "rb"), delimiter=",")
	except:
		print("QUEST results csv ({}) not found.".format(filename))
		return None, None, None
	df_results = pd.DataFrame(columns=["id", "overall_offset","overall_gaussw","overall_zero_phase","overall_lorentz_width","num_pat_Choline(Cho)_area","num_pat_Glutamate(Glu)_area","num_pat_Lactate(Lac)_area","num_pat_Glutamine_noNH2(Gln)_area","num_pat_creatine_3T_TE135_BW1200_47ppm_area","num_pat_Aspartate(Asp)_area","num_pat_NAcetylAspartate(NAA)_area","num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area","num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area"])
	new_list = []
	new_params = []

	#clean QUEST
	std_cho = np.std(lista[:,2]-parameters[dict_names[2]].values)
	std_naa = np.std(lista[:,0]-parameters[dict_names[0]].values)
	std_cr = np.std(lista[:,3]-parameters[dict_names[3]].values)
	outlier_thr = 4
	counter = 0
	lista_clean = []
	for i in range(0, len(lista)):
		if np.abs(lista[i, 2]-parameters[dict_names[2]].values[i])>outlier_thr*std_cho or np.abs(lista[i, 3]-parameters[dict_names[3]].values[i])>outlier_thr*std_cr or np.abs(lista[i, 0]-parameters[dict_names[0]].values[i])>outlier_thr*std_naa:
			counter = counter + 1
			lista_clean.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
		else:
			lista_clean.append(lista[i,:])
	lista_clean = np.array(lista_clean)

	print("dropped: ",counter)
	for m in range(0,9):
		err = 100*np.nanmean(np.abs((lista_clean[:,m] - parameters[dict_names[m]].values)/(parameters[dict_names[m]].values)))
		#set TensorFit parameter to 0 if negative, to have the same behaviour as QUEST.
		paramiter = np.where(netfit[dict_names[m]].values < 0, 0, netfit[dict_names[m]].values)
		err_tensorfit = 100*np.nanmean(np.abs(paramiter - parameters[dict_names[m]].values)/(parameters[dict_names[m]].values))
		print(dict_names[m])
		print(str(err).replace(".",","))
		print(str(err_tensorfit).replace(".",","))

	return dict_names, err_quest, err_tensorfit