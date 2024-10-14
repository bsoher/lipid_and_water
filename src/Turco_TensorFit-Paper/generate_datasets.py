import pandas as pd
from ReadBaselines import *
import json
from utils import flip_spectra
from tqdm import tqdm
from copy import deepcopy
import matplotlib.pyplot as plt

def read_dict(use_gauss):
	file = "distribution_3T.json"
	f = open(file)
	data = json.load(f)
	if use_gauss:
		data.update({"width_gaus":[1,7]})
	return data

def update_prior(prior, random_dict, original_prior, use_gauss):
	prior = deepcopy(original_prior)
	prior.update_param_value("area", "Choline(Cho)", random_dict["Choline(Cho)"])
	prior.update_param_value("area", "Glutamate(Glu)", random_dict["Glutamate(Glu)"])
	prior.update_param_value("area", "Lactate(Lac)", random_dict["Lactate(Lac)"])
	prior.update_param_value("area", "Glutamine_noNH2(Gln)", random_dict["Glutamine_noNH2(Gln)"])
	prior.update_param_value("area", "creatine_3T_TE135_BW1200_47ppm", random_dict["Creatine_3T_TE135_BW1200_47ppm"])
	prior.update_param_value("area", "Aspartate(Asp)", random_dict["Aspartate(Asp)"])
	prior.update_param_value("area", "NAcetylAspartate(NAA)", random_dict["NAcetylAspartate(NAA)"])
	prior.update_param_value("area", "myo_inositol_3T_TE135_BW1200_47ppm", random_dict["myo_inositol_3T_TE135_BW1200_47ppm"])
	prior.update_param_value("area", "NAcetylAspartate(NAA)-2_6ppmmultiplet", random_dict["NAcetylAspartate(NAA)-2_6ppmmultiplet"])
	if use_gauss:
		prior.update_commons(random_dict["offset"], random_dict["phase"], random_dict["width_loren"], random_dict["width_gaus"])
	else:
		prior.update_commons(random_dict["offset"], random_dict["phase"], random_dict["width_loren"], 0)

	return prior


def get_random_params(data, prior, original_prior, use_gauss=False):
	random_dict = {}
	for item, limits in data.items():
		random_dict[item] = np.random.uniform(limits[0],limits[1])

	prior = update_prior(prior, random_dict, original_prior, use_gauss)

	return prior, random_dict


def get_noisy_spectra(base_spectra, spectra_params, noise):
	spectra = np.concatenate([base_spectra.real,base_spectra.imag]) + noise*np.random.normal(size=(2*len(base_spectra)))
	spectra = np.concatenate([spectra_params, spectra])
	return spectra


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Generates simulated spectra and export it for the rest of codes to use')
	parser.add_argument('--include_gauss', action='store_true', help='Include Gaussian decay in simulation apart from the Lorentzian decay.')
	parser.add_argument('-o', '--output', default="Dataset\\Simulated_dataset.csv", help='Output filename for the simulation. Default: Dataset\\Simulated_dataset.csv')
	parser.add_argument('--size', default="1000", help="Amount of spectra to simulated. default: 10000", type=int)
	args = parser.parse_args()

	use_gauss = False
	if args.include_gauss:
		use_gauss = True

	model = "semiLaser135BW1200_3T"

	number_of_spectra = int(args.size)

	dictio = read_dict(use_gauss)

	file_manager = tdfdfit_parser(model_name = model)

	prior = file_manager.read_file()

	original_prior = deepcopy(prior)
	SNR = [50, 25, 10, 5]

	names = ["overall_offset","overall_gaussw","overall_zero_phase","overall_lorentz_width","num_pat_Choline(Cho)_area","num_pat_Glutamate(Glu)_area","num_pat_Lactate(Lac)_area","num_pat_Glutamine_noNH2(Gln)_area",
			 "num_pat_creatine_3T_TE135_BW1200_47ppm_area","num_pat_Aspartate(Asp)_area","num_pat_NAcetylAspartate(NAA)_area","num_pat_myo_inositol_3T_TE135_BW1200_47ppm_area","num_pat_NAcetylAspartate(NAA)-2_6ppmmultiplet_area"]
	for i in range(0,1024):
		name = "absFD[{}]".format(i)
		names.append(name)
	for i in range(0,1024):
		name = "dispFD[{}]".format(i)
		names.append(name)

	snr_abs = prior.get_snr_base(SNR)
	orig_spectra = prior.show_model(ret=True, snr=0)
	data = []
	for i in tqdm(range(0, number_of_spectra)):
		prior, random_dict = get_random_params(dictio, prior, original_prior, use_gauss)
		if use_gauss:
			spectra_params = [random_dict["offset"], random_dict["width_gaus"], random_dict["phase"], random_dict["width_loren"], random_dict["Choline(Cho)"], random_dict["Glutamate(Glu)"], random_dict["Lactate(Lac)"], random_dict["Glutamine_noNH2(Gln)"],
							  random_dict["Creatine_3T_TE135_BW1200_47ppm"], random_dict["Aspartate(Asp)"], random_dict["NAcetylAspartate(NAA)"], random_dict["myo_inositol_3T_TE135_BW1200_47ppm"], random_dict["NAcetylAspartate(NAA)-2_6ppmmultiplet"]]
		else:
			spectra_params = [random_dict["offset"], 0, random_dict["phase"], random_dict["width_loren"], random_dict["Choline(Cho)"], random_dict["Glutamate(Glu)"], random_dict["Lactate(Lac)"], random_dict["Glutamine_noNH2(Gln)"],
							  random_dict["Creatine_3T_TE135_BW1200_47ppm"], random_dict["Aspartate(Asp)"], random_dict["NAcetylAspartate(NAA)"], random_dict["myo_inositol_3T_TE135_BW1200_47ppm"], random_dict["NAcetylAspartate(NAA)-2_6ppmmultiplet"]]

		spectra = prior.show_model(ret=True, snr=0)
		noise = np.random.choice(snr_abs)
		new_spectra = get_noisy_spectra(spectra, spectra_params, noise)
		#prior.show_basis(show_all=True)
		#plt.plot(new_spectra[12:12+1024])
		#plt.show()
		data.append(new_spectra)
	df_total = pd.DataFrame(columns = names, data=np.squeeze(data))


	folder_path = os.path.dirname(args.output)
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)

	filename, file_extension = os.path.splitext(args.output)
	if file_extension.lower() != '.csv':
		args.output = args.output + ".csv"

	df_total.to_csv(args.output, index=False)