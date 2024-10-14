import re
import numpy as np
from pathlib import Path
import os
import math
import struct
import utils
import json
import torch
import matplotlib.pyplot as plt
import argparse
win_sep = os.sep

class prior_knowledge_parse():
	def __init__(self, metabolite, prior_knowledge_string, optim):
		self.nro_basis = len(metabolite)
		self.metabolites = metabolite
		self.ndp = int(metabolite.get_ndp())
		self.delimiter = "tdfdfit_prior_knowledge"
		self.prior = []
		self.scale = 1
		if optim:
			self.scale = 1000
		self.parse_prior_knowledge(prior_knowledge_string)
		self.min_area = 0
		self.max_area = 0
		self.offset_var = 0
		self.width_loren_var = 0
		self.phase_var = 0
		self.width_gaus_var = 0
		self.optim = optim
		self.metabolite_shape_final = []

	def __len__(self):
		return self.nro_basis

	def __getitem__(self, index):
		if index >= 0 and index < self.nro_basis:
			return self.prior[index]
		else:
			print(f"Index out of range. Should be in the range [0, {self.nro_basis-1}]")
			return None


	def parse_prior_knowledge(self, prior_knowledge_string):
		if prior_knowledge_string.find(self.delimiter)==-1:
			print("#### ERROR!: prior knowledge is not correctly formated, could not find the delimiter: ", self.delimiter)

		prior_knowledge_string = prior_knowledge_string.replace(self.delimiter, "")
		prior_knowledge_list = prior_knowledge_string.split("\n")
		prior_knowledge_list = [inform.strip() for inform in prior_knowledge_list if inform!=""]
		for i in range(0, len(self.metabolites)):
			self.prior.append(self.parse_metab_info(prior_knowledge_list[8*i:8*(i+1)], i))

		self.update_deltas()

	def get_initial_params(self):
		ret = []
		ret.append(self.width_loren_common)
		ret.append(self.offset_common)
		ret.append(self.phase_common)
		for i in range(0, len(self.prior)):
			ret.append(self.prior[i]["area"]["original_value"])

		if self.gauss_common!=0.0:
			ret.append(self.phase_common)

		return torch.from_numpy(np.array(ret))[:,None]

	def get_fixed_params(self, x=None, numpy=False, get_base=False):
		ret = []
		if x is None:
			len_spectra = self.prior[0]["metabolite_shape"].shape[0]
			x = np.linspace(0, len_spectra-1, len_spectra)
			x = x[:,None]
		if numpy:
			x_aux = torch.Tensor(x)
		else:
			x_aux = x.to("cpu")
		ndp = x_aux.shape[0]

		for i in range(0, len(self.prior)):
			if get_base:
				expon_im = (2*math.pi*self.prior[i]["offset"]["original_value"]*x_aux[:,0]/ndp + self.prior[i]["phase"]["delta_value"]*math.pi/180)
			else:
				expon_im = (2*math.pi*self.prior[i]["offset"]["delta_value"]*x_aux[:,0]/ndp + self.prior[i]["phase"]["delta_value"]*math.pi/180)
	
			expon_real = -self.prior[i]["width_loren"]["delta_value"]*math.pi*x_aux[:,0]/ndp
			
			exponential = torch.exp(expon_real + 1j*expon_im)
			self.metabolite_shape_final.append(torch.view_as_real((exponential.t()*self.prior[i]["metabolite_shape"])[None,:]).cpu().detach().numpy())

		if numpy:
			self.metabolite_shape_final = np.squeeze(np.array(self.metabolite_shape_final))
			return self.metabolite_shape_final[:,:,0] + 1j*self.metabolite_shape_final[:,:,1]
		else:
			return torch.view_as_complex(torch.Tensor(self.metabolite_shape_final))

	def get_absolute_params(self, x, numpy=False):
		ret = []
		if numpy:
			x_aux = torch.Tensor(x)
		else:
			x_aux = x.to("cpu")
		ndp = x_aux.shape[0]

		for i in range(0, len(self.prior)):
			if get_base:
				expon_im = (2*math.pi*self.prior[i]["offset"]["original_value"]*x_aux[:,0]/ndp + self.prior[i]["phase"]["original_value"]*math.pi/180)
				expon_real = -self.prior[i]["width_loren"]["original_value"]*math.pi*x_aux[:,0]/ndp
			else:
				expon_im = (2*math.pi*self.prior[i]["offset"]["original_value"]*x_aux[:,0]/ndp)
				
			exponential = torch.exp(expon_real + 1j*expon_im)
			self.metabolite_shape_final.append(torch.view_as_real((exponential.t()*self.prior[i]["metabolite_shape"])[None,:]).cpu().detach().numpy())

		if numpy:
			self.metabolite_shape_final = np.squeeze(np.array(self.metabolite_shape_final))
			return self.metabolite_shape_final[:,:,0] + 1j*self.metabolite_shape_final[:,:,1]
		else:
			return torch.view_as_complex(torch.Tensor(self.metabolite_shape_final))



	def update_commons(self, offset_common, phase_common, width_loren_common, gauss_common):
		for i in range(0, len(self.prior)):
			if self.prior[i]["is_common"]:
				self.prior[i]["offset"]["original_value"] = offset_common
				self.prior[i]["phase"]["original_value"] = phase_common
				self.prior[i]["width_loren"]["original_value"] = width_loren_common
				self.prior[i]["width_gaus"]["original_value"] = gauss_common
				self.offset_common = self.prior[i]["offset"]["original_value"]
				self.phase_common = self.prior[i]["phase"]["original_value"]
				self.width_loren_common = self.prior[i]["width_loren"]["original_value"]
				self.gauss_common = self.prior[i]["width_gaus"]["original_value"]


	def update_deltas(self):
		for info in self.prior:
			if info["is_common"]:
				self.offset_common = info["offset"]["original_value"]
				self.phase_common = info["phase"]["original_value"]
				self.width_loren_common = info["width_loren"]["original_value"]
				self.gauss_common = info["width_gaus"]["original_value"]

		for i in range(0, len(self.prior)):
			self.prior[i]["offset"]["delta_value"] = self.prior[i]["offset"]["original_value"] - self.offset_common
			self.prior[i]["phase"]["delta_value"] = self.prior[i]["phase"]["original_value"] - self.phase_common
			self.prior[i]["width_loren"]["delta_value"] = self.prior[i]["width_loren"]["original_value"] - self.width_loren_common
			self.prior[i]["width_gaus"]["delta_value"] = self.prior[i]["width_gaus"]["original_value"] - self.gauss_common

	def set_random_limits(min_area, max_area, offset_var, phase_var, width_gaus_var, width_loren_var):
		self.min_area = min_area
		self.max_area = max_area
		self.offset_var = offset_var
		self.phase_var = phase_var
		self.width_gaus_var = width_gaus_var
		self.width_loren_var = width_loren_var


	def random_params(self):
		parameters = {}
		parameters["Areas"] = {}
		parameters["Common"] = {}
		for i in range(0, nro_basis):
			parameters["Areas"][self.metabolites.get_metaboliteName(i)] = self.prior[i]["Area"]["original_value"]*np.random.uniform(self.min_area, self.max_area)

		parameters["Common"]["offset"] = np.random.uniform(self.offset_common-self.offset_var/2, self.offset_common+self.offset_var/2)
		parameters["Common"]["phase"] = np.random.uniform(self.phase_common-self.phase_var/2, self.phase_common+self.phase_var/2)
		parameters["Common"]["width_loren"] = np.random.uniform(self.width_loren_var[0], self.width_loren_var[1])
		parameters["Common"]["width_gaus"] = np.random.uniform(self.width_gaus_var[0], self.width_gaus_var[1])


	def __str__(self):
		for info in self.prior:
			print('-'*50,info["metabolite_name"])
			print("Area:", info["area"]["original_value"]/self.scale)
			print("Delta Offset (Hz):", info["offset"]["delta_value"]*1.243858)
			print("Delta Phase:", info["phase"]["delta_value"])
			print("Delta Loren:", info["width_loren"]["delta_value"])

		print("Offset:", self.offset_common*1.243858)
		print("Phase:", self.phase_common)
		print("Lorentz:", self.width_loren_common)
		print("Gauss:", self.gauss_common)
		return "Done"

	def parse_metab_info(self, prior_knowledge, i):
		#Here a substract 1 because in the model_verbose_master, the metabolites start in 1.
		#Info[0]: metabolite_nro +1
		#Info[1]: line shape fir, ignore
		#Info[2]: line type: 4 = parametrized, only option
		#Info[3]: Area: [reference, value, min, max]
		#Info[4]: Offset: [reference, value, min, max]
		#Info[5]: Width_gaus: [reference, value, min, max]
		#Info[6]: Phase: [reference, value, min, max]
		#Info[7]: Width_loren: [reference, value, min, max]
		fitting_type = int(prior_knowledge[2])
		if fitting_type != 4:
			print("#### ERROR!: the expected fitting method is not supported, only pattern fitting can be used (i.e. 4), type given:", fitting_type)
			return
		info_string = [params.split("\t") for params in prior_knowledge]
		#need to sabe the original reference, because I need to search the reference with that number, since somethimes it could be the metabolite nro 1.
		metabolite_ref = int(prior_knowledge[0])
		metabolite_name = self.metabolites.get_metaboliteName(i)
		metabolite_shape = self.metabolites.get_metaboliteTD(i)
		#print(info_string)
		is_common = (info_string[4][0]==info_string[6][0]==info_string[7][0]=='0')
		def get_dict(info):
			return {"relative_to": int(info[0]),"original_value":float(info[1]),"delta_value":0.0}

		area_info = {"relative_to":int(info_string[3][0]), "original_value": float(info_string[3][1])/self.scale}
		offset_info = get_dict(info_string[4])
		width_gaus_info = get_dict(info_string[5])
		phase_info = get_dict(info_string[6])
		width_loren_info= get_dict(info_string[7])

		meta_info = {	"metabolite_id": i,
						"metabolite_ref": metabolite_ref,
						"metabolite_name": metabolite_name,
						"metabolite_shape": metabolite_shape,
						"is_common": is_common,
						"area": area_info,
						"offset": offset_info,
						"width_gaus": width_gaus_info,
						"phase": phase_info,
						"width_loren": width_loren_info
					}
		'''
		{metabolite_id,
		 metabolite_name,
		 metabolite_shape,
		 is_common, #This will be True only for the NAA; or whatever the reference to everything is.
		 area:{ relative_to,
				original_value}
		 offset:{relative_to, #If 0, trainable
				 original_value,
				 delta_value	#Update in a next iteration, when I know which one is the reference, final value: original_value - original_value_ref}
		 width_gaus:{relative_to,
					 original_value,
					 delta_value}
		 phase:{relative_to,
				original_value,
				delta_value}
		 width_loren:{relative_to,
					  original_value,
					  delta_value}}
		'''
		return meta_info

	def update_param_value(self, param, metabolite_name, value):
		meta_id = self.metabolites.get_metaboliteID(metabolite_name)
		self.prior[meta_id][param]["original_value"] = value

	def generate_weighted_TD(self, i):
		area = self.prior[i]["area"]["original_value"]
		offset = self.prior[i]["offset"]["delta_value"]
		phase = self.prior[i]["phase"]["delta_value"]
		width_loren = self.prior[i]["width_loren"]["delta_value"]
		width_gaus = self.prior[i]["width_gaus"]["delta_value"]
		metabolite_spectra = self.prior[i]["metabolite_shape"]
		ndp = self.ndp
		pi = math.pi
		alpha = pi
		beta = 2*pi
		gamma = pi/180.0
		delta = pi*pi/(1.67*1.67)
		x = np.linspace(0, ndp-1, ndp)
		t = x/ndp
		expon_real = - alpha*width_loren*t - width_gaus*width_gaus*(t*t)*delta
		expon_imag = beta*offset*t + gamma*phase
		final = np.exp(expon_real + 1j*expon_imag)
		final = area * metabolite_spectra * np.exp(expon_real + 1j*expon_imag)

		return final

	def generate_common(self, use_phase=True, use_gauss=True, use_offset=True, use_loren=True):
		loren_common = 0
		phase_common = 0
		offset_common = 0
		gauss_common = 0
		if use_loren:
			loren_common = self.width_loren_common
		if use_phase:
			phase_common = self.phase_common
		if use_offset:
			offset_common = self.offset_common
		if use_gauss:
			gauss_common = self.gauss_common
		
		ndp = self.ndp
		pi = math.pi
		alpha = pi 
		beta = 2*pi
		gamma = pi/180.0
		delta = pi*pi/(1.67*1.67)
		x = np.linspace(0, ndp-1, ndp)
		t = x/ndp
		expon_real = - alpha*loren_common*t - gauss_common*gauss_common*(t*t)*delta
		expon_imag = beta*offset_common*t + phase_common*gamma
		return np.exp(expon_real + 1j*expon_imag)

	def get_snr_base(self, snr):
		base_model = self.show_model(ret=True)
		return base_model.real.max()/snr

	def show_basis(self, show_all=True, show_total=False):
		total = None
		common = self.generate_common(use_loren=True)
		for i in range(0, self.nro_basis):
			basis_i = self.generate_weighted_TD(i)
			if show_all:
				basis_i_fd = np.fft.fft(basis_i*common)
				plt.plot((basis_i_fd - np.mean(basis_i_fd[0:200])).real, label=self.get_metaboliteName(i))
			if total is None:
				total = basis_i
			else:
				total = total + basis_i
		total = total*common
		total[0] = total[0]/2
		if show_total:
			total_fd = np.fft.fft(total)
			plt.plot(total_fd - np.mean(total_fd[0:200]), label="Simulated Spectrum")
		plt.xlabel('points')
		plt.legend()
		plt.show()

	def show_model(self, meta_index=None, ret=False, snr=None, include_gaus=False):
		import matplotlib.pyplot as plt
		model_final = 0 + 1j*0
		for i in range(0, self.nro_basis):
			model = self.generate_weighted_TD(i)
			if i == meta_index:
				plt.plot(model.real)
				plt.show()
			model_final += model

		model_final = model_final*self.generate_common()
		model_final[0] = model_final[0]/2
		model_final = np.fft.fft(model_final)
		if snr is not None:
			noise = snr*np.random.normal(size=(len(model_final), 2))
			model_final.real = model_final.real + noise[:,0]
			model_final.imag = model_final.imag + noise[:,1]
		if ret==False:
			plt.plot(model_final.real, label="modelV2")
			plt.legend()
			plt.show()
		else:
			return model_final

	def get_metaboliteName(self, i):
		return self.prior[i]["metabolite_name"]

	def get_parameters(self, i):
		# area, delta_offset, delta_loren, delta_phase, shape, name 
		return self.prior[i]["area"]["original_value"], self.prior[i]["offset"]["delta_value"], self.prior[i]["width_loren"]["delta_value"], self.prior[i]["phase"]["delta_value"], self.prior[i]["metabolite_shape"], self.prior[i]["metabolite_name"]

	def get_parameters_common(self):
		return self.offset_common, self.width_loren_common, self.phase_common, self.gauss_common


class metabolite_info():
	def __init__(self, metabolite_list, nro_basis, ndp):
		list_aux = metabolite_list.split("\n")
		self.metaboliteList = [file for file in list_aux if file.find(".bin")>0]
		self.metaboliteName = [os.path.basename(file).replace(".bin","") for file in self.metaboliteList]
		self.nro_basis = nro_basis
		self.metaboliteFD = {}
		self.ndp = ndp

	def __len__(self):
		return self.nro_basis

	def get_ndp(self):
		return self.ndp

	def get_nro(self):
		return self.nro_basis

	def get_metaboliteList(self):
		return self.metaboliteList

	def get_metaboliteName(self, index=None):
		if index is None:
			return self.metaboliteName
		else:
			if index<0 or index>=self.nro_basis:
				print("#### ERROR!: {} is out of limits.".format(index))
			else:
				return self.metaboliteName[index]

	def load_metabolites(self):
		if len(self.metaboliteList) != self.nro_basis:
			print("#### ERROR!: there is more metabolites declared than .bin files")
		for meta, meta_name in zip(self.metaboliteList, self.metaboliteName):
			if Path(meta).stat().st_size != 8*self.ndp:
				print("\n#### ERROR!: The metabolite file: {} is incorrectly formated, should have {} (2*4*{}) but it has {}\n".format(meta, 8*self.ndp, self.ndp, Path(meta).stat().st_size))
				return

			file = open(meta, "rb")
			#print("Metabolite: {} .. OK".format(meta_name))

			imag = np.zeros(self.ndp)
			real = np.zeros(self.ndp)
			for j in range(0, self.ndp):
				test = file.read(4)
				(num_re,) = struct.unpack("f", test)
				real[j] = num_re
				test = file.read(4)
				(num_im,) = struct.unpack("f", test)
				imag[j] = num_im

			complex_list = real + 1j*imag
			complex_list = utils.flip_spectra(complex_list)
			self.metaboliteFD[meta_name] = {"FD": complex_list, "TD": np.fft.ifft(complex_list)}

	def get_metaboliteID(self, name):
		return self.metaboliteName.index(name)

	def get_metaboliteFD(self, index):
		if isinstance(index, int) == True:
			return self.metaboliteFD[self.metaboliteName[index]]["FD"]
		elif isinstance(index, str) == True:
			return self.metaboliteFD[index]["FD"]
		else:
			print("#### ERROR!: Method not implemented")

	def get_metaboliteTD(self, index):
		if isinstance(index, int) == True:
			return self.metaboliteFD[self.metaboliteName[index]]["TD"]
		elif isinstance(index, str) == True:
			return self.metaboliteFD[index]["TD"]
		else:
			print("#### ERROR!: Method not implemented")


class tdfdfit_parser():
	def __init__(self, filename=None, model_name=None, optim=False):
		self.model_master_file = None
		self.optim = optim
		if filename is not None:
			self.model_master_file = filename
			if model_name is not None:
				self.model_name = model_name
		elif model_name is not None:
			self.model_name = model_name
			self.base_path = "Models\\"
			if self.model_name == "2HG-7T-EpsiModel":
				self.epsi = True
			self.model_path = self.base_path + win_sep + self.model_name + win_sep
			self.MasterModelPath = self.model_path + "MasterTDFDFitData" + win_sep
			self.model_master_file = self.MasterModelPath + self.model_name + "_verbose_master-NET.par"
		else:
			print("ERROR!! Need to define which is the model.par to load..")
			exit()
		self.delimiters = ["tdfdfit_begin_of_file","tdfdfit_spectrum_information", "tdfdfit_pattern_references", "tdfdfit_fit_control", "tdfdfit_starting_values", "tdfdfit_prior_knowledge", "tdfdfit_line_shape_function", "tdfdfit_end_of_file"]


	def spectral_info(self, info):
		info_list = info.split("\n")
		info_list = [inform.strip() for inform in info_list]
		return int(info_list[1]), int(info_list[10]), float(info_list[6])


	def read_file(self):
		file = open(self.model_master_file).read()
		segments_list = self.separate_file(file)

		#Parse importantInfo
		ndp, nro_basis, n_shift = self.spectral_info(segments_list[1])
		metabolites = metabolite_info(segments_list[2], nro_basis, ndp)
		metabolites.load_metabolites()
		#Not used now, could be usefull in the future, in case that the convergence is not good
		fittingOrder = segments_list[3]
		#Should not be used in the future, and is not needed now neither since the same info is in the proor_knowledge
		initialValues = segments_list[4]

		prior = prior_knowledge_parse(metabolites, segments_list[5], self.optim)
		return prior


	def separate_file(self, file):
		file = re.sub(re.compile(r"#.*?\n") ,r"" , re.sub(re.compile(r"\n+"), r"\n", file))
		segment_list = []
		for i in range(0,len(self.delimiters)-1):
			segment = r"{}".format(file[file.find(self.delimiters[i]):file.find(self.delimiters[i+1])])
			segment_list.append(segment)
		segment = file[file.find(self.delimiters[i+1]):]
		segment_list.append(segment)
		return segment_list



def main():
	parser = argparse.ArgumentParser(description='Argument parser')
	parser.add_argument('-m', '--model', metavar='model_name', help='Specify model name')
	parser.add_argument('--show', action='store_true', help='Shows the full model with the corresponding initial values')
	parser.add_argument('--detail', action='store_true', help='Show ')
	parser.add_argument('--showBasis', action='store_true', help='Show basis')

	args = parser.parse_args()

	if not args.model:
		print("ERROR: Model is MANDATORY to run each options of this code.")
		exit()

	parser = tdfdfit_parser(model_name=args.model)
	prior = parser.read_file()

	print(f"Model name: {args.model}")

	if args.show:
		prior.show_model()
	if args.detail:
		print(prior)
	if args.showBasis:
		prior.show_basis(show_all=False, show_total=True)



if __name__=="__main__":
	main()




