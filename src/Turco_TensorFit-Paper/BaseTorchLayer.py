import torch
import math
import numpy as np

def prepare_metabolite(prior_knowledge, metabolite_id, offset_common, number_of_spectra, device, x, truncate, freedom, type="torch", ub=math.inf, lb=-math.inf, limit_area=False):
    ampl, delta_offset, delta_loren_width, delta_phase, metabolite_shape, metabolite_name = prior_knowledge.get_parameters(metabolite_id)
    metabolite_shape = torch.from_numpy(metabolite_shape).to(device).to(torch.complex64)
    return Metabolite_torch(metabolite_name, number_of_spectra, ampl, delta_offset, delta_loren_width, delta_phase, metabolite_shape, metabolite_shape.shape[0], device, x, truncate, freedom, ub, lb, limit_area)


class Spectra_model(torch.nn.Module):
    def __init__(self, prior_knowledge, number_of_spectra, device, x, input_size, initial_offset=None, use_td=False, limits="None", freedom="None", truncate=1024, add_gaussian=False, limit_area=False, ub=None, lb=None):
        super(Spectra_model, self).__init__()
        offset_common, width_common, phase_common, gaus_width_common = prior_knowledge.get_parameters_common()
        self.prior = prior_knowledge
        self.use_td = use_td
        self.nro_bases = len(prior_knowledge)
        self.limit_area = limit_area
        if limit_area == False:
            ub = torch.ones(self.nro_bases)
            lb = -torch.ones(self.nro_bases)

        if freedom not in ["None", "met_shift"]:
            print("The freedom parameter: {freedom} is not supported, only ['None','met_shift']. Using default instead (i.e. without freedom in the metabolite shift)")
            freedom = "None"
        if limits not in ["None","Sigmoid","Clamp"]:
            print("The border condition ({}) is not supported, using nothing instead".format(limits))
            limits = "None"
        self.truncate = truncate
        if truncate < input_size and truncate>0:
            self.zeroes = torch.zeros(number_of_spectra,input_size-truncate, device=device)

        self.common = Common_torch(number_of_spectra, offset_common, width_common, phase_common, gaus_width_common, device, input_size, x, truncate, initial_offset=initial_offset, limit_type=limits, add_gaussian=add_gaussian)
        self.metabolites = torch.nn.ModuleList([prepare_metabolite(prior_knowledge, i, offset_common, number_of_spectra, device, x, truncate, freedom=freedom, ub=ub[i], lb=lb[i], limit_area=limit_area) for i in range(0, self.nro_bases)])

    def get_metabolites(self, formatt="torch"):
        if formatt == "torch":
            return [self.metabolites[i].get_metabolite_shape() for i in range(0, self.nro_bases)]
        else:
            #No need for detach here, but... you never know..
            return [self.metabolites[i].get_metabolite_shape().detach().cpu().numpy() for i in range(0, self.nro_bases)]

    def forward(self, x, trunc=False):
        metas = self.metabolites[0](x, trunc=trunc)
        for i in range(1, self.nro_bases):
            metas.add_(self.metabolites[i](x, trunc=trunc))

        common_mult = self.common(x, trunc=trunc)
        #final = torch.einsum('ik,ik->ik', metas, common_mult)
        final = torch.mul(metas, common_mult)
        final[:,0] *= 0#.5
        if not self.use_td:
            if trunc:
                final = torch.cat((final, self.zeroes), dim=1)
            final = torch.fft.fftn(final, dim=1)
        return final

    def set_params(self, dictio):
        for meta in self.metabolites:
            meta.set_params(dictio)
        self.common.set_params(dictio)

    def get_tensor_params(self):
        final = {}
        for meta in self.metabolites:
            final.update(meta.get_tensor_params())
        final.update(self.common.get_tensor_params())
        return final

    def get_params(self):
        final = {}
        for meta in self.metabolites:
            final.update(meta.get_params())
        final.update(self.common.get_params())
        return final

class Metabolite_torch(torch.nn.Module):
    def __init__(self, metabolite_name, number_of_spectra, ampl, delta_offset, delta_loren_width, delta_phase, metabolite_shape, input_size, device, x, trunc, freedom, ub, lb, limit_area):
        super(Metabolite_torch, self).__init__()
        self.limit_area = limit_area
        self.metabolite_name = metabolite_name
        aux_ampl = torch.ones([number_of_spectra,1])
        torch.nn.init.constant_(aux_ampl, ampl)
        aux_offset = torch.zeros([number_of_spectra,1])
        self.Dense_ampl = torch.nn.parameter.Parameter(aux_ampl)
        self.trunc = trunc
        self.freedom = freedom
        if freedom == "met_shift":
            self.Dense_offset = torch.nn.parameter.Parameter(aux_offset)
        self.f2p = 0.854044
        self.ndp = input_size
        self.ub = ub
        self.lb = lb
        expon_im = (2*math.pi*delta_offset*x.t()[:,0]/self.ndp + delta_phase*math.pi/180)
        expon_real = -delta_loren_width*math.pi*x.t()[:,0]/self.ndp

        exponential = torch.exp(expon_real + 1j*expon_im)
        self.metabolite_shape = (exponential.t()*metabolite_shape)[None,:]
        self.const_3 = 2*math.pi/self.ndp

    def get_metabolite_shape(self):
        return self.metabolite_shape

    def limits(self, x, min_val, max_val):
        if self.limit_area:
            return torch.clamp(x, min=min_val, max=max_val)
        else:
            return x
    def forward(self, x, trunc):
        if self.freedom == "met_shift":
            exp_tmp = torch.exp(1J*(self.limits(self.Dense_offset, -20.0*self.f2p, 20.0*self.f2p)@x*self.const_3))
            return (self.limits(self.Dense_ampl, self.lb, self.ub).to(torch.complex64)@self.metabolite_shape)*exp_tmp
        else:
            if trunc:
                return self.limits(self.Dense_ampl, self.lb, self.ub).to(torch.complex64)@self.metabolite_shape[:, 0:self.trunc]
            return self.limits(self.Dense_ampl, self.lb, self.ub).to(torch.complex64)@self.metabolite_shape

    def set_params(self, dictio):
        self.Dense_ampl = torch.nn.parameter.Parameter(dictio[self.metabolite_name], requires_grad=True)

    def get_tensor_params(self):
        return {self.metabolite_name: self.Dense_ampl.to(torch.float32).detach()}

    def get_params(self):
        return {"{}".format(self.metabolite_name): self.Dense_ampl.to(torch.float32).cpu().detach().numpy().squeeze()}


class Common_torch(torch.nn.Module):
    def __init__(self, number_of_spectra, common_offset, common_width, common_phase, gaus_width_common, device, input_size, x, truncate, initial_offset=None, limit_type="None", add_gaussian=False):
        super(Common_torch, self).__init__()
        self.limit_type = limit_type
        self.number_of_spectra = number_of_spectra
        aux_offset = torch.ones([number_of_spectra,1])
        torch.nn.init.constant_(aux_offset, common_offset)
        if initial_offset is not None:
            aux_offset = aux_offset + initial_offset[:,None]
        self.offset = aux_offset
        self.trunc = truncate
        self.Dense_offset = torch.nn.parameter.Parameter(aux_offset)
        aux_width = torch.ones([number_of_spectra,1])
        torch.nn.init.constant_(aux_width, common_width)
        aux_phase = torch.ones([number_of_spectra,1])
        torch.nn.init.constant_(aux_phase, common_phase)
        self.Dense_width = torch.nn.parameter.Parameter(aux_width)
        self.Dense_phase = torch.nn.parameter.Parameter(aux_phase)
        if gaus_width_common!=0 or add_gaussian:
            aux_gaus_width = torch.ones([number_of_spectra, 1])
            torch.nn.init.constant_(aux_gaus_width, gaus_width_common)
            self.Dense_gaus_width = torch.nn.parameter.Parameter(aux_gaus_width)
        else:
            self.Dense_gaus_width = None

        ndp = x.shape[1]
        x = x[0:1,:]
        self.gaus_aux = (x*math.pi/(ndp*1.67))*(x*math.pi/(ndp*1.67))
        self.const = math.pi/ndp
        self.const_2 = math.pi/180
        self.const_3 = 2*self.const
        self.const_4 = self.const/1.67

    def limits(self, x, min_val, max_val):
        if self.limit_type == "Clamp":
            return torch.clamp(x, min=min_val, max=max_val)
        else:
            return x

    def forward(self, x, trunc):
        if not trunc:
            self.trunc = x.shape[1]

        if self.Dense_gaus_width is None:
            expon_re = -(self.const*(torch.abs(self.Dense_width@x[:,0:self.trunc])))
        else:
            expon_re = -(self.const*(torch.abs(self.Dense_width@x[:,0:self.trunc])) + torch.abs(self.Dense_gaus_width + (self.Dense_gaus_width==0).float()*10E-6)@self.gaus_aux[:, 0:self.trunc])

        expon_im = (self.const_2*self.Dense_phase + self.const_3*(self.Dense_offset@x[:,0:self.trunc]))
        return torch.exp(torch.complex(expon_re, expon_im))

    #get dictionary with format: {"offset":Tensor, "width":Tensor, "phase":Tensor, "gausw":Tensor}
    def set_params(self, dictio):
        self.Dense_offset = torch.nn.parameter.Parameter(dictio["offset"], requires_grad=True)
        self.Dense_width = torch.nn.parameter.Parameter(dictio["width"], requires_grad=True)
        self.Dense_phase = torch.nn.parameter.Parameter(dictio["phase"], requires_grad=True)
        if self.Dense_gaus_width is not None:
            self.Dense_gaus_width = torch.nn.parameter.Parameter(dictio["gausw"], requires_grad=True)

    def get_tensor_params(self):
        if self.Dense_gaus_width is None:
            return {"offset": self.Dense_offset.detach(),
                    "width": self.Dense_width.detach(),
                    "phase": self.Dense_phase.detach(),
                    "gausw": torch.zeros_like(self.Dense_width.detach())}
        else:
            return {"offset": self.Dense_offset.detach(),
                    "width": self.Dense_width.detach(),
                    "phase": self.Dense_phase.detach(),
                    "gausw": torch.sqrt(torch.abs(self.Dense_gaus_width.detach()))}

    def get_params(self):
        if self.Dense_gaus_width is None:
            return {"offset": self.Dense_offset.cpu().detach().numpy().squeeze(),
                    "width": self.Dense_width.cpu().detach().numpy().squeeze(),
                    "phase": self.Dense_phase.cpu().detach().numpy().squeeze(),
                    "gausw": np.zeros_like(self.Dense_width.cpu().detach().numpy().squeeze())}
        else:
            return {"offset": self.Dense_offset.cpu().detach().numpy().squeeze(),
                    "width": self.Dense_width.cpu().detach().numpy().squeeze(),
                    "phase": self.Dense_phase.cpu().detach().numpy().squeeze(),
                    "gausw": np.sqrt(np.abs(self.Dense_gaus_width.cpu().detach().numpy().squeeze()))}
