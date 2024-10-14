function [sol,out] = energy_tensor(T,R,Nmalise,wei,init)

size_tens = size(T);

if(isempty(init))
    model.variables.A = rand(size_tens(1),R);
    model.variables.C = rand(size_tens(3),R);
else
    model.variables.A = init.src;
    model.variables.C = init.hbd;
end
% model.variables.C = rand(R,size_tens(3));
if(isempty(Nmalise))
    model.factors.A = {'A',@struct_nonneg};
else
    if(Nmalise == 1)
%         model.factors.A = {'A',@struct_nonneg,@struct_normalize};
        model.factors.A = {'A',@struct_normalize};
    else
        model.factors.A = {'A',@struct_nonneg};
    end
end
model.factors.C = {'C',@struct_nonneg};
% model.factors.C = {'C'};
model.factorizations.mybtd.data = T;
model.factorizations.mybtd.cpd = {'A','A','C'};

if(~isempty(wei))
    model.factorizations.myreg1.regL1 = {'C'};
    options.Weights = wei;
end

options.Display = 10;
options.TolFun = 1e-9; 
% options.TolX = 1e-08;
[sol,out] = sdf_nls(model,options);


% size_tens = size(T);
% model.variables.A1 = rand(size_tens(1),1);
% model.variables.A2 = rand(size_tens(1),1);
% 
% model.variables.C = rand(size_tens(3),R);
% % model.variables.C = rand(R,size_tens(3));
% 
% model.factors.A = {'A1','A2',@struct_nonneg,@struct_normalize};
% model.factors.C = {'C',@struct_nonneg};
% % model.factors.C = {'C'};
% model.factorizations.mybtd.data = T;
% model.factorizations.mybtd.cpd = {'A','A','C'};
% 
% % model.factorizations.myreg1.regL1 = {'C'};
% 
% options.Display = 5;
% sol = sdf_nls(model,options);

