%  function [water_rmd_sig11,water_rmd_sig11_withpoly,water_rmd_sig11_withpoly_replaced] = water_supression_loewner(T,spectra,R,ppm_axis_seg,ppmaxis,ppmrange_HSVD,pylt_ord,freqindex)
 function [ret1,water_rmd_sig11_withpoly,ret2] = water_supression_loewner(T,spectra,R,ppm_axis_seg,ppmaxis,ppmrange_HSVD,pylt_ord,freqindex)


% options.Display = 10;
% % options.Compression = false;
% % U0{1} = randn(size(T,1),R);
% % U0{2} = randn(size(T,2),R);
% % U0{3} = randn(size(T,3),R);
% U0 = cpd_rnd(size(T),R);
% T_cpd = cpd_als(T,U0,options);
%  % [T_cpd,output] = cpd(T,U0);

options.Refinement = false;
options.MaxIter = 200;
options.Display = false;
options.ExploitStructure = false;
[T_cpd,output] = cpd(T,R,options);

for i = 1:R
    [~,rec_sig(:,i),para(i,:)] =  direct_lorentz_modelling(T_cpd{1}(:,i),T_cpd{2}(:,i),ppm_axis_seg',ppmaxis'); 
end
%%
% find the peaks that are outside the region of intrest (4.2 - 0.25)
[~,locs] = max(abs(rec_sig));
ppm_srcs = ppmaxis(locs);
water_resign_src = find((ppm_srcs >= ppmrange_HSVD(2)) | (ppm_srcs <= ppmrange_HSVD(1)));
water_resign_src = water_resign_src(:)';
S_new = rec_sig;

% water_resign_src = [water_resign_src,3]; 
H_new = S_new\spectra; % obtain miximg matrix
k2 = S_new(:,water_resign_src)*H_new(water_resign_src,:); % reconstruct water sources and baseline
water_rmd_sig11 = spectra - k2; %  remove water sources and baseline
% pylt_ord = 5;
x_ser = ppmaxis';
nw_mat = ones(size(spectra,1),pylt_ord+1);
for i = 1:pylt_ord
    nw_mat(:,i+1) = x_ser.^i;
end
nw_mat = normc(nw_mat);
S_new_withpoly = [rec_sig,nw_mat];
water_resign_src_withpoly = [water_resign_src,R+1:R+pylt_ord+1]; % select only the water sources and polynomial
H_new_withpoly = S_new_withpoly\spectra; % obtain miximg matrix
k2_withpoly = S_new_withpoly(:,water_resign_src_withpoly)*H_new_withpoly(water_resign_src_withpoly,:); % reconstruct water sources and baseline

water_rmd_sig11_withpoly = spectra - k2_withpoly;
water_rmd_sig11_withpoly_replaced = water_rmd_sig11_withpoly;
water_rmd_sig11_withpoly_replaced(freqindex(1):freqindex(2),:) = water_rmd_sig11(freqindex(1):freqindex(2),:);
% k3 = S_new_withpoly*H_new_withpoly;

ret1 = S_new_withpoly(:,water_resign_src_withpoly);
ret2 = H_new_withpoly(water_resign_src_withpoly,:);
