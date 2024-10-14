function [water_rmd_sig11_withpoly,ret1,ret2] = Loewner_water_supression(spectra,R,ppmaxis,freqindex,ppmrange_HSVD) 


ppm_axis_seg = ppmaxis(freqindex(1):freqindex(2));
spectra_croped = spectra(freqindex(1):freqindex(2),:);

N = length(ppm_axis_seg);
set1 = 1:2:N;
set2 = 2:2:N;
I = length(set1);
J = length(set2);
K = size(spectra_croped,2);
T = zeros(I,J,K);
for i= 1:K
%     T(:,:,i) = loewnerize(spectra_croped(:,i));
    T(:,:,i) = loewnerize(spectra_croped(:,i),'T', ppm_axis_seg);
%     T(:,:,i) = loewnerize(spectra_croped(:,i),'T',2*pi*freq_axis_seg);
end


water_region_ana = [4.2, 5.4];
water_region_ana_index = [find(ppmaxis>min(water_region_ana),1,'first'), find(ppmaxis<max(water_region_ana),1,'last')];
water_region_len = round((water_region_ana_index(2)- water_region_ana_index(1)+1)/2);
pylt_ord = 4;
loops = 0;
while(1)
 [ret1,water_rmd_sig11_withpoly,ret2] = water_supression_loewner(T,spectra,R,ppm_axis_seg,ppmaxis,ppmrange_HSVD,pylt_ord,freqindex);

    for i = 1:K
%         temp1 = norm(signal_without_water_fd(:,i),'fro');
%         err_tns(i) = norm(water_rmd_sig11(:,i) - signal_without_water_fd(:,i),'fro')/temp1; 
%         err_mat(i) = norm(signalfil_fd(:,i) - signal_without_water_fd(:,i),'fro')/temp1; 
%         err_tns_withpoly(i) = norm(water_rmd_sig11_withpoly(:,i) - signal_without_water_fd(:,i),'fro')/temp1; 

        water_region_part_ten = [real(water_rmd_sig11_withpoly(water_region_ana_index(1):water_region_ana_index(2),i))',imag(water_rmd_sig11_withpoly(water_region_ana_index(1):water_region_ana_index(2),i))'];
        noise_region_part_ten = [real(spectra(1:water_region_len,i))',real(spectra(end-water_region_len:end,i))',imag(spectra(1:water_region_len,i))',imag(spectra(end-water_region_len:end,i))'];
        noise_region_part_ten1 = [real(water_rmd_sig11_withpoly(1:water_region_len,i))',real(water_rmd_sig11_withpoly(end-water_region_len:end,i))',imag(water_rmd_sig11_withpoly(1:water_region_len,i))',imag(water_rmd_sig11_withpoly(end-water_region_len:end,i))'];
        norm_wat_reg(i) = var(water_region_part_ten);
        norm_noireg(i) = var(noise_region_part_ten);
        norm_noireg1(i) = var(noise_region_part_ten1);

    end
    mean(norm_wat_reg)
    mean(norm_noireg1)
    if(mean(norm_wat_reg)<4.5*mean(norm_noireg1))
        break
    end
    loops = loops + 1;
end
loops