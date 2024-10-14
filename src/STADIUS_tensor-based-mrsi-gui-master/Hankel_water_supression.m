function [water_rmd_sig11_TD] = Hankel_water_supression(with_water_signal,Kest,N,ppmrange_HSVD,step,frequency) 
ndp = size(with_water_signal,1);
with_water_signal_truncated = with_water_signal(1:N,:);
H_T = hankelize(with_water_signal_truncated,'Full', true);
[U,~] = mlsvd(H_T,[Kest,Kest,min([Kest,size(with_water_signal,2)])]);
u1 = U{1};
L = size(u1,1);
upt=u1(1:L-1,:);           % top part of matrix [u1] without bottum line 
upb=u1(2:L,:);             % bottom part of matrix [u1] without top line
q=pinv(upt)*upb;           % LS solution by computing the Pseudoinverse
z=eig(q).';

for i = 1:Kest
    sources(:,i) = ((z(i)*ones(1,ndp)).^(0:ndp-1)); 
    sources_fd(:,i) = fftshift(fft(sources(:,i)));
end
% ppmrange_HSVD = [0.25, 4.2];
% [~,locs] = max(abs(sources_fd));
% ppm_srcs = ppmaxis(locs);
ppm_srcs = kHz2ppm(((imag(log(z))/(2*pi*step))*1000), frequency*1000);

water_resign_src = find((ppm_srcs >= ppmrange_HSVD(2)) | (ppm_srcs <= ppmrange_HSVD(1)));
length(water_resign_src)
H_new = sources\with_water_signal;
k2 = sources(:,water_resign_src)*H_new(water_resign_src,:);
% H_new = sources(:,water_resign_src)\with_water_signal;
% k2 = sources(:,water_resign_src)*H_new;

water_rmd_sig11_TD = with_water_signal - k2; %


