function [signalfil,signalfil_fd] = water_removal_HSVD(signal,Kest,step,frequency,ppmrange)
% function [signalfil] = water_removal_HSVD(signal,Kest,step,frequency,ppmrange)
fsampl = 1/step*1000; %in Hz
xmin = (ppmrange(1)-4.7)*frequency/10^6;%in kHz
xmax = (ppmrange(2)-4.7)*frequency/10^6;%in kHz
freqrange1 = [-(fsampl*1000/2-1),xmin*1000]; %in Hz
freqrange2 = [xmax*1000,(fsampl*1000/2-1)]; %in Hz
freqrange = [freqrange1;freqrange2];
[nos,ndp] = size(signal);
t = [0:step:(ndp-1)*step]; 
MM = ndp/2;
for signum=1:nos

   [signalfil(signum,:),freq,damp,ampl,phas,SV] = HLSVDPROrec(signal(signum,:),Kest,fsampl,t/1000,MM,freqrange);
%     if kmodel== 1
%         boundL = ceil(ppmaxis(1)); %in ppm
%         boundH = floor(ppmaxis(end)); % in ppm
%         signalfil = basOffCor(signalfil(signum,:),step,frequency,ndp,begin,boundL,boundH);
%     end
    signalfil_fd(:,signum) = fftshift(fft(signalfil(signum,:))).';

end %i=1:nos
% damp/10
% kHz2ppm(freq,frequency*1000)

