function signal = align_spectra(signal,frequency,step,ndp)

% Creation of ppm axis:
max_shift = 15; % this is the maximum displacement (number of point on the frequency axis) we allow the spectra 

freq_axis = [-1/(2*step/1000):1/((ndp-1)*step/1000):1/(2*step/1000)];
% freq_axis = linspace(-(1/2/step)*1000,(1/2/step)*1000,ndp);
ppm_axis = kHz2ppm(freq_axis,frequency*1000);
ppm_low = find(ppm_axis>=0.25,1,'first');
ppm_high = find(ppm_axis<=4.2,1,'last');
ppm_axis = ppm_axis(ppm_low:ppm_high);

% Calculation of the windowed fft spectra, to get more accurate alignment: 

n = [0:ndp-1];
N = ndp/2;
expwin_vect = exp(-n/N);
expwin_mat = repmat(expwin_vect,size(signal,1),1);
signal_win=signal.*expwin_mat;
fftmap_win=fftshift(fft(signal_win,[],2),2);
fftmap_win = real(fftmap_win);

% Creation of reference spectrum used for correlation. Exact chemical shift
% for each metabolite is obtained from paper Govindaraju et al.
[NAA_diff,NAA_ppm] = min(abs(ppm_axis-2.008));
[Cre_diff,Cre_ppm] = min(abs(ppm_axis-3.027));
[Cho_diff,Cho_ppm] = min(abs(ppm_axis-3.208));
[Lip_diff,Lip_ppm] = min(abs(ppm_axis-1.28));

ref_peak = triang(7);
ref_spect = zeros(length(ppm_axis),1);

ref_spect(NAA_ppm-3:NAA_ppm+3) = ref_peak;
ref_spect(Cre_ppm-3:Cre_ppm+3) = ref_peak;
ref_spect(Cho_ppm-3:Cho_ppm+3) = ref_peak;
ref_spect(Lip_ppm-3:Lip_ppm+3) = ref_peak;


% finding the highest correltaion by shifting the spectra with respect to
% the reference spectrum:
cor_max_shift = zeros(size(fftmap_win,1),1);
%cor_max = zeros(size(fftmap_win,1),1);
max_reached = zeros(size(fftmap_win,1),1); % this vector registers monitors in which signals the maximum allowed shift has been reached

for i=1:size(fftmap_win,1)
    cor_max = 0;
    for j=-max_shift:max_shift
        cor = fftmap_win(i,ppm_low+j:ppm_high+j)*ref_spect/(norm(fftmap_win(i,ppm_low+j:ppm_high+j))*norm(ref_spect));
        if (cor > cor_max)
            cor_max = cor;
            cor_max_shift(i) = j;
        end
    end
    if (cor_max_shift(i) == -max_shift || cor_max_shift(i) == max_shift)
        max_reached(i) = 1;
    end
end
% max_shift_reached = sum(max_reached)
% if(max_shift_reached ~= 0)
%     uiwait(warndlg(['The maximum allowed spectral shift has been reached in ' num2str(max_shift_reached) ' voxels'],'Caution'));
% end

% Optimal shift has to be converted to frequency shift:
ppm_step = ppm_axis(2)-ppm_axis(1);
ppm_shift = ppm_step*cor_max_shift;
freq_shift = ppm_shift*frequency*1000/10^6;

t = [0:step/1000:(step/1000*(ndp-1))];
% signal = signal.*exp(-sqrt(-1)*freq_shift*2*pi*t);
signal = signal.*exp(-1i*freq_shift*2*pi*t);

% button = questdlg('Do you want to save the aligned time signals?','Saving option','No');
%         if (strcmp(button,'Yes'))
%              [filename, pathname] = uiputfile('*.mat', 'Save aligned time signals');
%              save(fullfile(pathname,filename),'signal','step','frequency','ndp','-v7.3');
%         end
