function ppm = kHz2ppm(kHz,frequency,ppmref)
if nargin<3
    ppmref = 4.7;
end
ppm =kHz/frequency*10^6+ppmref;