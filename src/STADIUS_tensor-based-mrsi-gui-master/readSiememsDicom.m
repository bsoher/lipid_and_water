% special function to disentagle MRSI FIDs read by the simpledicom
% function
function [signal, ndp, step, begin, nfids, volume, frequency, position] = readSiememsDicom(file2read)

h = simpledicom(file2read);

ndp = h.f.hdr.sSpecPara.lVectorSize;

rows   = h.f.hdr.sSpecPara.lFinalMatrixSizePhase;
cols   = h.f.hdr.sSpecPara.lFinalMatrixSizeRead;
slices = h.f.hdr.sSpecPara.lFinalMatrixSizeSlice;

frequency = h.f.SF/1e3;

for i=1:256, 
    fidmrsi(i,1:512)=h.fid((i-1)*512+1:i*512);
end