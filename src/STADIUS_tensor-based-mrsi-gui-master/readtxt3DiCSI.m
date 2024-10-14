function [signal,begin,ndp,frequency,step,signalnames,position] = readtxt3DiCSI(filename)
% This function reads a text file from jMRUI3.0 and gives the signal parameters
% as output
 [signalnames] = textread(filename,'CSI set name:%s',1);

 [nos] = textread(filename,'Number of voxels:%f',1,'headerlines',1);
 
 step = 0.5; % for 3T data (in ms)
 frequency = 127764.45; % for 3T data (in kHz)
 
 [ndp,val1,val2] = textread(filename,'Number of points per voxel:%d (%d,%d)',1,'headerlines',2);
 begin = val1*step;

 position = zeros(nos,3);
for i = 1:nos
    [pos1,pos2] = textread(filename,'%d%d',1,'headerlines',5+i-1);
    position(i,:) = [pos1+1 pos2+1 i];
%     [sigreal sigimag]= textread(filename,'%f%f',ndp,'headerlines',6+nos+(i-1)*(ndp+1));
%     signal(i,:) = (sigreal - sqrt(-1)*sigimag)';
end

fid=fopen(filename,'r');
st = fread(fid,'*char')';
fclose(fid);
k = findstr(st,'Imaginary');
st(1:k+length('Imaginary')) = [];
data=sscanf(st,'%e',[2,inf])';
data = data(:,1) + sqrt(-1)*data(:,2);
signal = reshape(data,ndp,nos).';
