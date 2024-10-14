function [signal,begin,ndp,frequency,step,signalnames,position] = readtxt(filename)
% This function reads a text file from jMRUI3.0 and gives the signal parameters
% as output
 firstword = textread(filename,'%s',1);
 if strcmp(firstword,'CSI')
     [signal,begin,ndp,frequency,step,signalnames,position] = readtxt3DiCSI(filename);
 else
     [signal,begin,ndp,frequency,step,signalnames] = readtxtjMRUI(filename);
     position = [];
 end
 