function [signal,begin,ndp,frequency,step,signalnames] = readtxtjMRUI(filename)
% This function reads a text file from jMRUI3.0 and gives the signal parameters
% as output
 [var,signalnames] = textread(filename,'%s%s',1,'headerlines',2);

 [var,val] = textread(filename,'%s%f',6,'headerlines',4);
 ndp = val(1);
 nos = val(2);
 step = val(3);
 begin = val(5);
 frequency = val(6)/1000; % in kHz
% [var,val] = textread(filename,'%s%s',2,'headerlines',16);
% def = strcmp(char(var(1)),'Signal');
 fil = fopen(filename, 'r');
 text = textscan(fil, '%s',50, 'delimiter', '\n');
 searching = strfind([text{1}],'Signal 1 out of');
 last_header_line = find(not(cellfun('isempty', searching)));
 
 if last_header_line > 0
     for i = 1:nos
%        [sigreal sigimag fftreal fftimag]= textread(filename,'%f%f%f%f',ndp,'headerlines',21+(i-1)*(ndp+1));
        [sigreal sigimag fftreal fftimag]= textread(filename,'%f%f%f%f',ndp,'headerlines',last_header_line+(i-1)*(ndp+1));
        signal(i,:) = (sigreal - sqrt(-1)*sigimag)'; 
    end
 else
  %% for signals saved from jmrui siemens format or svsviewer format
     
     [var,val] = textread(filename,'%s%s',2,'headerlines',418);
     def2 = strcmp(char(var(1)),'Signal');
     if def2==1
        for i = 1:nos
            [sigreal sigimag]= textread(filename,'%f%f',ndp,'headerlines',422+(i-1)*(ndp+1));
            signal(i,:) = (sigreal - sqrt(-1)*sigimag)'; 
        end
     else 
        for i = 1:nos
            [sigreal sigimag]= textread(filename,'%f%f',ndp,'headerlines',426+(i-1)*(ndp+1));
            signal(i,:) = (sigreal - sqrt(-1)*sigimag)'; 
        end
     end
  end