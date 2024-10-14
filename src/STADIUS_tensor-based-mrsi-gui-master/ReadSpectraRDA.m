function [Data_time_matrix, Data_time_reshaped, Data_fft_matrix, ...
    Data_fft_reshaped, Labeling, Parameters] = ReadSpectraRDA(filename)
% [Data_time_matrix, Data_time_reshaped, Data_fft_matrix, ...
%     Data_fft_reshaped, Labeling, Parameters] = ReadSpectraRDA(filename)
%
% Function to read metabolite data from *.RDA (Siemens) files.
% Parameters are read from the ASCII part and the data from the binary part 
% of the file.
% (Parameters from the binary part are stored in the array Parameters)
%
%
% INPUT:
% - filename            : Filename of the input data in .RDA format
%
% OUTPUT:
% - Data_time_matrix    : FID in a 4-D matrix notation (in time domain)
% - Data_time_reshaped  : FID listed in a 3-D matrix notation (in time domain)
%                         Each row corresponds with the fid of one voxel
% - Data_fft_matrix     : Spectra in a 4-D matrix notation (in Fourier domain)
% - Data_fft_reshaped   : Spectra listed in a 3-D matrix notation (in Fourier domain)
%                         Each row corresponds with the spectrum of one voxel
% - Labeling            : Labeling of the spectra
%                         Corresponds to the 3-D matrix notation:
%                         [voxel number ; row nr ; column nr ; slice nr]
% - Parameters          : Parameters, extracted from the binary part of the
%                         data file
%
%
%
%  Author:    P.W.T. Krooshof
%             Radboud University Nijmegen
%             Institute for Molecules and Materials
%             Analytical Chemistry
%             Nijmegen, The Netherlands
%  Created:   June 28, 2007
%  Updated:   June 03, 2009 Bob Hamans, comments marked BCH
%             Jan 08, 2013 Diana Sima, KU Leuven

% NOTE: Assumed is that the data is stored as Little Endian

%%% VARIABLES %%%
nr_points_complex = 2;      % The data is stored as complex double values
nr_bytes_datapoint = 8;     % Each data value is 8 bytes
%%% VARIABLES %%%

% Open the .rda file in matlab. Give the read format little-endian (ieee-le)
[FileId, file_open_error] = fopen(filename, 'r', 'ieee-le');
if ~isempty(file_open_error)
    disp(['File opening error: ' file_open_error]);
end

% Load the data in the array data, using unsigned character; 8 bits.
fseek(FileId, 0, 'bof');
filecontent = fread(FileId, 2^22); %BCH read first 4MB only
charfilecontent = char(filecontent');
% file_size = length(filecontent); %BCH commented out
fileinfo = dir(filename); %BCH
file_size = fileinfo.bytes; %BCH

% Read various parameters from the data
[Parameters.VectorSize] = FindAndGet(charfilecontent,'VectorSize');
% [Parameters.VectorSizeComplex] = 2 * [Parameters.VectorSize];
[Parameters.CSIMatrixSizeX] = FindAndGet(charfilecontent,'CSIMatrixSize[0]');
[Parameters.CSIMatrixSizeY] = FindAndGet(charfilecontent,'CSIMatrixSize[1]');
[Parameters.CSIMatrixSizeZ] = FindAndGet(charfilecontent,'CSIMatrixSize[2]');
[Parameters.frequency] = FindAndGet(charfilecontent,'MRFrequency')*1e3; % get spectrometer frequency in kHz (DS)
[Parameters.step] = FindAndGet(charfilecontent,'DwellTime')/1000; % get sampling frequency in ms (DS)
[Parameters.aver] = FindAndGet(charfilecontent,'NumberOfAverages');

[Parameters.TotalSize] = [Parameters.VectorSize] * ...
    [Parameters.CSIMatrixSizeX] * [Parameters.CSIMatrixSizeY] * ...
    [Parameters.CSIMatrixSizeZ];
[Parameters.NumberOfBytes] = nr_points_complex * nr_bytes_datapoint * ...
    [Parameters.TotalSize];
[Parameters.NumberOfPoints] = nr_points_complex * [Parameters.TotalSize];

startpos = findstr(charfilecontent, '>>> End of header <<<') + 21;
endpos = file_size;

% Load the data in the array data, using floating point; 64 bits.
fseek(FileId, double([file_size] - [Parameters.NumberOfBytes]), 'bof');
[temp_data, count] = fread(FileId, inf, 'float64');
temp_data = temp_data';
fclose(FileId);

% Extract the data of each voxel
nr_columns = [Parameters.CSIMatrixSizeX];
nr_rows    = [Parameters.CSIMatrixSizeY];
nr_slices  = [Parameters.CSIMatrixSizeZ];
nr_vars    = [Parameters.VectorSize];

Data_time_matrix   = zeros(nr_rows, nr_columns, nr_slices, nr_vars);
Data_time_reshaped = zeros(nr_rows*nr_columns, nr_vars, nr_slices);
Data_fft_matrix   = zeros(nr_rows, nr_columns, nr_slices, nr_vars);
Data_fft_reshaped = zeros(nr_rows*nr_columns, nr_vars, nr_slices);

pos = 1;
vect_size = nr_points_complex * [Parameters.VectorSize];
for slice_nr = 1:nr_slices
    voxel_nr = 1;
    for row_nr = 1:nr_rows
        for column_nr = 1:nr_columns
            % Split the data in real and imaginary datapoints
            data_real = ...
                temp_data([((pos-1) * vect_size) + 1]:2:[pos * vect_size]);
            data_imag = ...
                temp_data([((pos-1) * vect_size) + 2]:2:[pos * vect_size]);
            data_complex = data_real + i*data_imag;
            
            Data_time_matrix(row_nr, column_nr, slice_nr, :) = data_complex;
            Data_time_reshaped(voxel_nr,:,slice_nr) = data_complex;
            
            Data_fft_matrix(row_nr, column_nr, slice_nr, :) = ...
                real(fftshift(fft(data_complex)));
            Data_fft_reshaped(voxel_nr,:,slice_nr) = ...
                real(fftshift(fft(data_complex)));
            
            Labeling(pos,:) = [voxel_nr row_nr column_nr slice_nr];
            
            voxel_nr = voxel_nr + 1;
            pos = pos + 1;
        end
    end
end

return
 
% =========================================================================
function [value] = FindAndGet(charfilecontent, Paramstr)
% [value] = FindAndGet(charfilecontent, Paramstr)
%
% Function to extract file contents from the ASCII part of the file
%

ascii_param_pos = findstr(charfilecontent, '>>> Begin of header <<<') + 23;
ascii_param_pos = ascii_param_pos(1); %BCH
tpos = findstr(charfilecontent([ascii_param_pos]:[end]), Paramstr);
tstring = sscanf(charfilecontent([tpos + ascii_param_pos - 1]:[end]), '%s', 2);
t2pos = findstr(tstring, ':');
tstringl = length(tstring);
charval = tstring([t2pos(1,1) + 1]:tstringl);
value = str2num(charval);

return
