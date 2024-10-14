% function [signal,step,frequency,ndp,begin,signalnames,procres,classopt,simsig,position,averages] = loadsignals(fn,p)
function [signal,step,frequency,ndp,begin,signalnames,procres,classopt,simsig,position] = loadsignals(fn,p)

%fn = file(s)
%p = path(s)

[pathstr,name,ext] = fileparts(fn);
signal = [];
frequency = [];
signaltemp = [];
signalnames = [];
classopt= [];
procres = [];
simsig = [];
position = [];
iscfn = iscell(fn);
if iscfn
    loadf = fullfile(p,char(fn(1)));
else
    loadf = fullfile(p,fn);
end
astr=[];
switch ext
    case '.mat' %matlab data
        step=0;
        astr = load(loadf);
        names = fieldnames(astr);
        % if ~(isfield(astr,'ndp')&isfield(astr,'step')&isfield(astr,'begin'))
        %     set(handles.preprocess,'Enable','off');
        %     set(handles.process,'Enable','off');
        %     set(handles.plotting,'Enable','off');
        %     set(handles.analysis,'Enable','off');
        %     errordlg('The variables signal, ndp, begin and step must be defined');
        %     return
        % end

        if ~(isfield(astr,'ndp') && isfield(astr,'step') && isfield(astr,'begin'))
            %     set(handles.preprocess,'Enable','off');
            %     set(handles.process,'Enable','off');
            %     set(handles.plotting,'Enable','off');
            %     set(handles.analysis,'Enable','off');
            errordlg('The variables signal, ndp, begin and step must be defined');
            return
        end

        ndp = astr.ndp;
        begin = astr.begin;
        step = astr.step;
        signalnames ={};
        csn = 1;
        for i=1:length(names)
            if (astr.ndp == size(astr.(names{i}),2))
                signalnames(csn) = names(i);
                csn = csn + 1;
                signaltemp = [signaltemp;astr.(names{i})];
            else
                v = genvarname(names{i});
                eval([v ' = astr.(names{i});'])
            end
        end
        % if ~exist('signal')
        %     set(handles.preprocess,'Enable','off');
        %     set(handles.process,'Enable','off');
        %     set(handles.plotting,'Enable','off');
        %     set(handles.analysis,'Enable','off');
        %     errordlg('The variables signal, ndp, begin and step must be defined');
        %     return
        % end
        if ~exist('signal')
            %     set(handles.preprocess,'Enable','off');
            %     set(handles.process,'Enable','off');
            %     set(handles.plotting,'Enable','off');
            %     set(handles.analysis,'Enable','off');
            errordlg('The variables signal, ndp, begin and step must be defined');
            return
        end
        signal = [signal;signaltemp];
        % checking the hyperparameters (begin, step, ndp, frequency)
        if iscfn
            for i = 2:length(fn)
                signaltemp = [];
                loadf = fullfile(p,char(fn(i)));
                astrtemp = load(loadf);
                names = fieldnames(astrtemp);
                if (ndp ~= astrtemp.ndp)
                    errordlg('The signals should have the same length, i.e. same ''ndp''.');
                end
                if (begin ~= astrtemp.begin)
                    errordlg('The signals should have the same begin time.');
                end

                if (step ~= astrtemp.step)
                    errordlg('The signals should have the same time sampling, i.e. same ''step''.');
                end
                for i=1:length(names)
                    if (astrtemp.ndp == size(astrtemp.(names{i}),2))
                        signalnames(csn) = names(i);
                        csn = csn + 1;
                        signaltemp = [signaltemp;astrtemp.(names{i})];
                    else
                        v = genvarname(names{i});
                        eval([v ' = astrtemp.(names{i});'])
                    end
                end
                signal = [signal;signaltemp];
            end
        end
    case '.txt' %jMRUI or 3DiCSI data
        [signal,begin,ndp,frequency,step,signalnames,position] = readtxt(loadf);
        classopt = 0;
        procres = 0;
        
    case {'.spar','.sdat','.SPAR','.SDAT'} %philips data
        [signal, ndp, step, begin, nfids, volume, frequency, position]= read_philips_file(loadf(1:end-5));
        signalnames = name;
        %frequency = 63898; %this supposes a 1H MRS at 1.5 T => not OK if 3T or 13C or any other types of acquisition protocol
        classopt = 0;
        procres = 0;
    case {'.rda','.RDA'} %siemens data
        [Data_time_matrix, Data_time_reshaped, Data_fft_matrix, Data_fft_reshaped, Labeling, Parameters] = ReadSpectraRDA(loadf);
        ndp = Parameters.VectorSize;
        step = Parameters.step;
        frequency = Parameters.frequency;
        begin = 0;
        slices = Parameters.CSIMatrixSizeZ;
        averages = Parameters.aver;
        signal = [];
        for l = 1:slices
            signal = [signal; squeeze(Data_time_reshaped(:,:,l)).''];
        end
        position = Labeling(:,[1 3 2 4]); %Labeling has voxel_nr row_nr column_nr slice_nr, while position has voxel_nr column_nr row_nr slice_nr
    case {'.ima','.IMA','.dcm'} %siemens DICOM data
        h = simpledicom(loadf);
        %ndp = h.f.hdr.sSpecPara.lVectorSize;
        [text,r]  = strtok(h.f.hdr,'=');
        [t,value] = strtok(r,' ');
        ndp = str2double(readfieldparam({text value},'sSpecPara.lVectorSize'));
        step = 2*h.f.DW*1e3;
        begin = 0;
        nrows  = str2double(readfieldparam({text value},'sSpecPara.lFinalMatrixSizePhase'));
        ncols  = str2double(readfieldparam({text value},'sSpecPara.lFinalMatrixSizeRead'));
        slices = str2double(readfieldparam({text value},'sSpecPara.lFinalMatrixSizeSlice'));
        frequency = h.f.SF/1e3;

        nfids = nrows*ncols*slices;
        if isnan(nfids)
            nfids = 1;
        end
        signal = zeros(nfids,ndp);
        for i = 1:nfids
            signal(i,:) = (h.fid((i-1)*ndp+1:i*ndp)).'';
        end

        if (nfids > 1)
            pos1 =[]; pos2 = []; pos3 = [];
            for l = 1:slices
                for i = 1:nrows
                    pos1 = [pos1; [1:ncols]'];
                    pos2 = [pos2; i*ones(1,ncols)'];
                end
                pos3 = [pos3; l*ones(nrows*ncols,1)];
            end
            position = [[1:nfids]' pos1 pos2 pos3];
        else
            position = [1 1 1 1];
        end
    otherwise
        error('No signal has been loaded. File format unknown.')
        return
end
if isempty(frequency)
    frequency = 63830; %frequency at 1.5T for H
    warndlg('The spectrometer frequency has been set at 63.83 MHz. Correct (add ''frequency'' variable) if this is not right.')
end

if ~isempty(astr) 
    for i=1:length(names)
        if (strcmp(names{i},{'signal','step','frequency','ndp','begin','signalnames','procres','classopt','position'}))
            simsig.(names{i}) = astr.names{i};
        end
    end
end