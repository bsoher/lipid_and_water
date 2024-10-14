function varargout = MRSI_tissue(varargin)
% MRSI_TISSUE MATLAB code for MRSI_tissue.fig
%      MRSI_TISSUE, by itself, creates a new MRSI_TISSUE or raises the existing
%      singleton*.
%
%      H = MRSI_TISSUE returns the handle to a new MRSI_TISSUE or the handle to
%      the existing singleton*.
%
%      MRSI_TISSUE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRSI_TISSUE.M with the given input arguments.
%
%      MRSI_TISSUE('Property','Value',...) creates a new MRSI_TISSUE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MRSI_tissue_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MRSI_tissue_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MRSI_tissue

% Last Modified by GUIDE v2.5 20-Dec-2018 18:35:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRSI_tissue_OpeningFcn, ...
                   'gui_OutputFcn',  @MRSI_tissue_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before MRSI_tissue is made visible.
function MRSI_tissue_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MRSI_tissue (see VARARGIN)

% Choose default command line output for MRSI_tissue
handles.output = hObject;
handles.dialogcut_handle=[];
% Update handles structure

set(handles.pushbutton22,'enable','off');
set(handles.slider1,'enable','off');
set(handles.slider2,'enable','off');
set(handles.pushbutton21,'enable','off');
set(handles.HSVDbutton,'enable','off');
set(handles.Hankelbutton,'enable','off');
set(handles.Loewnerbutton,'enable','off');
set(handles.cutbutton,'enable','off');
set(handles.uncutbutton,'enable','off')

set(handles.unfiltbutton,'enable','off')
set(handles.msbutton,'enable','off')
set(handles.hzbutton,'enable','off')
set(handles.ppmbutton,'enable','off')
set(handles.realbutton,'enable','off')
set(handles.imagbutton,'enable','off')
set(handles.absbutton,'enable','off')
set(handles.NNTFbutton,'enable','off')
set(handles.savebutton,'enable','off')
set(handles.slidervalue,'enable','off')
set(handles.edit3,'enable','off')
set(handles.edit5,'enable','off')
set(handles.edit6,'enable','off')
set(handles.edit8,'enable','off')
set(handles.edit9,'enable','off')
set(handles.edit10,'enable','off')

guidata(hObject, handles);
    

% UIWAIT makes MRSI_tissue wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MRSI_tissue_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [fn p]=uigetfile({'*.*','MAT-files (*.mat)'; ...
    '*.txt','jMRUI files(*.txt)'; ...
    '*.sdat','Philips files(*.sdat)'; ...
    '*.spar','Philips files(*.spar)'; ...
    '*.rda','Siemens files(*.rda)'; ...
    '*.*',  'All Files (*.*)'},'Load file from','MultiSelect', 'on');
    if isequal(fn,0) || isequal(p,0)
        return
    end
    [~,~,ext] = fileparts(fn);
    filepath = strcat(p,fn);
    % if saved .mat file is loaded
%     if((length(ext) ==  4) && (ext == '.mat'))
    if(strcmp(ext,'.mat'))
        loadsig = load(filepath); 
        signal = loadsig.signal;
        step = loadsig.step;
        frequency = loadsig.frequency;
        ndp = loadsig.ndp;
        begin = loadsig.begin;
        signalnames = loadsig.signalnames;
        procres = loadsig.procres;
        classopt = loadsig.classopt;
        simsig = loadsig.simsig;
        position = loadsig.position;
%         averages = loadsig.averages;
%         process_info = 1
        if(isfield(loadsig,'process_info'))
            water_supressed_algo = loadsig.process_info.water_supressed_algo; %0--not supressed, 1--HSVD, 2--Hankel_tensor, 3-- Loewner tensor
            tissue_diff_applied  = loadsig.process_info.tissue_diff_applied; %0-- not applied, 1-- applied;
            if(tissue_diff_applied == 1)
                sources = loadsig.process_info.sources;
                abudances = loadsig.process_info.abundances;
                freqregion = loadsig.process_info.freqregion;
                source_list = loadsig.process_info.source_list;
            end
        else
            water_supressed_algo = 0;
            tissue_diff_applied = 0;
        end
    else % if raw data is loaded
        [signal,step,frequency,ndp,begin,signalnames,procres,classopt,simsig,position] = loadsignals(fn,p);
        water_supressed_algo = 0;
        tissue_diff_applied = 0;
    end
    
    handles.f = -1/(2*step/1000):1/((ndp-1)*step/1000):1/(2*step/1000);
    handles.f2 = kHz2ppm(handles.f,frequency*1000);
    handles.nos = size(signal,1);
    handles.col = numel(unique((position(:,2))));
    handles.row = numel(unique((position(:,3))));
    spectra = zeros(handles.nos,ndp);             
    for i=1:size(signal,1)
        spectra(i,:) = fftshift(fft(signal(i,:)));
    end
            
    handles.sec = [1:ndp]*step; % convert the ndp to ms
    handles.signal = signal;
    handles.step = step;
    handles.frequency = frequency; 
    handles.ndp = ndp; 
    handles.begin = begin;
    handles.signalnames = signalnames; 
    handles.procres = procres;
    handles.classopt = classopt;
    handles.simsig = simsig;
    handles.position = position;
    handles.new_position = position;
    handles.spectra= spectra;
    handles.current_signal = 1;
    
    handles.location_map = zeros(handles.row,handles.col);
    handles.cropped = 0; %0- not croped, 1- croped
    handles.crop = [0,0,0,0];
    %variables to be plotted
    handles.signal_croped = signal;
    handles.signal_croped_filtered = handles.signal_croped;
    handles.spectra_croped=spectra;
    handles.spectra_croped_filtered=handles.spectra_croped;
%     handles.plot_ms_hn_ppm = 1; % 1 -ms, 2- hz, 3- ppm                  
%     handles.real_imag_abs = 1; %1- real, 2-imag, 3 abs
    set(handles.edit3, 'String', handles.nos);
    
    toplot=signal;
    handles.toplot=toplot;
    newtoplot=real(toplot);
    handles.newtoplot=newtoplot;

    % optimal slider step (.mat data vs. raw data)
%     if (length(ext) ==  4) && (ext == '.mat')
    if(strcmp(ext,'.mat'))
        set(handles.slider1,'Min',1,'Max',handles.nos,'SliderStep',[1/(handles.nos-1) 2/(handles.nos-1)]);
    else
        set(handles.slider1,'Min',1,'Max',handles.nos)
    end

    handles.scrollbarValue = 1;
    set(handles.slidervalue, 'String', handles.scrollbarValue);  
    set(handles.msbutton,'Value',1);
    set(handles.hzbutton,'Value',0);
    set(handles.ppmbutton,'Value',0);
    
    set(gcf,'CurrentAxes',handles.axes1);    
    if(get(handles.msbutton,'Value') == 1)
        plot(handles.axes1,handles.sec,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xlim',[0 handles.ndp*handles.step]);
        set(get(handles.axes1, 'xlabel'), 'string', 'ms'); axis tight;
    end
    if(get(handles.hzbutton,'Value') == 1)
        plot(handles.axes1,handles.f,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xdir','reverse');
        set(get(handles.axes1, 'xlabel'), 'string', 'Hz');axis tight;
    end
    if(get(handles.ppmbutton,'Value') == 1)
        plot(handles.axes1,handles.f2,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xdir','reverse');
        set(get(handles.axes1, 'xlabel'), 'string', 'ppm');axis tight;
    end
    
    temp = handles.location_map;
    temp(1) = 1;
    rgb_temp = cat(3, temp, zeros(size(temp)), zeros(size(temp))); 
    set(gcf,'CurrentAxes',handles.axes5)
    imagesc(rgb_temp)
    handles.img_loaded  = 0;
    
    handles.row1=handles.row;
    handles.col1=handles.col;
    
    handles.water_supressed_algo = water_supressed_algo;
    handles.tissue_diff_applied = tissue_diff_applied;
    if(tissue_diff_applied == 1)
        handles.sources = sources;
        handles.abudances = abudances;
        handles.freqregion = freqregion;
        

        freq_axis = -1/(2*step/1000):1/((ndp-1)*step/1000):1/(2*step/1000);
        ppmaxis = kHz2ppm(freq_axis,frequency*1000);

        freqindex = [find(ppmaxis>min(freqregion),1,'first'), find(ppmaxis<max(freqregion),1,'last')];
        ppm_axis_seg = ppmaxis(freqindex(1):freqindex(2));
        
        
        set(gcf,'CurrentAxes',handles.axes4)
        plot(handles.axes4,ppm_axis_seg,real(handles.sources(:,1)));set(handles.axes4,'xdir','reverse');axis tight
        set(get(handles.axes4, 'xlabel'), 'string', 'ppm')
        R = size(sources,2);
        set(handles.slider2,'Min',1,'Max',R,'Value',1,'SliderStep',[1/(R-1) 2/(R-1)]); 
        temp_img = reshape(handles.abudances(1,:)',handles.row1,handles.col1);
        set(gcf,'CurrentAxes',handles.axes2)
        imagesc(temp_img');

        for i = 1:R
            handles.source_list{i} = ['source ',num2str(i)];
        end
        set(handles.edit5, 'String', handles.source_list{1});
        set(handles.edit5,'enable','on');         
        set(handles.slider2,'enable','on'); 
        set(handles.savebutton,'enable','on'); 
    else
        cla(handles.axes4);
        cla(handles.axes2);
        set(handles.edit5, 'String', [],'enable','off');   
        set(handles.slider2,'enable','off'); 
        set(handles.savebutton,'enable','off');
    end


    set(handles.pushbutton22,'enable','on'); %load img
    set(handles.slider1,'enable','on');
    set(handles.pushbutton21,'enable','on'); %save 
    set(handles.HSVDbutton,'enable','on');
    set(handles.Hankelbutton,'enable','on');
    set(handles.Loewnerbutton,'enable','on');
    set(handles.cutbutton,'enable','on');
    set(handles.msbutton,'enable','on')
    set(handles.hzbutton,'enable','on')
    set(handles.ppmbutton,'enable','on')
    set(handles.realbutton,'enable','on')
    set(handles.imagbutton,'enable','on')
    set(handles.absbutton,'enable','on')
    set(handles.NNTFbutton,'enable','on')
    temp_str = ['GRID:',num2str(handles.row),'X',num2str(handles.row),'; Crop: L-0; R-0; T-0; B-0'];
    set(handles.edit6, 'String', temp_str);

    
    if(strcmp(ext,'.mat'))
    if(isfield(loadsig,'process_info'))
        handles.img_loaded = loadsig.process_info.img_loaded;
        if(handles.img_loaded == 1)
            handles.back_img = loadsig.process_info.back_img;            
            if(loadsig.process_info.cropped == 0)
                out_num = handles.scrollbarValue;
                image_x = size(handles.back_img,1);
                image_y = size(handles.back_img,2);
                aa = round(linspace(0,image_x,handles.row+1));
                aa1 = round(linspace(0,image_y,handles.col+1));
                handles.image_block = zeros(size(handles.back_img));
                image_block1 = zeros(size(handles.back_img));
                count = 1;
            %     handles.col
            %     handles.row
                for j = 1:handles.col
                    for i= 1:handles.row
                        if(handles.location_map(i,j) == 2)
                            image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
                        end
                        if(out_num == count)
                            handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
                        end
                        count = count +1;
                    end
                end
                set(gcf,'CurrentAxes',handles.axes5)
                rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
                imshow(rgbImage,'parent',handles.axes5);
                hold on
                H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
                h = imshow(H_out,'parent',handles.axes5);
                hold off
                alpha_data = 0.5*(handles.image_block + image_block1);
                set(h, 'AlphaData', alpha_data);
            end            
        end
        if(loadsig.process_info.water_supressed_algo ~= 0)
            set(handles.HSVDbutton,'enable','off');
            set(handles.Hankelbutton,'enable','off');
            set(handles.Loewnerbutton,'enable','off');
            set(handles.cutbutton,'enable','off');
            set(handles.uncutbutton,'enable','off')
            set(handles.unfiltbutton,'enable','off')
        end
        if(loadsig.process_info.cropped == 1)            
            set(handles.cutbutton,'enable','off');
            set(handles.uncutbutton,'enable','off')

            handles.row = loadsig.process_info.row;
            handles.col = loadsig.process_info.col;            
            handles.crop = loadsig.process_info.crop;
            handles.location_map = zeros(handles.row,handles.col);
            handles.location_map(:,1:handles.crop(1)) = 2;
            handles.location_map(:,handles.col-handles.crop(2)+1:handles.col) = 2;
            handles.location_map(1:handles.crop(3),:) = 2;
            handles.location_map(handles.row-handles.crop(4)+1:handles.row,:) = 2;
            handles.cropped = 1;
            temp_str = ['GRID:',num2str(handles.row),'X',num2str(handles.row),'; Crop: L-',num2str(handles.crop(1)),'; R-',num2str(handles.crop(2)),'; T-',...
                num2str(handles.crop(3)),'; B-',num2str(handles.crop(4))];
            set(handles.edit6, 'String', temp_str);
            if(handles.cropped == 1)
                col_i = ceil(handles.scrollbarValue/handles.row1);
                row_i = mod(handles.scrollbarValue,handles.row1);
                if(row_i == 0)
                    row_i = handles.row1;
                end
                col_o = col_i + handles.crop(1);
                row_o = row_i + handles.crop(3);
                out_num = (col_o-1)*handles.row + row_o;
            else
                out_num = handles.scrollbarValue;
            end
            if(handles.img_loaded  == 1)    
                image_x = size(handles.back_img,1);
                image_y = size(handles.back_img,2);
                aa = round(linspace(0,image_x,handles.row+1));
                aa1 = round(linspace(0,image_y,handles.col+1));
                handles.image_block = zeros(size(handles.back_img));
                image_block1 = zeros(size(handles.back_img));
                count = 1;
            %     handles.col
            %     handles.row
                for j = 1:handles.col
                    for i= 1:handles.row
                        if(handles.location_map(i,j) == 2)
                            image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
                        end
                        if(out_num == count)
                            handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
                        end
                        count = count +1;
                    end
                end
                set(gcf,'CurrentAxes',handles.axes5)
                rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
                imshow(rgbImage,'parent',handles.axes5);
                hold on
                H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
                h = imshow(H_out,'parent',handles.axes5);
                hold off
                alpha_data = 0.5*(handles.image_block + image_block1);
                set(h, 'AlphaData', alpha_data);
            else
                temp = zeros(size(handles.location_map));
                temp(out_num) = 1;
                H_out = cat(3, temp, zeros(size(handles.location_map)),handles.location_map == 2);   
            %     set(gcf,'CurrentAxes',handles.axes5)
                imshow(H_out,'parent',handles.axes5);    
            end            
        end    
    end
    end
    
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ScrollBarMoved(handles)
% global kuts
handles.scrollbarValue = round(get(handles.slider1, 'Value'));
set(handles.slidervalue, 'String', handles.scrollbarValue);
%plot(handles.axes1,handles.newtoplot(handles.scrollbarValue,:))
% set(handles.cutsum,'String',sum(kuts))
set(gcf,'CurrentAxes',handles.axes1)
if get(handles.msbutton,'Value') == 1
    plot(handles.axes1,handles.sec,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xlim',[0 handles.ndp*handles.step]);
    set(get(handles.axes1, 'xlabel'), 'string', 'ms'); axis tight;
end
if get(handles.hzbutton,'Value') == 1
    plot(handles.axes1,handles.f,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xdir','reverse');
    set(get(handles.axes1, 'xlabel'), 'string', 'Hz'); axis tight;
end
if get(handles.ppmbutton,'Value') == 1
    plot(handles.axes1,handles.f2,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xdir','reverse');%set(gca,'XTick',-1024:129:1024);set(gca,'XTickLabel',-3.15:1:12.5);set(gca,'xlim',[-1024 1024])
    set(get(handles.axes1, 'xlabel'), 'string', 'ppm'); axis tight;
end

if(handles.cropped == 1)
    col_i = ceil(handles.scrollbarValue/handles.row1);
    row_i = mod(handles.scrollbarValue,handles.row1);
    if(row_i == 0)
        row_i = handles.row1;
    end
    col_o = col_i + handles.crop(1);
    row_o = row_i + handles.crop(3);
    out_num = (col_o-1)*handles.row + row_o;
else
    out_num = handles.scrollbarValue;
end

if(handles.img_loaded  == 1)
    image_x = size(handles.back_img,1);
    image_y = size(handles.back_img,2);
    aa = round(linspace(0,image_x,handles.row+1));
    aa1 = round(linspace(0,image_y,handles.col+1));
    handles.image_block = zeros(size(handles.back_img));
    set(gcf,'CurrentAxes',handles.axes5)  
    image_block1 = zeros(size(handles.back_img));
    count = 1;
    for j = 1:handles.col
        for i= 1:handles.row
            if(handles.location_map(i,j) == 2)
               image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
            end
            if(out_num == count)
                handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
            end
            count = count +1;
        end
    end
    rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
    imshow(rgbImage,'parent',handles.axes5);
    hold on
    H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
    h = imshow(H_out,'parent',handles.axes5);
    hold off
    alpha_data = 0.5*(handles.image_block + image_block1);
    set(h, 'AlphaData', alpha_data);
else
    temp = zeros(size(handles.location_map));
    temp(out_num) = 1;
    H_out = cat(3, temp, zeros(size(handles.location_map)),handles.location_map == 2);   
%     set(gcf,'CurrentAxes',handles.axes5)
    imshow(H_out,'parent',handles.axes5);
end   
guidata(hObject, handles);    
    
    
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slidervalue_Callback(hObject, eventdata, handles)
% hObject    handle to slidervalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slidervalue as text
%        str2double(get(hObject,'String')) returns contents of slidervalue as a double

% --- Executes during object creation, after setting all properties.
function slidervalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidervalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in msbutton.
function msbutton_Callback(hObject, eventdata, handles)
% hObject    handle to msbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.hzbutton,'Value',0)
% set(handles.ppmbutton,'Value',0)

% Hint: get(hObject,'Value') returns toggle state of msbutton


% --- Executes on button press in hzbutton.
function hzbutton_Callback(hObject, eventdata, handles)
% hObject    handle to hzbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.msbutton,'Value',0)
% set(handles.ppmbutton,'Value',0) 

% plot(handles.axes1,handles.toplot(round(handles.scrollbarValue),:))

  
% Hint: get(hObject,'Value') returns toggle state of hzbutton


% --- Executes on button press in ppmbutton.
function ppmbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ppmbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.msbutton,'Value',0)
% set(handles.hzbutton,'Value',0) 
% get(hObject,'Value')
% plot(handles.axes1,handles.toplot(round(handles.scrollbarValue),:))

% Hint: get(hObject,'Value') returns toggle state of ppmbutton


% --- Executes when selected object is changed in plotgroup.
function plotgroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in plotgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selob=get(handles.plotgroup,'SelectedObject');
select=get(selob,'Tag');
switch select % Get Tag of selected object.
	case 'msbutton'
		% Code for when radiobutton1 is selected.
         handles.toplot=handles.signal_croped_filtered;
	case 'hzbutton'
		% Code for when radiobutton2 is selected.
        handles.toplot=handles.spectra_croped_filtered;
    case 'ppmbutton'
        % Code for when radiobutton2 is selected.
        handles.toplot=handles.spectra_croped_filtered;   

end
guidata(hObject, handles);
showgroup_SelectionChangedFcn(hObject, eventdata, handles)



% --- Executes when selected object is changed in showgroup.

function showgroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in showgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selobShow=get(handles.showgroup,'SelectedObject');
selectShow=get(selobShow,'Tag');

switch selectShow % Get Tag of selected object.
	case 'realbutton'
		% Code for when radiobutton1 is selected.
        handles.newtoplot=real(handles.toplot);      
	case 'imagbutton'
		% Code for when radiobutton2 is selected.
        handles.newtoplot=imag(handles.toplot);   
    case 'absbutton'
        % Code for when radiobutton2 is selected.
        handles.newtoplot=abs(handles.toplot);   
end
%plot(handles.axes1,handles.newtoplot(round(handles.scrollbarValue),:))
guidata(hObject, handles);

set(gcf,'CurrentAxes',handles.axes1)
if get(handles.msbutton,'Value') == 1
    plot(handles.axes1,handles.sec,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xlim',[0 handles.ndp*handles.step]);
    set(get(handles.axes1, 'xlabel'), 'string', 'ms'); axis tight;
end
if get(handles.hzbutton,'Value') == 1
    plot(handles.axes1,handles.f,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xdir','reverse');
    set(get(handles.axes1, 'xlabel'), 'string', 'Hz'); axis tight;
end
if get(handles.ppmbutton,'Value') == 1
    plot(handles.axes1,handles.f2,handles.newtoplot(round(handles.scrollbarValue),:));set(gca,'xdir','reverse');
    set(get(handles.axes1, 'xlabel'), 'string', 'ppm'); axis tight;
end

% --- Executes on button press in cutbutton.
function cutbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.dialogcut_handle )
   % create GUI-3  
%    handles.dialogcut_handle  = dialogcut; 
prompt = {'left','right','top','bottom'};
title = 'Select #voxels to remove from each side';
dims = [1 35];
definput = {'0','0','0','0'};
answer = inputdlg(prompt,title,dims,definput);
crop = str2double(answer);
if(isempty(crop))
    return
end
% crop

if(crop(1) >= handles.col)
   uiwait(errordlg('Left crop value is out of limit','Error','modal'));
   return
end
if(crop(2) >= handles.col)
    uiwait(errordlg('Right crop value is out of limit','Error','modal'));
   return
end
if(crop(3) >= handles.row)
    uiwait(errordlg('Top crop value is out of limit','Error','modal'));
   return
end
if(crop(4) >= handles.row)
    uiwait(errordlg('Bottom crop value is out of limit','Error','modal'));
   return
end
if((crop(1)+crop(2)) >= handles.col)
    uiwait(errordlg('Left and right crop values are out of limit','Error','modal'));
   return
end
if((crop(3)+crop(4)) >= handles.row)
    uiwait(errordlg('Top and Bottom crop values are out of limit','Error','modal'));
   return
end

%Make values cut the right parts of signal

signal=handles.signal;
position=handles.position;
spectra=handles.spectra;
% signalppm=handles.signalppm;
tot_nos = size(signal,1);
col = numel(unique((position(:,2))));
row = numel(unique((position(:,3))));
row1 = row-(crop(2)+ crop(1));
col1 = col-(crop(3)+ crop(4));

handles.row1=row1;
handles.col1=col1;

if size(position,2)<5
    position(:,5) = 1;
end
if ~isempty(crop)
%     handles.crop = crop;
    matr = reshape(1:tot_nos,col,row);
    position(matr(1:crop(1),:),5) = 0;
    position(matr(:,1:crop(3)),5) = 0;
    position(matr(end-crop(2)+1:end,:),5) = 0;
    position(matr(:,end-crop(4)+1:end),5) = 0;
end

handles.signal_croped = signal(position(:,5)==1,:); 
handles.spectra_croped=spectra(position(:,5)==1,:);

handles.signal_croped_filtered=handles.signal_croped;
handles.spectra_croped_filtered=handles.spectra_croped;


% position calculation of new positions
nfids = row1*col1;
handles.nfids = nfids;
pos1 =[]; 
pos2 = []; 
pos3 = [];
for i = 1:row1
        pos1 = [pos1; [1:col1]'];
        pos2 = [pos2; i*ones(1,col1)'];
end
pos3 = [pos3; ones(row1*col1,1)];
new_position = [[1:nfids]' pos1 pos2 pos3];
handles.new_position = new_position;
handles.position = position;

nos = size(handles.signal_croped,1);
%set slider to 1 to prevent out of range
set(handles.slider1, 'Value',1)
set(handles.slidervalue, 'String', 1);
handles.scrollbarValue = 1;
%update plots
%plotgroup_SelectionChangedFcn(hObject, eventdata, handles)

set(handles.edit3, 'String', nos);
set(handles.slider1,'Min',1,'Max',nos,'SliderStep',[1/(nos-1) 2/(nos-1)]);
% guidata(hObject, handles);
handles.crop = crop;
handles.location_map(:,1:crop(1)) = 2;
handles.location_map(:,handles.col-crop(2)+1:handles.col) = 2;
handles.location_map(1:crop(3),:) = 2;
handles.location_map(handles.row-crop(4)+1:handles.row,:) = 2;
handles.cropped = 1;
temp_str = ['GRID:',num2str(handles.row),'X',num2str(handles.row),'; Crop: L-',num2str(handles.crop(1)),'; R-',num2str(handles.crop(2)),'; T-',...
        num2str(handles.crop(3)),'; B-',num2str(handles.crop(4))];
set(handles.edit6, 'String', temp_str);
            
    if(handles.cropped == 1)
        col_i = ceil(handles.scrollbarValue/handles.row1);
        row_i = mod(handles.scrollbarValue,handles.row1);
        if(row_i == 0)
            row_i = handles.row1;
        end
        col_o = col_i + handles.crop(1);
        row_o = row_i + handles.crop(3);
        out_num = (col_o-1)*handles.row + row_o;
    else
        out_num = handles.scrollbarValue;
    end
set(gcf,'CurrentAxes',handles.axes5)      
if(handles.img_loaded  == 1)    
    image_x = size(handles.back_img,1);
    image_y = size(handles.back_img,2);
    aa = round(linspace(0,image_x,handles.row+1));
    aa1 = round(linspace(0,image_y,handles.col+1));
    handles.image_block = zeros(size(handles.back_img));
    image_block1 = zeros(size(handles.back_img));
    count = 1;
%     handles.col
%     handles.row
    for j = 1:handles.col
        for i= 1:handles.row
            if(handles.location_map(i,j) == 2)
                image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
            end
            if(out_num == count)
                handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
            end
            count = count +1;
        end
    end
    
    rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
    imshow(rgbImage,'parent',handles.axes5);
    hold on
    H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
    h = imshow(H_out,'parent',handles.axes5);
    hold off
    alpha_data = 0.5*(handles.image_block + image_block1);
    set(h, 'AlphaData', alpha_data);
else
    temp = zeros(size(handles.location_map));
    temp(out_num) = 1;
    H_out = cat(3, temp, zeros(size(handles.location_map)),handles.location_map == 2);   
%     set(gcf,'CurrentAxes',handles.axes5)
    imshow(H_out,'parent',handles.axes5);    
end
    set(handles.cutbutton,'enable','off');
    set(handles.uncutbutton,'enable','on')
    
    
    
guidata(hObject, handles); 

% hahahhaha
end

plotgroup_SelectionChangedFcn(hObject, eventdata, handles)


% --- Executes on button press in preprocessingbutton.
function preprocessingbutton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessingbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.preprocessing_handle = preprocessingdialog; 


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
  
guidata(hObject, handles);
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in HSVDbutton.
function HSVDbutton_Callback(hObject, eventdata, handles)
% hObject    handle to HSVDbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Model order','ppm min','ppm max'};
title = 'Filter specifications';
dims = [1 35];
definput = {'25','0.25', '4.2'};
answerhsvd = inputdlg(prompt,title,dims,definput);

if(isempty(answerhsvd))
    return
end

handles.Morder=str2double(answerhsvd(1));
freq_pts = [-1/(2*handles.step/1000),1/(2*handles.step/1000)];
ppm_pts = kHz2ppm(freq_pts,handles.frequency*1000);
ppmrange=[str2double(answerhsvd(2)) str2double(answerhsvd(3))];
if((ppmrange(1)>=ppm_pts(2))||(ppmrange(1) < ppm_pts(1)))   
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return
end
if((ppmrange(2)>=ppm_pts(2))||(ppmrange(2) < ppm_pts(1)))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end
if(ppmrange(2)<=ppmrange(1))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end


Kest=handles.Morder;

waitfunc = waitbar(0,'Please wait...');
[handles.signal_croped_filtered,handles.spectra_croped_filtered] = water_removal_HSVD(handles.signal_croped,Kest,handles.step,handles.frequency,ppmrange);
    
[sigrow,sigcol]=size(handles.signal_croped_filtered);   
handles.spectra_croped_filtered=zeros(sigrow,sigcol);             
for i=1:size(handles.signal_croped_filtered,1)
    %spectra(i,:) = fliplr(fftshift(fft(signal(i,:))));  % from time domain to real part of the freq domain
    handles.spectra_croped_filtered(i,:) = fftshift(fft(handles.signal_croped_filtered(i,:)));
    waitbar(i/size(handles.signal_croped,1),waitfunc)
end          
% handles.signalppm_croped_filtered=handles.spectra_croped_filtered;
handles.water_supressed_algo = 1;
delete(waitfunc)
set(handles.HSVDbutton,'enable','off');
set(handles.Hankelbutton,'enable','off');
set(handles.Loewnerbutton,'enable','off');
set(handles.unfiltbutton,'enable','on')

guidata(hObject, handles);
plotgroup_SelectionChangedFcn(hObject, eventdata, handles)

% --- Executes on button press in Hankelbutton.
function Hankelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Hankelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



prompt = {'Model order','ppm min','ppm max','length'};
title = 'Filter specifications';
dims = [1 35];
definput = {'100','0.25', '4.2','1024'};
answerhankel = inputdlg(prompt,title,dims,definput);

if(isempty(answerhankel))
    return
end

handles.Morder=str2double(answerhankel(1));
freq_pts = [-1/(2*handles.step/1000),1/(2*handles.step/1000)];
ppm_pts = kHz2ppm(freq_pts,handles.frequency*1000);
ppmrange=[str2double(answerhankel(2)) str2double(answerhankel(3))];
if((ppmrange(1)>=ppm_pts(2))||(ppmrange(1) < ppm_pts(1)))   
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return
end
if((ppmrange(2)>=ppm_pts(2))||(ppmrange(2) < ppm_pts(1)))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end
if(ppmrange(2)<=ppmrange(1))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end
Kest=handles.Morder;
length=str2double(answerhankel(4));
waitfunc = waitbar(0,'Please wait...');

size(handles.signal_croped_filtered)
% for i=1:size(handles.signal_croped,1)
    hankelresult = Hankel_water_supression(handles.signal_croped.',Kest,length,ppmrange,handles.step,handles.frequency);
    handles.signal_croped_filtered = hankelresult.';
% end

for i=1:size(handles.signal_croped_filtered,1)
%     handles.spectra_croped_filtered(i,:) = fliplr(fftshift(fft(handles.signal_croped_filtered(i,:))));
    handles.spectra_croped_filtered(i,:) = fftshift(fft(hankelresult(:,i))).';
    waitbar(i/size(handles.signal_croped,1),waitfunc)
end
% handles.signalppm_croped_filtered=handles.spectra_croped_filtered;
handles.water_supressed_algo = 2;
delete(waitfunc);

set(handles.HSVDbutton,'enable','off');
set(handles.Hankelbutton,'enable','off');
set(handles.Loewnerbutton,'enable','off');
set(handles.unfiltbutton,'enable','on')

guidata(hObject, handles);
plotgroup_SelectionChangedFcn(hObject, eventdata, handles)

% --- Executes on button press in Loewnerbutton.
function Loewnerbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Loewnerbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prompt = {'Model order','ppm min','ppm max'};
title = 'Filter specifications';
dims = [1 35];
definput = {'50','0.25', '4.2'};
answerloewner = inputdlg(prompt,title,dims,definput);

if(isempty(answerloewner))
    return
end

handles.Morder=str2double(answerloewner(1));
freq_pts = [-1/(2*handles.step/1000),1/(2*handles.step/1000)];
ppm_pts = kHz2ppm(freq_pts,handles.frequency*1000);
ppmrange=[str2double(answerloewner(2)) str2double(answerloewner(3))];
if((ppmrange(1)>=ppm_pts(2))||(ppmrange(1) < ppm_pts(1)))   
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return
end
if((ppmrange(2)>=ppm_pts(2))||(ppmrange(2) < ppm_pts(1)))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end
if(ppmrange(2)<=ppmrange(1))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end

R=handles.Morder;

freqregion = [6 0];
spectra=handles.spectra_croped;
freq_axis = -1/(2*handles.step/1000):1/((handles.ndp-1)*handles.step/1000):1/(2*handles.step/1000);
ppmaxis = kHz2ppm(freq_axis,handles.frequency*1000);
freqindex = [find(ppmaxis>min(freqregion),1,'first'), find(ppmaxis<max(freqregion),1,'last')];


% [rowspec,colspec]=size(spectra);   
% spectra_loew = zeros(rowspec,colspec); 
% for i=1:size(spectra,1)
%     spectra_loew(i,:) = ifftshift(spectra(i,:));
% end
            


waitfunc = waitbar(0,'Please wait...');
[handles.spectra_croped_filtered,source_mat,H_mat] = Loewner_water_supression(spectra.',R,ppmaxis,freqindex,ppmrange);

handles.spectra_croped_filtered = handles.spectra_croped_filtered.';
    for i=1:size(handles.spectra_croped_filtered,1)
        handles.signal_croped_filtered(i,:) = ifft(ifftshift(handles.spectra_croped_filtered(i,:))); 
    end
% handles.signalppm_croped_filtered=handles.spectra_croped_filtered;

handles.water_supressed_algo = 3;
delete(waitfunc);
set(handles.HSVDbutton,'enable','off');
set(handles.Hankelbutton,'enable','off');
set(handles.Loewnerbutton,'enable','off');
set(handles.unfiltbutton,'enable','on')

guidata(hObject, handles);
plotgroup_SelectionChangedFcn(hObject, eventdata, handles);

% handles.spectra_croped_filtered=signalfil_Lowner_FD;
% handles.spectra_croped_filtered=signalfil_Lowner_FD;
% handles.signal_croped_filtered=[]
% for i=1:size(handles.spectra_croped_filtered,1)
%         handles.signal_croped_filtered(i,:) = ifft(fftshift(fliplr(handles.signal_croped_filtered(i,:))));  
%     end


% --- Executes on button press in unfiltbutton.
function unfiltbutton_Callback(hObject, eventdata, handles)
% hObject    handle to unfiltbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Loewnerbutton,'enable','on');
set(handles.Hankelbutton,'enable','on');
set(handles.HSVDbutton,'enable','on');
set(handles.unfiltbutton,'enable','off')

if(handles.water_supressed_algo~=0)
    handles.tissue_diff_applied = 0;
    handles.signal_croped_filtered = handles.signal_croped;
    handles.spectra_croped_filtered = handles.spectra_croped;
%     handles.signalppm_croped_filtered = handles.signalppm_croped;
    guidata(hObject, handles);
end
plotgroup_SelectionChangedFcn(hObject, eventdata, handles)



% --- Executes on button press in uncutbutton.
function uncutbutton_Callback(hObject, eventdata, handles)  
% hObject    handle to uncutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.signal_croped = handles.signal;
handles.spectra_croped= handles.spectra;
% handles.signalppm_croped= handles.signalppm;

handles.signal_croped_filtered=handles.signal_croped;
handles.spectra_croped_filtered=handles.spectra_croped;
% handles.signalppm_croped_filtered=handles.signalppm_croped;


handles.water_supressed_algo = 0;
handles.tissue_diff_applied = 0;


%set slider to 1 to prevent out of range
set(handles.slider1, 'Value',1)
set(handles.slidervalue, 'String', 1);
handles.scrollbarValue = 1;

nos = size(handles.signal_croped,1);
set(handles.edit3, 'String', nos);
set(handles.slider1,'Min',1,'Max',nos,'SliderStep',[1/(nos-1) 2/(nos-1)]);

handles.cropped = 0; %0- not croped, 1- croped
handles.crop = [0,0,0,0];
handles.location_map = zeros(handles.row,handles.col);
set(gcf,'CurrentAxes',handles.axes5)  
temp_str = ['GRID:',num2str(handles.row),'X',num2str(handles.row),'; Crop: L-',num2str(handles.crop(1)),'; R-',num2str(handles.crop(2)),'; T-',...
        num2str(handles.crop(3)),'; B-',num2str(handles.crop(4))];
set(handles.edit6, 'String', temp_str);

if(handles.img_loaded  == 1)    
    image_x = size(handles.back_img,1);
    image_y = size(handles.back_img,2);
    aa = round(linspace(0,image_x,handles.row+1));
    aa1 = round(linspace(0,image_y,handles.col+1));
    handles.image_block = zeros(size(handles.back_img));
    image_block1 = zeros(size(handles.back_img));
    count = 1;
%     out_num = 1;
    handles.image_block(aa(1)+1:aa(2),aa1(1)+1:aa1(2)) = 1;
%     handles.col
%     handles.row
%     for j = 1:handles.col
%         for i= 1:handles.row
%             if(handles.location_map(i,j) == 2)
%                 image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
%             end
%             if(out_num == count)
%                 handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
%             end
%             count = count +1;
%         end
%     end
    
    rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
    imshow(rgbImage,'parent',handles.axes5);
    hold on
    H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
    h = imshow(H_out,'parent',handles.axes5);
    hold off
    alpha_data = 0.5*(handles.image_block + image_block1);
    set(h, 'AlphaData', alpha_data);
else
    temp = zeros(size(handles.location_map));
    temp(1) = 1;
    H_out = cat(3, temp, zeros(size(handles.location_map)),handles.location_map == 2);   
%     set(gcf,'CurrentAxes',handles.axes5)
    imshow(H_out,'parent',handles.axes5);    
end  

    cla(handles.axes4);
    cla(handles.axes2);    
    set(handles.uncutbutton,'enable','off');
    set(handles.cutbutton,'enable','on')
    set(handles.slider2,'enable','off'); 
    set(handles.HSVDbutton,'enable','on');
    set(handles.Hankelbutton,'enable','on');
    set(handles.Loewnerbutton,'enable','on');
    set(handles.unfiltbutton,'enable','off')
    set(handles.savebutton,'enable','off')
    set(handles.edit5,'enable','off')

guidata(hObject, handles);
plotgroup_SelectionChangedFcn(hObject, eventdata, handles)

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NNTFbutton.
function NNTFbutton_Callback(hObject, eventdata, handles)
% hObject    handle to NNTFbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nmalise = 1;
negative_constrain = 5;
freqregion = [0.25, 4.2];
wei = [1 0.05]; % default weight for L1 regularizarion

title = 'NCPD specifications';
prompt = {'ppm min','ppm max','Rank','lambda'};
% prompt = {'ppm min','ppm max','lambda'};

dims = [1 35];
definput = {'0.25', '4.2','0','0.05'}; % rank 0 means detect automaticallt max number of sources is 9
% definput = {'0.25', '4.2','0.05'};

answernntf = inputdlg(prompt,title,dims,definput);

if(isempty(answernntf))
    return
end
freqregion=[str2double(answernntf(1)) str2double(answernntf(2))];

freq_pts = [-1/(2*handles.step/1000),1/(2*handles.step/1000)];
ppm_pts = kHz2ppm(freq_pts,handles.frequency*1000);
if((freqregion(1)>=ppm_pts(2))||(freqregion(1) < ppm_pts(1)))   
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return
end
if((freqregion(2)>=ppm_pts(2))||(freqregion(2) < ppm_pts(1)))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end
if(freqregion(2)<=freqregion(1))
   uiwait(errordlg('Selected ppm range is out of limit','Error','modal'));
   return    
end


rank = str2double(answernntf(3));
wei=[1, str2double(answernntf(4))];

handles.freqregion = freqregion;

signal=handles.signal_croped_filtered;
frequency=handles.frequency;
step=handles.step;
ndp=handles.ndp;
row1=handles.row1;
col1=handles.col1;
signal = align_spectra(signal,frequency,step,ndp);
freq_axis = -1/(2*step/1000):1/((ndp-1)*step/1000):1/(2*step/1000);
ppmaxis = kHz2ppm(freq_axis,frequency*1000);

freqindex = [find(ppmaxis>min(freqregion),1,'first'), find(ppmaxis<max(freqregion),1,'last')];
ppm_axis_seg = ppmaxis(freqindex(1):freqindex(2));
nos = size(signal,1);

for i=1:nos
      spectra_uncmp(i,:) = fftshift(fft(signal(i,:)));  % from time domain to real part of the freq domain
end

spectra_uncmp = spectra_uncmp(:,freqindex(1):freqindex(2)).';  % select frequency region of interest
no_of_voxels = size(spectra_uncmp,2);

% normalize spectra
spectra = spectra_uncmp;
norm_fac = sqrt(sum(spectra.*conj(spectra)));
norm_fac = repmat(norm_fac,[size(spectra,1),1]);
nor_spectra = spectra./norm_fac;

win_len = 25; % change this to change window length
junp_val = 5; % change this to change window slide length

[U1,~,~] = svd(nor_spectra);
max_mno_src = 8;
init_sp1 = U1(:,1:max_mno_src);

N = size(spectra,1); % for smoothed spectra
for i = 1:floor(N/junp_val)
    stp = (i-1)*junp_val + 1;
    if(stp+win_len>N)
        break;
    end
    sig_mat_nor(:,i,:) = nor_spectra(stp:stp+win_len,:);
    sig_init(:,i,:) = init_sp1(stp:stp+win_len,:);
end

for i = 1:no_of_voxels  % for tensoer construction 
    temp = sum(sig_mat_nor(:,:,i).*conj(sig_mat_nor(:,:,i)));
    energy_both_nor(:,i) = temp';
    T(:,:,i) = temp'*temp;
end

for i = 1:max_mno_src
    temp = sum(sig_init(:,:,i).*conj(sig_init(:,:,i)));
    energy_both_init(:,i) = temp';
end
energy_both_init = normc(energy_both_init);
cov_mat = cov(energy_both_nor); % for estimating rank.
eig_val = sort(eig(cov_mat),'descend');
R_cut = min(size(energy_both_nor));
while(1)
    for j = 1:R_cut
         if(sum(eig_val(1:j))/sum(eig_val(1:R_cut))>= 0.99)
             break;
         end
    end
    if(j<=8)
        break;
    else
        R_cut = 9;
    end
end
R = j; 
init.src = (energy_both_init(:,1:R));
init.hbd = ((init.src\energy_both_nor)');
if(rank>0)
    R = rank;
end
[sol_ten, out2] = energy_tensor(T,R,1,wei,init); % tensor decomposition
if(exist('l1_ls_nonneg') == 2)
    nnls_type = 1;
else
    nnls_type = 2;
end
[avg_src,H1] = weighted_source(sol_ten.factors.C,nor_spectra,spectra,3,nnls_type,200,3,norm_fac(1:R,:)'); % part output in paper

% assignment = tisue_type_identification(real(avg_src),ppm_axis_seg); % automatic assignment of tissue type 
handles.tissue_diff_applied = 1;

handles.sources = avg_src;
handles.abudances = H1;

set(gcf,'CurrentAxes',handles.axes4)
plot(handles.axes4,ppm_axis_seg,real(handles.sources(:,1)));set(handles.axes4,'xdir','reverse');axis tight
set(get(handles.axes4, 'xlabel'), 'string', 'ppm')

set(handles.slider2,'Min',1,'Max',R,'Value',1,'SliderStep',[1/(R-1) 2/(R-1)]); 
temp_img = reshape(handles.abudances(1,:)',row1,col1);
set(gcf,'CurrentAxes',handles.axes2)
imagesc(temp_img');

for i = 1:R
    handles.source_list{i} = ['source ',num2str(i)];
end
set(handles.edit5, 'String', handles.source_list{1});

% set(handles.edit3, 'String', handles.nos);

set(handles.slider2,'enable','on'); 
set(handles.savebutton,'enable','on')
set(handles.edit5, 'enable','on');    

guidata(hObject, handles);

 
    
% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savebutton

% save all relevant variables in a .mat file


startingFolder = pwd;
% Put in the name of the mat file that the user wants to save.
% defaultFileName = fullfile(startingFolder, '*.mat')
[file,path] = uiputfile('figure','Save file name');
if file == 0
    % User clicked the Cancel button.
    return;
end

[~,only_name] = fileparts(file);
fullFileName = fullfile(path, only_name);

freqregion = handles.freqregion;
frequency = handles.frequency;
step=handles.step;
ndp=handles.ndp;
row1=handles.row1;
col1=handles.col1;
freq_axis = -1/(2*step/1000):1/((ndp-1)*step/1000):1/(2*step/1000);
ppmaxis = kHz2ppm(freq_axis,frequency*1000);

freqindex = [find(ppmaxis>min(freqregion),1,'first'), find(ppmaxis<max(freqregion),1,'last')];
ppm_axis_seg = ppmaxis(freqindex(1):freqindex(2));


for i = 1:size(handles.sources,2)
    hdl = figure();
    subplot(1,2,1);plot(ppm_axis_seg,real(handles.sources(:,i)));set(gca,'xdir','reverse');axis tight
    xlabel('ppm');
    temp_img = reshape(handles.abudances(i,:)',row1,col1);
    subplot(1,2,2);imagesc(temp_img');
    saveas(hdl,[fullFileName,'_',num2str(i)],'png')
    close(hdl);
end


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    scrollbarValue = round(get(handles.slider2, 'Value'));
    freqregion = handles.freqregion;
    frequency = handles.frequency;
    step=handles.step;
    ndp=handles.ndp;
    row1=handles.row1;
    col1=handles.col1;
    freq_axis = -1/(2*step/1000):1/((ndp-1)*step/1000):1/(2*step/1000);
    ppmaxis = kHz2ppm(freq_axis,frequency*1000);
    
    freqindex = [find(ppmaxis>min(freqregion),1,'first'), find(ppmaxis<max(freqregion),1,'last')];
    ppm_axis_seg = ppmaxis(freqindex(1):freqindex(2));

    set(gcf,'CurrentAxes',handles.axes4)    
    plot(handles.axes4,ppm_axis_seg,real(handles.sources(:,scrollbarValue)));set(handles.axes4,'xdir','reverse');axis tight
    set(get(handles.axes4, 'xlabel'), 'string', 'ppm')
    
    temp_img = reshape(handles.abudances(scrollbarValue,:)',row1,col1);
    set(gcf,'CurrentAxes',handles.axes2)
    imagesc(temp_img');
    
    set(handles.edit5, 'String', handles.source_list{scrollbarValue});

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
scrollbarValue = round(get(handles.slider2, 'Value'));
handles.source_list{scrollbarValue} = get(handles.edit5, 'String');
guidata(hObject, handles);   


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles structure with handles and user data (see GUIDATA)

signal = handles.signal_croped_filtered;            
step = handles.step;  
frequency = handles.frequency;  
ndp = handles.ndp;  
begin = handles.begin; 
signalnames = handles.signalnames;
procres = handles.procres;  
classopt = handles.classopt;  
simsig = handles.simsig;
position = handles.new_position;
process_info.water_supressed_algo = handles.water_supressed_algo;
process_info.tissue_diff_applied = handles.tissue_diff_applied;
process_info.img_loaded = handles.img_loaded;
if(process_info.tissue_diff_applied == 1)
    process_info.freqregion = handles.freqregion;
    process_info.source_list = handles.source_list;
    process_info.sources = handles.sources;
    process_info.abundances = handles.abudances;
end

if(process_info.img_loaded  == 1)
    process_info.back_img = handles.back_img;
end
process_info.cropped = handles.cropped;
if(process_info.cropped ==1)
    process_info.row = handles.row;
    process_info.col = handles.col;
    process_info.crop = handles.crop;
end


% row = handles.row1;
% col = handles.col1;
% row1 = handles.row1;
% col1 = handles.col1;
        
        


startingFolder = pwd;
% Put in the name of the mat file that the user wants to save.
% defaultFileName = fullfile(startingFolder, '*.mat')
[file,path] = uiputfile('signal.mat','Save file name');
if file == 0
    % User clicked the Cancel button.
    return;
end
fullFileName = fullfile(path, file);
% processed_info.
save(fullFileName,'signal','step','frequency','ndp','begin','signalnames','procres','classopt','simsig','position','process_info')
% [signal,step,frequency,ndp,begin,signalnames,procres,classopt,simsig,position]


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [fn p]=uigetfile({'*.*','MAT-files (*.mat)'; ...
    '*.jpg','JPEG images(*.jpg)'; ...
    '*.png','PNG images(*.png)';});

%     [fn p]=uigetfile({'*.*','MAT-files (*.mat)','*.png','*.jpg'});
    if isequal(fn,0) | isequal(p,0)
        return
    end
    [~,~,ext] = fileparts(fn);
    filepath = strcat(p,fn);
    if(strcmp(ext,'.mat'))
        temp_img = load(filepath);
        handles.back_img = rgb2gray(temp_img.img);
    else
        temp_img = imread(filepath);
        handles.back_img = rgb2gray(temp_img);
    end
    
    if(handles.cropped == 1)
        col_i = ceil(handles.scrollbarValue/handles.row1);
        row_i = mod(handles.scrollbarValue,handles.row1);
    if(row_i == 0)
        row_i = handles.row1;
    end
        col_o = col_i + handles.crop(1);
        row_o = row_i + handles.crop(3);
        out_num = (col_o-1)*handles.row + row_o;
    else
        out_num = handles.scrollbarValue;
    end
    
    handles.img_loaded  = 1;
%     image_x = size(handles.back_img,1);
%     image_y = size(handles.back_img,2);
%     aa = round(linspace(0,image_x,handles.row+1));
%     aa1 = round(linspace(0,image_y,handles.col+1));
%     handles.image_block = zeros(size(handles.back_img));
%     image_block1 = zeros(size(handles.back_img));
%     count = 1;
%     for j = 1:handles.col
%         for i= 1:handles.row
%             if(handles.location_map(i,j) == 2)
% %                image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
%                image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;                
%             end
%             if(out_num == count)
%                 handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
%             end
%             count = count +1;
%         end
%     end
%     set(gcf,'CurrentAxes',handles.axes5);
%     rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
% %     size(rgbImage)
%     imshow(rgbImage,'parent',handles.axes5);
%     hold on
%     H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
%     h = imshow(H_out,'parent',handles.axes5);
%     hold off
%     alpha_data = 0.5*handles.image_block;
%     set(h, 'AlphaData', alpha_data);
%     guidata(hObject, handles);
    image_x = size(handles.back_img,1);
    image_y = size(handles.back_img,2);
    aa = round(linspace(0,image_x,handles.row+1));
    aa1 = round(linspace(0,image_y,handles.col+1));
    handles.image_block = zeros(size(handles.back_img));
    set(gcf,'CurrentAxes',handles.axes5)  
    image_block1 = zeros(size(handles.back_img));
    count = 1;
    for j = 1:handles.col
        for i= 1:handles.row
            if(handles.location_map(i,j) == 2)
               image_block1(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1; 
            end
            if(out_num == count)
                handles.image_block(aa(i)+1:aa(i+1),aa1(j)+1:aa1(j+1)) = 1;
            end
            count = count +1;
        end
    end
    rgbImage = cat(3, handles.back_img, handles.back_img, handles.back_img);   
    imshow(rgbImage,'parent',handles.axes5);
    hold on
    H_out = cat(3, handles.image_block, zeros(size(handles.back_img)),image_block1);
    h = imshow(H_out,'parent',handles.axes5);
    hold off
    alpha_data = 0.5*(handles.image_block + image_block1);
    set(h, 'AlphaData', alpha_data);
    guidata(hObject, handles);

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
