%% Nuisance signa removal using a union-of-subspaces model
% Chao Ma, Dec. 8th, 2014
%
clear all;
close all;
clc;

%% Load data
% CSI data
load('./data/exp2_no_ovs_csi.mat','ktCsiAll');
xtCsiAll                         = mifft2d(ktCsiAll);
xtCsi                            = combineData(xtCsiAll);
csiInfo                          = [];
csiInfo.Ny                       = size(xtCsiAll,1);
csiInfo.Nx                       = size(xtCsiAll,2);
csiInfo.Nt                       = size(xtCsiAll,3);
csiInfo.dt                       = 1/2000;    % sampling bandwidth
csiInfo.NtEcho                   = 0;

% EPSI data
load('./data/exp2_no_ovs_epsi.mat','ktEpsiAll','sampleMaskEpsi','tEpsi');
epsiInfo                         = [];
epsiInfo.echospace               = 2.06e-3/2; % echo spacing
epsiInfo.Ny                      = size(ktEpsiAll,1);
epsiInfo.Nx                      = size(ktEpsiAll,2);
epsiInfo.NechoBeforeEcho         = 0;

% mask
load('./data/exp2_no_ovs_mask.mat','ImMask','fatMask','waterMask');

% B0 map
load('./data/exp2_no_ovs_B0.mat','B0Map');

%% Water and fat removal for EPSI data
algParamsWFrEpsi                 = struct('B0Map',B0Map,'sampleMask',sampleMaskEpsi,...
                                  'ImMask', ImMask,'waterMask', waterMask, 'noMeta', 0,'NosD1',2,...
                                  'NosEpsi', 1, 'tEpsi', tEpsi,'RW', 15, 'RL', 9, 'RM',3, ...
                                  'verbose', 0, 'W', 0.01, 'lambdaWat', 1e-4, 'lambdaMet', 1e-4, ...
                                  'lambdaLip',1e-4,'lambdaDfWat',0,'lambdaDfLip',0, ...
                                  'lambdaDfMet', 0, 'tolCG', 1e-3, 'MaxIterCG', 100, 'Nc', 4, ...
                                  'csiInfo', csiInfo, 'epsiInfo',epsiInfo,'debug_on',1); 
wfRmEpsi_out                     = wfRmEpsi_run(ktEpsiAll,algParamsWFrEpsi,xtCsi,[],[]);


%% water removal for CSI data
optWat                           = struct('fSel',[-180, 200],'fSel2',[-180, -120, 200],'maxT2',[10,1e6]);
optOther                         = struct('fSel',[-180, 200],'fSel2',[-180, 200],'maxT2',[1e6]);
% water removal
algParamsWrCsi                   = struct('B0Map',B0Map,'waterMask',waterMask,...
                                   'dt',csiInfo.dt,'NtEcho',csiInfo.NtEcho,'optWat',optWat,...
                                   'optOther',optOther);
[~, xtCsiAllWr]                  = wRmCsi_run(xtCsiAll, algParamsWrCsi);

%% fat removal for CSI data
NosCsi                           = epsiInfo.Ny/csiInfo.Nx; 
algParamsFrCsiEpsi               = struct('NosCsi', NosCsi, 'csiInfo', csiInfo, ...
                                   'epsiInfo', epsiInfo, 'ImMask', ImMask, ...
                                   'waterMask', waterMask, 'RM', 4, 'W', 0.05, ...
                                   'B0Map', B0Map, 'NtTr', 100, 'GsReg', 1e-4, 'maxIter', 10, ...
                                   'tol', 2e-6, 'Nc', 4);
xtCsiAllWrLr                     = fRmCsiEpsi_run(xtCsiAllWr, wfRmEpsi_out.VLD1, ...
                                    wfRmEpsi_out.xtL, algParamsFrCsiEpsi);
xtCsiWr                          = combineData(xtCsiAllWr);
xtCsiWrLr                        = combineData(xtCsiAllWrLr);

%% display EPSI
xfEpsi                           = fftshift(fft(combineData(mifft2d(ktEpsiAll)),[],3),3);
xfEpsiW                          = fftshift(fft(combineData(mifft2d(wfRmEpsi_out.usktW)),[],3),3);
xfEpsiL                          = fftshift(fft(combineData(mifft2d(wfRmEpsi_out.usktL)),[],3),3);
xfEpsiWrLr                       = fftshift(fft(combineData(mifft2d(wfRmEpsi_out.usktWrLr)),[],3),3);
scale                            = 20;
figure;
imagesc(sum(abs(xfEpsi),3)*scale,[0,1]);axis image off; colormap(gray);
figure;
imagesc(sum(abs(xfEpsiW),3)*scale,[0,1]);axis image off; colormap(gray);
figure;
imagesc(sum(abs(xfEpsiL),3)*scale,[0,1]);axis image off; colormap(gray);
figure;
imagesc(sum(abs(xfEpsiWrLr),3)*scale,[0,0.2]);axis image off; colormap(gray);

%% display CSI
BW                               = 1/csiInfo.dt;
hzpppm                           = 2.8936*42.5786;
ppm                              = linspace(-BW/2,BW/2,csiInfo.Nt)/hzpppm; 
ppm                              = 4.7 + ppm - 0.034;

NosCsi                           = 4;
tImMask                          = logical(imresize(logical(ImMask),...
                                      round([csiInfo.Ny*NosCsi,csiInfo.Nx*NosCsi])));
twaterMask                       = logical(imresize(logical(waterMask),...
                                      round([csiInfo.Ny*NosCsi,csiInfo.Nx*NosCsi])));
tImMask                          = bwmorph(tImMask,'dilate');
twaterMask                       = bwmorph(twaterMask,'erode');
tlipMask                         = tImMask - twaterMask;
water_edge                       = edge(bwmorph(twaterMask,'dilate'));

% B0 map
tB0Map                           = imresize(B0Map,round([csiInfo.Ny*NosCsi,csiInfo.Nx*NosCsi]));
tCsi                             = (0:csiInfo.Nt-1)*csiInfo.dt;
temp                             = reshape(kron(tB0Map(:),tCsi),[],csiInfo.Nt);
Psi                              = exp(1i*2*pi*temp);
Psi                              = reshape(Psi,round([csiInfo.Ny*NosCsi,csiInfo.Nx*NosCsi,csiInfo.Nt]));

xtCsiWr_Zp                       = mcprec2d(xtCsiWr,Psi,NosCsi,0);
xfCsiWr_Zp                       = fftshift(fft(xtCsiWr_Zp,[],3),3);

xtCsiWrLr_Zp                     = mcprec2d(xtCsiWrLr,Psi,NosCsi,0);
xfCsiWrLr_Zp                     = fftshift(fft(xtCsiWrLr_Zp,[],3),3);

naa_ind                          = find(ppm>1.9 & ppm<2.1);
scale                            = 0.5e3;
temp                             = [0.05,0.9];
figure;
imagesc(sum(abs(xfCsiWr_Zp(:,:,naa_ind)),3)*scale,temp);axis image off;colormap(gray);

figure;
imagesc(sum(abs(xfCsiWrLr_Zp(:,:,naa_ind)),3)*scale+water_edge,temp);axis image off;colormap(gray);

xpts                             = [27,42,26];
ypts                             = [15,28,32];
scale                            = 0.5*1e4;
x_range                          = [0.5,4.2];
y_range                          = [0,2];
tfname                           = 'spec';
for ind = 1:length(xpts)
    figure;
    plot(ppm,squeeze(abs(xfCsiWr_Zp(ypts(ind),xpts(ind),:)))*scale,'k','LineWidth',2);
    set(gca,'FontSize',22,'XDir','Reverse');
    xlabel('ppm');
    xlim(x_range);

    figure;
    plot(ppm,squeeze(abs(xfCsiWrLr_Zp(ypts(ind),xpts(ind),:)))*scale,'k','LineWidth',2);
    set(gca,'FontSize',22,'XDir','Reverse');
    xlabel('ppm');
    xlim(x_range);
    ylim(y_range);
end