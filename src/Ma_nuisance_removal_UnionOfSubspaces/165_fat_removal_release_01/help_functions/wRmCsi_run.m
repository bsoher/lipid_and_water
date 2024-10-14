%% water removal for CSI data

function [xtCsiAllW, xtCsiAllWr] = wRmCsi_run(xtCsiAll, params)
    disp('Water removal for the CSI data ...');
    % default parameters
    selParams                    = [];
    selParams.dt                 = params.dt;
    selParams.NtEcho             = params.NtEcho;
    selParams.signalType         = 'water';

    [optWat,optLip,optMeta]      = setDefaultSpectInfo;
    params                       = setDefaultField(params,'optWat',optWat);
    params                       = setDefaultField(params,'optLip',optLip);
    params                       = setDefaultField(params,'optMeta',optMeta);
    params                       = setDefaultField(params,'optOther',[]);
    
    % water mask  
    Ny                           = size(xtCsiAll,1);
    Nx                           = size(xtCsiAll,2);
    Nt                           = size(xtCsiAll,3);
    Nc                           = size(xtCsiAll,4);
    ImMask                       = logical(imresize(logical(params.waterMask),[Ny,Nx])); 

    % water removal
    xtCsiAllW                    = zeros([Ny,Nx,Nt,Nc]);
    xtCsiAllWr                   = zeros([Ny,Nx,Nt,Nc]);
    tic;
    for indc = 1:Nc
        [xtCsiAllW(:,:,:,indc), xtCsiAllWr(:,:,:,indc)] ...
                                 = nsRm(xtCsiAll(:,:,:,indc), ImMask, params.B0Map, selParams, ...
                                   params.optWat, params.optLip, params.optMeta,params.optOther);
    end
    temp                         = toc;
    disp(['Water removal for Csi data finished in ',num2str(temp), ' seconds.']);
end