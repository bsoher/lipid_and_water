%% water/fat removal for EPSI data
function output        = wfRmEpsi_run(ktEpsi, params, xtCsi, xtCsiW, xtCsiWr)
    disp('Water removal for the EPSI data ...');

    % prepare other basic parameters
    % dimension and timing
    params.ND1         = [params.csiInfo.Ny,params.csiInfo.Nx,(params.csiInfo.Nt)];
    params.ND2         = round([params.epsiInfo.Ny*params.NosEpsi,params.epsiInfo.Nx*params.NosEpsi]);
    params.tD1         = reshape((0:params.csiInfo.Nt-1)*params.csiInfo.dt,1,[]) - ...
                         params.csiInfo.NtEcho*params.csiInfo.dt;
    params.tD2         = reshape(params.tEpsi(:),1,[]);
    params.dt          = params.csiInfo.dt;
    params             = nsRmSparseInitParams(params);

    % this where to make additional changes
    tmask              = logical(params.waterMaskD2);
    temp               = params.W*ones(size(params.waterMaskD2));
    temp(tmask)        = 1;
    params.waterMaskD2 = temp;

    tmask              = logical(params.lipMaskD2);
    temp               = params.W*ones(size(params.lipMaskD2));
    temp(tmask)        = 1;
    params.lipMaskD2   = temp;

    % spectral mask
    tfCsi              = linspace(-1/params.csiInfo.dt/2, 1/params.csiInfo.dt/2, params.csiInfo.Nt-1);
    tfCsi_ind          = find((tfCsi>-50)|(tfCsi<-330));
    params.specMask(:,tfCsi_ind) = 0;
    
    % water and fat removal
    uskt               = ktEpsi(:,:,:,1:params.Nc);
    [temp, output]     = wfRmEpsi(uskt,params,1, xtCsi, xtCsiW, xtCsiWr);
    output.usktWrLr    = temp;
end