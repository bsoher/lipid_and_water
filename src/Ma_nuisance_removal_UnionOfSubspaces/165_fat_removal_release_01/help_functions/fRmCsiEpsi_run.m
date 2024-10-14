%% fat removal for CSI data
% 
function [xtCsiAllWrLr, xtCsiAllWrM, xtCsiAllWrL, xtLGs] = fRmCsiEpsi_run(xtCsiAllWr, VtL, xtL, params)
    % set up algParamsFrCsi
    algParamsFrCsi.NosD1         = params.NosCsi;
    
    if params.epsiInfo.NechoBeforeEcho == 0 && params.csiInfo.NtEcho ~= 0
        algParamsFrCsi.ND1       = round([params.csiInfo.Ny,params.csiInfo.Nx,...
                                      params.csiInfo.Nt-params.csiInfo.NtEcho]);
        algParamsFrCsi.ND2       = round([params.csiInfo.Ny*params.NosCsi,params.csiInfo.Nx*params.NosCsi,...
                                      params.csiInfo.Nt-params.csiInfo.NtEcho]);
        algParamsFrCsi.tD1       = (0:params.csiInfo.Nt-1-params.csiInfo.NtEcho)*params.csiInfo.dt;
    else
        algParamsFrCsi.ND1       = round([params.csiInfo.Ny,params.csiInfo.Nx,params.csiInfo.Nt]);
        algParamsFrCsi.ND2       = round([params.csiInfo.Ny*params.NosCsi,params.csiInfo.Nx*params.NosCsi,...
                                      params.csiInfo.Nt]);
        algParamsFrCsi.tD1       = (0:params.csiInfo.Nt-1)*params.csiInfo.dt - ...
                                   params.csiInfo.NtEcho*params.csiInfo.dt;
    end
    algParamsFrCsi.dt            = params.csiInfo.dt;
    algParamsFrCsi.ImMask        = params.ImMask;
    algParamsFrCsi.waterMask     = params.waterMask;
    algParamsFrCsi.RW            = 0;
    algParamsFrCsi.RL            = 0;
    algParamsFrCsi.RM            = params.RM;
    algParamsFrCsi.B0Map         = params.B0Map;
    if params.epsiInfo.NechoBeforeEcho == 0 && params.csiInfo.NtEcho ~= 0
        temp                     = ones([params.csiInfo.Ny,params.csiInfo.Nx+2,...
                                         params.csiInfo.Nt-params.csiInfo.NtEcho]);
    else
        temp                     = ones([params.csiInfo.Ny,params.csiInfo.Nx+2,params.csiInfo.Nt]);
    end
    sampleMaskCsi                = mzp2d(temp,[params.csiInfo.Ny,params.csiInfo.Nx]*params.NosCsi);
    tind                         = (params.csiInfo.Nx*params.NosCsi/2+1-params.csiInfo.Nx/2-1)...
                                   :(params.csiInfo.Nx*params.NosCsi/2+1-params.csiInfo.Nx/2);
    sampleMaskCsi(:,tind,:)      = 0;
    algParamsFrCsi.sampleMask    = sampleMaskCsi;
    algParamsFrCsi.csiInfo       = params.csiInfo;
    algParamsFrCsi               = setFieldSafe(algParamsFrCsi,'lambdaWat',0);
    algParamsFrCsi               = setFieldSafe(algParamsFrCsi,'lambdaLip',0);
    algParamsFrCsi               = setFieldSafe(algParamsFrCsi,'lambdaMet',1e-2);
    algParamsFrCsi               = setFieldSafe(algParamsFrCsi,'lambdaDfWat',0);
    algParamsFrCsi               = setFieldSafe(algParamsFrCsi,'lambdaDfLip',0);
    algParamsFrCsi               = setFieldSafe(algParamsFrCsi,'lambdaDfMeta',0);
    algParamsFrCsi               = nsRmSparseInitParams(algParamsFrCsi);
    algParamsFrCsi.optWat        = [];
    % this where to make additional changes
    tmask                        = logical(algParamsFrCsi.metaMaskD2);
    temp                         = params.W*ones(size(algParamsFrCsi.metaMaskD2));
    temp(tmask)                  = 1;
    algParamsFrCsi.metaMaskD2    = temp;
    algParamsFrCsi.verbose       = 0;

    % algParamsGs
    algParamsGs                  = [];
    algParamsGs.NtTr             = params.NtTr;
    algParamsGs.VtL              = [];
    algParamsGs.VtM              = [];
    algParamsGs.B0Map            = params.B0Map;
    if params.epsiInfo.NechoBeforeEcho == 0 && params.csiInfo.NtEcho ~= 0
        algParamsGs.t            = reshape((0:params.csiInfo.Nt-params.csiInfo.NtEcho-1)*params.csiInfo.dt,1,[]);
    else
        algParamsGs.t            = reshape((0:params.csiInfo.Nt-1)*params.csiInfo.dt...
                                   -params.csiInfo.NtEcho*params.csiInfo.dt,1,[]);
    end
    algParamsGs.reg              = params.GsReg;

    % fat removal
    algParams.maxIter            = params.maxIter;
    algParams.tol                = params.tol;
    if params.epsiInfo.NechoBeforeEcho == 0 && params.csiInfo.NtEcho ~= 0
        NtEcho                   = params.csiInfo.NtEcho;
    else
        NtEcho                   = 0;
    end
    [xtCsiAllWrLr, xtCsiAllWrM, xtCsiAllWrL, xtLGs] ...
                                 = fRmCsiEpsi(xtCsiAllWr(:,:,:,1:params.Nc), xtL(:,:,:,1:params.Nc), ...
                                   VtL, algParams, algParamsGs, algParamsFrCsi,NtEcho);
end
