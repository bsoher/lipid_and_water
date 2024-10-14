% fRmCsiEpsi.m
% Remove fat signal from CSI data with aid of EPSI datas
% Created by Chao Ma
%
function [xtCsiWrLr, xtCsiWrM, xtCsiWrL, xtLGs] = fRmCsiEpsi(xtCsiWr, xtL, VtL, algParams, ...
                                                algParamsGs, algParamsFrCsi, NtEcho)
    % initialization
    if ~exist('NtEcho','var')
        NtEcho                   = 0;
    end
    
    algParams                    = prepInit(algParams);

    if NtEcho == 0 % single side signal
        Nt                       = size(xtCsiWr,3);
        [xtCsiWrLr, xtCsiWrM, xtCsiWrL, xtLGs] ...
                                 = main_fRmCsiEpsi(xtCsiWr, xtL(:,:,1:Nt,:), VtL(:,1:Nt), algParams, ...
                                   algParamsGs, algParamsFrCsi);
    else % echo signal
        xtCsiWr_R                = xtCsiWr(:,:,NtEcho+1:end,:);
        xtCsiWr_L                = xtCsiWr(:,:,1:NtEcho,:);
        Nt                       = size(xtCsiWr,3);
        algParamsFrCsi.csiInfo.NtEcho = 0;

        % process right hand signal
        [xtCsiWrLr_R, xtCsiWrM_R, xtCsiWrL_R, xtLGs_R] ...
                                 = main_fRmCsiEpsi(xtCsiWr_R, xtL(:,:,1:Nt-NtEcho,:), ...
                                   VtL(:,1:Nt-NtEcho), algParams, ...
                                   algParamsGs, algParamsFrCsi);

        % process left hand signal
        txtCsiWrL                = conj(xtCsiWr_L(:,:,end:-1:1,:));
        algParamsFrCsi.ND1(3)    = NtEcho;
        algParamsFrCsi.ND2(3)    = NtEcho;
        algParamsFrCsi.tD1       = algParamsFrCsi.tD1(1:NtEcho);
        algParamsFrCsi.sampleMask= algParamsFrCsi.sampleMask(:,:,1:NtEcho);
        algParamsGs.NtTr         = NtEcho;
        algParamsGs.t            = algParamsGs.t(1:NtEcho);
        [txtCsiWrLr_L, txtCsiWrM_L, txtCsiWrL_L, txtLGs_L] ...
                                 = main_fRmCsiEpsi(txtCsiWrL, xtL(:,:,1:NtEcho,:), ...
                                   VtL(:,1:NtEcho), algParams, ...
                                   algParamsGs, algParamsFrCsi);
        xtCsiWrLr_L              = conj(txtCsiWrLr_L(:,:,end:-1:1,:));
        xtCsiWrM_L               = conj(txtCsiWrM_L(:,:,end:-1:1,:));
        xtCsiWrL_L               = conj(txtCsiWrL_L(:,:,end:-1:1,:));
        xtLGs_L                  = conj(txtLGs_L(:,:,end:-1:1,:));

        % assemble data
        xtCsiWrLr                = cat(3,xtCsiWrLr_L,xtCsiWrLr_R);
        xtCsiWrM                 = cat(3,xtCsiWrM_L,xtCsiWrM_R);
        xtCsiWrL                 = cat(3,xtCsiWrL_L,xtCsiWrL_R);
        xtLGs                    = cat(3,xtLGs_L,xtLGs_R);
    end
end

function algParams = prepInit(algParams)
    if ~isfield(algParams,'maxIter') || isempty(algParams.maxIter)
        algParams.maxIter        = 10;
    end

    if ~isfield(algParams,'tol') || isempty(algParams.tol)
        algParams.tol            = 1e-5;
    end
end

function [xtCsiWrLr, xtCsiWrM, xtCsiWrL, xtLGs] = main_fRmCsiEpsi(xtCsiWr, xtL, VtL, algParams, ...
                                                  algParamsGs, algParamsFrCsi)
    Nc                           = size(xtCsiWr,4);

    % initial guess of xtLGsDs
    temp                         = algParamsGs.reg;
    algParamsGs.reg              = 1e-1;
    txtCsiWrL                    = xtCsiWr(:,:,:,1:Nc);
    algParamsGs.VtL              = VtL(1:4,:);
    [~,xtLGsDs]                  = gsRecon(xtL(:,:,:,1:Nc),txtCsiWrL,algParamsGs);
    
    % start iteration
    algParamsGs.reg              = temp;
    algParamsGs.VtL              = [];
    iter                         = 0;
    
    while (iter<algParams.maxIter)    
        txtCsiWrLr               = xtCsiWr(:,:,:,1:Nc) - xtLGsDs(:,:,:,1:Nc);
        
        % fat removal
        [~,~,xtCsiWrM,outputCsi] ...
                                 = fRmCsi(txtCsiWrLr(:,:,:,1:Nc), algParamsFrCsi, xtCsiWr, []);
        algParamsFrCsi.VMD1      = outputCsi{1}.VMD1;
        algParamsFrCsi.UMD1      = outputCsi{1}.UMD1;
        algParamsFrCsi.SMD1      = outputCsi{1}.SMD1;
        
        % GS reconstruction to use both CSI and EPSI data
        txtCsiWrL                = xtCsiWr(:,:,:,1:Nc) - xtCsiWrM(:,:,:,1:Nc);
        [~,xtLGsDs_new]          = gsRecon(xtL(:,:,:,1:Nc),txtCsiWrL,algParamsGs);

        % calculate difference
        if (norm(xtLGsDs_new(:)-xtLGsDs(:)) > algParams.tol) && (iter < (algParams.maxIter-1))
            xtLGsDs              = xtLGsDs_new;
            iter                 = iter + 1;
        else 
            algParamsGs.VtL      = VtL;
            [xtLGs,xtLGsDs]      = gsRecon(xtL(:,:,:,1:Nc),txtCsiWrL,algParamsGs);
            xtCsiWrLr            = xtCsiWr(:,:,:,1:Nc) - xtLGsDs(:,:,:,1:Nc);
            break;
        end
    end

    xtCsiWrL                     = xtLGsDs;
end