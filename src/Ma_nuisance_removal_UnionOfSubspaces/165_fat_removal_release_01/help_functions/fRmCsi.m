%% Prepare in vivo CSI data
% Created by Chao Ma
% Last modified on 01-02-2013
%

function [xtCsiWrL,xtCsiWrLr,xtCsiWrM,outputArray] = fRmCsi(xtCsiWr, params, xtCsiWr2, sampleMask)
    Ny                                 = size(xtCsiWr,1);
    Nx                                 = size(xtCsiWr,2);
    Nt                                 = size(xtCsiWr,3);
    Nc                                 = size(xtCsiWr,4);
    Na                                 = size(xtCsiWr,5);

    if ~exist('sampleMask','var') || isempty(sampleMask)
        sampleMask                     = true([Ny,Nx,Nt]);
    end

    if ~exist('xtCsiWr2','var')
        xtCsiWr2                       = xtCsiWr;
    end
    
    % fat removal
    xtCsiWrL                           = zeros([Ny,Nx,Nt,Nc,Na]);
    xtCsiWrLr                          = zeros([Ny,Nx,Nt,Nc,Na]);
    xtCsiWrM                           = zeros([Ny,Nx,Nt,Nc,Na]);
    sampleIndex                        = find(params.sampleMask==1);
    sampleMaskShift                    = ifftshift(ifftshift(params.sampleMask,1),2);
    sampleIndexShift                   = find(sampleMaskShift==1);
    params.sampleIndex                 = sampleIndex;
    params.sampleIndexShift            = sampleIndexShift;
    for inda = 1:Na
        for indc = 1:Nc
            params.xtD1W               = [];
            if size(xtCsiWr2,4)        == Nc
    		  params.xtD1WrL           = squeeze(xtCsiWr2(:,:,:,indc,inda));
    		  params.xtD1WrLr          = squeeze(xtCsiWr2(:,:,:,indc,inda));
            else
              params.xtD1WrL           = squeeze(xtCsiWr2);
              params.xtD1WrLr          = squeeze(xtCsiWr2);
            end
            tdWr                       = mfft2d(squeeze(xtCsiWr(:,:,:,indc,inda)));
            [dWL, dWrLr, output]       = nsRmSparse( tdWr(sampleMask), params );
            temp                       = zeros(size(tdWr)); 
            temp(sampleMask)           = dWL;                
            xtCsiWrL(:,:,:,indc,inda)  = mifft2d(temp);
            temp                       = zeros(size(tdWr)); 
            temp(sampleMask)           = dWrLr;                
            xtCsiWrLr(:,:,:,indc,inda) = mifft2d(temp);   
            temp                       = zeros(size(tdWr)); 
            temp(sampleMask)           = output.dM;                
            xtCsiWrM(:,:,:,indc,inda)  = mifft2d(temp);    
            outputArray{indc,inda}     = output;
        end
    end
end
