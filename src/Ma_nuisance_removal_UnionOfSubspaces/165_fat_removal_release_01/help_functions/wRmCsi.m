%% Prepare in vivo CSI data
% function [xtCsiW,xtCsiWr] = wRmCsi(xtCsi, waterMaskCsi, optWat, optOther)
%     Input:
%     xtCsi: CSI data
%     output:
%     xtCsiW: water signal
%     xtCsiWr: water removed CSI data
% Created by Chao Ma
% Last modified on 07-16-2014
% History
% 07-16-2014: New approach to handle echo signal
% 04-20-2014: Field inhomogeneity
% 03-29-2014: Write it as a fucntion

function [xtCsiW,xtCsiWr]       = wRmCsi(xtCsi, params)
    if ~isfield(params,'NtEcho')
        params.NtEcho           = 0;
    end
    
    % water removal coil by coil
    Ny                          = size(xtCsi,1);
    Nx                          = size(xtCsi,2);
    Nt                          = size(xtCsi,3);
    Nc                          = size(xtCsi,4);

    xtCsiW                      = zeros([Ny,Nx,Nt,Nc]);
    xtCsiWr                     = zeros([Ny,Nx,Nt,Nc]);
    for indc = 1:Nc
        [xtCsiW(:,:,:,indc,inda), xtCsiWr(:,:,:,indc,inda)] ...
                                = nsRmNew(xtCsi(:,:,:,indc), ImMask, B0Map, params, ...
                                  optWat, optLip, optMeta);
    end
end
