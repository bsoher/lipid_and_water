function [xtGs, xtGsDs, coef, coefx] = gsRecon(xtEpsi, xtCsi, params)
% Perform general series reconstruction to compensate mismatches between EPSI and CSI data
% Created by Chao Ma
% Last modified 04-21-2014
%

    % data size
    ktEpsi                         = mfft2d(xtEpsi);
    NyEpsi                         = size(ktEpsi,1);
    NxEpsi                         = size(ktEpsi,2);
    ktCsi                          = mfft2d(xtCsi);
    NyCsi                          = size(ktCsi,1);
    NxCsi                          = size(ktCsi,2);
    Nt                             = size(ktCsi,3);
    Nc                             = size(ktCsi,4);

    % some default values
    if ~isfield(params,'reg')  || isempty(params.reg)
        params.reg                 = 1e-5;
    end
    if ~isfield(params,'NtTr') || isempty(params.NtTr) % if not given, estimate params.NtTr
        [~,ind]                    = max(reshape(xtCsi(:,:,1,1),[],1));
        [indy,indx]                = ind2sub([NyCsi,NxCsi],ind);
        temp                       = xtCsi(indy,indx,:,1);
        figure;
        plot(abs(temp(:)));
        params.NtTr                = input('NtTr = ');
    end
    if ~isfield(params,'dorig') || isempty(params.dorig)
        params.dorig               = [NyCsi/2+1,NxCsi/2];
    end
    if ~isfield(params,'rorig') || isempty(params.rorig)
        params.rorig               = [NyEpsi/2+1,NxEpsi/2+1];
    end
    if ~isfield(params,'pmode')
        params.pmode               = 0;
    end
    if ~isfield(params,'gs') || isempty(params.gs)
        gs                         = 1;
    end

    % roughly remove metabolite signals through projection
    if isfield(params,'VtM') && ~isempty(params.VtM)
        % projection matrix
        A                          = params.VtM.';
        Ainv                       = pinv(A);
        % B0 correction using CP
        tB0Map                     = imresize(params.B0Map,round([NyCsi,NxCsi]*params.NosCsi));
        temp                       = reshape(kron(tB0Map(:),params.t),[],Nt);
        Psi                        = exp(1i*2*pi*temp);
        Psi                        = reshape(Psi,[NyCsi*params.NosCsi,NxCsi*params.NosCsi,Nt]);
        for indc                   = 1:Nc
            xtCsiBc(:,:,:,indc)    = mcprec2d(xtCsi(:,:,:,indc),Psi,params.NosCsi,0);
        end
        % projection onto VtL and VtM
        xtCsiProjBc                = zeros([NyCsi*params.NosCsi,NxCsi*params.NosCsi,Nt]);
        for indc                   = 1:Nc
            for indy               = 1:NyCsi*params.NosCsi
                for indx           = 1:NxCsi*params.NosCsi
                    txt            = xtCsiBc(indy,indx,:,indc);
                    c              = Ainv*reshape(txt(:),[],1);
                    xtCsiProjBc(indy,indx,:,indc) = txt(:) - A*c;
                end
            end
        end
        % add field inhomogenetiy back
        for indc                   = 1:Nc
            xtCsiProj(:,:,:,indc)  = mcprec2d(xtCsiProjBc(:,:,:,indc),conj(Psi),1,1);
        end
        ktCsiProj                  = mtr2d(mfft2d(xtCsiProj),[NyCsi,NxCsi]);
    else
        ktCsiProj                  = ktCsi;
    end

    % GS reconstruction in the time domain
    ktGs                           = zeros(size(ktEpsi));
    coef                           = zeros([NyCsi,NxCsi,params.NtTr,Nc]);
    coefx                          = zeros([NyEpsi,NxEpsi,params.NtTr,Nc]);
    for indc = 1:Nc
        for indt = 1:params.NtTr
            kdyn                   = squeeze(ktCsiProj(:,:,indt,indc));
            kref                   = squeeze(ktEpsi(:,:,indt,indc));
            [r,temp]               = mgs2d(kdyn, kref, params.dorig, params.rorig, params.reg, 0, params.pmode);
            kr                     = Pishft2(fft2(Pishft2(r)));
            ktGs(:,:,indt,indc)    = kr;
            coef(:,:,indt,indc)    = temp.coef;
            coefx(:,:,indt,indc)   = temp.coefx;
        end
    end
    xtGs                           = mifft2d(ktGs);

    % projection onto VtL, if VtL exist
    if isfield(params,'VtL') && ~isempty(params.VtL)
        % projection matrix
        A                          = params.VtL.';
        Ainv                       = pinv(A);
        % B0 correction using CP
        tB0Map                     = imresize(params.B0Map,[NyEpsi,NxEpsi]);
        temp                       = reshape(kron(tB0Map(:),params.t),[],Nt);
        Psi                        = exp(1i*2*pi*temp);
        for indc                   = 1:Nc
            xtGsBc(:,:,:,indc)     = reshape(reshape(xtGs(:,:,:,indc),NyEpsi*NxEpsi,Nt).*conj(Psi),[NyEpsi,NxEpsi,Nt]);
        end
        % projection onto VtL
        xtGsProjBc                 = zeros([NyEpsi,NxEpsi,Nt]);
        for indc                   = 1:Nc
            for indy               = 1:NyEpsi
                for indx           = 1:NxEpsi
                    txt            = xtGsBc(indy,indx,:,indc);
                    c              = Ainv*reshape(txt(:),[],1);
                    xtGsProjBc(indy,indx,:,indc) = A*c;
                end
            end
        end
        % add field inhomogenetiy back
        for indc                   = 1:Nc
            xtGsProj(:,:,:,indc)   = reshape(reshape(xtGsProjBc(:,:,:,indc),NyEpsi*NxEpsi,Nt).*Psi,[NyEpsi,NxEpsi,Nt]);
        end
    else
        xtGsProj                   = xtGs;
    end
    xtGs                           = xtGsProj;

    % truncation
    ktGsProj                       = mfft2d(xtGsProj);
    ktGsTr                         = ktGsProj((NyEpsi/2+1-NyCsi/2):(NyEpsi/2+NyCsi/2), ...
                                              (NxEpsi/2+1-NxCsi/2+1):(NxEpsi/2+NxCsi/2+1),:,:);
    xtGsDs                         = mifft2d(ktGsTr);
end

