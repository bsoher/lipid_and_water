function xtpocs    = mpocs(xt,mask,opt)
%% POCS (projection on to convex set) reconstruction
% Created by Chao Ma
% Last modified on 04-23-2013
%

% some default values
if ~exist('opt','var') || isempty(opt)
	MaxIter        = 5;
	tolFunc        = 1e-5;
    debug_on       = 0;
else
	MaxIter        = opt.MaxIter;
	tolFunc        = opt.tolFunc;
    debug_on       = opt.debug_on;
end

%% preparation
Ny                 = size(xt,1);
Nx                 = size(xt,2);
Nt                 = size(xt,3);
Nc                 = size(xt,4);
NyInterp           = size(mask,1);
NxInterp           = size(mask,2);

%% pocs iteration
kt                 = mfft2d(xt);
iter               = 0;
xtpocs             = zeros(NyInterp,NxInterp,Nt,Nc);
tmask              = repmat(mask,[1,1,Nt,Nc]);
while iter < MaxIter
    ktpocs         = mfft2d(xtpocs,1);
    
    % calculate difference
    ktpocs_center  = mtr2d(ktpocs,[Ny,Nx]);
    k_diff         = norm(ktpocs_center(:)-kt(:));
    if debug_on    == 1
        iter, k_diff,
    end
    if (k_diff<tolFunc)
        break;
    else
        iter       = iter + 1;
    end

    % replace the center of the kspace data
    ktpocsNew      = mreplace2d(kt,ktpocs);
    
    % transform back to the image domain
    xtpocs         = mifft2d(ktpocsNew,1);
    
    % apply mask
    xtpocs         = xtpocs.*tmask; 
end