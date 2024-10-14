% function output = subASpatialBasisPR( a, transp_option,A,L,signalmask3D,lambda )
%
%
% Parameters:
%   - a: [P Nf x 1] complex amplitude vector
%        (P is the number of high-res voxels possibly containing fat)
%        (Nf is number of points of the spectra)
%   - transp_option: 'notransp' or 'transp'
%   - A: "spatial" basis for subspace of signals that may come from lipid region
%   - L: resolution of the output low-res dataset
%   - signalmask3D: [Lx x Ly x P] indicator function for the potential presence
%   of metabolites
%   - lambda: regularization parameter
%
% Returns:
%   - output: [LxLyP x 1] the low spatial frequency k-space data samples
%
% Author: Chao Ma, modified from Diego's code
% Created: July 15, 2012
% Author: Diego Hernando
% Created: October 22, 2006
% Last modified: November 17, 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output              = subSpatBasisSpecRoiMat(a, transp_option, A, L, signalmask3D, lambda, watermask3D)

if strcmp('notransp',transp_option)   
    Nf                       = size(signalmask3D,3);
    P                        = size(A,2);
    
    % x
    a_array                  = reshape(a,[P Nf]);
    % Vx
    xfcsilow                 = reshape(A*a_array,[L(1) L(2) Nf]);
    % Wvx
    xfcsilow2                = xfcsilow.*signalmask3D;
    
    % Penalize roughness in spectral domain, lambda*Df*rhou
    % $$$   for m=1:L(1)
    % $$$     for n=1:L(2)
    % $$$       diff_xfcsilow(m,n,:) = lambda*D*reshape(xfcsilow(m,n,:),[],1);
    % $$$     end
    % $$$   end
    diff_xfcsilow            = lambda*diff(xfcsilow.*watermask3D,1,3);
    
    % Reshape into vector, [WVx,lambda*Df*x]
    output                   = [reshape(xfcsilow2,[],1); reshape(diff_xfcsilow,[],1)];
else % Transpose 
    Nf                       = size(signalmask3D,3);
    
    % Reshape into low-res CSI data
    xfcsilow                 = reshape(a(1:L(1)*L(2)*Nf), [L(1), L(2), Nf]);    
    xfcsilow                 = xfcsilow.*signalmask3D;
    
    % 
    diff_xfcsilow            = reshape( a(L(1)*L(2)*Nf+1:end), [L(1) L(2) (Nf-1)]);   
    xfcsilow2                = zeros([L(1), L(2), Nf]);
    xfcsilow2(:,:,1)         = -diff_xfcsilow(:,:,1);
    xfcsilow2(:,:,2:(end-1)) = -diff(diff_xfcsilow,1,3);
    xfcsilow2(:,:,end)       = diff_xfcsilow(:,:,end); 
    xfcsilow2                = lambda*(watermask3D.*xfcsilow2);
    % $$$   for m=1:L(1)
    % $$$     for n=1:L(2)
    % $$$       xfcsilow2(m,n,:) = lambda*D'*reshape(diff_xfcsilow(m,n,:),[],1);
    % $$$     end
    % $$$   end
    
    % 
    xfcsilow                 = xfcsilow+xfcsilow2;
    a2                       = A'*reshape(xfcsilow,L(1)*L(2),Nf);
    
    % Reshape into vector
    output                   = reshape(a2,[],1);
end