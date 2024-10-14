% function output = subSpatBasisSpecBasisMat(x,Aspat,Aspec,transOpt)
%
%
% Parameters:
%   - x: [P, L] complex amplitude vector
%        (P is the number of high-res voxels possibly containing fat)
%
%   - transp_option: 'notransp' or 'transp'
%   - Aspat: [Nkx*Nky,P]  matrix, "spatial" basis for subspace of signals that may come from lipid region
%   - Aspec: [P, Nt, L] matrix, "spec" basis for subspace of signals
%
% Returns:
%   - output: 
%
% Author: Chao Ma,
% Created: August 6, 2012
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = subSpatBasisSpecBasisMat(x, Aspat, Aspec,transpOpt,W)
NLr = size(Aspat,1);
P = size(Aspat,2);
Nt = size(Aspec,2);
L = size(Aspec,3);

if (nargin < 4) || (isempty(transpOpt ))
    transpOpt = 'notransp';
end

if (nargin < 5) || (isempty(W))
    W = ones(NLr,Nt);
end

if strcmp('notransp',transpOpt)
    x = reshape(x,P,L);
    output = zeros(NLr,Nt);
    for indL = 1:L
        output = output + subSpatBasisSpecFreqBasisMat( x(:,indL),Aspat,Aspec(:,:,indL), transpOpt);
    end
    output = output.*W;
else
    output = zeros(P,L);
    x = x.*W;
    for indL = 1:L
        output(:,indL) = subSpatBasisSpecFreqBasisMat( x,Aspat,Aspec(:,:,indL), transpOpt);
    end
    output = output(:);
end