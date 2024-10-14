% function output = subSpatBasisSpecFreqBasisMat( x, transp_option,Ax,Phif )
%
%
% Parameters:
%   - x: [P, 1] complex amplitude vector
%        (P is the number of high-res voxels possibly containing fat)
%
%   - transp_option: 'notransp' or 'transp'
%   - Ax: [Nkx*Nky,P]  matrix, "spatial" basis for subspace of signals that may come from lipid region
%   - Phi: [P,Nf] spectral waveform
%
% Returns:
%   - output: [P x 1] vector
%
% Author: Chao Ma,
% Created: August 6, 2012
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = subSpatBasisSpecFreqBasisMat( x,Ax,Phif, transpOpt )
NLr = size(Ax,1);
P = size(Ax,2);
Nf = size(Phif,2);
if (nargin < 4) || (isempty(transpOpt ))
    transpOpt = 'notransp';
end

Nseg = floor(P/1000);       % number of segment
if strcmp('notransp',transpOpt)
    output = zeros(NLr,Nf);
    %     Ax(:,1:N)*diag(x)*phi(1:N,:);
    if Nseg > 0
        for ind = 1:Nseg
            indP = ((ind-1)*1000 + 1):(ind*1000);
            output = Ax(:,indP)*diag(x(indP))*Phif(indP,:) + output;
        end
    end
    indP = (Nseg*1000 + 1):P;
    output = Ax(:,indP)*diag(x(indP))*Phif(indP,:) + output;
else % Transpose
    output = zeros(P,1);
    %     a2 = diag(Ax(:,1:N)'*x*phi(1:N,:)');
    if Nseg>0
        for ind = 1:Nseg
            indP = ((ind-1)*1000 + 1):(ind*1000);
            output(indP) = diag(Ax(:,indP)'*x*Phif(indP,:)');
        end
    end
    indP = (Nseg*1000 + 1):P;
    output(indP) = diag(Ax(:,indP)'*x*Phif(indP,:)');
end
