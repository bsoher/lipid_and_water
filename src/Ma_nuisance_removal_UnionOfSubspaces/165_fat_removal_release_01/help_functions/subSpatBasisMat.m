% subSpatialBasis( lx,ly,fatmask )
%
% Arguments:
%   - lx: number of pixels in the x dimension for the low-res dataset
%   - ly: number of pixels in the y dimension for the low-res dataset
%   - fatmask: fat mask of size HX x HY (indicator function of fat
%   locations)
% 
% Returns:
%   - A: matrix contains the basis functions
%
% Author: Chao Ma, modified from Diego's code
% Date created: July 16, 2012
%
% Author: Diego Hernando 
% Date created: June 26, 2006
% Date last update: June 27, 2006
%--------------------------------------------------------------------------

function A              = subSpatBasisMat( lx,ly,fatmask )

hx                      = size(fatmask,1);
hy                      = size(fatmask,2);

% Let us create a basis for the signal subspace (brute-force for
% now)
indFat                  = find(fatmask==1);
A                       = zeros(lx*ly,length(indFat));
for ind                 = 1:length(indFat)
      temp              = zeros(size(fatmask));
      temp(indFat(ind)) = 1e8;
      ftemp             = fftshift(fft2(temp));
      ftemp2            = ftemp((hx/2-lx/2+1):(hx/2+lx/2),(hy/2-ly/2+1+1):(hy/2+ly/2+1));
      % normalization
      temp3             = reshape(ifft2(fftshift(ftemp2)),lx*ly,1);
      temp3             = temp3./norm(temp3);
      A(:,ind)          = temp3;
end