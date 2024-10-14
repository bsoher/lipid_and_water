% ---------------------------------------------------------------------
%  GS2D:  Two-dimensional generalized series reconstruction
%
%  Accepts two arrays, *dyn and *ref, corresponding to dynamic 
%  and reference data, respectively.  The generalized series 
%  reconstruction of the dynamic data using the given constraint is
%  returned in the array *rec.
%
%    dorig - index of the k=0 signal in the dynamic data set 
%    rorig - index of the k=0 signal in the reference data set
%    reg   - regularization parameter (between 0 and 1) (optional)
% 
%  Implementation notes:
%
%  Created 11/18/2001; Zhi-Pei Liang/Stanford University (Sabbatical)
%
% ---------------------------------------------------------------------

function [r,output] = gs2d(dyn, ref, dorig, rorig, reg, filter, pmode)


if (ndims(dyn) ~= ndims(ref) | ndims(ref) ~= 2)
   error ('Invalid data sets');
end
if(nargin < 5 | reg < 0 | reg > 1) reg = 0; end
if(nargin < 6) filter = 0; end
if(nargin < 7) pmode = 0; end
 
if pmode == 1
    ref_ima = Pishft2(ifft2(Pishft2(ref)));
    ref     = Pishft2(fft2(Pishft2(ref_ima)));
else
    ref_ima = Pishft2(ifft2(Pishft2(ref)));
    ref     = Pishft2(fft2(Pishft2(abs(ref_ima))));
end

%
% Create a block-Toeplitz coefficient Matrix:
%
[Mx My] = size(ref);
[Nx Ny] = size(dyn);
npt     = Nx*Ny;
H       = zeros([npt npt]);
ref(rorig(1),rorig(2)) = ref(rorig(1),rorig(2))*(1+reg);

for n1=1:npt
  for n2=1:npt
      m1 = mod(n1-1,Nx)-mod(n2-1,Nx)+rorig(1);
      m2 = floor((n1-1)/Nx)-floor((n2-1)/Nx)+rorig(2);
      if (m1>0 & m1 <= Mx & m2> 0 & m2 <= My) 
         H(n1,n2) = ref(m1,m2);
      end
  end
end

% coef = pinv(H)*dyn(:);
coef = H\dyn(:);
coef = reshape(coef,[Nx,Ny]);
output = [];
output.coef = coef;

if (filter == 1)
   h = ones([Nx,1]);
   for nn=1:Nx
       h(nn)=0.7+0.3*cos(2*pi*(nn-dorig(1))/Nx); 
   end
   for nn=1:Ny
       coef(:,nn) = coef(:,nn).*h; 
   end
   clear h;
   h = ones([1,Ny]);
   for nn=1:Ny
       h(nn)=0.7+0.3*cos(2*pi*(nn-dorig(2))/Ny); 
   end
   for nn=1:Nx
       coef(nn,:) = coef(nn,:).*h; 
   end
end

r  = zeros(size(ref));
ist1 = rorig(1)-dorig(1)+1; iend1 = ist1 + Nx -1;
ist2 = rorig(2)-dorig(2)+1; iend2 = ist2 + Ny -1;
r(ist1:iend1,ist2:iend2) = coef;

r = Pishft2(ifft2(Pishft2(r)))*Mx*My;
output.coefx = r;

if (pmode == 0)
    r = abs(ref_ima).*r;
elseif pmode == 1
    r = ref_ima.*r;
else
    ref = Pishft2(fft2(Pishft2(ref_ima.*abs(r))));
    ref(ist1:iend1,ist2:iend2) = dyn;
    r = Pishft2(ifft2(Pishft2(ref)));
end

return;



