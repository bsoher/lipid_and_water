% ---------------------------------------------------------------------
%  GS2D_F:  Fast algorithm for 
%           two-dimensional generalized series reconstruction
%
%  Accepts two arrays, *dyn and *ref, corresponding to dynamic 
%  and constraint data (k-space data; recontruction along the readout 
%  direction is alreday performed), respectively.  The generalized series 
%  reconstruction of the dynamic data using the given constraint is
%  returned in the array *r. 
%
%    dorig(1), dorig(2) - indices of the k=0 signal in the dynamic data set 
%    rorig(1), rorig(2) - indices of the k=0 signal in the reference data set
%	  reg - regularization parameter (suggestion: in the range 0.001 to 0.1)
%	  filter: 0 (no filter), 1 (filter activated)
%    pmode: normally set to 0.
% 
%  Implementation notes:
%
%  (C) Zhi-Pei Liang
%  Last revised 11/15/2001 (Zhi-Pei Liang)/Stanford University (sabbatical)
% -------------------------------------------------------------------------

function [r,output] = gs2d_f(dyn, ref, dorig, rorig, reg, filter,pmode)

if(nargin < 5 | reg < 0 | reg > 1) reg = 0; end
if(nargin < 6) filter = 1; end
if(nargin < 7) pmode = 0; end

ist1  = rorig(1)-dorig(1)+1;
iend1 = size(dyn,1)+rorig(1)-dorig(1);
ist2  = rorig(2)-dorig(2)+1;
iend2 = size(dyn,2)+rorig(2)-dorig(2);

H_ref_ima = Pishft2(ifft2(Pishft2(ref)));
if (pmode == 0) ref = Pishft2(fft2(Pishft2(abs(H_ref_ima)))); end

[Nx,Ny] = size(dyn);
if (filter == 1)
   h = ones([Nx,1]);
   for nn=1:Nx
       h(nn)=0.7+0.3*cos(2*pi*(nn-dorig(1))/Nx); 
   end
   for nn=1:Ny
       dyn(:,nn) = dyn(:,nn).*h; 
       ref(ist1:iend1,ist2+nn-1) = ref(ist1:iend1,ist2+nn-1).*h;;
   end
   clear h;
   h = ones([1,Ny]);
   for nn=1:Ny
       h(nn)=0.8+0.2*cos(2*pi*(nn-dorig(2))/Ny); 
   end
   for nn=1:Nx
       dyn(nn,:) = dyn(nn,:).*h; 
       ref(ist1+nn-1,ist2:iend2) = ref(ist1+nn-1,ist2:iend2).*h;;
   end
end

L_ref_ima = zeros(size(ref));
L_ref_ima(ist1:iend1,ist2:iend2) = ref(ist1:iend1,ist2:iend2);
L_ref_ima = Pishft2(ifft2(Pishft2(L_ref_ima)));

L_dyn_ima = zeros(size(ref));
L_dyn_ima(ist1:iend1,ist2:iend2) = dyn;
L_dyn_ima = Pishft2(ifft2(Pishft2(L_dyn_ima)));
  
mean_value = mean(abs(L_ref_ima(:)));
L_ref_ima  = (abs(L_ref_ima) + reg*mean_value).*exp(i*angle(L_ref_ima));

r = zeros(size(ref));
if (mean_value > eps) r = H_ref_ima.*abs(L_dyn_ima)./abs(L_ref_ima); end
output = [];
output.coefx = abs(L_dyn_ima)./abs(L_ref_ima);
return;

