%  ktcsilow = fftshift(fftshift(fft(fft(xtfatremoved_lr_PG_Mask2,[],1),[],2),1),2);

%  data = mrsiWithKnownLineshapesAndFieldMap( reshape(ktcsilow,[],1), 'transp',fimatrix, lineshapes,[64 64],[16 16],watermask3D );
%  [meta_estim,FLAG,RELRES,ITER] = pcg(@(a) mrsiWithKnownLineshapesAndFieldMapAndTranspose( a, fimatrix, lineshapes,[64 64],[16 16],watermask3D ),data,1e-6,10 );
%
%
%
%


function x = subSpatBasisSpecRoiProjMat(a, A, L, signalmask3D, lambda, watermask3D)

x1         = subSpatBasisSpecRoiMat(a, 'notransp', A, L, signalmask3D, lambda, watermask3D);
x          = subSpatBasisSpecRoiMat(x1, 'transp', A, L, signalmask3D, lambda, watermask3D);




