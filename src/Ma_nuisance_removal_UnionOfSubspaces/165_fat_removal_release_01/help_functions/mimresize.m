function ImRs = mimresize( Im, baseSize, Nos )
% resize
%

if nargin < 3
	Nos 	= 1;
end
if length(Nos) == 1
	Nos = [Nos,Nos];
end
Ny 		 	= size(Im,1);
Nx 	     	= size(Im,2);
NyBase   	= baseSize(1);
NxBase   	= baseSize(2);

% set up coordinates based on the base size
x 		 	= (-Nx/2:(Nx/2-1))*(1/Nx) - (1/Nx)/2;
y 		 	= (-Ny/2:(Ny/2-1))*(1/Ny) - (1/Ny)/2;
[X,Y] 	 	= meshgrid(x,y);
xNew 	 	= -(1/2-(1/NxBase/2)) + (0:NxBase*Nos(2)-1)*1/NxBase/Nos(2);
yNew 	 	= -(1/2+(1/NyBase/2)) + (0:NyBase*Nos(1)-1)*1/NyBase/Nos(1);
[XNew,YNew] = meshgrid(xNew,yNew);

% interpolation
ImRs 		= interp2(X, Y, Im, XNew,YNew,'cubic',0);

end