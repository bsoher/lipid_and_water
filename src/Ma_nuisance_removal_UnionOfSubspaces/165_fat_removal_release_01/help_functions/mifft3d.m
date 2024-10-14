function x = mifft3d(k,norm_on)
	% perform inverse 3d fft with correct normalization and shift
	% Chao Ma
	%

Ny = size(k,1);
Nx = size(k,2);
Nz = size(k,3);

% by default, we do normalization
if nargin<2
	norm_on = 1;
end

if Nz > 1
	% ifft with correct shift
	x = ifftshift(ifft(fftshift(...
	    ifftshift(ifft(fftshift(...
	    ifftshift(ifft(fftshift(...
	    k,...
	    1),[],1),1),...
	    2),[],2),2),...
	    3),[],3),3); 
else
	% ifft with correct shift
	x = ifftshift(ifft(fftshift(...
	    ifftshift(ifft(fftshift(...
	    k,...
	    1),[],1),1),...
	    2),[],2),2); 
end

% ifft with normalization
if norm_on == 1
	x = x.*sqrt(Ny).*sqrt(Nx).*sqrt(Nz);
end