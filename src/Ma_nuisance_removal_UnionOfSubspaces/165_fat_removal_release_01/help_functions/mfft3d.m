function k = mfft3d(x,norm_on)
% perform 3d fft with correct normalization and shift
% Chao Ma
% 

Ny = size(x,1);
Nx = size(x,2);
Nz = size(x,3);

% by default, we do normalization
if nargin<2
	norm_on = 1;
end

if Nz > 1
	% fft with correct shift
	k = fftshift(fft(ifftshift(...
	    fftshift(fft(ifftshift(...
	    fftshift(fft(ifftshift(...
	    x,...
	    1),[],1),1),...
	    2),[],2),2),...
	    3),[],3),3);

else
	% fft with correct shift
	k = fftshift(fft(ifftshift(...
	    fftshift(fft(ifftshift(...
	    x,...
	    1),[],1),1),...
	    2),[],2),2);
end

% fft with normalization
if norm_on == 1
	k = k./sqrt(Ny)./sqrt(Nx)./sqrt(Nz);
end