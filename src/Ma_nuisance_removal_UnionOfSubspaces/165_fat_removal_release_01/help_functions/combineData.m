function xt = combineData( xtAll )
% combine phased array data using SVD
% Created by Chao Ma
% Last modified 12-20-2013
%

[Ny,Nx,Nt,Nc] = size(xtAll);
xt 			  = zeros(Ny,Nx,Nt); 

% combine data	
for indy = 1:Ny
    for indx = 1:Nx
        [~,S,V] = svd(reshape(xtAll(indy,indx,:,:),[Nt,Nc]).'); % num_coils*num_echoes
        xt(indy,indx,:) = S(1,1)*V(:,1)'; % take conjugate here to flip the spectrum, if not taking conjugate, the peaks at right side have lower ppms
    end
end

disp('Coil successfully combined');

end
