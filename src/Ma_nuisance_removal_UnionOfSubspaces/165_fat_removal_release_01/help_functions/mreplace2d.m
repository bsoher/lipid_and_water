function kReplaced = mreplace2d(kCenter, kOrig, Shifted)
% replace the center part of data kOrig by kCenter
% Chao Ma
%

if nargin<3
    Shifted                                          = 1;
end
temp                                                 = size(kCenter);
Ny                                                   = temp(1);
Nx                                                   = temp(2);
temp                                                 = size(kOrig);
Ny2                                                  = temp(1);
Nx2                                                  = temp(2);

kReplaced = kOrig;
if Shifted
    NyFloor                                          = Ny2/2 + 1 - Ny/2;
    NyCeil                                           = Ny2/2 + Ny/2;
    NxFloor                                          = Nx2/2 + 1 - Nx/2;
    NxCeil                                           = Nx2/2 + Nx/2;
    kReplaced(NyFloor:NyCeil,NxFloor:NxCeil,:,:)     = kCenter;
else
    % top left corner
    kReplaced(1:Ny/2,1:Nx/2,:,:)                     = kCenter(1:Ny/2,1:Nx/2,:,:);
    % top right corner
    kReplaced(1:Ny/2,(Nx2-Nx/2+1):Nx2,:,:)           = kCenter(1:Ny/2,Nx/2+1:Nx,:,:);
    % bottom left corner
    kReplaced((Ny2-Ny/2+1):Ny2,1:Nx/2,:,:)           = kCenter(Ny/2+1:Ny,1:Nx/2,:,:);
    % bottom right corner
    kReplaced((Ny2-Ny/2+1):Ny2,(Nx2-Nx/2+1):Nx2,:,:) = kCenter(Ny/2+1:Ny,Nx/2+1:Nx,:,:);
end