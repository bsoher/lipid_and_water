% function [specProtectRange] = subSpecProtectRange( b0map,signalmask,L,...
%     protectfactor,ppm,metaPos,metaLineWidth ) 
% Estimate protect range in the spectral domain from a given B0 inhomogeneity map
% Input:
%     b0map: field inhomogeneity map in ppm
%     singalmask: high resolution support mask Nx*Ny
%     L: the resolution of the low resolution CSI data, [Lx,Ly]
%     protectfactor: not used, 1 by default
%     ppm: 
%     metaPos: position of the metabolites in ppm
%     metaLineWidth: line width of the metabolites
%  output:
%     specProtectRange: [Lx,Ly,Nf]
% Modified from Diego's Code
% Chao Ma July 23 2012

function specProtectRange          = subSpecProtectRange(b0map, signalmask, L, ppm, metaPos, metaLineWidth) 

protectfactor                      = 1;
N                                  =size(b0map);
W                                  = protectfactor*N./L;

% estimate B0max and B0min for each spatial point of the low resolution csi
% data
pr                                 = zeros([L(1) L(2) 2]);
for m = 1:L(1)
  for n = 1:L(2)    
    % Find aliasing pattern
    Alr                            = zeros(L);
    Alr(m,n)                       =1;
    fAlr                           = fftshift(fft2(ifftshift(Alr)));
    fAhr                           = zeros(N);
    fAhr(N(1)/2+1-L(1)/2:N(1)/2+L(1)/2,N(2)/2+1-L(2)/2:N(2)/2+L(2)/2)...
                                   =fAlr;
    Ahr                            =fftshift(ifft2(ifftshift(fAhr)));
    
    % Find max point of aliasing pattern
    [temp,mptemp]                  = max(abs(Ahr));
    [~,mp(2)]                      = max(temp);
    mp(1)                          = mptemp(mp(2));
    
    % Find region of influence (around max point)
    lowx                           = mp(1)-W(1);
    lowy                           = mp(2)-W(2);
    highx                          = mp(1)+W(1);
    highy                          = mp(2)+W(2);
    
    AVOID_WRAPAROUND               = 1;
    if AVOID_WRAPAROUND
      lowx                         = max(lowx,1);
      lowy                         = max(lowy,1);
      highx                        = min(highx,N(1));
      highy                        = min(highy,N(2));
    end
    
    pointsx                        = floor(mod([lowx:highx]-1,N(1))+1);
    pointsy                        = floor(mod([lowy:highy]-1,N(2))+1);
    maskInfluence                  = zeros(N);

    maskInfluence(pointsx,pointsy) = 1;
    
    % Find composite mask
    mask                           = maskInfluence>0 & signalmask>0;
    summask(m,n)                   = sum(sum(mask));
    
    % Apply composite mask
    if(sum(sum(mask))>0)
      pr(m,n,1)                    = min(b0map(mask));
      pr(m,n,2)                    = max(b0map(mask));
    else
      pr(m,n,:)                    = zeros(2,1);
    end
    
  end
end

% estimate the protect r
T                                  = length(ppm);
specProtectRange                   = zeros([L(1) L(2) T]);
for m=1:L(1)
  for n=1:L(2)
    curpr                          = zeros(T,1);
    for p=1:length(metaPos)
      if summask(m,n) > 0
        [~,i1]                     = min(abs(ppm-(metaPos(p)+pr(m,n,2))));
        [~,i2]                     = min(abs(ppm-(metaPos(p)+pr(m,n,1))));
        curpr(i1)                  = 1;
        curpr(i2)                  = 1;
        curpr(ppm<=(metaPos(p)+metaLineWidth(p)/2+pr(m,n,2)) ...
            & ppm>=(metaPos(p)-metaLineWidth(p)/2+pr(m,n,1))) = 1;
      end
    end
    specProtectRange(m,n,:)        = curpr;
  end
end