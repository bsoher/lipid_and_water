%% water and fat removal from undersampled data
% Created by Chao Ma
% Last modified: 03-29-2013
%

% water and fat removal from undersampled data
function [usktWrLr,output,params]       = wfRmEpsi(uskt, params, echoMode, xtCsi, xtCsiW, xtCsiWr)

    % default parameters
    if ~exist('xtCsiW','var') || isempty(xtCsiW)
        xtCsiW                          = xtCsi;
    end
    if ~exist('xtCsiWr','var') || isempty(xtCsiWr)
        xtCsiWrL                        = xtCsi;
        xtCsiWrLr                       = xtCsi;
    else
        xtCsiWrL                        = xtCsiWr;
        xtCsiWrLr                       = xtCsiWr;
    end
    Nc                                  = size(uskt,4);
     
    % water and fat removal using both even and odd echoes 
    if echoMode                         == 1
        NyEpsi                          = size(uskt,1);
        NxEpsi                          = size(uskt,2);
        NEchoEpsi                       = size(uskt,3);

        % flip the second dimension
        uskt(:,:,2:2:end,:,:)           = flipdim(uskt(:,:,2:2:end,:,:),2);
        uskt                            = reshape(uskt,NyEpsi*NxEpsi*NEchoEpsi,Nc);
        usktWL                          = zeros(size(uskt));
        usktWrLr                        = zeros(size(uskt));
        usktW                           = zeros(size(uskt));
        usktL                           = zeros(size(uskt));

        sampleIndex                     = find(params.sampleMask==1);
        sampleMaskShift                 = ifftshift(ifftshift(params.sampleMask,1),2);
        sampleIndexShift                = find(sampleMaskShift==1);
        params.sampleIndex              = sampleIndex;
        params.sampleIndexShift         = sampleIndexShift;

        xtL                             = [];
        xtW                             = [];
        xSol                            = [];
        tind                            = find(uskt(:,1)~=0);
        for indc = 1:Nc
        	params.xtD1W                = squeeze(xtCsiW);
        	params.xtD1WrL              = squeeze(xtCsiWrL);
            params.xtD1WrLr             = squeeze(xtCsiWrLr);
        	[dWL,dWrLr,output]          = nsRmSparse(uskt(tind,indc), params);
        	if indc == 1
        		xtL 			        = output.xtL;
        		xtW 			        = output.xtW;
                xSol                    = output.xSol;
                VLD1                    = output.VLD1;
                VWD1                    = output.VWD1;
                VMD1                    = output.VMD1;
                params.RW               = size(VWD1,1);
                params.RL               = size(VLD1,1);
                params.RM               = size(VMD1,1);
                params.verbose          = 0;
        	else
        		xtL 			        = cat(4,xtL,output.xtL);
        		xtW 			        = cat(4,xtW,output.xtW);
                xSol                    = cat(3,xSol,output.xSol);
        	end
            usktWL(tind,indc)           = dWL;
            usktWrLr(tind,indc)         = dWrLr;
            usktW(tind,indc)            = output.dW;
            usktL(tind,indc)            = output.dL;
        end

        % reshape data
        uskt                            = reshape(uskt,[NyEpsi,NxEpsi,NEchoEpsi,Nc]);
        usktWrLr                        = reshape(usktWrLr,[NyEpsi,NxEpsi,NEchoEpsi,Nc]);
        usktWL                          = reshape(usktWL,[NyEpsi,NxEpsi,NEchoEpsi,Nc]);
        usktW                           = reshape(usktW,[NyEpsi,NxEpsi,NEchoEpsi,Nc]);
        usktL                           = reshape(usktL,[NyEpsi,NxEpsi,NEchoEpsi,Nc]);
        
        uskt(:,:,2:2:end,:)             = flipdim(uskt(:,:,2:2:end,:),2);
        usktWrLr(:,:,2:2:end,:)         = flipdim(usktWrLr(:,:,2:2:end,:),2);
        usktWL(:,:,2:2:end,:)           = flipdim(usktWL(:,:,2:2:end,:),2);
        usktW(:,:,2:2:end,:)            = flipdim(usktW(:,:,2:2:end,:),2);
        usktL(:,:,2:2:end,:)            = flipdim(usktL(:,:,2:2:end,:),2);

        output                          = [];
        output.usktWL                   = usktWL;
        output.usktW                    = usktW;
        output.usktL                    = usktL;
        output.xtL                      = xtL;
        output.xtW                      = xtW;
        output.xSol                     = xSol;
        output.VLD1                     = VLD1;
        output.VWD1                     = VWD1;
        output.VMD1                     = VMD1;
    end
end
