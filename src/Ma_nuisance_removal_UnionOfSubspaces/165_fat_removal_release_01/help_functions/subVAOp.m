function [y,yk] = subVAOp(x, opSign, designPara)

	% default values
    if isfield(designPara,'PsiWat')
        PsiWat                = designPara.PsiWat;
    else
        PsiWat                = 1;
    end
    if isfield(designPara,'PsiLip')
        PsiLip                = designPara.PsiLip;
    else
        PsiLip                = 1;
    end
    if isfield(designPara,'PsiMet')
        PsiMet                = designPara.PsiMet;
    else
        PsiMet                = 1;
    end
    if isfield(designPara,'UxMet')
        UxMet                 = designPara.UxMet;
        if isfield(designPara,'metaMask')
            metaMask          = designPara.metaMask;
        else
            metaMask          = 1;
        end
    else
        UxMet                 = [];
    end 
    if isfield(designPara,'UxWat')
        UxWat                 = designPara.UxWat;
        if isfield(designPara,'waterMask')
            waterMask         = designPara.waterMask;
        else
            waterMask         = 1;
        end
    else
        UxWat                 = [];
    end
    if isfield(designPara,'UxLip')
        UxLip                 = designPara.UxLip;
        if isfield(designPara,'fatMask')
            lipMask           = designPara.fatMask;
        elseif isfield(designPara,'lipMask')
            lipMask           = designPara.lipMask;
        else
            lipMask           = 1;
        end
    else
        UxLip                 = [];
    end
	
	sampleIndex               = designPara.sampleIndex;
	NHr                       = designPara.NHr;
	NyHr                      = NHr(1);
	NxHr                      = NHr(2);
	RWat                      = size(UxWat,2);
	RLip                      = size(UxLip,2);
	RMet                      = size(UxMet,2);
    Nt                        = designPara.Nt;

	if strcmp(opSign,'notransp')
		y                     = AOp(x);
	elseif strcmp(opSign,'transp')
		y                     = AcOp(x);
	elseif strcmp(opSign,'recon')
		[y,yk]                = AOp(x);
	elseif strcmp(opSign,'reconWL') || strcmp(opSign,'reconW') ...
		   || strcmp(opSign,'reconL') || strcmp(opSign,'reconM') ...
		   || strcmp(opSign,'reconWM')
		[y,yk]                = AOpOut(x);
	else
		warning('opSign is either notransp or transp!')
	end

	% generate output
	function [yWLM,reconWLM]  = AOpOut(x)
		% obtain UsLip and UsMet from x
		x                     = reshape(x,RWat+RLip+RMet,Nt);
		if RWat               ~= 0
			VtWat             = x(1:RWat,:);
        else
        	UxWat             = 0; 
			VtWat             = 0;
        end
        if RLip               ~= 0
            VtLip             = x((RWat+1):(RWat+RLip),:);
        else
        	UxLip             = 0;
			VtLip             = 0;
        end
        if RMet ~= 0     
			VtMet             = x((RWat+RLip+1):end,:);
        else
        	UxMet             = 0;
			VtMet             = 0;
        end

        if strcmp(opSign,'reconW')
			if RWat           ~= 0
				temp          = reshape(PsiWat.*(UxWat*VtWat),[NyHr,NxHr,Nt]);
				reconWLM      = temp;
				% x to k
				tktHr 	      = fft(fft(temp,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
				% truncation
				yWLM 	      = tktHr(sampleIndex);
			else
				reconWLM      = 0;
				yWLM          = 0;
			end
		end
        if strcmp(opSign,'reconL')
			if RLip           ~= 0
				temp          = reshape(PsiLip.*(UxLip*VtLip),[NyHr,NxHr,Nt]);
				reconWLM      = temp;
				% x to k
				tktHr 	      = fft(fft(temp,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
				% truncation
				yWLM 	      = tktHr(sampleIndex);
		else
				reconWLM      = 0;
				yWLM          = 0;
		end
		end
		if strcmp(opSign,'reconM')
		  if RMet             ~= 0
		  	temp              = reshape(PsiMet.*(UxMet*VtMet),[NyHr,NxHr,Nt]);
		  	reconWLM          = temp;
		  	% x to k
		  	tktHr 	          = fft(fft(temp,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
		  	% truncation
		  	yWLM 	          = tktHr(sampleIndex);
		  else
		  	reconWLM          = 0;
		  	yWLM              = 0;
		  end
		end
        if strcmp(opSign,'reconWL')
			temp              = reshape(PsiWat.*(UxWat*VtWat)+PsiLip.*(UxLip*VtLip),...
				                [NyHr,NxHr,Nt]);
			reconWLM          = temp;
			% x to k
			tktHr 	          = fft(fft(temp,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
			% truncation
			yWLM 	          = tktHr(sampleIndex);
		end
		if strcmp(opSign,'reconWLM')
			temp              = reshape(PsiWat.*(UxWat*VtWat)+PsiLip.*(UxLip*VtLip)...
				                        +PsiMet.*(UxMet*VtMet), [NyHr,NxHr,Nt]);
			reconWLM          = temp;
			% x to k
			tktHr 	          = fft(fft(temp,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
			% truncation
			yWLM 	          = tktHr(sampleIndex);
		end
	end

	% forward operator for updating V with given U
	function y                = AOp(x)
		% obtain Vt from x
		x                     = reshape(x,RWat+RLip+RMet,Nt);
        ImtHr                 = 0;
		if RWat               ~= 0
			VtWat             = x(1:RWat,:);
			ImtHr             = ImtHr + PsiWat.*(UxWat*VtWat);
        end
        if RLip               ~= 0
            VtLip             = x((RWat+1):(RWat+RLip),:);
            ImtHr             = ImtHr + PsiLip.*(UxLip*VtLip);
        end
        if RMet               ~= 0     
			VtMet             = x((RWat+RLip+1):end,:);
            ImtHr             = ImtHr + PsiMet.*(UxMet*VtMet); 
        end

		ImtHr                 = reshape(ImtHr,[NyHr,NxHr,Nt]);
		% downsampling
		% x to k
		ktHr                  = fft(fft(ImtHr,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
		% truncation
		y                     = ktHr(sampleIndex);
	end % end of AOp

	% Adjoint of the forward operator for updating U with given V
	function y                = AcOp(x)
		% zeropadding
		ktHr                  = zeros(NyHr,NxHr,Nt);
		ktHr(sampleIndex)     = x;

		ImtHr                 = ifft(ifft(ktHr,[],1),[],2)*sqrt(NyHr)*sqrt(NxHr); 
		ImtHr                 = reshape(ImtHr,[NyHr*NxHr,Nt]);

		% estimate UsLip and UsMet
		y                     = [];
        if RWat               ~= 0		
			VtWat             = UxWat'*(conj(PsiWat).*ImtHr);
        	y                 = [y;VtWat(:)];
        end
        if RLip               ~= 0		
			VtLip             = UxLip'*(conj(PsiLip).*ImtHr);
        	y                 = [y;VtLip(:)];
        end
		if RMet               ~= 0		
			VtMet             = UxMet'*(conj(PsiMet).*ImtHr);
            y                 = [y;VtMet(:)];
        end
	end % end of AcOp

end % end of subVAOp 
