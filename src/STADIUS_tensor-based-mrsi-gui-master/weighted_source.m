function [avg_src,H1] = weighted_source(distr,nor_spectra,spectra,src_est_type,nnls_type,l1_wei,H_type,norn_fac )

% no_src = size(distr,2);
% init = norn_fac.*distr;
init = distr;

no_of_voxels = size(nor_spectra,2);
switch src_est_type
    case 1
        avg_src = nor_spectra*distr;
    case 2
        [~,idx] = max(distr,[],2);
        weights = zeros(size(distr));
        R = size(distr,2);
        for i = 1:R
            temp_pts = idx == i;
            weights(temp_pts,i) = distr(temp_pts,i);
        %     weights(:,i) = weights(:,i)./max(weights(:,i));
        end
        avg_src = nor_spectra*weights;
    case 3
        avg_src = (distr\nor_spectra')';
%         A_B = [distr,nor_spectra'];
%         [U,si,V] = svd(A_B);
%         temp1 = 
        
    case 4
        [~,idx] = max(distr,[],2);
        weights = zeros(size(distr));
        R = size(distr,2);
        for i = 1:R
            temp_pts = idx == i;
            weights(temp_pts,i) = distr(temp_pts,i);
        %     weights(:,i) = weights(:,i)./max(weights(:,i));
        end
        avg_src = (weights\nor_spectra')';      
    case 5 
        A_B = [distr,nor_spectra'];
        [~,~,V] = svd(A_B);
        n = size(distr,2);
%         d = size(nor_spectra,1);
        V12 = V(1:n,n+1:end);
        V22 = V(n+1:end,n+1:end);
        avg_src = (-V12*inv(V22))';
    case 6
        avg_src = distr;
end


norm_fac1 = sqrt(sum(avg_src.*conj(avg_src)));
norm_fac1 = repmat(norm_fac1,[size(avg_src,1),1]);
avg_src = avg_src./norm_fac1;


% H1 = distr';
if(H_type == 1)
    H1 = distr';
elseif(H_type == 2)    
    total_src = avg_src;
    for i = 1:no_of_voxels
        if(nnls_type == 1)
%             H1(:,i) = l1_ls_nonneg(total_src, spectra(:,i), l1_wei,[],true,[],[],init(i,:),[]);
            H1(:,i) = l1_ls_nonneg(total_src, spectra(:,i), l1_wei);
        else
            H1(:,i) = lsqnonneg(total_src, spectra(:,i));
        end
    end
elseif(H_type == 3)
    total_src = [real(avg_src);imag(avg_src)];
    for i = 1:no_of_voxels
        if(nnls_type == 1)
%             H1(:,i) = l1_ls_nonneg(total_src, [real(spectra(:,i));imag(spectra(:,i))], l1_wei,[],true,[],[],init(i,:)',[]);
            H1(:,i) = l1_ls_nonneg(total_src, [real(spectra(:,i));imag(spectra(:,i))], l1_wei);
        else
            H1(:,i) = lsqnonneg(total_src, [real(spectra(:,i));imag(spectra(:,i))]);
        end
    end
elseif(H_type == 4)
    total_src = abs(avg_src);
    for i = 1:no_of_voxels
        if(nnls_type == 1)
%             H1(:,i) = l1_ls_nonneg(total_src, abs(spectra(:,i)), l1_wei,[],true,[],[],init(i,:),[]);
            H1(:,i) = l1_ls_nonneg(total_src, abs(spectra(:,i)), l1_wei);
        else
            H1(:,i) = lsqnonneg(total_src, abs(spectra(:,i)));
        end
    end
end

% for i = 1:no_of_voxels
%     if(nnls_type == 1)
%         H1(:,i) = l1_ls_nonneg(total_src, [real(nor_spectra(:,i));imag(nor_spectra(:,i))], l1_wei);
%     else
%         H1(:,i) = lsqnonneg(total_src, [real(nor_spectra(:,i));imag(nor_spectra(:,i))]);
%     end
% end

% H1 = total_src\[real(nor_spectra);imag(nor_spectra)];

% plot_results_subplot(real(avg_src),H1',row-(crop(2)+ crop(1)),col-(crop(3)+ crop(4)),1);