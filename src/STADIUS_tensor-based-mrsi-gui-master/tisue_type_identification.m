function assignment = tisue_type_identification(W,ppmaxis)
% assign spectrum to a certain tissue pattern according to prior knowledge
% priorknow can currently be 'GBM' or 'lowgradeglioma'

% use default prior knowledge for brain tumors

% ppmaxis = ppm_axis_seg;
% W = real(avg_src_6);

ppmrange = [0.75 1.45;     % Lips, marker for necrosis
    1.95 2.10;... % NAA, marker for normal
    2.92 3.08;... % Creatine marker 
    3.13 3.30];% Cho, marker for active tumor

freqindex = [];
for i = 1:size(ppmrange,1)
    freqindex = [freqindex; [find(ppmaxis>min(ppmrange(i,:)),1,'first'), find(ppmaxis<max(ppmrange(i,:)),1,'last')]];
end
no_of_sig  = size(W,2);
for j = 1:size(freqindex,1)
    maxpeak(j,:) = max(W(freqindex(j,1):freqindex(j,2),:));
end

for i = 1:no_of_sig
    [~,max_ind] = max(W(:,i));
%     loc = max_ind>freqindex(:,1).*max_ind>freqindex(:,1);
    if(max_ind>freqindex(1,1) && max_ind<freqindex(1,2))
        if(maxpeak(1,i)/sum(maxpeak(:,i)) > .6) % Lips/(lips+NAA+Cho)>60% => infer necrosis
            assignment(i) = 2; %'N'
        elseif(0.8*maxpeak(2,i) > maxpeak(4,i))
            assignment(i) = 3; %'C'
        elseif(maxpeak(4,i) > maxpeak(2,i))
            assignment(i) = 1; %'T'
        else
            assignment(i) = 5; %'U'
        end
%         disp('i was here 1')
    elseif(max_ind>freqindex(2,1) && max_ind<freqindex(2,2))
        if(maxpeak(2,i)/sum(maxpeak(:,i)) > .6)
            assignment(i) = 4; %'A'
        elseif(0.8*maxpeak(2,i) > maxpeak(4,i))
            assignment(i) = 3; %'C'
        elseif(0.8*maxpeak(2,i) < maxpeak(4,i))
            assignment(i) = 1; %'T'
        else
            assignment(i) = 5; %'U'
        end
%         disp('i was here 2')
    elseif(max_ind>freqindex(3,1) && max_ind<freqindex(3,2))
        if(maxpeak(3,i)/sum(maxpeak(:,i)) > .6)
            assignment(i) = 4; %'A'
        elseif(0.8*maxpeak(2,i) > maxpeak(4,i))
            assignment(i) = 3; %'C'
        elseif(0.8*maxpeak(2,i) < maxpeak(4,i))
            assignment(i) = 1; %'T'
        else
            assignment(i) = 5; %'U'
        end
%         disp('i was here 3')
    elseif(max_ind>freqindex(4,1) && max_ind<freqindex(4,2))
        if(maxpeak(4,i)/sum(maxpeak(:,i)) > .6)
            assignment(i) = 4; %'A'
        else
            assignment(i) = 1; %'T'
        end
%         disp('i was here 4')
    else
        assignment(i) = 5; %'U'
%         disp('i was here 5')
    end
end
%% 

% x = [freqindex(:,1),freqindex(:,1);freqindex(:,2),freqindex(:,2)];
% y1 = [min(real(avg_src_2(:,1)))*ones(8,1),max(real(avg_src_2(:,1)))*ones(8,1)];
% y2 = [min(real(avg_src_2(:,2)))*ones(8,1),max(real(avg_src_2(:,2)))*ones(8,1),];
% x_ppm = [ppmaxis(freqindex(:,1))',ppmaxis(freqindex(:,1))';ppmaxis(freqindex(:,2))',ppmaxis(freqindex(:,2))'];
% figure();plot(ppm_axis_seg,real(avg_src_2(:,2)));line(x_ppm',y2');axis tight







% % extract max value in each frequency region
% maxpeak = zeros(size(freqindex,1),2);
% for i = 1:size(freqindex,1)
%     maxpeak(i,:) = max(W(freqindex(i,1):freqindex(i,2),:));
% end
% 
% if (nargin<3) || isempty(priorknow)
%     % infer whether necrosis might be present
%     if (maxpeak(3,1)/sum(maxpeak(:,1)) > .6) || (maxpeak(3,2)/sum(maxpeak(:,2)) > .6) % Lips/(lips+NAA+Cho)>60% => infer necrosis
%         assignment = assign_2tissues1(W,ppmaxis,'GBM');
%     else
%         assignment = assign_2tissues1(W,ppmaxis,'lowgradeglioma');
%     end
% end
% 
% if strcmp(priorknow,'GBM') == 1
%     % assign necrosis vs tumor
%     [m,ind] = max(maxpeak(3,:)./maxpeak(2,:)); % max Lip/Cho is necrosis
%     
%     assignment{ind} = 'necrosis';
%     assignment{3 - ind} = 'tumor';
% elseif strcmp(priorknow,'lowgradeglioma') == 1
%     % assign normal vs abnormal
%     [m,ind] = max(maxpeak(1,:)./maxpeak(2,:)); % max NAA/Cho is normal
% 
%     assignment{ind} = 'normal';
%     assignment{3 - ind} = 'tumor';
% end

