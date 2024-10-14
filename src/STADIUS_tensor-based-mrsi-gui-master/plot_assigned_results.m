function plot_assigned_results(src,H,row,col,dist_max,ppm_axis_seg,no_tp,assignment,ratio)

no_src = size(H,2);

if(~isempty(src))
    plot_assigned_sources(src,ppm_axis_seg,assignment);
end




[~,idx] = max(H,[],2);

titl{1} = 'Tumor';
titl{2} = 'Necrosis';
titl{3} = 'Control';
titl{4} = 'Bad';
titl{5} = 'Unidentified';
titl{6} = 'Nornal/Tumor';
titl{7} = 'Tumor/Necrosis';
titl{8} = 'Nornal/Necrosis';
titl{9} = 'Nornal/Tumor/Necrosis';

 
switch(no_src)
    case 1
        figure();
        if(no_tp == 1)
            imagesc(reshape(H,row,col));
            title(titl{assignment});
        else
            imagesc(reshape(H,row,col)');title(titl{assignment});
        end
    case 2
        X = 1;Y = 2;
    case 3
        X = 2;Y = 2;
    case 4        
        X = 2;Y = 2;
    case 5
        X = 2;Y = 3;
    case 6
        X = 2;Y = 3;
    case 7
%         X = 2;Y = 2;
%         X1 = 2;Y1 = 2;
        X = 3;Y = 3;
        X1 = 2;Y1 = 3;
    case 8        
        X = 3;Y = 3;
        X1 = 2;Y1 = 2;
    case 9
        X = 3;Y = 3;
        X1 = 2;Y1 = 3;
    case 10        
        X = 4;Y = 3;
        X1 = 2;Y1 = 3;
end

if(no_src>1)
    %if(no_src<7)
        figure();
        for i = 1:no_src
            if(dist_max == 1)
                subplot(X,Y,i);
                if(no_tp == 1)
                    imagesc(reshape(H(:,i),row,col));
                else
                    imagesc(reshape(H(:,i),row,col)');
                end
                if(isempty(ratio))
                    title(titl{assignment(i)}); 
                else
                    title(strcat(num2str(ratio(i)), titl{assignment(i)})); 
                end
            else
                temp = idx ==i;
                subplot(X,Y,i);
                if(no_tp == 1)
                    imagesc(reshape(temp,row,col));
                else
                    imagesc(reshape(temp,row,col)');
                end
                if(isempty(ratio))
                    title(titl{assignment(i)}); 
                else
                    title(strcat(num2str(ratio(i)), titl{assignment(i)})); 
                end
            end
        end
        
%     else
%         figure();
%         for i = 1:4
%             if(dist_max == 1)
%                 subplot(X,Y,i);
%                 if(no_tp == 1)
%                     imagesc(reshape(H(:,i),row,col));
%                 else
%                     imagesc(reshape(H(:,i),row,col)');
%                 end
%                 if(isempty(ratio))
%                     title(titl{assignment(i)}); 
%                 else
%                     title(strcat(num2str(ratio(i)), titl{assignment(i)})); 
%                 end
%             else
%                 temp = idx ==i;
%                 subplot(X,Y,i);
%                 if(no_tp == 1)
%                     imagesc(reshape(temp,row,col));
%                 else
%                     imagesc(reshape(temp,row,col)');
%                 end
%                 if(isempty(ratio))
%                     title(titl{assignment(i)}); 
%                 else
%                     title(strcat(num2str(ratio(i)), titl{assignment(i)})); 
%                 end
%             end
%         end
%         figure();
%         for i = 5:no_src
%             if(dist_max == 1)
%                 subplot(X1,Y1,i-4);
%                 if(no_tp == 1)
%                     imagesc(reshape(H(:,i),row,col));
%                 else
%                     imagesc(reshape(H(:,i),row,col)');
%                 end
%                 if(isempty(ratio))
%                     title(titl{assignment(i)}); 
%                 else
%                     title(strcat(num2str(ratio(i)), titl{assignment(i)})); 
%                 end
%             else
%                 temp = idx ==i;
%                 subplot(X1,Y1,i-4);
%                 if(no_tp == 1)
%                     imagesc(reshape(temp,row,col));
%                 else
%                     imagesc(reshape(temp,row,col)');
%                 end
%                 if(isempty(ratio))
%                     title(titl{assignment(i)}); 
%                 else
%                     title(strcat(num2str(ratio(i)), titl{assignment(i)})); 
%                 end
%             end
%         end
    end
end
% colormap hot