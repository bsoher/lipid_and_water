 function [] = plot_assigned_sources(src,ppm_axis_seg,assignment)

no_src = size(src,2);
% for i = 1:no_src
%     figure();plot(src(:,i));
%     set(gca, 'xdir','reverse')
%     axis tight
% end

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
        if(isempty(ppm_axis_seg))
            figure();plot(src(:,1));
        else
            figure();plot(ppm_axis_seg,src(:,1));
            xlabel('ppm');
        end
        
        set(gca, 'xdir','reverse')
        ylabel('a.u');
%         title(titl{assignment}); 
        axis tight
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
           subplot(X,Y,i);
           if(isempty(ppm_axis_seg))
               plot(src(:,i));
           else
            plot(ppm_axis_seg,src(:,i));
            xlabel('ppm');
           end
           set(gca, 'xdir','reverse')
           ylabel('a.u');
           axis tight
%            title(titl((i-1)*3+1:(i-1)*3+3));
           title(titl{assignment(i)}); 
        end
        
%     else
%         figure();
%         for i = 1:4
%             subplot(X,Y,i);
%             if(isempty(ppm_axis_seg))
%                 plot(src(:,i));
%             else
%                 plot(ppm_axis_seg,src(:,i));
%                 xlabel('ppm');
%             end
%             set(gca, 'xdir','reverse')
%             ylabel('a.u');
%             axis tight
%             title(titl{assignment(i)}); 
%         end
%         figure();
%         for i = 5:no_src
%             subplot(X1,Y1,i-4);
%             if(isempty(ppm_axis_seg))
%                 plot(src(:,i));
%             else
%                 plot(ppm_axis_seg,src(:,i));
%                 xlabel('ppm');
%             end
%             set(gca, 'xdir','reverse')
%             ylabel('a.u');
%             axis tight
%             title(titl{assignment(i)}); 
%         end
    end
end