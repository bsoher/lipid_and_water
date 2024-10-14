function [rec_sig,rec_sig1,para] = direct_lorentz_modelling(X_data,Y_data,t,teval) 

% t_len = length(t);
% inp = [ones(t_len,1),-sig_appr];
% outt = 1i*sig_appr.*t;
% x_0 = inp\outt;
if(isempty(t))
    N = length(X_data) + length(Y_data);
    X_pnt = 1:2:N;
    Y_pnt = 2:2:N;
else
    N = length(t);
    X_pnt = t(1:2:N);
    Y_pnt = t(2:2:N);
end

funp = @(x,xdata)((x(1))./(1i*xdata + x(2)));

% fitting 1
t_len = length(X_pnt);
inp = [ones(t_len,1),-X_data];
outt = 1i*X_data.*X_pnt;
para1 = inp\outt;
% para1 = lsqcurvefit(funp,para1,X_pnt,X_data);

% fitting 2
t_len = length(Y_pnt);
inp = [ones(t_len,1),-Y_data];
outt = 1i*Y_data.*Y_pnt;
para2 = inp\outt;
% para2 = lsqcurvefit(funp,para2,Y_pnt,Y_data);

% fitting 3
% funp1 = @(x,xdata,ydata)((-1i*x(1))./((1i*xdata + x(2)).*(1i*ydata + x(2))));
% para = lsqcurvefit(@lowner_fitting,x_0,t,L(:));

para(1) = -1i*para1(1)*para2(1);
para(2) = (para1(2) + para2(2))/2;
rec_sig = funp([1/(2*pi),para(2)],teval);
rec_sig1 = funp(para,teval);