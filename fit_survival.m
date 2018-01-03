function [estimates,eval] = fit_survival(xdata, ydata, x)
% options = optimset('MaxFunEvals',5e6,'MaxIter',5e6,'TolFun',1e-4);
estimates = lsqcurvefit(@testfun, x,xdata,ydata);%,options);
    function F = testfun(x,xdata)
        F=(x(1)-x(4))*(1-0.5*(1+erf((xdata-x(2))/(x(3)*sqrt(2)))))+x(4);
        %        F=0.5*x(1)*(1-erf((xdata-x(2))/(x(3)*sqrt(2))))+x(4);
    end
eval=(estimates(1)-estimates(4)).*(1-0.5.*(1+erf((xdata-estimates(2))./(estimates(3).*sqrt(2)))))+estimates(4);
end