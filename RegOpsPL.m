tic;
xdata = T.Heights;
ydata = flip(D.Shear(:,1)');

% fun = @(x)sseval(x,xdata,ydata);
fun = @(x)sum((ydata - x(1).*(xdata/43).^x(2)).^2);
x0 = [ydata(1),-0.2];
options = optimset('display','final');
bestx = fminsearch(fun,x0,options)
toc


uref = bestx(1);
a = bestx(2);
yfit = uref.*(xdata/43).^a;
% plot(xdata,ydata,'*');
plot(ydata,xdata,'*');
hold on
% plot(xdata,yfit,'r');
plot(yfit,xdata,'r');
xlabel('tdata')
ylabel('Response Data and Curve')
title('Data and Best Fitting Exponential Curve')
legend('Data','Fitted Curve')
hold off