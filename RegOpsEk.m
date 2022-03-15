tic;
xdata = T.Heights;
ydata = flip(D.Shear(:,161)');

fc = 2 * 7.2921159e-5 * sind(T.Lat);

fun = @(x)sum((ydata - sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2)).^2);


% x0 = [ydata(12),1e-6];
x0 = [ydata(12),0.04];
options = optimset('display','final');
bestx = fminsearch(fun,x0,options)
toc

yfit = sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2);

1-sum((ydata - yfit).^2)/sum((ydata - mean(ydata)).^2)

x = bestx;
% yfit = sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2)
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