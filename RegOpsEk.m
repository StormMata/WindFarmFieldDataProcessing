tic;
xdata = T.Heights;
ydata = flip(D.Shear(:,1)');

fc = 2 * 7.2921159e-5 * sind(T.Lat);

% fun = @(x)sum((ydata - sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2)).^2);

% fun = @(x)sqrt(mean((ydata - sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2)).^2));

% sqrt((G * (1 - exp(-x * sqrt(f/(2 * k)))*cos(x * sqrt(f/(2 * k)))))^2 + (G * (exp(-x * sqrt(f/(2 * k))) * sin(x * sqrt(f/(2 * k)))))^2)

% x0 = [ydata(12),1e-6];
% x0 = [ydata(12),0.2];

        options = optimset('MaxIter',2500,'MaxFunEvals',2500,'display','off');

        fun = @(x)sum((ydata - sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2)).^2);
        
        kops = [1e-6 0.1];

        for i = 1:2

            x0 = [ydata(12),kops(i)];
            
            fits(:,i) = fminsearch(fun,x0,options)';
            
            yfit = sqrt((fits(1,i) .* (1 - exp(-xdata .* sqrt(fc/(2 * fits(2,i)))).*cos(xdata .* sqrt(fc/(2 * fits(2,i)))))).^2 + (fits(1,i) .* (exp(-xdata .* sqrt(fc/(2 * fits(2,i)))) .* sin(xdata .* sqrt(fc/(2 * fits(2,i)))))).^2);
            
            Rtemp(i) = 1-sum((ydata - yfit).^2)/sum((ydata - mean(ydata)).^2);

        end

[v,i] = max(Rtemp)




% options = optimset('display','final','Maxiter',7000,'maxfunevals',7000);
% bestx = fminsearch(fun,x0,options)
% toc
% 
% yfit = sqrt((bestx(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * bestx(2)))).*cos(xdata .* sqrt(fc/(2 * bestx(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * bestx(2)))) .* sin(xdata .* sqrt(fc/(2 * bestx(2)))))).^2);
% 
% 1-sum((ydata - yfit).^2)/sum((ydata - mean(ydata)).^2)

% sqrt(mean((ydata - sqrt((bestx(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * bestx(2)))).*cos(xdata .* sqrt(fc/(2 * bestx(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * bestx(2)))) .* sin(xdata .* sqrt(fc/(2 * bestx(2)))))).^2)).^2))

% 1-sum((ydata - sqrt((bestx(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * bestx(2)))).*cos(xdata .* sqrt(fc/(2 * bestx(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * bestx(2)))) .* sin(xdata .* sqrt(fc/(2 * bestx(2)))))).^2)).^2)/sum((ydata - mean(ydata)).^2)

% x = bestx;
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

% x = 0:500;
% G = 10;
% k = 6;
% f = 2 * 7.2921159e-5 * sind(T.Lat);
% u = sqrt((G .* (1 - exp(-x .* sqrt(f/(2 * k))).*cos(x .* sqrt(f/(2 * k))))).^2 + (G .* (exp(-x .* sqrt(f/(2 * k))) .* sin(x .* sqrt(f/(2 * k))))).^2);
% hold on
% plot(u,x)
% 
% x = 0:200;G = Ekmantest2.G(1);k = Ekmantest2.K(1);f = 2 * 7.2921159e-5 * sind(T.Lat); u = sqrt((G .* (1 - exp(-x .* sqrt(f/(2 * k))).*cos(x .* sqrt(f/(2 * k))))).^2 + (G .* (exp(-x .* sqrt(f/(2 * k))) .* sin(x .* sqrt(f/(2 * k))))).^2);