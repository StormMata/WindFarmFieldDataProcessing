function [Ekman] = EkmanFit(Shear,T)
%EkmanFit Calculates the eddy diffusivity and geostrophic wind values for
%the given shear profile.
%   [A] = EkmanFit(B,C)
%           A = Structure containing the optimum fit variables
%               A.K       = Eddy diffusivity
%               A.G       = Geostrophic wind
%               A.rsquare = R^2 value of fit
%               A.rmse    = Root Mean Square Error of fit
%           B = Vector containing shear profile
%           C = Vector containing heights of wind speed measurements 

TimeStart = tic;

fprintf('\n------------------------')
fprintf('\n-------Ekman Fits-------')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

fc = 2 * 7.2921159e-5 * sind(T.Lat);                                        % Coriolis term

fits  = zeros(2,2);                                                         % Preallocate matrix holding fit parameters
Rtemp = zeros(1,2);                                                         % Preallocate vector holding goodness-of-fit descriptor
RMSE  = Rtemp;                                                              % Preallocate vector holding goodness-of-fit descriptor

for i = 1:size(Shear,2)

% ----------------- Ekman Fit -----------------
    [xdata, ydata] = prepareCurveData(flip(T.Heights)',Shear(:,i));         % Set x and y-data for fit

%     x0 = Shear(1,i);                                                        % Initial guess for geostrophic wind, G     [m/s]
% 
%     model = strcat('sqrt((G * (1 - exp(-x * sqrt(',fc,...                   % Define Ekman model
%         '/(2 * K)))*cos(x * sqrt(',fc,...
%         '/(2 * K)))))^2 + (G * (exp(-x * sqrt(',fc,...
%         '/(2 * K))) * sin(x * sqrt(',fc,'/(2 * K)))))^2)');
% 
%     ft = fittype(model,'dependent','y','independent','x',...                % Classify variables
%         'coefficients',{'G','K'});
% 
%     opts = fitoptions('Method', 'NonlinearLeastSquares');                   % Regression method
%     opts.Display    = 'Off';
%     opts.StartPoint = [x0, 0.01];                                           % [geostrophic wind, eddy diffusivity]

        options = optimset('MaxIter',2500,'MaxFunEvals',2500, ...           % Set fminsearch options
                           'display','off');

        Norm = mean(Shear(:,i));
   try

        fun = @(x)sum((ydata - sqrt((x(1) .* (1 - exp(-xdata .* sqrt(fc/(2 * x(2)))).*cos(xdata .* sqrt(fc/(2 * x(2)))))).^2 + (x(1) .* (exp(-xdata .* sqrt(fc/(2 * x(2)))) .* sin(xdata .* sqrt(fc/(2 * x(2)))))).^2)).^2);
        
        kops = [1e-6 0.2];

        for j = 1:2

            x0 = [ydata(12),kops(j)];
            
            fits(:,j) = fminsearch(fun,x0,options)';

            if fits(2,j) < 0 
                fits(2,j) = eps(1);
            end
            
            yfit = sqrt((fits(1,j) .* (1 - exp(-xdata .* sqrt(fc/(2 * fits(2,j)))).*cos(xdata .* sqrt(fc/(2 * fits(2,j)))))).^2 + (fits(1,j) .* (exp(-xdata .* sqrt(fc/(2 * fits(2,j)))) .* sin(xdata .* sqrt(fc/(2 * fits(2,j)))))).^2);

            Rtemp(j) = 1-sum((ydata - yfit).^2)/sum((ydata - mean(ydata)).^2);

            RMSE(j) = sqrt(mean((ydata - yfit).^2))/mean(Shear(:,i));

        end

%         if fits(2,1) == fits(2,2)
%             v    = Rtemp(1);
%             maxi = 1;
%         else
            [v,maxi] = max(Rtemp);
%         end

        Ekman.K(i)       = fits(2,maxi);                                     % Store eddy diffusivity for each profile
        Ekman.G(i)       = fits(1,maxi);                                     % Store geostrophic wind for each profile
        Ekman.R(i)       = v;
        Ekman.NRMSE(i)   = RMSE(maxi);                                      % Store normalized Root Mean Square Error for each fit
    
    catch
    
        Ekman.K(i)       = NaN;                                             % If curve fit fails, store NaN for eddy diffusivity
        Ekman.G(i)       = NaN;                                             % If curve fit fails, store NaN for geostrophic wind
        Ekman.R(i)       = NaN;                                             % If curve fit fails, store NaN for R^2
        Ekman.NRMSE(i)   = NaN;                                             % If curve fit fails, store NaN for RMSE
    
    end

    if mod(i,10)==0
        p = i/size(Shear,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % print progress to screen
    end

end

TimeEnd = toc(TimeStart);

fprintf('\n\nTime: %10d minutes \n%16.0f seconds\n', floor(TimeEnd/60), rem(TimeEnd,60));

end