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
fc = num2str(fc);                                                           % Convert Coriolis to string

for i = 1:size(Shear,2)

% ----------------- Ekman Fit -----------------
    [xdata, ydata] = prepareCurveData(flip(T.Heights)',Shear(:,i));         % Set x and y-data for fit

    x0 = Shear(1,i);                                                        % Initial guess for geostrophic wind, G     [m/s]

    model = strcat('sqrt((G * (1 - exp(-x * sqrt(',fc,...                   % Define Ekman model
        '/(2 * K)))*cos(x * sqrt(',fc,...
        '/(2 * K)))))^2 + (G * (exp(-x * sqrt(',fc,...
        '/(2 * K))) * sin(x * sqrt(',fc,'/(2 * K)))))^2)');

    ft = fittype(model,'dependent','y','independent','x',...                % Classify variables
        'coefficients',{'G','K'});

    opts = fitoptions('Method', 'NonlinearLeastSquares');                   % Regression method
    opts.Display    = 'Off';
    opts.StartPoint = [x0, 0.1];                                            % [geostrophic wind, eddy diffusivity]

   try
        
        [fitresult, gof] = fit(xdata, ydata, ft, opts);                     % Perform fit
    
        Ekman.K(i)       = fitresult.K;                                     % Store eddy diffusivity for each profile
        Ekman.G(i)       = fitresult.G;                                     % Store geostrophic wind for each profile
        Ekman.rsquare(i) = gof.rsquare;                                     % Store R^2 value for each fit
        Ekman.RMSE(i)    = gof.rmse;                                        % Store Root Mean Square Error for each fit
    
    catch
    
        Ekman.K(i)       = NaN;                                             % If curve fit fails, store NaN for eddy diffusivity
        Ekman.G(i)       = NaN;                                             % If curve fit fails, store NaN for geostrophic wind
        Ekman.R(i)       = NaN;                                             % If curve fit fails, store NaN for R^2
        Ekman.RMSE(i)    = NaN;                                             % If curve fit fails, store NaN for RMSE
    
    end

    if mod(i,10)==0
        p = i/size(Shear,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % print progress to screen
    end

end

TimeEnd = toc(TimeStart);

fprintf('\n\nTime: %10d minutes \n%16.0f seconds\n', floor(TimeEnd/60), rem(TimeEnd,60));

end