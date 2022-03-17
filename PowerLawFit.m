function [PowLaw] = PowerLawFit(Shear,T,FitProfile)
%PowerLawFit Calculates the optimum alpha exponent for the power law fit
%for the given shear profile.
%   [A] = PowerLawFit(B,C,D)
%           A = Structure containing the optimum fit variables
%               A.alpha   = alpha exponent
%               A.R       = R^2 value of fit
%               A.RMSE    = Root Mean Square Error of fit
%           B = Vector containing shear profile
%           C = Vector containing heights of wind speed measurements
%           D = String indicating shear profile to fit curve to
%               'Full'   = Full profile
%               'Inflec' = Up to inflection point (negative shear)

TimeStart = tic;

if strcmp(FitProfile,'Full') == 1

    fprintf('\n------------------------')
    fprintf('\n-----Power Law Fits-----')
    fprintf('\n------------------------\n')
    fprintf('\nComplete:             0')

else

    fprintf('\n************************')
    fprintf('\n*Modified Power Law Fit*')
    fprintf('\n************************\n')
    fprintf('\nComplete:             0')

end

warning('off')                                                              % Turn off warnings from prepareCurveData()

% ----------------- Power Law Fit -----------------

    for i = 1:size(Shear,2)

        Heights = flip(T.Heights)';

        PowLaw.InflecHeight(i)  = NaN;                                      % Store Inflection height       [m]

        if strcmp(FitProfile,'Inflec') == 1

            dudz  = gradient(Shear(:,i))./gradient(Heights);                % Calculate sign of shear profile
            index = find(dudz < 0, 1, 'last');                              % Find first point of negative shear
            Shear(1:index,i) = NaN;                                         % Set all measurements above that point to NaN
            Heights(1:index) = NaN;                                         % Do same for heights

            if index ~= 1 & index ~= 12
    
                PowLaw.InflecHeight(i) = Heights(index+1);                  % Store Inflection height       [m]
    
            end

        end
    
        [xdata, ydata] = prepareCurveData(Heights,Shear(:,i));              % Prepare data

        options = optimset('MaxIter',2500,'MaxFunEvals',2500, ...           % Set fminsearch options
                           'display','off');

        try

            fun = @(x)sum((ydata - x(1).*(xdata/43).^x(2)).^2);             % Define objective function (SSE)

            x0 = [ydata(1),0];                                              % Initial guesses for Uref and a

            fits = fminsearch(fun,x0,options);                              % Optimization routine

            yfit = fits(1).*(xdata/43).^fits(2);                            % Calculate fit values          [m/s]

            PowLaw.Uref(i)  = fits(1);                                      % Store Uref for each profile   [m/s]
            PowLaw.alpha(i) = fits(2);                                      % Store alpha for each profile
            PowLaw.R(i)     = 1-sum((ydata - yfit).^2)/sum((ydata - ...     % Store R^2 value for each fit
                              mean(ydata)).^2);                                  
            PowLaw.NRMSE(i) = sqrt(mean((ydata - yfit).^2))/mean(Shear(:,i));% Store normalized Root Mean Square Error for each fit
    
        catch
        
            PowLaw.alpha(i) = NaN;                                          % If curve fit fails, store NaN for alpha
            PowLaw.R(i)     = NaN;                                          % If curve fit fails, store NaN for R^2
            PowLaw.NRMSE(i) = NaN;                                          % If curve fit fails, store NaN for Root Mean Square Error
        
        end
    
        if mod(i,10)==0
            p = i/size(Shear,2)*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end
    
    end

TimeEnd = toc(TimeStart);

fprintf('\n\nTime: %10d minutes \n%16.0f seconds\n', floor(TimeEnd/60), rem(TimeEnd,60));

end