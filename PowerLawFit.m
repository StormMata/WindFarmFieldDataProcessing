function [PowLaw] = PowerLawFit(Shear,T,FitProfile)
%PowerLawFit Calculates the optimum alpha exponent for the power law fit
%for the given shear profile.
%   [A] = PowerLawFit(B,C,D)
%           A = Structure containing the optimum fit variables
%               A.alpha   = alpha exponent
%               A.rsquare = R^2 value of fit
%               A.rmse    = Root Mean Square Error of fit
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

        if strcmp(FitProfile,'Inflec') == 1

            dudz  = gradient(Shear(:,i))./gradient(Heights);                % Calculate sign of shear profile
            index = find(dudz < 0, 1, 'last');                              % Find first point of negative shear
            Shear(1:index,i) = NaN;                                         % Set all measurements above that point to NaN
            Heights(1:index) = NaN;                                         % Do same for heights

        end
    
        [xdata, ydata] = prepareCurveData(Heights,Shear(:,i));              % Recommended for curve-fitting toolbox
    
        u = num2str(Shear(end,i));                                          % Reference wind speed          [m/s]
        z = num2str(Heights(end));                                          % Reference height              [m]
        
        model = strcat(u,'*(x/',z,')^a');                                   % Define power law model
        
        ft = fittype(model, 'independent', 'x', 'dependent', 'y');          % Classify variables
    
        opts = fitoptions('Method', 'NonlinearLeastSquares');               % Regression method
        opts.Display    = 'Off';
        opts.StartPoint = 1/7;                                              % Required for curve-fitting toolbox
    
        PowLaw.InflecHeight(i)  = NaN;                                      % Store Inflection height       [m]

        try
        
            [fitresult, gof] = fit(xdata, ydata, ft, opts);                 % Perform fit
            
            PowLaw.alpha(i) = fitresult.a;                                  % Store alpha for each profile
            PowLaw.R(i)     = gof.rsquare;                                  % Store R^2 value for each fit
            PowLaw.RMSE(i)  = gof.rmse;                                     % Store Root Mean Square Error for each fit
    
        catch
        
            PowLaw.alpha(i)         = NaN;                                  % If curve fit fails, store NaN for alpha
            PowLaw.R(i)             = NaN;                                  % If curve fit fails, store NaN for R^2
            PowLaw.RMSE(i)          = NaN;                                  % If curve fit fails, store NaN for Root Mean Square Error
        
        end

        if ~isnan(PowLaw.alpha(i)) & index ~= 1 & index ~= 12

            PowLaw.InflecHeight(i) = Heights(index+1);                      % Store Inflection height       [m]

        end
    
        if mod(i,10)==0
            p = i/size(Shear,2)*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end
    
    end

TimeEnd = toc(TimeStart);

fprintf('\n\nTime: %10d minutes \n%16.0f seconds\n', floor(TimeEnd/60), rem(TimeEnd,60));

end