function [PowLaw] = PowerLawFit(Shear,T)
%PowerLawFit Calculates the optimum alpha exponent for the power law fit
%for the given shear profile.
%   [A] = PowerLawFit(B,C)
%           A = Structure containing the optimum fit variables
%               A.alpha   = alpha exponent
%               A.rsquare = R^2 value of fit
%               A.rmse    = Root Mean Square Error of fit
%           B = Vector containing shear profile
%           C = Vector containing heights of wind speed measurements

warning('off')

fprintf('\n------------------------')
fprintf('\n-----Power Law Fits-----')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

for i = 1:size(Shear,2)

% ----------------- Power Law Fit -----------------

    Filter = ~isnan(Shear(:,i));
    Heights = flip(T.Heights)'.*Filter;
    Heights(Heights == 0) = NaN;

    [xdata, ydata] = prepareCurveData(Heights,Shear(:,i));                  % Required for curve-fitting toolbox

    u = num2str(Shear(end,i));                                              % Reference wind speed          [m/s]
    z = num2str(T.Heights(1));                                              % Reference height              [m]
    
    model = strcat(u,'*(x/',z,')^a');                                       % Define power law model
    
    ft = fittype(model, 'independent', 'x', 'dependent', 'y');              % Classify variables

    opts = fitoptions('Method', 'NonlinearLeastSquares');                   % Regression method
    opts.Display    = 'Off';
    opts.StartPoint = 1/7;                                                  % Required for curve-fitting toolbox

try

    [fitresult, gof] = fit(xdata, ydata, ft, opts);                         % Perform fit
    
    PowLaw.alpha(i) = fitresult.a;                                          % Store alpha for each profile
    PowLaw.R(i)     = gof.rsquare;                                          % Store R^2 value for each fit
    PowLaw.RMSE(i)  = gof.rmse;                                             % Store Root Mean Square Error for each fit

catch

    PowLaw.alpha(i) = NaN;                                                  % If curve fit fails, store Nan for alpha
    PowLaw.R(i)     = NaN;                                                  % If curve fit fails, store NaN for R^2
    PowLaw.RMSE(i)  = NaN;                                                  % If curve fit fails, store NaN for Root Mean Square Error

end

    if mod(i,10)==0
        p = i/size(Shear,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % print progress to screen
    end

end

fprintf('\n')

end