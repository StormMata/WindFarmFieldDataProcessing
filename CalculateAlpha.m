function [PowLaw] = CalculateAlpha(Shear,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n-----Power Law Fits-----')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

for i = 1:size(Shear,2)

% ----------------- Power Law Fit -----------------
    [xData, yData] = prepareCurveData( flip(T.Heights)',Shear(:,i));        % Required for curve-fitting toolbox
    
    u = num2str(Shear(12,i));                                               % Reference wind speed          [m/s]
    z = num2str(T.Heights(1));                                              % Reference height              [m]
    
    func = strcat(u,'*(x/',z,')^a');                                        % Generate the fit equation per power law
    
    % Set up fittype and options.
    ft = fittype( func, 'independent', 'x', 'dependent', 'y' );             % Required for curve-fitting toolbox
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );                 % Required for curve-fitting toolbox
    opts.Display = 'Off';                                                   % Required for curve-fitting toolbox
    opts.StartPoint = 0.63235924622541;                                     % Required for curve-fitting toolbox
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );                         % Find alpha
    
    PowLaw.alpha(i) = fitresult.a;                                                 % Store alpha
    PowLaw.R(i)     = gof.rsquare;
    PowLaw.RMSE(i)  = gof.rmse;

    if mod(i,10)==0
        p = i/size(Shear,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % Let me know it's working
    end

end

fprintf('\n')

end