function [PLInflec] = InflectionPL(Shear, T)
%InflectionPL Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n************************')
fprintf('\n*Modified Power Law Fit*')
fprintf('\n************************\n')

ShearProfile = zeros(size(Shear,1),size(Shear,2));

% ----------------- Determine Inflection Points -----------------

for i = 1:size(Shear,2)
    
    dUdz = gradient(Shear(:,i))./gradient(flip(T.Heights)');                % Calculate slope of shear profile
    
    Filter = (dUdz>0);                                                      % Construct a filter for only positive shear
    
    ShearFiltered = Shear(:,i).*Filter;                                     % Apply filter to shear profile
    ShearFiltered(ShearFiltered == 0) = NaN;                                % Set all undesired points to NaN

    [indices] = find(isnan(ShearFiltered));                                 % Find indices of NaN values

    InflecPoint = max(indices);                                             % Find the inflection point of the shear profile

    ShearFiltered(1:InflecPoint) = NaN;                                     % Set all values above inflection to NaN
    
    ShearProfile(:,i) = ShearFiltered;

end

% --------------- Fit Power Law to Modified Shear ---------------

[PLInflec] = PowerLawFit(ShearProfile,T);

end