function [NormalWind] = NormalInflow(D)
%NormalInflow Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n-----Normal Inflow------')
fprintf('\n------------------------\n')

    NacelleMatrix = D.Nacelle .* ones(size(D.Veer,1),size(D.Veer,2));       % Preallocate array
    
    Diff = rad2deg(angdiff(deg2rad(NacelleMatrix),deg2rad(D.Veer)));        % Calculate smallest angle      [deg]
    
    NormalPart = cosd(Diff);                                                % Calculate cosine of angle
    
    NormalWind = D.Shear .* NormalPart;                                     % Cosine projection of velocity [m/s]

fprintf('\n\nComplete.\n\n')

end