function [NormalWind] = NormalInflow(Nacelle,Veer,Shear)
%NormalInflow Summary of this function goes here
%   Detailed explanation goes here

    NacelleMatrix = Nacelle .* ones(size(Veer,1),size(Veer,2));             % Preallocate array
    
    Diff = rad2deg(angdiff(deg2rad(NacelleMatrix),deg2rad(Veer)));          % Calculate smallest angle      [deg]
    
    NormalPart = cosd(Diff);                                                % Calculate cosine of angle
    
    NormalWind = Shear .* NormalPart;                                       % Cosine projection of velocity [m/s]

end