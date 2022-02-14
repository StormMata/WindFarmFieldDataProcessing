function [] = PlotShearProfile(T,IncidentWind)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% 1. Meshgrid of coordinates

    [ycoords,zcoords] = meshgrid(-63:1:63,169:-1:43);                           % Grid has 1-meter spacing

% 2. Interpolate sheer profile to match grid spacing    

    SP = interp1(T.Heights,IncidentWind,43:1:169);                              % Interpolate to 1-meter spacing (linear)

% 3. Create grid with sheer profile surf(ycoords,zcoords,SPgrid.*RotorFilter)
   
    SPgrid = SP' .* ones(size(ycoords));                                        % Wind speed is assumed to be constant in y

% 4. Logic matrix to extract values over rotor area

    RotorFilterO = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) >= T.HubR;            % Only capture values within the rotor area

    RotorFilterI = abs(sqrt(ycoords.^2 + (zcoords-HT.ub).^2)) <= T.R;               % Only capture values within the rotor area

% 5. Compute dy-dz integral over rotor area
    
    AreaInt = trapz(1:length(ycoords),trapz(1:length(zcoords),SPgrid.*RotorFilterO.*RotorFilterI,2));

% 6. Divide by area of rotor
        
    AreaMean = AreaInt/(pi*T.R^2);                                                % Return area average wind      [m/s]

    %AreaMean = mean(nonzeros(SPgrid.*RotorFilter));

figure
surf(ycoords,zcoords,SPgrid.*RotorFilterO.*RotorFilterI)
colorbar

end