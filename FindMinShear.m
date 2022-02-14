function [ShearIndices,ShearInt] = FindMinShear(T,NormalWind)

fprintf('\n------------------------')
fprintf('\n----Shear Magnitude-----')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

ShearInt = zeros(1,size(NormalWind,2));                                             % Preallocate array

% 1. Meshgrid of coordinates

    [ycoords,zcoords] = meshgrid(-63:1:63,169:-1:43);                               % Grid has 1-meter spacing

% 2. Logic matrix to extract values over rotor area

    RotorFilterOutter = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) >= T.HubR;       % Exclude hub area

    RotorFilterInner = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) <= T.R;           % Only capture values within the outer rotor radius

for i = 1:size(NormalWind,2)

    % 3. Interpolate sheer profile to match grid spacing    
    
        SP = interp1(T.Heights,NormalWind(:,i),43:1:169);                           % Interpolate to 1-meter spacing (linear)

        SP = SP - min(SP);                                                          % only want the shear magnitude
    
    % 4. Create grid with sheer profile
    
        SPgrid = SP' .* ones(size(ycoords));                                        % Wind speed is assumed to be constant in y
    
    % 5. Compute dy-dz integral over rotor area
        
        ShearInt(i) = trapz(1:length(ycoords),trapz(1:length(zcoords),SPgrid.*RotorFilterOutter.*RotorFilterInner,2));

        p = i/size(NormalWind,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % Let me know it's working
end

ShearIndices = 1:length(ShearInt);

[~,SortIndices] = sort(ShearInt,'ascend');

ShearIndices = ShearIndices(SortIndices);

fprintf('\n')

% figure
% surf(ycoords,zcoords,SPgrid.*RotorFilterOutter.*RotorFilterInner)
% colorbar
% title('Shear Profile Over Rotor Area')
% xlabel('y (m)')
% ylabel('z (m)')
% ylabel(colorbar('eastoutside'),'Wind Speed (m/s)')

end