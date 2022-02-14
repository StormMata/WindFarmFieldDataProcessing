function [VeerIndices] = FindMinVeer(T,Diff)

fprintf('\n------------------------')
fprintf('\n--Direction Magnitude---')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

VeerInt = zeros(1,size(Diff,2));                                                    % Preallocate array

% 1. Meshgrid of coordinates

    [ycoords,zcoords] = meshgrid(-63:1:63,169:-1:43);                               % Grid has 1-meter spacing

% 2. Logic matrix to extract values over rotor area

    RotorFilterOutter = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) >= T.HubR;       % Exclude hub area

    RotorFilterInner = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) <= T.R;           % Only capture values within the outer rotor radius

for i = 1:size(Diff,2)

    % 3. Interpolate sheer profile to match grid spacing    
    
        VP = interp1(T.Heights,Diff(:,i),43:1:169);                                 % Interpolate to 1-meter spacing (linear)

        VP = abs(VP);
    
    % 4. Create grid with sheer profile
    
        SPgrid = VP' .* ones(size(ycoords));                                        % Wind speed is assumed to be constant in y
    
    % 5. Compute dy-dz integral over rotor area
        
        VeerInt(i) = trapz(1:length(ycoords),trapz(1:length(zcoords),SPgrid.*RotorFilterOutter.*RotorFilterInner,2));

        p = i/size(Diff,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                          % Let me know it's working
end

VeerIndices = 1:length(VeerInt);

[~,SortIndices] = sort(VeerInt,'ascend');

VeerIndices = VeerIndices(SortIndices);

fprintf('\n')

% figure
% surf(ycoords,zcoords,SPgrid.*RotorFilterOutter.*RotorFilterInner)
% colorbar
% title('Shear Profile Over Rotor Area')
% xlabel('y (m)')
% ylabel('z (m)')
% ylabel(colorbar('eastoutside'),'Wind Direction (deg)')

end