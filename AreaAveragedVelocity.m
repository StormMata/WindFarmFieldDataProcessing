function [AreaMean,AvgComp] = AreaAveragedVelocity(T,Shear,NormalWind)
%AreaAveragedVelocity Calculates the average incoming windspeed over the
%rotor using linear interpolation between measurements points and the
%assumption of horizontal homogeneity.
%   [A,B] = AreaAveragedVelocity(C,D,E)
%           A = Vector containing the average windspeeds for each time
%           B = Structure containing the logic vectors indicating when the
%               area-averaged windspeed is less or greater than the hub
%               height velocity
%               B.Greater = When the area average is greater than hub
%                           height
%               B.Less    = When the area average is less than hub height
%           C = Vector containing heights of wind speed measurements
%           D = Structure containing the optimum fit variables
%               A.alpha   = alpha exponent
%               A.rsquare = R^2 value of fit
%               A.rmse    = Root Mean Square Error of fit
%           E = Matrix containing cosine projection of the incoming wind
%               speeds for each measurement in time and space

TimeStart = tic;

fprintf('\n------------------------')
fprintf('\n------Area Average------')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

AreaMean = zeros(1,size(NormalWind,2));                                     % Preallocate array

% 1. Meshgrid of coordinates

    [ycoords,zcoords] = meshgrid(-63:1:63,169:-1:43);                       % Grid has 1-meter spacing

% 2. Logic matrix to extract values over rotor area

    RotorFilterOutter = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) >= T.HubR;   % Exclude hub area

    RotorFilterInner = abs(sqrt(ycoords.^2 + (zcoords-T.Hub).^2)) <= T.R;   % Only capture values within the outer rotor radius

for i = 1:size(NormalWind,2)

    % 3. Interpolate sheer profile to match grid spacing    
    
        SP = interp1(T.Heights,NormalWind(:,i),43:1:169);                   % Interpolate to 1-meter spacing (linear)
    
    % 4. Create grid with sheer profile
    
        SPgrid = SP' .* ones(size(ycoords));                                % Wind speed is assumed to be constant in y
    
    % 5. Compute dy-dz integral over rotor area
        
        AreaInt = trapz(1:length(ycoords),trapz(1:length(zcoords),...
            SPgrid.*RotorFilterOutter.*RotorFilterInner,2));
    
    % 6. Divide by area of rotor
        
        AreaMean(i) = AreaInt/(pi*T.R^2);                                   % Return area average wind      [m/s]
    
        %AreaMean = mean(nonzeros(SPgrid.*RotorFilterOutter.*RotorFilterInner));

        p = i/size(NormalWind,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % Let me know it's working
end

TimeEnd = toc(TimeStart);

fprintf('\n\nTime: %10d minutes \n%16.0f seconds\n', floor(TimeEnd/60), rem(TimeEnd,60));

AvgComp.Greater = AreaMean > Shear(T.HubRow,:);                             % Area average greater than hub height
AvgComp.Less    = AreaMean < Shear(T.HubRow,:);                             % Area average less than hub height

%% Optional Plot 

% figure
% surf(ycoords,zcoords,SPgrid.*RotorFilterOutter.*RotorFilterInner)
% colorbar
% title('Shear Profile Over Rotor Area')
% xlabel('y (m)')
% ylabel('z (m)')
% ylabel(colorbar('eastoutside'),'Wind Speed (m/s)')

end