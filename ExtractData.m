function [Shear,Veer,Power,I,Time,Nacelle,OrigIndices] = ExtractData(T,data,Indices)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%OrigIndices = 1;

%% Shear Data

    Shear = zeros(length(T.Heights),length(data.DateTime));                                 % Preallocate array             [-]
    
    for i = 1:length(T.Heights)                                                             % Retrieve wind shear           [m/s]
        Shear(i,:) = data.Lidar.(strcat('H',num2str(T.Heights(i),'%02.f'),'m')).WndSpd;
    end

    Shear = Shear(:,Indices);                                                               % Apply initial filter          [m/s]

    Shear(isnan(Shear)) = 0;                                                                % Convert NaN to 0              [m/s]

    Shear(Shear<0)      = 0;                                                                % Convert negative values to 0  [m/s]

    ShearIndices        = all(Shear);                                                       % Identify rows that have zeros [-]

%% Veer Data

    Veer = zeros(length(T.Heights),length(data.DateTime));                                  % Preallocate array             [-]
    
    for i = 1:length(T.Heights)                                                             % Retrieve wind veer            [deg]
        Veer(i,:) = data.Lidar.(strcat('H',num2str(T.Heights(i),'%02.f'),'m')).WndDir;
    end

    Veer = Veer(:,Indices);                                                                 % Apply initial filter          [deg]

    Veer(isnan(Veer))   = 0;                                                                % Convert NaN to 0              [deg]

    Veer(Veer<0)        = 0;                                                                % Convert negative values to 0  [deg]

    VeerIndices         = all(Veer);                                                        % Identify rows that have zeros [-]

%% Power

    Power = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).Power;                              % Retrieve TOI power            [kW]

    Power = Power(Indices);                                                                 % Apply initial filter          [kW]

    Power(isnan(Power)) = 0;                                                                % Convert NaN to 0              [kW]

    Power(Power<0)      = 0;                                                                % Convert negative values to 0  [kW]

    PowerIndices        = (Power ~= 0);                                                     % Identify rows that have zeros [-]

%% Nacelle Heading

    Nacelle = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).NacellePosition_corrected;        % Retrieve nacelle heaing       [deg]

    Nacelle = Nacelle(Indices);                                                             % Apply initial filter          [deg]

%% Turbulence Intensity

    I = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).TurbulenceIntensity;                    % Retrieve turbulence intensity [%]

    I = I(Indices);                                                                         % Apply initial filter          [%]

    I(isnan(I))         = 0;                                                                % Convert NaN to 0              [%]

    I(I<0)              = 0;                                                                % Convert negative values to 0  [%]

    TurbIndices         = (I ~= 0);                                                         % Identify rows that have zeros [-]

%% Datetime

    Time = data.DateTime;                                                                   % Retrieve date and timestamps  [-]

    Time = Time(Indices);                                                                   % Apply initial filter          [-]

%% Create new indices vector from deletions

    OrigIndices = Indices .* ShearIndices .* VeerIndices .* PowerIndices .* TurbIndices;    % Preserve the indices of the original data

    OrigIndices = nonzeros(OrigIndices)';

    NewIndices  = 1:length(Indices);

    Indices     = ShearIndices .* VeerIndices .* PowerIndices .* TurbIndices .* NewIndices; % Update indices                [-]

    Indices     = nonzeros(Indices)';                                                       % Remove deleted elements       [-]

    Shear       = Shear(:,Indices);                                                         % Filter shear array            [m/s]

    Veer        = Veer(:,Indices);                                                          % Filter veer array             [deg]

    Power       = Power(Indices);                                                           % Filter power array            [kW]

    I           = I(Indices);                                                               % Fitler turbulence intensity   [%]

    Nacelle     = Nacelle(Indices);                                                         % Filter power array            [deg]

    Time        = Time(Indices);                                                            % Filter time array             [day/time]
    
    Shear       = flip(Shear);                                                              % Flip array so that top row is 200m and bottom is 43m

    Veer        = flip(Veer);                                                               % Flip array so that top row is 200m and bottom is 43m


%% Validation

% The original data.XXX(OrigIndices) values should match 100% of the new
% filtered data sets.

Elevations = flip(T.Heights);

% SHEAR
fprintf('\n------------------------')
fprintf('\n-------Validation-------')
fprintf('\n------------------------\n')
fprintf('\nSHEAR DATA\n\n')
for i = 1:length(Elevations)
    Equals = sum(Shear(i,:) == data.Lidar.(strcat('H',num2str(Elevations(i),'%02.f'),'m')).WndSpd(OrigIndices));

    if Equals == size(Shear,2)
        fprintf('%3.0f m  -  MATCH\n',Elevations(i))
    else
        fprintf('%3.0f m  -  NO MATCH\n',Elevations(i))
    end
end

% VEER
fprintf('\nVEER DATA\n\n')
for i = 1:length(Elevations)
    Equals = sum(Veer(i,:) == data.Lidar.(strcat('H',num2str(Elevations(i),'%02.f'),'m')).WndDir(OrigIndices));

    if Equals == size(Veer,2)
        fprintf('%3.0f m  -  MATCH\n',Elevations(i))
    else
        fprintf('%3.0f m  -  NO MATCH\n',Elevations(i))
    end
end

% Power
fprintf('\nPOWER DATA\n\n')
    Equals = sum(Power == data.(strcat('BHR',num2str(T.TOI,'%02.f'))).Power(OrigIndices));

    if Equals == length(Power)
        fprintf('          MATCH\n')
    else
        fprintf('          NO MATCH\n')
    end

% Nacelle
fprintf('\nNACELLE HEADING DATA\n\n')
    Equals = sum(Nacelle == data.(strcat('BHR',num2str(T.TOI,'%02.f'))).NacellePosition_corrected(OrigIndices));

    if Equals == length(Nacelle)
        fprintf('          MATCH\n')
    else
        fprintf('          NO MATCH\n')
    end

% Turbulence Intensity
fprintf('\nTURBULENCE INTENSITY DATA\n\n')
    Equals = sum(I == data.(strcat('BHR',num2str(T.TOI,'%02.f'))).TurbulenceIntensity(OrigIndices));

    if Equals == length(I)
        fprintf('          MATCH\n')
    else
        fprintf('          NO MATCH\n')
    end

% Date
fprintf('\nDATETIME DATA\n\n')
    Equals = sum(Time == data.DateTime(OrigIndices));

    if Equals == length(Time)
        fprintf('          MATCH\n')
    else
        fprintf('          NO MATCH\n')
    end

%% Apply Veer Offset

    Veer = Veer + T.Offset;                                                 % Apply LiDAR offset            [deg]
    
    Veer(Veer > 360) = Veer(Veer > 360) - 360;                              % Correct values > 360 degrees  [deg]
    
    Veer(Veer == 360) = 0;  

end