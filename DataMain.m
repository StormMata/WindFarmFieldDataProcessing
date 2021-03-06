%% ---------------------------------------------------------------
%   Author: Storm Mata
%   Functions:
%       1. Preliminary physical values
%       2. Load full windfarm data
%       3. Construct initial filter
%       4. Extract and filter sheer data
% ----------------------------------------------------------------

% Glossary -------------------------------------------------------
% T - structure containing constant values pertaining to turbine
% LiCo, TuCo, InCo - structures containing initial filter values
% Shear   - matrix containing wind speed measurements where m dimension
%           contains measurements in space, n dimension contains
%           measurements in time
% Veer    - matrix containing wind direction measurements. Same
%           orientation as shear
% VeerOff - matrix containing adjusted lidar measurements with offset 8.5
% Power   - vector containing turbine power
% I       - vector containing turbulence intensity
% Time    - vector of datetime values corresponding to measurements
% Nacelle - vector containing turbine headings
% Indices - vector showing last round of filtered value positions
% Diff    - matrix of smallest angle bewteen turbine heading and incoming
%           wind at each height
% Normal  - The cosine of the Diff angles
% Part
% Normal  - The cosine projection of the wind at each height
% Wind


clearvars -except Ekman PLFull PLInflec WndFamMap
%close all
clc

%% 1 - Physical constants and global variables

% Windfarm
    T.TOI      = 62; %[59 62 58];                                                            % Turbine(s) of interest           [-]
    T.Lat      = 23;                                                            % Latitude of wind farm         [deg]

% LiDAR
    T.Heights  = [43 55 67 80 91 104 117 128 141 153 165 200];                  % LiDAR measurement heights     [m]
    T.Offset   = 8.5;                                                           % LiDAR measurement offset      [deg]

% Turbine parameters
    T.Hub      = 106;                                                           % Turbine hub height            [m]
    T.R        = 57;                                                            % Rotor radius                  [m]
    T.HubR     = 5;                                                             % Hub radius                    [m]
    T.GR       = 128.5;                                                         % Generator gear ratio          [-]

% Operational paramters
    T.CutIn    = 2.5;                                                           % Turbine cut-in speed          [m/s]
    T.CutOut   = 25;                                                            % Turbine cut-out speed         [m/s]

% Reference curve
    T.RefWind  = linspace(3, 20, 18);                                           % Reference x axis              [m/s]
    T.RefCurve = [  33,  146,  342,  621, 1008, 1501, 1909, 2076, 2099,...      % Reference power curve         [kW]
                  2100, 2100, 2100, 2100, 2100, 2100, 2100, 2100, 2100];

% Data binning parameters
    T.HubRow   = find(flip(T.Heights)==104);                                    % Find row with hub-height measurements (104m)
    T.numbins  = 40;                                                            % Number of bins                [-]         

% Random seed for sampling
    T.Seed     = 4096;                                                          

%% 2 - Load Data

    if exist('data','var') == 0
        addpath(['/Users/stormmata/Library/Mobile Documents/com~apple~' ...
                 'CloudDocs/Courses/Research/Conditional Averages/' ...
                 'WindFarmFieldDataProcessing']);
        [data] = LoadFullData();
    end

%% 3 - Construct Initial Filter

% LiDAR conditions
    C.LiAv      = 1;                                                            % Lidar.Available               [0-1]

% Turbine conditions
    C.TAv       = 1;                                                            % Turbine.Available             [0-1]
    C.OpSt      = 5;                                                            % Turbine.OpSt                  [1-5]
    C.WFCStrat  = 3;                                                            % Turbine.WFCStrategy           [0-3]
    C.OffsetAp  = 0;                                                            % Turbine.YawOffsetApplied      [0-10]
    C.NacPos    = 1;                                                            % Turbine.NacellePosition       [0,1]
    C.OffsetStp = NaN;                                                          % Turbine.YawOffsetSetPoint     [0-10]
    C.WFCErrCo  = NaN;                                                          % Turbine.WFCErrorCode          [0-120]

% Inflow conditions
    C.E         = 45;                                                           % Max eastward inflow angle     [deg]
    C.W         = 330;                                                          % Max westward inflow angle     [deg]

%% 4 - Extract data and apply initial filter

    [D,Indices] = ExtractDataM(C,T,data);                                       % Performs initial filtering

%% 5 - Find adjusted wind angle

    [NormalWind] = NormalInflow(D);                                             % Find cosine projection of freestream wind [m/s]

%% 6 - Compute Area Average of Incident Wind

    [AAWind,AvgComp] = AreaAveragedVelocity(T,D.Shear,NormalWind);              % Area-averaged wind speed      [m/s]

%% 7 - Fit Power Law to Shear

    [PLFull]   = PowerLawFit(D.Shear,T,'Full');                                 % Fit power law to shear profiles

    [PLInflec] = PowerLawFit(D.Shear,T,'Inflec');                               % Fit power law to shear profiles up to inflection point

%% 8 - Fit Power Law to Shear

    [Ekman] = EkmanFit(D.Shear,T);                                              % Fit power law to shear profiles

%% 9 - Conditional Averages and Stats

    [WindBins,Mean,Num,Dist,STD,Sig] = CondAvgs(T,D,AvgComp);                   % Calculate descriptive and inferential stats for data

%% 10 - Time Series of Power Error

    [Ep] = ErrorTimeSeries(D,WindBins,Mean,T);                                  % Calculate and plot error time series for data

%% 11 - Direction Shear and Probability Plot

    [PDFs] = ShearCharacterization(D,T,PLFull,PLInflec,Ekman);                  % Calculate and plot PDF for direction and speed shear

%% 12 - Alpha/Beta Relationship

    N.SS = [-0.1 0.8];                                                          % Speed shear range [low high]
    N.DS = [-0.1 0.6];                                                          % Direction shear range [low high]
    N.s  = 0.1;                                                                 % Increment size

    [AlphaBeta,AB] = AlphaBetaRel(D,PLFull,PDFs,WindBins,T,N);                  % Calculate and plot the relationship between shear and power

%% 13 - Alpha/Beta U Relationship

    N.SS = [-0.1 0.8];                                                          % Speed shear range [low high]
    N.DS = [-0.1 0.6];                                                          % Direction shear range [low high]
    N.s  = 0.1;                                                                 % Increment size

    [AlphaBetaU,~] = AlphaBetaURel(AAWind,D,PLFull,PDFs,WindBins,T,N);  

%     [AAWind,AvgComp] = AreaAveragedVelocity(T,NormalWind,NormalWind);              % Area-averaged wind speed      [m/s]

%     [AlphaBetaU,~] = AlphaBetaURel(AAWind,NormalWind,PLFull,PDFs,WindBins,T,N);     
%% 14 - Km/Beta Relationship

    N.Km = [0 1];                                                               % Speed shear range [low high]
    N.DS = [-0.1 0.6];                                                          % Direction shear range [low high]
    N.s  = 0.1;                                                                 % Increment size

    [KmBeta,KB] = KmBetaRel(D,PLFull,PDFs,WindBins,T,N);                        % Calculate and plot the relationship between shear and power

%% 15 - Wind Farm Map

    [WndFamMap] = WndFamMap();

%% 16 - Plot Selections

% --- Normalized Wind Speed -----------------------------------------------
P.ShearProfMedians  = 0;    % Median hourly speed shear normalized by U(43)

% --- Turbine Drawing -----------------------------------------------------
P.TurbineSchema     = 0;    % Wind turbine diagram
P.WindSchema        = 0;    % Wind turbine with power law and ekman profiles

% --- Installed Capacity --------------------------------------------------
P.Capacity          = 0;    % Yearly global capacity
P.CapMap            = 0;    % Global capacity heat map
P.USCapMap          = 0;    % US capacity heat map

% --- Map -----------------------------------------------------------------
P.WndFmMap          = 0;    % Wind farm contour map

% --- Polar Plots ---------------------------------------------------------
P.WindRoseFull      = 0;    % Full LiDAR wind data
P.WindRoseAnalysis  = 0;    % Only filtered LiDAR data

% --- Time ----------------------------------------------------------------
P.WdSpHH            = 0;    % Unfiltered Hub height wind speed
P.TI                = 0;    % Unfiltered Turbulence intensity
P.UnFilPow          = 0;    % Unfiltered Power

% --- Time of Day ---------------------------------------------------------
P.EpTODHisto        = 0;    % Day and night hourly power error histograms
P.EpTODMag          = 0;    % Power Error
P.EpTODDiff         = 0;    % Power Error
P.EpTODAvg          = 0;
P.SSTOD             = 0;    % Speed and Direction shear evolution
P.SSTODAVG          = 0;    % Speed and Direction shear evolution, 10-min averages
P.EkTOD             = 0;    % Ekman parameter evolution
P.EkTODAVG          = 0;    % Ekman parameter evolution, 10-min AVERAGES
P.EkTODMED          = 0;    % Ekman parameter evolution, 10-min MEDIANS
P.InflecTODHisto    = 0;    % Inflection height hourly histograms
P.InflecTOD         = 0;
P.FitErrorTOD       = 0;
P.ShearMedError     = 1;    % NRMSE of Ekman and PL by time of day

% --- Counts --------------------------------------------------------------
P.AAPD              = 0;    % Vertical hist of power by wind speed bins, all average
P.GPD               = 0;    % Vertical hist of power by wind speed bins, AA > Hub
P.LPD               = 0;    % Vertical hist of power by wind speed bins, AA < Hub

% --- Wind Speed Bins -----------------------------------------------------
P.HistoPowerAA      = 0;    % Average power with histogram overlaid, all average
P.HistoPowerAALG    = 0;    % Average power with histogram overlaid, all average, AA < Hub, AA > Hub

% --- Speed Shear Alpha ---------------------------------------------------
P.AlphaBetaFull     = 0;    % Full color direction shear heatmap
P.AlphaBetaUFull    = 0;    % Full color direction shear heatmap
P.AlphaBetaMono     = 0;    % Monochromatic color direction shear heatmap
P.AlphaBetaLH       = 0;    % Two-color direction shear heatmap
P.AlphaBetaULH      = 0;    % Two-color direction shear heatmap
P.SSFull            = 0;    % Probability of occurence, Full profile fit
P.SSInflec          = 0;    % Probability of occurence, Partial profile fit

% --- Ekman Parameter -----------------------------------------------------
P.KmBetaFull        = 0;

% --- Direction Shear -----------------------------------------------------
P.DSprob            = 0;    % Probability of occurence

% --- Ekman Parameter -----------------------------------------------------
P.EkmanProb         = 0;    % Probability of occurence   

PlotSelections(P,data,Dist,D,T,WindBins,Mean,STD,Num,AB,AlphaBeta,KB,KmBeta,PDFs,PLFull,PLInflec,Ekman,Ep,WndFamMap)

