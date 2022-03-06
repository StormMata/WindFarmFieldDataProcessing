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


clearvars -except data Ekman PLFull PLInflec
%close all
clc

%% Section 1 - Physical constants and global variables

% Windfarm
    T.TOI      = 62;                                                            % Turbine of interest           [-]
    T.Lat      = 23;                                                            % Latitude of wind farm         [deg]

% LiDAR
    T.Heights  = [43 55 67 80 91 104 117 128 141 153 165 200];                  % LiDAR measurement heights     [m]
    T.Offset   = 8.5;                                                           % LiDAR measurement offset      [deg]

% Turbine parameters
    T.Hub      = 106;                                                           % Turbine hub height            [m]
    T.R        = 57;                                                            % Rotor radius                  [m]
    T.HubR     = 5;                                                             % Hub radius                    [m]

% Operational paramters
    T.CutIn    = 2.5;                                                           % Turbine cut-in speed          [m/s]
    T.CutOut   = 25;                                                            % Turbine cut-out speed         [m/s]

% Reference curve
    T.RefWind  = linspace(3, 20, 18);                                           % Reference x axis              [m/s]
    T.RefCurve = [33, 146, 342, 621, 1008, 1501, 1909, 2076, 2099,...           % Reference power curve         [kW]
                   2100, 2100, 2100, 2100, 2100, 2100, 2100, 2100, 2100];

% Data binning parameters
    T.HubRow   = find(flip(T.Heights)==104);                                    % Find row with hub-height measurements (104m)
    T.numbins  = 40;                                                            % Number of bins                [-]         

%% Section 2 - Load Data

if exist('data','var') == 0
    [data] = LoadFullData();
end

%% Section 3 - Construct Initial Filter

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

% Create filter
    Indices = IntFilter(C, T, data);                                            % Return indices of filtered data

clear C

%% Section 4 - Extract data and apply initial filter

    [D,Indices] = ExtractData(T,data,Indices);                                  % Performs initial filtering

%% Section 5 - Find adjusted wind angle

    [NormalWind] = NormalInflow(D);                                             % Find cosine projection of freestream wind [m/s]

%% Section 6 - Compute Area Average of Incident Wind

    [AAWind,AvgComp] = AreaAveragedVelocity(T,D.Shear,NormalWind);              % Area-averaged wind speed      [m/s]

%% Section 7 - Fit Power Law to Shear

    [PLFull]   = PowerLawFit(D.Shear,T,'Full');                                 % Fit power law to shear profiles

    [PLInflec] = PowerLawFit(D.Shear,T,'Inflec');                               % Fit power law to shear profiles up to inflection point

%% Section 8 - Fit Power Law to Shear

    [Ekman] = EkmanFit(D.Shear,T);                                              % Fit power law to shear profiles

%% Section 9 - Conditional Averages and Stats

    [WindBins,Mean,Num,Dist,STD,Sig] = CondAvgs(T,D,AvgComp);                   % Calculate descriptive and inferential stats for data

%% Section 10 - Time Series of Power Error

    [Ep] = ErrorTimeSeries(D,WindBins,Mean,T);                                  % Calculate and plot error time series for data

%% Section 11 - Direction Shear and Probability Plot

    [PDFs] = ShearCharacterization(D,T,PLFull,PLInflec,Ekman);                  % Calculate and plot PDF for direction and speed shear

%% Section 12 - Alpha/Beta Relationship

    N.SS = [-0.1 0.8];                                                          % Speed shear range [low high]
    N.DS = [-0.1 0.6];                                                          % Direction shear range [low high]
    N.s  = 0.1;                                                                 % Increment size

    [AlphaBeta,AB] = AlphaBetaRel(D,PLFull,PDFs,WindBins,T,N);                  % Calculate and plot the relationship between shear and power

%% 

angavg  = mean(reshape(data.BHR62.Pitch,10,[]));
wndbin  = data.Lidar.H104m.WndSpd(5:10:end);

angavg(isnan(angavg)) =  0;
wndbin(isnan(angavg)) =  0;

angavg(isnan(wndbin)) =  0;
wndbin(isnan(wndbin)) =  0;

bins = 0:0.5:floor(max(wndbin));

for i = 1:length(bins)
    MAD(i)  = mad((nonzeros((wndbin > bins(i) & wndbin <= (bins(i)+0.5)) .* angavg)));
end
%     MADhigh = angavg+4.5*MAD;
%     MADlow  = angavg-4.5*MAD;

scatter(wndbin,angavg)
scatter([D.Shear(7,5:10:end) 0 0 0 0 0 0 0],mean(reshape([D.Pitch zeros(1,67)],10,[])))
figure

    hold on
    scatter(wndbin(angavg>MADhigh),angavg(angavg>MADhigh),'red',Marker='.')
    scatter(wndbin(angavg<MADlow),angavg(angavg<MADlow),'red',Marker='.')
    scatter(wndbin(angavg<=MADhigh & angavg>=MADlow),angavg(angavg<=MADhigh & angavg>=MADlow),'blue',Marker='.')





%% Plot Selection

% --- Polar Plots ---------------------------------------------------------
P.WindRoseFull      = 0;    % Full LiDAR wind data
P.WindRoseAnalysis  = 0;    % Only filtered LiDAR data

% --- Time ----------------------------------------------------------------
P.WdSpHH            = 0;    % Hub height wind speed
P.TI                = 0;    % Turbulence intensity

% --- Time of Day ---------------------------------------------------------
P.EpTODHisto        = 0;    % Day and night hourly power error histograms
P.EpTOD             = 0;    % Power Error
P.EpTODAvg          = 0;
P.SSTOD             = 0;    % Speed and Direction shear evolution
P.SSTODAVG          = 0;    % Speed and Direction shear evolution, 10-min averages
P.EkTOD             = 0;    % Ekman parameter evolution
P.EkTODAVG          = 0;    % Ekman parameter evolution, 10-min AVERAGES
P.EkTODMED          = 0;    % Ekman parameter evolution, 10-min MEDIANS
P.InflecTODHisto    = 0;    % Inflection height hourly histograms
P.InflecTOD         = 0;
P.FitErrorTOD       = 1;


% --- Counts --------------------------------------------------------------
P.AAPD              = 0;    % Vertical hist of power by wind speed bins, all average
P.GPD               = 0;    % Vertical hist of power by wind speed bins, AA > Hub
P.LPD               = 0;    % Vertical hist of power by wind speed bins, AA < Hub

% --- Wind Speed Bins -----------------------------------------------------
P.HistoPowerAA      = 0;    % Average power with histogram overlaid, all average
P.HistoPowerAALG    = 0;    % Average power with histogram overlaid, all average, AA < Hub, AA > Hub

% --- Speed Shear Alpha ---------------------------------------------------
P.AlphaBetaFull     = 0;    % Full color direction shear heatmap
P.AlphaBetaMono     = 0;    % Monochromatic color direction shear heatmap
P.AlphaBetaLH       = 0;    % Two-color direction shear heatmap
P.SSFull            = 0;    % Probability of occurence, Full profile fit
P.SSInflec          = 0;    % Probability of occurence, Partial profile fit

% --- Direction Shear -----------------------------------------------------
P.DSprob            = 0;    % Probability of occurence

% --- Ekman Parameter -----------------------------------------------------
P.EkmanProb         = 0;    % Probability of occurence   




PlotSelections(P,data,Dist,D,T,WindBins,Mean,STD,Num,AB,AlphaBeta,PDFs,PLFull,PLInflec,Ekman,Ep)










%% Ep v. Alpha Full 

R = range(PLFull.alpha);
Bins = 40;
Rstep = R/Bins;
xaxis = linspace(min(PLFull.alpha),max(PLFull.alpha),Bins);
    
    for i = 1:length(xaxis)

        Indices = (PLFull.alpha <= min(PLFull.alpha)+(i*Rstep));
    
        Ep.EpvALPHA(i,:) = Indices .* Ep.abs;

    end
    
    for i = 1:size(Ep.EpvALPHA,1)
    
        Ep.EpvALPHAmean(i) = mean(nonzeros(Ep.EpvALPHA(i,:)));
    
    end

    figure
        plot(xaxis,Ep.EpvALPHAmean,'Color','k')
            ylabel('\epsilon_P (kW)')
            xlabel('\alpha (-)')
            title('\epsilon_P vs. Power Law Alpha (Full Profile)')

%% Ep v. Alpha Partial

R = range(PLInflec.alpha);
Bins = 40;
Rstep = R/Bins;
xaxis = linspace(min(PLInflec.alpha),max(PLInflec.alpha),Bins);
    
    for i = 1:length(xaxis)

        Indices = (PLInflec.alpha <= min(PLInflec.alpha)+(i*Rstep));
    
        Ep.EpvALPHAinflec(i,:) = Indices .* Ep.abs;

    end
    
    for i = 1:size(Ep.EpvALPHAinflec,1)
    
        Ep.EpvALPHAinflecmean(i) = mean(nonzeros(Ep.EpvALPHAinflec(i,:)));
    
    end

    figure
        plot(xaxis,Ep.EpvALPHAinflecmean,'Color','k')
            ylabel('\epsilon_P (kW)')
            xlabel('\alpha (-)')
            title('\epsilon_P vs. Power Law Alpha (Partial Profile)')

%% ------ Item 1 ------

% for i = 1:length(Power)
% 
% end
% 
% plot(Power)
% 
% dateaxis = datetime(Time,'ConvertFrom','datenum');
% scatter(dateaxis,Power,'.')

[ShearIndices,ShearInt] = FindMinShear(T,NormalWind);

figure
histogram(ShearInt,100)
hold on
xline(quantile(ShearInt,0.25),'Color','r','LineWidth',1)
xline(quantile(ShearInt,0.5),'Color','r','LineWidth',1)
xline(quantile(ShearInt,0.75),'Color','r','LineWidth',1)
title('Distribution of Shear Magnitude')
xlabel('Integral of Shear Profile Over Rotor Area (m^3/s)')
ylabel('Counts')
hold off

FirstQuart = ShearIndices(1:floor(0.25*length(ShearIndices)));

PowerFQ = Power(FirstQuart);

AAWindFQ = AAWind(FirstQuart);

FirstQuartShear = Shear(:,FirstQuart);

AvgCompFQ.Greater = AAWindFQ > FirstQuartShear(T.HubRow,:);                               % Area average greater than hub height
AvgCompFQ.Less    = AAWindFQ < FirstQuartShear(T.HubRow,:);                               % Area average less than hub height

[WindBinsFQ,MeanFQ,NumFQ,DistFQ,STDFQ,SigFQ] = CondAvgs(T,FirstQuartShear,PowerFQ,AvgCompFQ);

figure
plot(WindBins,Mean.All,'LineWidth',1.5)
hold on
plot(WindBinsFQ,MeanFQ.All,'LineWidth',1.5)
title('Effect of Shear on Power Output')
xlabel('Wind Speed (m/s)')
ylabel('Power (kW)')
legend('All Average','All Average with First Quartile Shear Magnitudes')
grid on
hold off

[VeerIndices] = FindMinVeer(T,Diff);













%% Section 10 - Plot

BetzWind  = T.CutIn:0.5:9;                                                      % Reference Betz limit plot
BetzPower = (0.593 * 0.5 * 1.162 * pi * T.R^2 .* BetzWind.^3)/1e3;              % Reference Betz limit plot

xi = linspace(min(T.RefWind), max(T.RefWind), 100);                             % Evenly-Spaced Interpolation Vector
yi = interp1(T.RefWind, T.RefCurve, xi, 'spline', 'extrap');

figure;
%    plot(BetzWind,BetzPower)%,'color','k','LineStyle',':','Marker','*')
    hold on
  plot(xi,yi,'LineWidth',0.85,'Color','k','LineStyle',':')
   plot(WindBins,Mean.All,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
   plot(WindBins,Mean.Larger,'LineWidth',1.5,'Color','#D95319','LineStyle','--')
   plot(WindBins,Mean.Less,'LineWidth',1.5,'Color','#A2142F','LineStyle','-.')
%      errorbar(WindBins, AllMean, AllSTD, '-')
%      errorbar(WindBins, LargerMean, LargerSTD, '-')
%      errorbar(WindBins, LessMean, LessSTD, '-')
     ylabel('Power (kW)')
    yyaxis right
    ylabel('Counts')
    bar(WindBins,Num.All,'FaceColor','k')
    bar(WindBins,Num.Less,'FaceColor','b')
    bar(WindBins,Num.Larger,'FaceColor','r')
    alpha(0.05)
    str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.All))),'\d{3}(?=\d)', '$0,')));
    str2 = sprintf('AA > Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Larger))),'\d{3}(?=\d)', '$0,')));
    str3 = sprintf('AA < Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Less))),'\d{3}(?=\d)', '$0,')));
    legend('Reference Curve', str1, str2, str3,'All Data','AA < Hub','AA > Hub')
%     legend(str1, str2, str3,'All Data','AA < Hub','AA > Hub')
    xlabel('u(z = z_h) (m/s)')
    grid on
    titstr = sprintf('Turbine %2.0f',T.TOI);
    title(titstr)
    axprop = gca;
    axprop.YAxis(1).Color = 'k';
    axprop.YAxis(2).Color = '#800000';
    xlim([0 16])
%    ylim([-25 2500])
hold off
% 
% 
% bar(AlphaBins,Mean.Alpha)
% title('Power Output by Power Law Exponent')
% xlabel('alpha (u = u_{ref} \cdot (z/z_{ref})^\alpha')
% ylabel('Power (kW)')

