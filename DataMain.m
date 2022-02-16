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

cd('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Conditional Averages/WindFarmFieldDataProcessing');

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

    [Shear,Veer,Power,I,Time,Nacelle,Indices] = ExtractData(T,data,Indices);    % Performs initial filtering

%% Section 5 - Find adjusted wind angle

    [NormalWind] = NormalInflow(Nacelle,Veer,Shear);                            % Find cosine projection of freestream wind [m/s]

%% Section 6 - Compute Area Average of Incident Wind

    [AAWind,AvgComp] = AreaAveragedVelocity(T,Shear,NormalWind);                % Area-averaged wind speed      [m/s]

%% Section 7 - Fit Power Law to Shear

    [PLFull]   = PowerLawFit(Shear,T,0);                                        % Fit power law to shear profiles

    [PLInflec] = PowerLawFit(Shear,T,1);                                        % Fit power law to shear profiles up to inflection point

%% Section 8 - Fit Power Law to Shear

    [Ekman] = EkmanFit(Shear,T);                                                % Fit power law to shear profiles

%% Section 9 - Conditional Averages and Stats

    [WindBins,Mean,Num,Dist,STD,Sig] = CondAvgs(T,Shear,Power,AvgComp);

%% New Stuff

[Ep] = ErrorTimeSeries(Shear,Power,WindBins,Mean,Time,T);

%% Alpa v. TOD

    [~,~,~,H,MN,~] = datevec(Time);
timeaxis = duration(minutes(0:1:1439),'Format','hh:mm');

    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);
        
            PL.TODvALPHA(index,:) = Indices .* PLFull.alpha;
    
            j     = j + 1;
    
            index = index + 1;
    
        end
    
        i = i + 1;
        j = 0;
    
    end
    
    for i = 1:size(Ep.TOD,1)
    
        PL.TODvALPHAmean(i) = mean(nonzeros(PL.TODvALPHA(i,:)));
    
    end

    figure
        plot(timeaxis,PL.TODvALPHAmean,'LineStyle','none','Marker','.','Color','k')
            xtickformat('hh:mm')
            xlim([timeaxis(1) timeaxis(end)+minutes(1)])
            ylabel('\alpha (-)')
            title('Average Time-of-Day vs. Power Law Alpha (Full Profile)')

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

% BetzWind  = T.CutIn:0.5:9;                                                      % Reference Betz limit plot
% BetzPower = (0.593 * 0.5 * 1.162 * pi * T.R^2 .* BetzWind.^3)/1e3;              % Reference Betz limit plot
% 
% xi = linspace(min(T.RefWind), max(T.RefWind), 100);                             % Evenly-Spaced Interpolation Vector
% yi = interp1(T.RefWind, T.RefCurve, xi, 'spline', 'extrap');
% 
% figure;
% %    plot(BetzWind,BetzPower)%,'color','k','LineStyle',':','Marker','*')
%     hold on
% %   plot(xi,yi,'LineWidth',0.85,'Color','k','LineStyle',':')
%    plot(WindBins,Mean.All,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
%    plot(WindBins,Mean.Larger,'LineWidth',1.5,'Color','#D95319','LineStyle','--')
%    plot(WindBins,Mean.Less,'LineWidth',1.5,'Color','#A2142F','LineStyle','-.')
% %      errorbar(WindBins, AllMean, AllSTD, '-')
% %      errorbar(WindBins, LargerMean, LargerSTD, '-')
% %      errorbar(WindBins, LessMean, LessSTD, '-')
%      ylabel('Power (kW)')
%     yyaxis right
%     ylabel('Counts')
%     bar(WindBins,Num.All,'FaceColor','k')
%     bar(WindBins,Num.Less,'FaceColor','b')
%     bar(WindBins,Num.Larger,'FaceColor','r')
%     alpha(0.05)
%     str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.All))),'\d{3}(?=\d)', '$0,')));
%     str2 = sprintf('AA > Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Larger))),'\d{3}(?=\d)', '$0,')));
%     str3 = sprintf('AA < Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Less))),'\d{3}(?=\d)', '$0,')));
%     legend('Reference Curve', str1, str2, str3,'All Data','AA < Hub','AA > Hub')
%     legend(str1, str2, str3,'All Data','AA < Hub','AA > Hub')
%     xlabel('u(z = z_h) (m/s)')
%     grid on
%     titstr = sprintf('Turbine %2.0f',T.TOI);
%     title(titstr)
%     axprop = gca;
%     axprop.YAxis(1).Color = 'k';
%     axprop.YAxis(2).Color = '#800000';
%     xlim([0 16])
% %    ylim([-25 2500])
% hold off
% 
% 
% bar(AlphaBins,Mean.Alpha)
% title('Power Output by Power Law Exponent')
% xlabel('alpha (u = u_{ref} \cdot (z/z_{ref})^\alpha')
% ylabel('Power (kW)')

