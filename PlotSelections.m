function [] = PlotSelections(P,data,Dist,D,T,WindBins,Mean,STD,Num,AB,AlphaBeta,PDFs,PLFull,PLInflec,Ekman,Ep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if P.WdSpHH == 1

    figure;
        plot(datetime(data.DateTime,'ConvertFrom','datenum'),data.Lidar.H104m.WndSpd,'linestyle','none','marker','.')
            title('LiDAR Hub Height Wind Speed Measurements')
            ylabel('WindSpeed (z = z_h) (m/s)')

end

if P.TI == 1

    figure;
        plot(datetime(data.DateTime,'ConvertFrom','datenum'),data.BHR62.TurbulenceIntensity,'linestyle','none','marker','.')
            title('Turbulence Intensity at Hub Height')
            ylabel('I (%)')
            ylim([0 100])

end

if P.AAPD == 1

    BinWidth = 25;

    for i = 7:19

        h(i-6) = subplot(1,13,i-6);
    
        [counts,edges,~] = histcounts(nonzeros(Dist.All(i,:)),'BinWidth',BinWidth);
    
        barh(edges(2:end),counts)
    
        str = sprintf('N = %3.0f',length(nonzeros(Dist.All(i,:))));
    
            if i < 16
                text(137.5,2100,str,'HorizontalAlignment','center','FontSize',12)
            else
                text(137.5,100,str,'HorizontalAlignment','center','FontSize',12)
            end
        
            xlim([0 275])
            ylim([0 2200])
            xtickangle(45)
        
            if i == 7
                ylabel('Power (kW)')
            end
        
            if i == 13
                xlabel('Counts')
                title('Power Distributions of Bins in Region 2 - All Average')
            end
        
            if i ~= 7
                set(gca,'ytick',[])
            end

    end
    
    for i = 1:13
        set(h(i), 'Position', [(1/14.5*(i)) 0.1 1/17 0.81])
    end

end

if P.GPD == 1

    BinWidth = 25;

    for i = 7:19
    
        h(i-6) = subplot(1,13,i-6);
        
        [counts,edges,~] = histcounts(nonzeros(Dist.Larger(i,:)),'BinWidth',BinWidth);
        
        barh(edges(2:end),counts)
        
        str = sprintf('N = %3.0f',length(nonzeros(Dist.Larger(i,:))));
        
            if i < 16
                text(25,2100,str,'HorizontalAlignment','center','FontSize',12)
            else
                text(25,100,str,'HorizontalAlignment','center','FontSize',12)
            end
            
            xlim([0 50])
            ylim([0 2200])
            xtickangle(45)
            
            if i == 7
                ylabel('Power (kW)')
            end
            
            if i == 13
                xlabel('Counts')
                title('Power Distributions of Bins in Region 2 Where Area-Average Velocity is Greater than Hub Height Velocity')
            end
            
            if i ~= 7
                set(gca,'ytick',[])
            end
    
    end

    for i = 1:13
        set(h(i), 'Position', [(1/14.5*(i)) 0.1 1/17 0.81])
    end

end

if P.LPD == 1

    BinWidth = 25;

    for i = 7:19
    
        h(i-6) = subplot(1,13,i-6);
        
        [counts,edges,~] = histcounts(nonzeros(Dist.Less(i,:)),'BinWidth',BinWidth);
        barh(edges(2:end),counts)
        
        str = sprintf('N = %3.0f',length(nonzeros(Dist.Less(i,:))));
        
            if i < 16
                text(112.5,2100,str,'HorizontalAlignment','center','FontSize',12)
            else
                text(112.5,100,str,'HorizontalAlignment','center','FontSize',12)
            end
            
            xlim([0 225])
            ylim([0 2200])
            xtickangle(45)
            
            if i == 7
                ylabel('Power (kW)')
            end
            
            if i == 13
                xlabel('Counts')
                title('Power Distributions of Bins in Region 2 Where Area-Average Velocity is Less than Hub Height Velocity')
            end
            
            if i ~= 7
                set(gca,'ytick',[])
            end
    
    end
    
    for i = 1:13
        set(h(i), 'Position', [(1/14.5*(i)) 0.1 1/17 0.81])
    end

end

if P.WindRoseFull == 1

    WindSpeed     = data.Lidar.H104m.WndSpd;
    WindDirection = data.Lidar.H104m.WndDir;
    
    figure;
    
    bins = fliplr([2:1.5:max(WindSpeed)+1]);                                   % Wind speed bins
    
    pax       = polaraxes;                                                      % Polar axes handle
    
    colors    = parula(length(bins)+2);                                         % Set colormap
    
        for i = 1:length(bins)-1
        
            LegendEntries = strcat(sprintf('%5.1f',bins(i+1)),' - ',...         % Make legend entry strings
                sprintf('%5.1f',bins(i)));
        
            polarhistogram(deg2rad(WindDirection(WindSpeed<bins(i))),...        % Add entries to legend
                deg2rad(0:10:360),'FaceColor',colors(i+1,:),'FaceAlpha',1,...
                'displayname',LegendEntries)
            hold on
        
        end
    
    pax.ThetaDir          = 'clockwise';                                        % Orientation
    pax.ThetaZeroLocation = 'top';                                              % North pointing
    
    Leg = legend('Show');
    title(Leg,'Wind speed [m s^{-1}]')                                          % Set legend title
    legend('AutoUpdate','off')                                                  % Prevent subsequent plot from being added to legend
    box(Leg,'off')                                                              % Remove outline from legend
    set(gca,'rticklabel',[])                                                    % Remove default radial markers
    set(gca,'FontSize',16)
    
    TotalPoints = sum(WindSpeed>0);
    axlims    = rlim;                                                           % Find radial limits
    positions = linspace(0,axlims(2),5);                                        % Define new radial marker spacing
    rticks(positions);                                                          % Set new radial marker positions
    angle     = deg2rad(95);                                                    % Angle of radial marker vector
    
    h = polarhistogram(deg2rad(WindDirection(WindDirection>=330 | ...           % Add outline to wind of interest
        WindDirection<=50)),deg2rad(0:10:360),'LineWidth',2.5);
    h.DisplayStyle = 'stairs';                                                  % Set outline style
    h.EdgeColor    = [0 0 0];                                                   % Set color to black                                                   
    
    text(angle, positions(2), sprintf('%2.0f%%',100*positions(2)/TotalPoints),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(3), sprintf('%2.0f%%',100*positions(3)/TotalPoints),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(4), sprintf('%2.0f%%',100*positions(4)/TotalPoints),'HorizontalAlignment','center','FontSize',12);
    
    hold off;

end

if P.WindRoseAnalysis == 1

    WindSpeed     = D.Shear(T.HubRow,:);
    WindDirection = D.Veer(T.HubRow,:);
    
    figure;
    
    bins = fliplr([2:2:max(WindSpeed)-1,20]);                                   % Wind speed bins
    
    pax       = polaraxes;                                                      % Polar axes handle
    
    colors    = parula(length(bins)+2);                                         % Set colormap
    
        for i = 1:length(bins)-1
        
            LegendEntries = strcat(sprintf('%5.1f',bins(i+1)),' - ',...         % Make legend entry strings
                sprintf('%5.1f',bins(i)));
        
            polarhistogram(deg2rad(WindDirection(WindSpeed<bins(i))),...        % Add entries to legend
                deg2rad(0:10:360),'FaceColor',colors(i+1,:),'FaceAlpha',1,...
                'displayname',LegendEntries)
            hold on
        
        end
    
    pax.ThetaDir          = 'clockwise';                                        % Orientation
    pax.ThetaZeroLocation = 'top';                                              % North pointing
    
    Leg = legend('Show');
    title(Leg,'Wind speed [m s^{-1}]')                                          % Set legend title
    legend('AutoUpdate','off')                                                  % Prevent subsequent plot from being added to legend
    box(Leg,'off')                                                              % Remove outline from legend
    set(gca,'rticklabel',[])                                                    % Remove default radial markers
    set(gca,'FontSize',16)
    
    TotalPoints = sum(WindSpeed>0);
    axlims    = rlim;                                                           % Find radial limits
    positions = linspace(0,axlims(2),5);                                        % Define new radial marker spacing
    rticks(positions);                                                          % Set new radial marker positions
    angle     = deg2rad(70);                                                    % Angle of radial marker vector                                                 
    
    text(angle, positions(2), sprintf('%2.0f%%',100*positions(2)/TotalPoints),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(3), sprintf('%2.0f%%',100*positions(3)/TotalPoints),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(4), sprintf('%2.0f%%',100*positions(4)/TotalPoints),'HorizontalAlignment','center','FontSize',12);
    
    hold off;

end

if P.HistoPowerAA == 1
    
    figure;
        hold on
        
        plot(WindBins,Mean.All,'LineWidth',1.5,'Color','k','LineStyle','-')
        
        M = Mean.All;
        S = STD.All;
        
        M(isnan(M)) = 0;
        S(isnan(S)) = 0;
        
        Top = M + S;
        Bot = M - S;
        
        patch([WindBins fliplr(WindBins)], [Top fliplr(Bot)],'k','facealpha',0.1,'edgecolor','None')
        
        ylim([0 1.2])
        xlim([0.5 14])
        ylabel('P/P_{rated}')
        yyaxis right
        ylabel('Counts')
        bar(WindBins,Num.All,'FaceColor','k')
        
        alpha(0.05)

        xlabel('u(z = z_h) (m/s)')
        
        titstr = sprintf('Turbine %2.0f',T.TOI);
        title(titstr)
        axprop = gca;
        axprop.YAxis(1).Color = 'k';
        axprop.YAxis(2).Color = '#800000';
        hold off

end

if P.HistoPowerAALG == 1
    
    figure;
        hold on
        plot(WindBins,Mean.All,'LineWidth',1.5,'Color','k','LineStyle','-')
        plot(WindBins,Mean.Larger,'LineWidth',1.5,'Color','k','LineStyle','--')
        plot(WindBins,Mean.Less,'LineWidth',1.5,'Color','k','LineStyle','-.')

        ylim([0 1.2])
        xlim([0.5 14])
        ylabel('P/P_{rated}')
        yyaxis right
        ylabel('Counts')
        bar(WindBins,Num.All,'FaceColor','k')

        alpha(0.05)
        str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.All))),'\d{3}(?=\d)', '$0,')));
        str2 = sprintf('AA > Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Larger))),'\d{3}(?=\d)', '$0,')));
        str3 = sprintf('AA < Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Less))),'\d{3}(?=\d)', '$0,')));

        legend(str1, str2, str3,'')
        xlabel('u(z = z_h) (m/s)')

        titstr = sprintf('Turbine %2.0f',T.TOI);
        title(titstr)
        axprop = gca;
        axprop.YAxis(1).Color = 'k';
        axprop.YAxis(2).Color = '#800000';
        xlim([0.5 14])

        legend({'Average Power Curve','$\int_A u(z) > u(z=z_h)$','$\int_A u(z) < u(z=z_h)$'},'Interpreter','latex');

    hold off

end

if P.AlphaBetaFull == 1

    figure;
        imagesc(AB.SSrange,AB.DSrange,AlphaBeta.Power);
        axis xy                                                                 % Flip image
        colormap([1 1 1; parula(256)])                                          % Bi-tone colomap
        colorbar
    
        xlabel('Speed Shear (\alpha)')
        ylabel('Direction Shear (\circ m^{-1})')
        ylabel(colorbar('eastoutside'),'Normalized Power (-)')
        xlim([AB.xlow AB.xhigh])
        caxis([0.75 1.1])
        hold on
    
        for i = 1:length(AB.SSrange)                                            % Vertical brid lines
            xline(AB.SSrange(i)-AB.s/2)
        end
    
        for i = 1:length(AB.DSrange)                                            % Horizontal brid lines
            yline(AB.DSrange(i)-AB.s/2)
        end
    
        x = [AB.SSrange 0.8];                                                   % Add threshold line
        y = 2/3*x - 0.1;
        plot(x,y,'color','r','LineWidth',2)
    
        hold off

end

if P.AlphaBetaMono == 1

    figure;
        imagesc(AB.SSrange,AB.DSrange,AlphaBeta.Power);
        axis xy                                                                 % Flip image
        colormap([1 1 1; gray(256)])                                            % Bi-tone colomap
        colorbar
    
        xlabel('Speed Shear (\alpha)')
        ylabel('Direction Shear (\circ m^{-1})')
        ylabel(colorbar('eastoutside'),'Normalized Power (-)')
        xlim([AB.xlow AB.xhigh])
        hold on
    
        for i = 1:length(AB.SSrange)                                            % Vertical brid lines
            xline(AB.SSrange(i)-AB.s/2)
        end
    
        for i = 1:length(AB.DSrange)                                            % Horizontal brid lines
            yline(AB.DSrange(i)-AB.s/2)
        end
    
        x = [AB.SSrange 0.8];                                                   % Add threshold line
        y = 2/3*x - 0.1;
        plot(x,y,'color','r','LineWidth',2)
    
        hold off

end

if P.AlphaBetaLH == 1

    low  = min(AlphaBeta.Power,[],'all');                                       % Find minimum
    high = max(AlphaBeta.Power,[],'all');                                       % Find maximum
    
    AlphaBeta.PowerLH(AlphaBeta.Power<1 & ~isnan(AlphaBeta.Power)) = low;       % Convert all values >1 to maximum
    AlphaBeta.PowerLH(AlphaBeta.Power>=1) = high;                               % Convert all values <1 to minimum
    
        if abs(high-1) > abs(1-low)                                             % Calibrate colormap values
            AlphaBeta.PowerLH(end,end) = 1-abs(high-1);
        else
            AlphaBeta.PowerLH(end,end) = 1+abs(low-1);
        end
    
    figure;
%         h = pcolor(AB.SSrange,AB.DSrange,AlphaBeta.PowerLH);
        h = imagesc(AB.SSrange,AB.DSrange,AlphaBeta.PowerLH);
        axis xy                                                                 % Flip image
        colormap([0.66 0.66 0.66; .33 .33 .33])                                 % Bi-tone colomap
        colorbar
        set(h,'alphadata',~isnan(AlphaBeta.PowerLH))                            % Set NaNs to white
    
        xlabel('Speed Shear (\alpha)')
        ylabel('Direction Shear (\circ m^{-1})')
        ylabel(colorbar('eastoutside'),'Normalized Power (-)')
        xlim([AB.xlow AB.xhigh])
        hold on
    
            for i = 1:length(AB.SSrange)                                        % Vertical brid lines
                xline(AB.SSrange(i)-AB.s/2)
            end
        
            for i = 1:length(AB.DSrange)                                        % Horizontal brid lines
                yline(AB.DSrange(i)-AB.s/2)
            end
    
        x = [AB.SSrange 0.8];                                                   % Add threshold line
        y = 2/3*x - 0.1;
        plot(x,y,'color','r','LineWidth',2)
    
        hold off

end

if P.DSprob == 1

    figure;                                                                 % 100 bins
        [PDFs.DSProb1,PDFs.DSEdges1] = histcounts(PDFs.DSrate,110,...
            'Normalization','probability','BinLimits',[-0.6 0.8]);
        plot(PDFs.DSEdges1(1:end-1),PDFs.DSProb1,'Color','k','LineWidth',1.5);
        title('Direction Shear Probability Distribution')
        xlabel('Direction Shear (\circ m^{-1})')
        ylabel('Probability of Occurrence')
        xlim([-0.6 0.8])

end

if P.SSFull == 1

    figure;                                                                 % Full profile fit
        [PDFs.SSProb1,PDFs.SSEdges1] = histcounts(PLFull.alpha,...
            'Normalization','probability','BinLimits',[-0.75 1.5]);
        plot(PDFs.SSEdges1(1:end-1),PDFs.SSProb1,'Color','k','LineWidth',1.5);
        title('Speed Shear Probability Distribution - Full Fit')
        xlabel('Speed Shear (\alpha)')
        ylabel('Probability of Occurrence')
        xlim([-0.75 1.5])

end

if P.SSInflec == 1

    figure;                                                                 % Inflection profile fit
        [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(PLInflec.alpha,...
            'Normalization','probability','BinLimits',[-0.5 1.5]);
        plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','k','LineWidth',1.5);
        title('Speed Shear Probability Distribution - Inflection Fit')
        xlabel('Speed Shear (\alpha)')
        ylabel('Probability of Occurrence')
        xlim([-0.35 1.4])

end

if P.EkmanProb == 1

    figure;                                                                 % Inflection profile fit
        [EkProb,EkEdges] = histcounts(Ekman.K,...
            'Normalization','probability','BinLimits',[-0.5 1.5]);
        plot(EkEdges(1:end-1),EkProb,'Color','k','LineWidth',1.5);
        title('Ekman Parameter Probability Distribution')
        xlabel('Eddy Diffusivity (m^2 s ^{-1})')
        ylabel('Probability of Occurrence')

end

if P.SSTOD == 1

    figure;
        plot(PDFs.TODAxis,PDFs.DSTOD,'LineStyle','none','Marker','.','Color','k')
        title('Direction Shear Evolution by Time of day')
        xlabel('Time of Day')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylabel('Direction Shear (\circ m^{-1})')

    figure;
        plot(PDFs.TODAxis,PDFs.SSTODFull,'LineStyle','none','Marker','.','Color','k')
        title('Speed Shear Evolution by Time of day - Full Profile Fit')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylabel('Speed Shear (\alpha)')

    figure;
        plot(PDFs.TODAxis,PDFs.SSTODInflec,'LineStyle','none','Marker','.','Color','k')
        title('Speed Shear Evolution by Time of day - Partial Profile Fit')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylabel('Speed Shear (\alpha)')

end


if P.SSTODAVG == 1

    figure;
        plot(PDFs.TODAvgAxis,PDFs.DSTODAvg,'LineStyle','-','LineWidth',1.5,'Color','k')
        hold on
        title('Speed and Direction Shear Evolution by Time of day')
        xlabel('Time of Day')
        ylabel('Direction Shear (\circ m^{-1})')
        ylim([-0.05 0.5])

        yyaxis right
        plot(PDFs.TODAvgAxis,PDFs.SSTODAvgFull,'LineStyle','-','LineWidth',1.5)
        plot(PDFs.TODAvgAxis,PDFs.SSTODAvgInflec,'LineStyle','--','LineWidth',1.5)
        ylabel('Speed Shear (\alpha)')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylim([-0.05 0.5])
        legend('Direction Shear','Speed Shear - Full Alpha','Speed Shear - Partial Alpha')
        hold off

end

if P.EkTOD == 1

    figure;
        plot(PDFs.TODAxis,PDFs.EkTOD,'LineStyle','none','Marker','.','Color','k')
        title('Ekman Parameter Evolution by Time of day')
        xlabel('Time of Day')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylabel('Eddy Diffusivity (m^2 s ^{-1})')
        ylim([0 0.5])

end

if P.EkTODAVG == 1

    figure;
        plot(PDFs.TODAvgAxis,PDFs.DSTODAvg,'LineStyle','-','LineWidth',1.5,'Color','k')
        hold on
        title('[Average] Ekman Parameter and Direction Shear Evolution by Time of day')
        xlabel('Time of Day')
        ylabel('Direction Shear (\circ m^{-1})')
        ylim([-0.05 0.5])

        yyaxis right
        plot(PDFs.TODAvgAxis,PDFs.EkTODAvg,'LineStyle','-','LineWidth',1.5)
        ylabel('Eddy Diffusivity (m^2 s ^{-1})')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
%         ylim([-0.05 0.5])
        hold off

    figure;
        plot(PDFs.TODAvgAxis,PDFs.SSTODAvgFull,'LineStyle','-','LineWidth',1.5,'Color','k')
        ylabel('Speed Shear (\alpha)')
        hold on
        title('[Average] Ekman Parameter and Speed Shear Evolution by Time of day (Full Alpha)')
        xlabel('Time of Day')
        ylabel('Speed Shear (\alpha)')
        ylim([-0.05 0.5])

        yyaxis right
        plot(PDFs.TODAvgAxis,PDFs.EkTODAvg,'LineStyle','-','LineWidth',1.5)
        ylabel('Eddy Diffusivity (m^2 s ^{-1})')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
%         ylim([-0.05 0.5])
        hold off

end

if P.EkTODMED == 1

    figure;
        plot(PDFs.TODAvgAxis,PDFs.DSTODAvg,'LineStyle','-','LineWidth',1.5,'Color','k')
        hold on
        title('[Median] Ekman Parameter and Direction Shear Evolution by Time of day')
        xlabel('Time of Day')
        ylabel('Direction Shear (\circ m^{-1})')
        ylim([-0.05 0.5])

        yyaxis right
        plot(PDFs.TODAvgAxis,PDFs.EkTODMed,'LineStyle','-','LineWidth',1.5)
        ylabel('Eddy Diffusivity (m^2 s ^{-1})')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylim([0 0.35])
        hold off

    figure;
        plot(PDFs.TODAvgAxis,PDFs.SSTODAvgFull,'LineStyle','-','LineWidth',1.5,'Color','k')
        ylabel('Speed Shear (\alpha)')
        hold on
        title('[Median] Ekman Parameter and Speed Shear Evolution by Time of day (Full Alpha)')
        xlabel('Time of Day')
        ylabel('Speed Shear (\alpha)')
        ylim([-0.05 0.5])

        yyaxis right
        plot(PDFs.TODAvgAxis,PDFs.EkTODMed,'LineStyle','-','LineWidth',1.5)
        ylabel('Eddy Diffusivity (m^2 s ^{-1})')
        xlim([PDFs.TODAxis(1) PDFs.TODAxis(end)+minutes(1)])
        ylim([-0.015 0.25])
        hold off

end

if P.EpTODMag == 1

    figure;
        scatter(Ep.TODAxis,Ep.TODMagMean,'Marker','.','MarkerEdgeColor','k')
        hold on;

        y = mean(reshape(Ep.TODMagMean,60,[]));
        x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');
    
        plot(x,y,'Color','r','LineWidth',1.5)

        title('mean $\left( \frac{P(t) - \overline{P}}{\overline{P}}\right)$','Interpreter','latex')
        xlabel('Time of Day')
        ylabel('\epsilon_P (-)')
        legend('1 - Minute Average', '1 - Hour Average')

end

if P.EpTODDiff == 1

    figure;
        scatter(Ep.TODAxis,Ep.TODDiffMean,'Marker','.','MarkerEdgeColor','k')
        hold on;

        y = mean(reshape(Ep.TODDiffMean,60,[]));
        x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');
    
        plot(x,y,'Color','r','LineWidth',1.5)

        yline(0,'LineWidth',1.5)

        title('$\frac{|P(t) - \overline{P}|}{\overline{P}}$','Interpreter','latex')
        xlabel('Time of Day')
        ylabel('\epsilon_P (-)')
        legend('1 - Minute Average', '1 - Hour Average')

end

if P.EpTODAvg == 1

    figure;
        y = mean(reshape(Ep.TODmean,10,[]));
        x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');
    
        plot(x,y,'Color','k','LineWidth',1.5)

end

if P.EpTODHisto == 1

% [1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24] - Ep.dist index
% [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1]  - Cooresponding hour in row

    % -------------- Hourly Histogram: Day Time --------------
    
        DayFig   = figure;
        BinWidth = 0.02;
    
            HourIndices = [6 7 8 9  10 11 12 13 14 15 16 17;
                           7 8 9 10 11 12 13 14 15 16 17 18];
        
            for i = 1:12
    
                h(i) = subplot(1,12,i);
    
                [counts,edges,~] = histcounts(nonzeros(Ep.dist(HourIndices(2,i),:)),'BinWidth',BinWidth);
                barh(edges(2:end),counts)
    
                str = sprintf('N = %3.0f',length(nonzeros(Ep.dist(HourIndices(2,i),:))));
                text(237.5,1.1,str,'HorizontalAlignment','center','FontSize',12)
                xlim([0 475])
                ylim([0 1.2])
                xtickangle(45)
                str = sprintf('%02.0f:00',HourIndices(1,i));
                title(str)
        
                    if i ~= 1
        
                        set(gca,'ytick',[])
        
                    end
            end
            
            AxHandle=axes(DayFig,'visible','off'); 
                AxHandle.Title.Visible   ='on';
                AxHandle.XLabel.Visible  ='on';
                AxHandle.XLabel.Position = [0.5 0];
                xlabel(AxHandle,'Counts');
                AxHandle.YLabel.Visible  ='on';
                AxHandle.YLabel.Position = [0 0.5];
                ylabel(AxHandle,'\epsilon_P (-)');
                title(AxHandle,'\epsilon_P Distribution by Hour');
            
            for i = 1:12
        
                set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])
        
            end
    
    % -------------- Hourly Histogram: Night Time --------------
    
        NightFig   = figure;
        BinWidth = 0.02;
    
        HourIndices = [18 19 20 21 22 23 0 1 2 3 4 5;
                       19 20 21 22 23 24 1 2 3 4 5 6];
    
        for i = 1:12
    
            h(i) = subplot(1,12,i);
    
            [counts,edges,~] = histcounts(nonzeros(Ep.dist(HourIndices(2,i),:)),'BinWidth',BinWidth);
            barh(edges(2:end),counts)
    
            str = sprintf('N = %3.0f',length(nonzeros(Ep.dist(HourIndices(2,i),:))));
            text(237.5,1.1,str,'HorizontalAlignment','center','FontSize',12)
            xlim([0 475])
            ylim([0 1.2])
            xtickangle(45)
            str = sprintf('%02.0f:00',HourIndices(1,i));
            title(str)
    
                if i ~= 1
    
                    set(gca,'ytick',[])
    
                end
    
        end
        
        AxHandle=axes(NightFig,'visible','off'); 
            AxHandle.Title.Visible   ='on';
            AxHandle.XLabel.Visible  ='on';
            AxHandle.XLabel.Position = [0.5 0];
            xlabel(AxHandle,'Counts');
            AxHandle.YLabel.Visible  ='on';
            AxHandle.YLabel.Position = [0 0.5];
            ylabel(AxHandle,'\epsilon_P (-)');
            title(AxHandle,'\epsilon_P Distribution by Hour');
        
        for i = 1:12
    
            set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])
    
        end
    
end

if P.InflecTODHisto == 1

    [~,~,~,H,~,~] = datevec(D.Time);

    HourlyDist = zeros(24,length(PLInflec.InflecHeight));

    for i = 1:24
    
        Indices = (H == i-1);
    
        HourlyDist(i,:) = Indices .* PLInflec.InflecHeight;
    
    end

    % -------------- Hourly Histogram: Night Time --------------

    NightFig = figure;

    HourIndices = [18 19 20 21 22 23 0 1 2 3 4 5;
                   19 20 21 22 23 24 1 2 3 4 5 6];

    for i = 1:12

        h(i) = subplot(1,12,i);

        [counts,edges,~] = histcounts(nonzeros(HourlyDist(HourIndices(2,i),:)),'binedges',T.Heights);
        barh(edges(2:end),counts)

        str = sprintf('N = %3.0f',length(nonzeros(HourlyDist(HourIndices(2,i),:))));
%         text(237.5,1.1,str,'HorizontalAlignment','center','FontSize',12)
        xlim([0 200])
        ylim([45 175])
        xtickangle(45)
        str = sprintf('%02.0f:00',HourIndices(1,i));
        title(str)

            if i ~= 1

                set(gca,'ytick',[])

            end

    end
    
    AxHandle=axes(NightFig,'visible','off'); 
        AxHandle.Title.Visible   ='on';
        AxHandle.XLabel.Visible  ='on';
        AxHandle.XLabel.Position = [0.5 0];
        xlabel(AxHandle,'Counts');
        AxHandle.YLabel.Visible  ='on';
        AxHandle.YLabel.Position = [0 0.5];
        ylabel(AxHandle,'z (m)');
        title(AxHandle,'Inflection Point Height Distribution by Hour');
    
    for i = 1:12

        set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])

    end

    % -------------- Hourly Histogram: Day Time --------------

    DayFig = figure;

    HourIndices = [6 7 8 9  10 11 12 13 14 15 16 17;
                   7 8 9 10 11 12 13 14 15 16 17 18];

    for i = 1:12

        h(i) = subplot(1,12,i);

        [counts,edges,~] = histcounts(nonzeros(HourlyDist(HourIndices(2,i),:)),'binedges',T.Heights);
        barh(edges(2:end),counts)

        str = sprintf('N = %3.0f',length(nonzeros(HourlyDist(HourIndices(2,i),:))));
%         text(237.5,1.1,str,'HorizontalAlignment','center','FontSize',12)
        xlim([0 200])
        ylim([45 175])
        xtickangle(45)
        str = sprintf('%02.0f:00',HourIndices(1,i));
        title(str)

            if i ~= 1

                set(gca,'ytick',[])

            end

    end
    
    AxHandle=axes(DayFig,'visible','off'); 
        AxHandle.Title.Visible   ='on';
        AxHandle.XLabel.Visible  ='on';
        AxHandle.XLabel.Position = [0.5 0];
        xlabel(AxHandle,'Counts');
        AxHandle.YLabel.Visible  ='on';
        AxHandle.YLabel.Position = [0 0.5];
        ylabel(AxHandle,'z (m)');
        title(AxHandle,'Inflection Point Height Distribution by Hour');
    
    for i = 1:12

        set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])

    end


end

if P.InflecTOD == 1

    Heights = PLInflec.InflecHeight;

    Heights(isnan(Heights)) = 0;

    [~,~,~,H,~,~] = datevec(D.Time);

    MinuteMean = zeros(1440,1);

    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);
        
            MinuteMean(index) = mean(nonzeros(Indices .* Heights));
    
            j     = j + 1;
    
            index = index + 1;
    
        end
    
        i = i + 1;
        j = 0;
    
    end

    x = duration(minutes(linspace(0,1439,length(MinuteMean))),'Format','hh:mm');

    figure;
        scatter(x,MinuteMean,'Marker','.','MarkerEdgeColor','k')
        hold on;

        y = mean(reshape(MinuteMean,60,[]));
        x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');
    
        plot(x,y,'Color','r','LineWidth',1.5)

        title('Inflection Height Evolution by Time of Day')
        xlabel('Time of Day')
        ylabel('z (m)')
        legend('1 - Minute Average', 'Hourly Average')

end

if P.FitErrorTOD == 1

    figure;
FullAlphaIndices   = ~isnan(PLFull.alpha);
InflecAlphaIndices = ~isnan(PLInflec.alpha);
EkmanIndices       = ~isnan(Ekman.K);

FullIndices = FullAlphaIndices .* InflecAlphaIndices .* EkmanIndices;

FullIndices(FullIndices == 0) = NaN;

y = FullIndices .* PLFull.R;
x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');

plot(x,y)
hold on

y = FullIndices .* PLInflec.R;
x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');

plot(x,y)

y = FullIndices .* Ekman.R;
x = duration(minutes(linspace(0,1439,length(y))),'Format','hh:mm');

plot(x,y)

% ylim([-60 1])
ylabel('Coefficient of Determination - R^2')
legend('Power Law - Full Profile','Power Law - Partial Profile','Ekman Fit')

end

figure;
    plot(datetime(D.Time,'ConvertFrom','datenum'),Ekman.R)

    figure;                                                                 % Inflection profile fit
        [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(Ekman.R,...
            'Normalization','probability','BinLimits',[-0.5 1.5]);
        plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','k','LineWidth',1.5);

        hold on

        [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(PLFull.R,...
            'Normalization','probability','BinLimits',[-0.5 1.5]);
        plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','r','LineWidth',1.5);


        [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(PLInflec.R,...
            'Normalization','probability','BinLimits',[-0.5 1.5]);
        plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','b','LineWidth',1.5);

        legend('Ekman','Full','Inflec')


        title('Speed Shear Probability Distribution - Inflection Fit')
        xlabel('Speed Shear (\alpha)')
        ylabel('Probability of Occurrence')
        xlim([-0.35 1.4])

% if P.HistoPower == 1
% 
%     BetzWind  = T.CutIn:0.5:9;                                                  % Reference Betz limit plot
%     BetzPower = (0.593 * 0.5 * 1.162 * pi * T.R^2 .* BetzWind.^3)/1e3;          % Reference Betz limit plot
%     
%     xi = linspace(min(T.RefWind), max(T.RefWind), 100);                         % Evenly-Spaced Interpolation Vector
%     yi = interp1(T.RefWind, T.RefCurve, xi, 'spline', 'extrap');
%     
%     figure;
% %         zWind,BetzPower)%,'color','k','LineStyle',':','Marker','*')
%         hold on
% %         plot(xi,yi,'LineWidth',0.85,'Color','k','LineStyle',':')
% %         plot(WindBins,Mean.All,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
% %         plot(WindBins,Mean.Larger,'LineWidth',1.5,'Color','#D95319','LineStyle','--')
% %         plot(WindBins,Mean.Less,'LineWidth',1.5,'Color','#A2142F','LineStyle','-.')
% %         color = jet(3);
% color = [0 0 0; 0 0 0; 0 0 0];
%         plot(WindBins,Mean.All,'LineWidth',1.5,'Color',color(1,:),'LineStyle','-')
%         plot(WindBins,Mean.Larger,'LineWidth',1.5,'Color',color(2,:),'LineStyle','--')
%         plot(WindBins,Mean.Less,'LineWidth',1.5,'Color',color(3,:),'LineStyle','-.')
% 
% M = Mean.All;
% S = STD.All;
% 
% M(isnan(M)) = 0;
% S(isnan(S)) = 0;
% 
% Top = M + S;
% Bot = M - S;
% 
% patch([WindBins fliplr(WindBins)], [Top fliplr(Bot)],'k','facealpha',0.1,'edgecolor','None')
% 
% % X=[WindBins,fliplr(WindBins)];                %#create continuous x value array for plotting
% % Y=[Mean.All - STD.All,fliplr(Mean.All + STD.All)];              %#create y values for out and then back
% % fill(X,Y,'blue'); 
% 
%     %      errorbar(WindBins, AllMean, AllSTD, '-')
%     %      errorbar(WindBins, LargerMean, LargerSTD, '-')
%     %      errorbar(WindBins, LessMean, LessSTD, '-')
%     ylim([0 1.2])
%     ylabel('P/P_{rated}')
%     yyaxis right
%     ylabel('Counts')
%     bar(WindBins,Num.All,'FaceColor',color(1,:))
% %         bar(WindBins,Num.Less,'FaceColor',color(2,:))
% %         bar(WindBins,Num.Larger,'FaceColor',color(3,:))
%         alpha(0.05)
%         str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.All))),'\d{3}(?=\d)', '$0,')));
%         str2 = sprintf('AA > Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Larger))),'\d{3}(?=\d)', '$0,')));
%         str3 = sprintf('AA < Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(Num.Less))),'\d{3}(?=\d)', '$0,')));
% %         legend('Reference Curve', str1, str2, str3,'All Data','AA < Hub','AA > Hub')
%         legend(str1, str2, str3,'')
%         xlabel('u(z = z_h) (m/s)')
% %         grid on
%         titstr = sprintf('Turbine %2.0f',T.TOI);
%         title(titstr)
%         axprop = gca;
%         axprop.YAxis(1).Color = 'k';
%         axprop.YAxis(2).Color = '#800000';
%         xlim([0.5 14])
%     %    ylim([-25 2500])
%     hold off
%     % 
%     % 
%     % bar(AlphaBins,Mean.Alpha)
%     % title('Power Output by Power Law Exponent')
%     % xlabel('alpha (u = u_{ref} \cdot (z/z_{ref})^\alpha')
%     % ylabel('Power (kW)')
% 
% end






end