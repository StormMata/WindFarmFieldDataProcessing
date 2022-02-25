function [] = PlotSelections(P,data,Dist,D,T,WindBins,Mean,STD,Num,AB,AlphaBeta)
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
    
    axlims    = rlim;                                                           % Find radial limits
    positions = linspace(0,axlims(2),5);                                        % Define new radial marker spacing
    rticks(positions);                                                          % Set new radial marker positions
    angle     = deg2rad(95);                                                    % Angle of radial marker vector
    
    h = polarhistogram(deg2rad(WindDirection(WindDirection>=330 | ...           % Add outline to wind of interest
        WindDirection<=50)),deg2rad(0:10:360),'LineWidth',2.5);
    h.DisplayStyle = 'stairs';                                                  % Set outline style
    h.EdgeColor    = [0 0 0];                                                   % Set color to black                                                   
    
    text(angle, positions(2), sprintf('%2.0f%%',100*positions(2)/(axlims(2))),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(3), sprintf('%2.0f%%',100*positions(3)/(axlims(2))),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(4), sprintf('%2.0f%%',100*positions(4)/(axlims(2))),'HorizontalAlignment','center','FontSize',12);
    
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
    
    axlims    = rlim;                                                           % Find radial limits
    positions = linspace(0,axlims(2),5);                                        % Define new radial marker spacing
    rticks(positions);                                                          % Set new radial marker positions
    angle     = deg2rad(70);                                                    % Angle of radial marker vector                                                 
    
    text(angle, positions(2), sprintf('%2.0f%%',100*positions(2)/(axlims(2))),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(3), sprintf('%2.0f%%',100*positions(3)/(axlims(2))),'HorizontalAlignment','center','FontSize',12);
    text(angle, positions(4), sprintf('%2.0f%%',100*positions(4)/(axlims(2))),'HorizontalAlignment','center','FontSize',12);
    
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