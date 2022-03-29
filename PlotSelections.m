function [] = PlotSelections(P,data,Dist,D,T,WindBins,Mean,STD,Num,AB,AlphaBeta,KB,KmBeta,PDFs,PLFull,PLInflec,Ekman,Ep,WndFamMap)
%PlotSelections produces the selected plots specified in DataMain
%   Plot Selections takes the selections and data from DataMain and
%   produces the specifed plots.

if P.TurbineSchema == 1

    % Load chord length
    
    addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/BEMT/BEMT');
    
    [~,~,~,~,c] = LoadGeometry('Airfoil_17.txt');
    
    % Plot
    
    figure;
    hold on
    
    % -------------------- Front View --------------------
    
        % Blade 1
        
            angle = 90;
            
            Ys = c';
            Xs = linspace(0,T.R,length(Ys));
            
            a_rad=((angle*pi)./180);
            [THETA,R] = cart2pol(Xs,Ys);        %Convert to polar coordinates
            THETA = THETA+a_rad;                %Add a_rad to theta
            [xr,yr] = pol2cart(THETA,R);        %Convert back to Cartesian coordinates
            
            plot(xr+T.HubR/2, T.Hub+yr+2,'k');
            plot([T.HubR/2 T.HubR/2+xr(end)],[T.Hub+2 T.Hub+yr(end)+2],'k')
        
        % Blade 2
        
            angle = 210;
            
            Ys = c';
            Xs = linspace(0,T.R,length(Ys));
            
            a_rad=((angle*pi)./180);
            [THETA,R] = cart2pol(Xs,Ys);        %Convert to polar coordinates
            THETA = THETA+a_rad;                %Add a_rad to theta
            [xr,yr] = pol2cart(THETA,R);        %Convert back to Cartesian coordinates
            
            plot(xr-T.HubR/2, T.Hub+yr,'k');
            plot([-T.HubR/2 -T.HubR/2+xr(end)],[T.Hub T.Hub+yr(end)],'k')
        
        % Blade 3
        
            angle = 330;
            
            Ys = c';
            Xs = linspace(0,T.R,length(Ys));
            
            a_rad=((angle*pi)./180);
            [THETA,R] = cart2pol(Xs,Ys);        %Convert to polar coordinates
            THETA = THETA+a_rad;                %Add a_rad to theta
            [xr,yr] = pol2cart(THETA,R);        %Convert back to Cartesian coordinates
            
            plot(xr, (T.Hub-T.HubR/2)+yr,'k');
            plot([0 xr(end)],[(T.Hub-T.HubR/2) (T.Hub-T.HubR/2)+yr(end)],'k')
        
        % Tower
        
            x = [-3.5 -2.5 2.5 3.5];
            y = [0 T.Hub-T.HubR/2 T.Hub-T.HubR/2 0];
            fill(x,y,'w')
        
        % Nacelle
        
            x = [-2.5 -2.5 2.5 2.5];
            y = [T.Hub-T.HubR/2 T.Hub+T.HubR/2 T.Hub+T.HubR/2 T.Hub-T.HubR/2];
            fill(x,y,'w')
        
        % Hub
        
            th = 0:pi/50:2*pi;
            xunit = T.HubR * cos(th) + 0;
            yunit = T.HubR * sin(th) + T.Hub;
            plot(xunit, yunit,'k');
            
            fill(xunit,yunit,'w')
        
        % Radius circle
        
            xunit = (T.R+2) * cos(th) + 0;
            yunit = (T.R+2) * sin(th) + T.Hub;
            plot(xunit, yunit,'--k');
        
            ha = annotation('doublearrow');                                     % store the arrow information in ha
            ha.Parent = gca;                                                    % associate the arrow the the current axes
            ha.X = [-T.R-2 T.R+2];                                              % the location in data units
            ha.Y = [T.Hub T.Hub]; 
            
            text(-T.R*2/3,T.Hub+3,'D = 114 m')
    
    % -------------------- Side View --------------------
    
        % Tower
        
            x = [-100-3.5 -100-2.5 -100+2.5 -100+3.5];
            y = [0 T.Hub-T.HubR/2 T.Hub-T.HubR/2 0];
            fill(x,y,'w')
        
        % Hub
        
            x = [-100-8.5 -100-8.5 -100+4.5 -100+4.5];
            y = [T.Hub-T.HubR/2 T.Hub+T.HubR/2 T.Hub+T.HubR/2 T.Hub-T.HubR/2];
            fill(x,y,'w')
        
            x = [-100+4.5 -100+4.5 -100+7.5 -100+7.5];
            y = [T.Hub-T.HubR T.Hub+T.HubR T.Hub+T.HubR T.Hub-T.HubR];
            fill(x,y,'w')
        
        % Hub cap
        
            th = linspace( pi/2, -pi/2, 100);
            xunit = T.HubR * cos(th) + -100+7.5;
            yunit = T.HubR * sin(th) + T.Hub;
            plot(xunit, yunit,'k');
            
            fill(xunit,yunit,'w')
        
        % Blade 1
        
            x = [-100+4.5 -100+5 -100+5 -100+7.5];
            y = [T.Hub-T.HubR yr(end)+T.Hub-1 yr(end)+T.Hub T.Hub-T.HubR];
            fill(x,y,'w')
        
        % Blade 2
        
            x = [-100+4.5 -100+7 -100+7 -100+7.5];
            y = [T.Hub-T.HubR yr(end)+T.Hub-4 yr(end)+T.Hub T.Hub-T.HubR];
            fill(x,y,'w')
        
        % Blade 3
        
            x = [-100+4.5 -100+5 -100+6 -100+7.5];
            y = [T.Hub+T.HubR T.Hub+T.R-20 T.Hub+T.R+2 T.Hub+T.HubR];
            fill(x,y,'w')
    
    axis equal
    ylabel('z (m)')
    ylim([0 170])
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[]);

end

if P.WindSchema == 1

    % Load chord length
    
    addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/BEMT/BEMT');
    
    [~,~,~,~,c] = LoadGeometry('Airfoil_17.txt');
    
    % Plot
    
    figure;
    hold on
    
    % -------------------- Side View --------------------

        % Blade 3
        
            angle = 330;
            
            Ys = c';
            Xs = linspace(0,T.R,length(Ys));
            
            a_rad=((angle*pi)./180);
            [THETA,R] = cart2pol(Xs,Ys);                                        %Convert to polar coordinates
            THETA = THETA+a_rad;                                                %Add a_rad to theta
            [~,yr] = pol2cart(THETA,R);                                         %Convert back to Cartesian coordinates
            
        % Tower
        
            x = [-100-3.5 -100-2.5 -100+2.5 -100+3.5];
            y = [0 T.Hub-T.HubR/2 T.Hub-T.HubR/2 0];
            fill(x,y,'w')
        
        % Hub
        
            x = [-100+8.5 -100+8.5 -100-4.5 -100-4.5];
            y = [T.Hub-T.HubR/2 T.Hub+T.HubR/2 T.Hub+T.HubR/2 T.Hub-T.HubR/2];
            fill(x,y,'w')
        
            x = [-100-4.5 -100-4.5 -100-7.5 -100-7.5];
            y = [T.Hub-T.HubR T.Hub+T.HubR T.Hub+T.HubR T.Hub-T.HubR];
            fill(x,y,'w')
        
        % Hub cap
        
            th = linspace( 3*pi/2, pi/2, 100);
            xunit = T.HubR * cos(th) + -100-7.5;
            yunit = T.HubR * sin(th) + T.Hub;
            plot(xunit, yunit,'k');
            
            fill(xunit,yunit,'w')
        
        % Blade 1
        
            x = [-100-4.5 -100-5 -100-5 -100-7.5];
            y = [T.Hub-T.HubR yr(end)+T.Hub-1 yr(end)+T.Hub T.Hub-T.HubR];
            fill(x,y,'w')
        
        % Blade 2
        
            x = [-100-4.5 -100-7 -100-7 -100-7.5];
            y = [T.Hub-T.HubR yr(end)+T.Hub-4 yr(end)+T.Hub T.Hub-T.HubR];
            fill(x,y,'w')
        
        % Blade 3
        
            x = [-100-4.5 -100-5 -100-6 -100-7.5];
            y = [T.Hub+T.HubR T.Hub+T.R-20 T.Hub+T.R+2 T.Hub+T.HubR];
            fill(x,y,'w')
    
    % -------------------- Wind Profiles --------------------
    
        Z = 0:0.5:300;

        k = 0.03;                                                               % Eddy diffusivity
        g = 4;                                                                  % Geostrophic wind

        EkmanProfile = sqrt((g .* (1 - exp(-sqrt((2*7.292e-5*sind(23))/ ...
            (2*k)).*Z) .* cos(sqrt((2*7.292e-5*sind(23))/(2*k)).*Z))).^2 ...
            +(g .* exp(-sqrt((2*7.292e-5*sind(23))/(2*k)).*Z) .* ...
            sin(sqrt((2*7.292e-5*sind(23))/(2*k)).*Z)).^2);

        PowerLaw = 300*(0.1*(Z./43).^(0.2))-190;

        EkmanProfile = (EkmanProfile * 10)-190;
        
        plot(PowerLaw,Z,'LineWidth',2)
        plot(EkmanProfile,Z,'LineWidth',2)
        xline(-165,'LineStyle','--')

        axis equal
        legend('','','','','','','','','Power Law','Ekman Layer','Location','northwest')
        ylabel('z (m)')
        xlabel('x (m)')
        ylim([0 170])
        xlim([-200 -75])
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[]);

end

if P.Capacity == 1

    figure;    
    subplot(2,1,1);
        Y = [6.5  0.1;7.3  0.2;8.1  0.3;8.2  0.1;11.5 0.1;14.7 0.1;20.3 0.3;    %Yearly new installations
             26.9 0.4;38.5 0.6;39.1 0.9;40.6 0.9;45   1.2;36   1.6;51.7 1.5;
             63.8 3.4;54.9 2.2;53.5 4.5;50.7 4.4;60.8 6.2;93   6.1];
    
        Years = 2001:2025;                                                      % Years in analysis
    
        hold on
        bar(Years(1:20),Y,'stacked')
    
        YPred = [87.487*.93 87.487*.07;                                         % Yearly new installation predictions
                 81.060*.93 81.060*.07;
                 90.520*.93 90.520*.07;
                 98.015*.93 98.015*.07;
                 112.224*.93 112.224*.07];
    
        b = bar(Years(21:25),YPred,'stacked','LineStyle',':','linewidth',2);
        
        b(1).FaceColor = '#0072BD';                                             % Face color of predictions
        b(2).FaceColor = '#D95319';
        b(1).EdgeColor = '#0072BD';                                             % Edge color of predictions
        b(2).EdgeColor = '#D95319';
        b(1).FaceAlpha = 0.5;                                                   % Face color opacity of predictions
        b(2).FaceAlpha = 0.5;
        ylabel('Yearly New Capacity (GW)','FontSize',10)
        legend('Onshore','Offshore','Location','northwest')
        xlim([Years(1)-0.5 Years(end)+0.5])
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[]);
        set(gca,'TickLength',[0 0])
        set(gca, 'Position', [0.1, 0.55, 0.85, 0.425]);
        text(2001,110,'(a)','FontSize',14)
        box on
    
    subplot(2,1,2); 
        Y = [24   0;31   0;39  -1;48  -1;59  -1;74  -1;94   1;121  1;159  2;    % Yearly installed capacity
             198  3;238  4;283  5;319  7;370  8;433 12;488 14;540 19;591 23;
             650 29;743 35];
    
        YPred2 = [824.3629 41.1241;                                             % Yearly installed capacity predictions
                  899.7487 46.7983;
                  983.9323 53.1347;
                  1.0751e3 59.9958;
                  1.1795e3 67.8515];
        hold on
        bar(Years(1:20),Y,'stacked')
    
        b = bar(Years(21:25),YPred2,'stacked','LineStyle',':','linewidth',2);
        
        b(1).FaceColor = '#0072BD';                                             % Face color of predictions
        b(2).FaceColor = '#D95319';
        b(1).EdgeColor = '#0072BD';                                             % Edge color of predictions
        b(2).EdgeColor = '#D95319';
        b(1).FaceAlpha = 0.5;                                                   % Face color opacity of predictions
        b(2).FaceAlpha = 0.5;
        ylabel('Total Capacity (GW)','FontSize',10)
        ylim([0 1300])
        xlim([Years(1)-0.5 Years(end)+0.5])
        set(gca,'TickLength',[0 0])
        text(2001,1175,'(b)','FontSize',14)
        box on
    
        set(gca, 'Position', [0.1, 0.1, 0.85, 0.425]);

end

if P.CapMap == 1

    addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Conditional Averages/borders');
    addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Conditional Averages');
        
    GT = readtable('GlobalTotals.csv');                                     % Import data
    
    onshore = table2array(GT(:,2));                                         % Extract capacity data
    countries = table2cell(GT(:,1));                                        % Extract country names
 
%     bins = [logspace(-6,2,10) 150];                                     % Capacity bins
%     bins = [0 10*bins(1) bins(3:7) 0.750 bins(8:end)];
    bins = [0 5/1e6 100/1e6 500/1e6 1e3 5/1e3 50/1e3 100/1e3 100/1e3 500/1e3 1 5 50 100 250];

    cmap = parula(length(bins));                                            % Construct colormap
    cmap = [[1 1 1]; cmap];

    hold on
    for i = [1 3 4 5 6 7 8 9 10 11 12 13 14 15]
        bar(i,0.1,'FaceColor',cmap(i,:))                                    % This is entirely so the legend works
    end

    borders('countries','k')                                                % Load global map
    for i = 1:length(onshore)
        colorI = find(bins < onshore(i), 1, 'last');                        % Find colormap index for each country
        borders(countries(i),'facecolor',cmap(colorI, :))                   % Apply color to each country
    end

    LegInfo = {'NO DATA';
               '0   kW - 100 kW';
               '100 kW - 500 kW';
               '500 kW - 1   MW';
               '1   MW - 5   MW';
               '5   MW - 50  MW';
               '50  MW - 100 MW';
               '100 MW - 500 MW';
               '500 MW - 1   GW';
               '1   GW - 5   GW';
               '5   GW - 50  GW';
               '50  GW - 100 GW';
               '100 GW - 250 GW';
               '> 250 GW'};

    Leg = legend(LegInfo);
    title(Leg,'Capacity') 
    set(Leg,'FontName','monospaced')

    axis tight
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    box on

end

if P.USCapMap == 1

    addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Conditional Averages/borders');
    addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Conditional Averages');
        
    GT = readtable('USTotals.csv');                                         % Import data
    
    onshore = table2array(GT(:,2));                                         % Extract capacity data
    countries = table2cell(GT(:,1));                                        % Extract country names
 
    bins = [1 50 100 250 500 1000 2500 5000 10000 15000 50000];                                     % Capacity bins
    cmap = parula(length(bins));                                            % Construct colormap
    cmap = [[1 1 1];cmap];

    hold on
    for i = 1:length(bins)
        bar(i,0.1,'FaceColor',cmap(i,:))                                    % This is entirely so the legend works
    end

    borders('continental US','k')                                           % Load global map
    for i = 1:length(onshore)
        colorI = find(bins < onshore(i), 1, 'last');                        % Find colormap index for each country
        borders(countries(i),'facecolor',cmap(colorI, :))                   % Apply color to each country
    end

    LegInfo = {'NO DATA';
               '1   MW - 50  MW';
               '50  MW - 100 MW';
               '100 MW - 250 MW';
               '250 MW - 500 MW';
               '500 MW - 1   GW';
               '1   GW - 2.5 GW';
               '2.5 GW - 5   GW';
               '5   GW - 10  GW';
               '10  GW - 15  GW'};

    Leg = legend(LegInfo);
    title(Leg,'Capacity') 
    set(Leg,'FontName','monospaced')

    axis tight
    xlim([-125 -64.5])
    ylim([23 52])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    box on

end

if P.WndFmMap == 1

    longs = [69.915731;69.921821;69.927336;69.930022;69.918603;69.931198;
             69.942441;69.941349;69.938305;69.928344;69.928476;69.931471;
             69.929320;69.933209;69.947703;69.950334;69.945089;69.951292;
             69.955339;69.971150;69.920886;69.948364;69.952242;69.912458;
             69.898884;69.953340;69.949652;69.980175;69.995010;70.002455;
             69.921132;69.914808;69.912355;69.905367;69.905153;69.904846;
             69.894236;69.887205;69.891898;69.884196;69.876398;69.868485;
             69.888151;69.883213;69.888244;69.898765;69.898859;69.890207;
             69.887086;69.979744;69.996783;69.943999;69.950159;69.949391;
             69.948115;69.946593;69.960869;69.959968;69.958938;69.964501;
             69.971110];

    lats  = [23.014843;23.012679;23.009665;23.006146;23.003806;23.003806;
             23.003399;22.998739;22.991458;22.987741;22.978922;22.969157;
             22.963916;22.957746;22.949595;22.956619;22.962667;23.007928;
             23.003962;23.006133;22.947109;22.917127;22.913670;22.921978;
             22.941076;22.997640;22.991812;23.006073;22.987714;22.984743;
             22.961927;22.963397;22.968327;22.962935;22.967161;22.971358;
             22.989367;22.992495;22.957835;22.957870;22.957405;22.956112;
             22.950510;22.947430;22.942954;22.949704;22.941063;22.936586;
             22.931929;22.949876;22.953555;23.026591;23.023910;23.027153;
             23.031698;23.038400;23.028809;23.032030;23.036633;23.038332;
             23.036943];

figure;

    % Contour

    imagesc(WndFamMap.long,WndFamMap.lat,WndFamMap.El);axis xy
    colormap summer;
    
    hold on; contour(WndFamMap.long,WndFamMap.lat,WndFamMap.El,'showtext','on','linecolor','#525252')

    % Turbine markers

    for i = 1:length(longs)

        plot(longs(i),lats(i),'marker','o','markersize',10,'color','k','markerfacecolor','#fc4903')

    end

    % Scale bar

    xoffset = 0.005;
    yoffset = 0.005;
    leftx   = WndFamMap.long(1) + xoffset;
    yheight = WndFamMap.lat(end) + yoffset;
%     leftx   = 69.87544;
%     yheight = 22.92045;

    barlength = 0.014565;
    barheight = 0.0007;

    rectangle('Position',[leftx,yheight,barlength,barheight],'FaceColor','#a7a7a7','EdgeColor','#a7a7a7','LineWidth',3)
    rectangle('Position',[leftx + barlength,yheight,barlength,barheight],'FaceColor','#c5c5c5','EdgeColor','#c5c5c5','LineWidth',3)
    text(leftx,yheight + 4*barheight,'0','HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
    text(leftx + barlength,yheight + 4*barheight,'1.5','HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
    text(leftx + 2 * barlength,yheight + 4*barheight,'3 km','HorizontalAlignment','center','FontSize',10,'FontWeight','bold')

    % Cardinal direction pointer

    centerx = 69.99;
    centery = 23.03;
    
        x = [centerx - 0.002, centerx, centerx];
        y = [centery - 0.002, centery, centery + 0.0032];
        c = [0 0 0];
        fill(x, y, c)
        
        x = [centerx, centerx, centerx + 0.002];
        y = [centery, centery + 0.0032, centery - 0.002];
        c = [1 1 1];
        fill(x, y, c)
    
    text(centerx,centery + 0.006,'N','HorizontalAlignment','center','FontSize',13.5,'FontWeight','bold')

    % LIDAR marker

    centerx = 69.9122;
    centery = 23.018;
    side = 0.00075;
    
        x = [centerx + side, centerx + side, centerx - side, centerx - side];
        y = [centery - side, centery + side, centery + side, centery - side];
        c = [0 0 0];
        fill(x, y, c)
    
%     centerx = 69.91875;
    centerx = 69.9122;
    centery = 23.018;
    side = 0.00075;
    
        x = [centerx, centerx - side, centerx + side];
        y = [centery + side, centery + 3*side, centery + 3*side];
        c = [0 0 0];
        fill(x, y, c)
    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    axis equal
    ylabel(colorbar('eastoutside'),'Elevation (m)')

    % Array marker

    rectangle('Position',[69.9142 23.0076 0.0153 0.0090],'linestyle','--','LineWidth',1.25)
%     rectangle('Position',[69.9126 22.9997 0.0225 0.0229],'LineWidth',1.25)

end

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
    title(Leg,'Wind Speed [m s^{-1}]')                                          % Set legend title
    set(Leg,'FontName','monospaced')                                            % Proper spacing
    legend('AutoUpdate','off')                                                  % Prevent subsequent plot from being added to legend
    box(Leg,'off')                                                              % Remove outline from legend
    set(gca,'rticklabel',[])                                                    % Remove default radial markers
    set(gca,'FontSize',18)
    
    TotalPoints = sum(WindSpeed>0);
    axlims    = rlim;                                                           % Find radial limits
    positions = linspace(0,axlims(2),5);                                        % Define new radial marker spacing
    rticks(positions);                                                          % Set new radial marker positions
    angle     = deg2rad(95);                                                    % Angle of radial marker vector
    
    h = polarhistogram(deg2rad(WindDirection(WindDirection>=330 | ...           % Add outline to wind of interest
        WindDirection<=50)),deg2rad(0:10:360),'LineWidth',2.5);
    h.DisplayStyle = 'stairs';                                                  % Set outline style
    h.EdgeColor    = [0 0 0];                                                   % Set color to black                                                   
    
    text(angle, positions(2), sprintf('%2.0f%%',100*positions(2)/TotalPoints),'HorizontalAlignment','center','FontSize',15);
    text(angle, positions(3), sprintf('%2.0f%%',100*positions(3)/TotalPoints),'HorizontalAlignment','center','FontSize',15);
    text(angle, positions(4), sprintf('%2.0f%%',100*positions(4)/TotalPoints),'HorizontalAlignment','center','FontSize',15);
    
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
    set(Leg,'FontName','monospaced')                                            % Proper spacing
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

if P.KmBetaFull == 1

    figure;
        imagesc(KB.Kmrange,KB.DSrange,KmBeta.Power);
        axis xy                                                                 % Flip image
        colormap([1 1 1; parula(256)])                                          % Bi-tone colomap
        colorbar
    
        xlabel('Ekman Parameter (K_m)')
        ylabel('Direction Shear (\circ m^{-1})')
        ylabel(colorbar('eastoutside'),'Normalized Power (-)')
        caxis([0.75 1.1])
        hold on
    
        for i = 1:length(KB.Kmrange)                                            % Vertical brid lines
            xline(KB.Kmrange(i)-KB.s/2)
        end
    
        for i = 1:length(KB.DSrange)                                            % Horizontal brid lines
            yline(KB.DSrange(i)-KB.s/2)
        end
    
        x = [-0.05 KB.Kmrange];                                                   % Add threshold line
        y = 2/3*x - 0.1;
        plot(x,y,'color','r','LineWidth',2)
        axis equal
        xlim([KB.xlow KB.xhigh])
%         xlim([-0.05 0.95])
        ylim([-0.15 0.65])
    
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

if P.ShearProfMedians == 1

    rng(T.Seed);
    
    [~,~,~,H,~,~] = datevec(D.Time);
    
    figure;
    hold on
    
    for i = 1:24
    
        h(i) = subplot(4,6,i);
        
        TOD = i-1;
        
        ind = 1:size(D.Shear,2);
        ind = ind(H == TOD);
        
        MedianU = median(D.Shear(:,ind)./D.Shear(end,ind),2);                           % Calculate median hourly profiles
        STDU = std(D.Shear(:,ind)./D.Shear(end,ind),[],2);                              % Calculate standard deviation of profiles
        indplot = datasample(ind,10);                                                   % Select 10 random profiles
        
            shadedErrorBar(MedianU, flip(T.Heights), STDU, 'lineProps',{'-', ...
                'Color', 'k', 'LineWidth', 2.5, 'MarkerSize', 2}, 'vertical', 0);
            plot(D.Shear(:,indplot)./D.Shear(end,indplot), flip(T.Heights), 'LineWidth', 1);
                yline(T.Hub,'k-','LineWidth',1.75)
                yline(T.Hub+T.R,'k--','LineWidth',0.75)
                yline(T.Hub-T.R,'k--','LineWidth',0.75)
                xlim([0.3,2])
                ylim([40,T.Heights(end)])
                titstr = sprintf('%2.0f:00',TOD);
                text(0.5,180,titstr,'HorizontalAlignment','center','FontSize',12)
                box on
        
            if i == 1 || i == 7 || i == 13 || i ==19
                ylabel('z (m)')
            end

            if i ~= 1 && i ~= 7 && i ~= 13 && i ~=19
                set(gca,'yticklabel',[])
            end

            if i <= 18
                set(gca,'xticklabel',[])
            end

            if i >= 19
                xlabel('U(z)/U_{43 m}')
            end
    
    end
    
    for i = 1:24
    
        ht = 0.23;
        s  = 0.01;
        b  = 0.04;
    
        if i >= 1 && i <= 6 
            set(h(i), 'Position', [0.03+(i-1)/6.2 s*4+3*ht+b 1/6.35 ht])
        elseif i >= 7 && i <= 12
            set(h(i), 'Position', [0.03+(i-7)/6.2 s*3+2*ht+b 1/6.35 ht])
        elseif i >= 13 && i <= 18
            set(h(i), 'Position', [0.03+(i-13)/6.2 s*2+ht+b 1/6.35 ht])
        elseif i >= 19
            set(h(i), 'Position', [0.03+(i-19)/6.2 s+b 1/6.35 ht])
        end
    
    end

end

if P.ShearMedError == 1

    [~,~,~,H,MN,~] = datevec(D.Time);
    
    h = 0:23;
    
    M.A6 = [0   4; 5   9; 10 14; 15 19; 20 24; 25 29; 30 34; 35 39; 40 44;
            45 49; 50 54; 55 59];
    
    M.A5 = [0   9; 10 19; 20 29; 30 39; 40 49; 50 59];
    
    M.A4 = [0  14; 15 29; 30 44; 45 59];
    
    M.A3 = [0  19; 20 39; 40 59];
    
    M.A2 = [0  29; 30 59];
    
    M.A1 = [0  59];
    
    tit  = [60 30 20 15 10 5];
    
    figure;
    hold on
    
    for j = 1:6
 
        clear m index taxis MedProfiles PLFullMedians EkmanMedians vector
        m           = M.(sprintf('A%1.0f',j));
        MedProfiles = zeros(12,size(m,1));
        index       = zeros(size(m,1),size(D.Shear,2));
    
        i = 1;
        
        for hour = 1:length(h)
            for minute = 1:size(m,1)
        
                index(i,:) = (H == h(hour) & MN >= m(minute,1) & MN <= m(minute,2));
        
                vector = index(i,:) .* D.Shear;
                vector = vector(:,any(vector));
        
                MedProfiles(:,i) = median(vector ,2);
                i = i+1;
            end
        end
        
        [PLFullMedians] = PowerLawFit(MedProfiles,T,'Full');
        [EkmanMedians]  = EkmanFit(MedProfiles,T);
        
        SB(j) = subplot(2,3,j);
            hold on
            taxis = duration(minutes(linspace(0,1439,size(MedProfiles,2))),'Format','hh:mm');
            plot(taxis,EkmanMedians.NRMSE,'LineWidth',0.75)
            plot(taxis,PLFullMedians.NRMSE,'LineWidth',0.75)
            xlim([taxis(1) taxis(end)])
            ylim([0 0.09])
            title(sprintf('%2.0f - Minute',tit(j)))
    
            if j == 2 || j == 3 || j == 5 || j == 6
                set(gca,'YTickLabel',[]);
            end
            
            if j == 1 || j == 2 || j == 3
                set(gca,'XTickLabel',[]);
            end
    
            if j == 1
                set(SB(j), 'Position', [0.05 0.52 1/3.3 0.45])
                ylabel('NRMSE')
                legend('Ekman','Power Law')
            elseif j == 2
                set(SB(j), 'Position', [0.06+1/3.3 0.52 1/3.3 0.45])
            elseif j == 3
                set(SB(j), 'Position', [0.07+2/3.3 0.52 1/3.3 0.45])
            elseif j == 4
                set(SB(j), 'Position', [0.05 0.04 1/3.3 0.45])
                ylabel('NRMSE')
            elseif j == 5
                set(SB(j), 'Position', [0.06+1/3.3 0.04 1/3.3 0.45])
            else
                set(SB(j), 'Position', [0.07+2/3.3 0.04 1/3.3 0.45])
            end

    end

end

% if P.PLFullRProb == 1
% 
%     figure;                                                                 % Inflection profile fit
%         [Probs,Edges] = histcounts(PLFull.R,250,...
%             'Normalization','probability','BinLimits',[-1.5 1.1]);
%         plot(Edges(1:end-1),Probs,'LineWidth',1.5);
% 
%         hold on
% 
%         [Probs,Edges] = histcounts(PLFullold.R,250,...
%             'Normalization','probability','BinLimits',[-1.5 1.1]);
%         plot(Edges(1:end-1),Probs,'LineWidth',1.5);
% 
%         xlim([-1.5 0.99])
% 
%         legend('Full Profile Optimization','Full Profile Toolbox')
% 
% end
% 
% if P.PLInflecRProb == 1
% 
%     figure;                                                                 % Inflection profile fit
%         [Probs,Edges] = histcounts(PLInflec.R,150,...
%             'Normalization','probability','BinLimits',[0 1.1]);
%         plot(Edges(1:end-1),Probs,'LineWidth',1.5);
% 
%         hold on
% 
%         [Probs,Edges] = histcounts(PLInflecold.R,150,...
%             'Normalization','probability','BinLimits',[0 1.1]);
%         plot(Edges(1:end-1),Probs,'LineWidth',1.5);
% 
%         xlim([0 0.99])
% 
%         legend('Partial Profile Optimization','Partial Profile Toolbox')
% 
% end

% figure;
%     plot(datetime(D.Time,'ConvertFrom','datenum'),Ekman.R)
% 
%     figure;                                                                 % Inflection profile fit
%         [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(Ekman.R,...
%             'Normalization','probability','BinLimits',[-0.5 1.5]);
%         plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','k','LineWidth',1.5);
% 
%         hold on
% 
%         [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(PLFull.R,...
%             'Normalization','probability','BinLimits',[-0.5 1.5]);
%         plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','r','LineWidth',1.5);
% 
% 
%         [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(PLInflec.R,...
%             'Normalization','probability','BinLimits',[-0.5 1.5]);
%         plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','b','LineWidth',1.5);
% 
%         legend('Ekman','Full','Inflec')
% 
% 
%         title('Speed Shear Probability Distribution - Inflection Fit')
%         xlabel('Speed Shear (\alpha)')
%         ylabel('Probability of Occurrence')
%         xlim([-0.35 1.4])

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