%% Additional Plots

%fig1
plot(datetime(data.DateTime,'ConvertFrom','datenum'),data.Lidar.H104m.WndSpd,'linestyle','none','marker','.')
ylabel('WindSpeed (z = z_h) (m/s)')
title('LiDAR Hub Height Wind Speed Measurements')

%fig3
plot(datetime(data.DateTime,'ConvertFrom','datenum'),data.BHR62.TurbulenceIntensity,'linestyle','none','marker','.')
ylabel('I (%)')
title('Turbulence Intensity at Hub Height')
ylim([0 100])
xlim([datetime(datenum('01-Apr-2020 00:00:00'),'ConvertFrom','datenum') datetime(datenum('02-Apr-2020 00:00:00'),'ConvertFrom','datenum')])

%fig2
plot(datetime(data.DateTime,'ConvertFrom','datenum'),data.BHR62.TurbulenceIntensity,'linestyle','none','marker','.')
ylabel('I (%)')
title('Turbulence Intensity at Hub Height')
ylim([0 100])

%fig4
plot(datetime(data.DateTime,'ConvertFrom','datenum'),data.Lidar.H104m.WndSpd,'linestyle','none','marker','.')
ylabel('WindSpeed (z = z_h) (m/s)')
title('LiDAR Hub Height Wind Speed Measurements')
xlim([datetime(datenum('29-Mar-2020 00:00:00'),'ConvertFrom','datenum') datetime(datenum('30-Mar-2020 00:00:00'),'ConvertFrom','datenum')])

%fig5
WindRose(gca, VeerOff(7,:), Shear(7,:), 5 : 10 : 365, 1 : 2 : 15, [5 10 15 20], 'Incoming Wind at Hub Height')
%WindRoseV(VeerOff(7,:),Shear(7,:),'anglenorth',0,'angleeast',90)

%fig6
plot(datetime(Time,'ConvertFrom','datenum'),Shear(7,:) - NormalWind(7,:),'linestyle','none','marker','.')
ylabel('WindSpeed (z = z_h) (m/s)')
title('Difference in Incoming Wind and Incident Wind at Hub Height')
legend('(U_{\infty} - U_{\perp})')

%fig7
xmin = 3;
xmax = 9.25;
RT = T.Hub + T.R;
RL = T.Hub - T.R;
subplot(1,5,1)
    plot(Shear(:,(Time==datenum('13-Feb-2020 12:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('13-Feb; 12:00')
    xlim([xmin xmax])
    ylabel('z (m)')
subplot(1,5,2)
    plot(Shear(:,(Time==datenum('16-Feb-2020 12:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('16-Feb; 12:00')
    xlim([xmin xmax])
subplot(1,5,3)
    plot(Shear(:,(Time==datenum('21-Feb-2020 12:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('21-Feb; 12:00')
    xlim([xmin xmax])
    xlabel('u (m/s)')
subplot(1,5,4)
    plot(Shear(:,(Time==datenum('22-Feb-2020 12:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('22-Feb; 12:00')
    xlim([xmin xmax])
subplot(1,5,5)
    plot(Shear(:,(Time==datenum('23-Feb-2020 12:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('23-Feb; 12:00')
    xlim([xmin xmax])

%fig8
xmin = 4;
xmax = 14;
RT = T.Hub + T.R;
RL = T.Hub - T.R;
subplot(1,5,1)
    plot(Shear(:,(Time==datenum('13-Feb-2020 05:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('13-Feb; 05:00')
    xlim([xmin xmax])
    ylabel('z (m)')
subplot(1,5,2)
    plot(Shear(:,(Time==datenum('16-Feb-2020 05:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('16-Feb; 05:00')
    xlim([xmin xmax])
subplot(1,5,3)
    plot(Shear(:,(Time==datenum('21-Feb-2020 05:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('21-Feb; 05:00')
    xlim([xmin xmax])
    xlabel('u (m/s)')
subplot(1,5,4)
    plot(Shear(:,(Time==datenum('22-Feb-2020 05:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('22-Feb; 05:00')
    xlim([xmin xmax])
subplot(1,5,5)
    plot(Shear(:,(Time==datenum('23-Feb-2020 05:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('23-Feb; 05:00')
    xlim([xmin xmax])

%fig9
xmin = 4;
xmax = 14;
RT = T.Hub + T.R;
RL = T.Hub - T.R;
subplot(1,5,1)
    plot(Shear(:,(Time==datenum('13-Feb-2020 02:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('13-Feb; 02:00')
    xlim([xmin xmax])
    ylabel('z (m)')
subplot(1,5,2)
    plot(Shear(:,(Time==datenum('16-Feb-2020 00:45:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('16-Feb; 00:45')
    xlim([xmin xmax])
subplot(1,5,3)
    plot(Shear(:,(Time==datenum('21-Feb-2020 02:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('21-Feb; 02:00')
    xlim([xmin xmax])
    xlabel('u (m/s)')
subplot(1,5,4)
    plot(Shear(:,(Time==datenum('22-Feb-2020 02:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('22-Feb; 02:00')
    xlim([xmin xmax])
subplot(1,5,5)
    plot(Shear(:,(Time==datenum('23-Feb-2020 02:00:00'))),flip(T.Heights'),'LineWidth',1.25)
    yline(RT,'LineStyle','--')
    yline(RL,'LineStyle','--')
    title('23-Feb; 02:00')
    xlim([xmin xmax])



%fig10
   subplot(2,2,[1 2])
   hold on
   plot(WindBins,AllMean,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
   plot(WindBins,LargerMean,'LineWidth',1.5,'Color','#D95319','LineStyle','--')
   plot(WindBins,LessMean,'LineWidth',1.5,'Color','#A2142F','LineStyle','-.')
   subplot(2,2,[3 4])
   bar(WindBins,Sig)
   xlim([0 15])


%fig11.1
NBINS = 40;
for i = 7:19
h(i-6) = subplot(1,13,i-6);
%histogram(nonzeros(LargerDist(i,:)),NBINS)
[counts,edges,bins] = histcounts(nonzeros(AllDist(i,:)),NBINS);
barh(edges(2:end),counts)
str = sprintf('N = %3.0f',length(nonzeros(AllDist(i,:))));
if i < 16
    text(150,2100,str,'HorizontalAlignment','center','FontSize',12)
else
    text(150,100,str,'HorizontalAlignment','center','FontSize',12)
end
xlim([0 300])
ylim([0 2200])
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
set(h(i), 'Position', [(1/14.5*(i)) .1 1/17 .81])
end

%fig11.2
NBINS = 40;
for i = 7:19
h(i-6) = subplot(1,13,i-6);
%histogram(nonzeros(LargerDist(i,:)),NBINS)
[counts,edges,bins] = histcounts(nonzeros(LargerDist(i,:)),NBINS);
barh(edges(2:end),counts)
str = sprintf('N = %3.0f',length(nonzeros(LargerDist(i,:))));
if i < 16
    text(55,2100,str,'HorizontalAlignment','center','FontSize',12)
else
    text(55,100,str,'HorizontalAlignment','center','FontSize',12)
end
xlim([0 110])
ylim([0 2200])
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
set(h(i), 'Position', [(1/14.5*(i)) .1 1/17 .81])
end


%fig12
NBINS = 40;
for i = 7:19
h(i-6) = subplot(1,13,i-6);
%histogram(nonzeros(LargerDist(i,:)),NBINS)
[counts,edges,bins] = histcounts(nonzeros(LessDist(i,:)),NBINS);
barh(edges(2:end),counts)
str = sprintf('N = %3.0f',length(nonzeros(LessDist(i,:))));
if i < 16
    text(112.5,2100,str,'HorizontalAlignment','center','FontSize',12)
else
    text(112.5,100,str,'HorizontalAlignment','center','FontSize',12)
end
xlim([0 225])
ylim([0 2200])
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
set(h(i), 'Position', [(1/14.5*(i)) .1 1/17 .81])
end

%fig13
h(1) = subplot(2,1,1);
    hold on
    plot(WindBins,AllMean,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
    plot(WindBins,LargerMean,'LineWidth',1.5,'Color','#D95319','LineStyle','--')
    ylabel('Power (kW)')
    yyaxis right
    ylabel('Counts')
    bar(WindBins,AllNum,'FaceColor','k')
    bar(WindBins,LargerNum,'FaceColor','r')
    alpha(0.05)
    str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(AllNum))),'\d{3}(?=\d)', '$0,')));
    str2 = sprintf('AA > Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(LargerNum))),'\d{3}(?=\d)', '$0,')));
    legend(str1, str2, 'All Data','AA > Hub')
    xlabel('u(z = z_h) (m/s)')
    grid on
    titstr = sprintf('T-Test for Turbine %2.0f Where AA > Hub',T.TOI);
    title(titstr)
    axprop = gca;
    axprop.YAxis(1).Color = 'k';
    axprop.YAxis(2).Color = '#800000';
    xlim([0 16])
    hold off
h(2) = subplot(2,1,2);
    bar(WindBins,SigLarge,'FaceColor','#308538')
    xlim([0 16])
    set(gca,'ytick',[])
    set(gca,'xtick',[])
set(h(1), 'Position', [.1 .25 .8 .71])
set(h(2), 'Position', [.1 .1 0.8 .05])


%fig14
h(1) = subplot(2,1,1);
    hold on
    plot(WindBins,AllMean,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
    plot(WindBins,LessMean,'LineWidth',1.5,'Color','#A2142F','LineStyle','-.')
    ylabel('Power (kW)')
    yyaxis right
    ylabel('Counts')
    bar(WindBins,AllNum,'FaceColor','k')
    bar(WindBins,LessNum,'FaceColor','b')
    alpha(0.05)
    str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(AllNum))),'\d{3}(?=\d)', '$0,')));
    str3 = sprintf('AA < Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(LessNum))),'\d{3}(?=\d)', '$0,')));
    legend(str1, str3, 'All Data','AA < Hub')
    xlabel('u(z = z_h) (m/s)')
    grid on
    titstr = sprintf('T-Test for Turbine %2.0f Where AA < Hub',T.TOI);
    title(titstr)
    axprop = gca;
    axprop.YAxis(1).Color = 'k';
    axprop.YAxis(2).Color = '#800000';
    xlim([0 16])
    hold off
h(2) = subplot(2,1,2);
    bar(WindBins,SigLess,'FaceColor','#308538')
    xlim([0 16])
    set(gca,'ytick',[])
    set(gca,'xtick',[])
set(h(1), 'Position', [.1 .25 .8 .71])
set(h(2), 'Position', [.1 .1 0.8 .05])

%fig15
h(1) = subplot(2,1,1);
    hold on
    plot(WindBins,AllMean,'LineWidth',1.5,'Color','#0072BD','LineStyle','-')
    plot(WindBins,LargerMean,'LineWidth',1.5,'Color','#D95319','LineStyle','--')
    plot(WindBins,LessMean,'LineWidth',1.5,'Color','#A2142F','LineStyle','-.')
    ylabel('Power (kW)')
    yyaxis right
    ylabel('Counts')
    bar(WindBins,AllNum,'FaceColor','k')
    bar(WindBins,LessNum,'FaceColor','b')
    bar(WindBins,LargerNum,'FaceColor','r')
    alpha(0.05)
    str1 = sprintf('All Average, N = %s',fliplr(regexprep(fliplr(num2str(sum(AllNum))),'\d{3}(?=\d)', '$0,')));
    str2 = sprintf('AA > Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(LargerNum))),'\d{3}(?=\d)', '$0,')));
    str3 = sprintf('AA < Hub, N = %s',fliplr(regexprep(fliplr(num2str(sum(LessNum))),'\d{3}(?=\d)', '$0,')));
    legend(str1, str2, str3,'All Data','AA < Hub','AA > Hub')
    xlabel('u(z = z_h) (m/s)')
    grid on
    titstr = sprintf('T-Test for Turbine %2.0f',T.TOI);
    title(titstr)
    axprop = gca;
    axprop.YAxis(1).Color = 'k';
    axprop.YAxis(2).Color = '#800000';
    xlim([0 16])
    hold off
h(2) = subplot(2,1,2);
    bar(WindBins,SigBoth,'FaceColor','#308538')
    xlim([0 16])
    set(gca,'ytick',[])
    set(gca,'xtick',[])
set(h(1), 'Position', [.1 .25 .8 .71])
set(h(2), 'Position', [.1 .1 0.8 .05])


%fig12
index = find(alpha>1/10 & alpha<1/5);
for i = 1:13
h(i) = subplot(1,13,i);
hold on
plot(LWP(:,index(i)),Z)
plot(Shear(4:12,index(i)),Zlin,'*','color','k')
%
%    text(6,150,str,'HorizontalAlignment','center','FontSize',12)
xlim([4 8])
ylim([50 150])
hold off

if i == 1
    ylabel('z (m)')
end
if i == 7
    xlabel('Wind Speed (m/s)')
    title('Select Log Wind Profile Fits')
end
if i ~= 1
    set(gca,'ytick',[])
end
end

for i = 1:13
set(h(i), 'Position', [(1/14.5*(i)) .1 1/17 .81])
end


%% Plots

%dateaxis = datetime(Time,'ConvertFrom','datenum');
% figure;
%     plot(Time,I,'LineStyle','none','Marker','.','Color','k')
%     titstr = sprintf('Turbulence Intensity at Hub Height for Turbine %2.0f',T.TOI);
%     title(titstr)
%     datetick('x','mm/yy','keepticks')
%     ylabel('I (%)')

% figure;
%     plot([0 2.9999 xi 25 25.0001 30],[0 0 yi 1600 0 0],'LineWidth',1,'Color','k')
%     hold on
%     xlabel('Average Wind Speed (m/s)')
%     ylabel('Power (kW)')
%     xline(3,'LineStyle','--')
%     xline(10.5,'LineStyle','--')
%     xline(T.CutOut,'LineStyle','--')
%     text(1.5,2350,'I','HorizontalAlignment','center','FontSize',16)
%     text(6.75,2350,'II','HorizontalAlignment','center','FontSize',16)
%     text(17.5,2350,'III','HorizontalAlignment','center','FontSize',16)
%     text(27.5,2350,'IV','HorizontalAlignment','center','FontSize',16)
%     hold off
%     ylim([0 2600])
%     grid on

% figure;
%     plot(Time,Veer(7,:),'LineStyle','none','Marker','.','Color','k')
%     titstr = sprintf('Wind Direction at Hub Height for Turbine %2.0f',T.TOI);
%     title(titstr)
%     datetick('x','mm/yy','keepticks')
%     ylabel('Wind Direction (deg)')

% figure;
%     plot(Time,Nacelle,'LineStyle','none','Marker','.','Color','k')
%     titstr = sprintf('Nacelle Position for Turbine %2.0f',T.TOI);
%     title(titstr)
%     datetick('x','mm/yy','keepticks')
%     ylabel('Heading (deg)')

% figure;
%     plot(dateaxis,Diff(7,:),'LineStyle','none','Marker','.','Color','k')
%     titstr = sprintf('Relative Angle between Nacelle and Incoming Wind for Turbine %2.0f',T.TOI);
%     title(titstr)
%     ylabel('Relative Bearing (deg)')
%     xlabel('Time')

% figure; hold on; scatter(LessAAWind,LessAAPower); plot(XMarks,AllPowerMean,'linewidth',2); title('AA<Hub'); hold off
% figure; hold on; scatter(GreaterAAWind,GreaterAAPower); plot(XMarks,AllPowerMean,'linewidth',2); title('AA>Hub'); hold off
% figure; hold on; scatter(SheerFltd(HubRow,:),PowerFltd); plot(XMarks,AllPowerMean,'linewidth',2); title('All Data'); hold off