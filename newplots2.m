%%

z = 43:200;
i= 150;
u = PLFull.Uref(i)*(z/43).^PLFull.alpha(i);

plot(D.Shear(:,i),flip(T.Heights))
hold on
plot(u,z)
yline(T.Hub-T.R,'LineStyle','--')
yline(T.Hub+T.R,'LineStyle','--')
legend('Shear Profile','Power Law Fit')
xlabel('Wind Speed (m/s)')
ylabel('z (m)')


%%

[~,~,~,H,~,~] = datevec(D.Time);

TOD = 12;

ind = 1:size(D.Shear,2);
ind = ind(H == TOD);

rng(T.Seed);

indplot = datasample(ind,10);

% ind_hub = find(T.Heights==104);
shadedErrorBar(median(D.Veer - D.Veer(T.HubRow,:),2), flip(T.Heights), std(D.Veer - D.Veer(T.HubRow,:),[],2), ...
    'lineProps',{'-', 'Color', 'k', 'LineWidth', 2.5, ...
    'MarkerSize', 2}, 'vertical', 0);
hold on
plot([-30,30], [104,104], 'k-', 'LineWidth',1.75);
plot(D.Veer(:,indplot) - D.Veer(T.HubRow,indplot), flip(T.Heights), 'LineWidth', 1);
xlim([-30,30]); 
ylim([-10,200]);


%% Normal Wind

shadedErrorBar(median(NormalWind - NormalWind(T.HubRow,:),2), flip(T.Heights), std(NormalWind - NormalWind(T.HubRow,:),[],2), ...
    'lineProps',{'-', 'Color', 'k', 'LineWidth', 2.5, ...
    'MarkerSize', 2}, 'vertical', 0);
hold on
plot([-30,30], [104,104], 'k-', 'LineWidth',1.75);
plot(NormalWind(:,indplot) - NormalWind(T.HubRow,indplot), flip(T.Heights), 'LineWidth', 1);
xlim([-5,5]); 





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
  plot(xi,yi,'LineWidth',0.85,'Color','k','LineStyle','--')
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

datetime('01 FEB 2020','format','dd MMM yyyy'):datetime('29 Mar 2020','format','dd MMM yyyy')

suntime = {'7:29';
           '7:28';
           '7:28';
           '7:27';
           };
