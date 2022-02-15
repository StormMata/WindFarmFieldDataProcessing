function [Ep] = ErrorTimeSeries(Shear,Power,WindBins,Mean,Time,T)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:size(Shear,2)

    BinNum = find(Shear(T.HubRow,i) <= WindBins,1,'first');

    Ep.abs(i)  = abs(Power(i) - Mean.All(BinNum));
    Ep.diff(i) = Power(i) - Mean.All(BinNum);

end

dateaxis = datetime(Time,'ConvertFrom','datenum');

figure
    plot(dateaxis,Ep.abs,'LineStyle','none','Marker','.','Color','k')
    title('$\epsilon_P = |\overline{P} - P(t)|$','interpreter','latex','FontSize',18)
    ylabel('$\epsilon_P\;(kW)$','interpreter','latex','FontSize',18)

figure
    plot(dateaxis,Ep.diff,'LineStyle','none','Marker','.','Color','k')
    title('$\epsilon_P = \overline{P} - P(t)$','interpreter','latex','FontSize',18)
    ylabel('$\epsilon_P\;(kW)$','interpreter','latex','FontSize',18)

for i = 1:5:size(Ep.abs,2)
    if i > floor(size(Ep.abs,2)/5)*5
        m(i) = mean(Ep.abs(i:end));
    else
        m(i) = mean(Ep.abs(i:i+4));
    end
end

m = nonzeros(m);

figure
bar(dateaxis(1:5:end),m)


[~,~,~,H,MN,~] = datevec(Time);

HourBins = 0:23;

for i = 1:length(HourBins)

    Indices = (H == i-1);

    Ep.dist(i,:) = Indices .* Ep.abs;
    Ep.mean(i)   = mean(nonzeros(Ep.dist(i,:)));

end

HourBins = 0:23;
MinBins  = 0:59;

for i = 1:length(HourBins)

    Indices = (H == i-1);

    Ep.dist(i,:) = Indices .* Ep.abs;
    Ep.mean(i)   = mean(nonzeros(Ep.dist(i,:)));

end

index = 1;
i     = 0;
j     = 0;

while i <= 23

    while j <= 59

        Indices = (H == i & MN == j);
    
        Ep.TOD(index,:) = Indices .* Ep.abs;

        j = j + 1;

        index = index + 1;

    end

    i = i + 1;
    j = 0;
end

for i = 1:size(Ep.TOD,1)
    TODmean(i) = mean(nonzeros(Ep.TOD(i,:)));
end

alltime=linspace(...
    datetime('00:00','InputFormat','HH:mm'),...
    datetime('23:59','InputFormat','HH:mm'),1440);

%alltime = timeofday(alltime);
plot(alltime,TODmean,'LineStyle','none','Marker','.','Color','k')
xtickformat('hh:mm')
xlim([alltime(1) alltime(end)+minutes(1)])
title('Time-of-Day Error in Power')
ylabel('\epsilon_P (kW)')

%% plot
% -------------- Day Time --------------
fig = figure;
    BinWidth = 20;

    for i = 13:24
        h(i-12) = subplot(1,12,i-12);
        [counts,edges,bins] = histcounts(nonzeros(Ep.dist(i,:)),'BinWidth',BinWidth);
        barh(edges(2:end),counts)
        str = sprintf('N = %3.0f',length(nonzeros(Ep.dist(i,:))));
        text(57.5,1250,str,'HorizontalAlignment','center','FontSize',12)
        xlim([0 115])
        ylim([0 1300])
        xtickangle(45)
        str = sprintf('%02.0f:00',i-1);
        title(str)
            if i ~= 13
                set(gca,'ytick',[])
            end
    end
    
    AxHandle=axes(fig,'visible','off'); 
    AxHandle.Title.Visible='on';
    AxHandle.XLabel.Visible='on';
    AxHandle.XLabel.Position = [0.5 0];
    xlabel(AxHandle,'Counts');
    AxHandle.YLabel.Visible='on';
    AxHandle.YLabel.Position = [0 0.5];
    ylabel(AxHandle,'\epsilon_P (kW)');
    title(AxHandle,'\epsilon_P Distribution by Hour');
    
    for i = 1:12
        set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])
    end

%% plot
% -------------- Night Time --------------
fig = figure;
    BinWidth = 20;
    for i = 1:12
        h(i) = subplot(1,12,i);
        [counts,edges,bins] = histcounts(nonzeros(Ep.dist(i,:)),'BinWidth',BinWidth);
        barh(edges(2:end),counts)
        str = sprintf('N = %3.0f',length(nonzeros(Ep.dist(i,:))));
        text(150,1250,str,'HorizontalAlignment','center','FontSize',12)
        xlim([0 300])
        ylim([0 1300])
        xtickangle(45)
        str = sprintf('%02.0f:00',i-1);
        title(str)
            if i ~= 1
                set(gca,'ytick',[])
            end
    end
    
    AxHandle=axes(fig,'visible','off'); 
    AxHandle.Title.Visible='on';
    AxHandle.XLabel.Visible='on';
    AxHandle.XLabel.Position = [0.5 0];
    xlabel(AxHandle,'Counts');
    AxHandle.YLabel.Visible='on';
    AxHandle.YLabel.Position = [0 0.5];
    ylabel(AxHandle,'\epsilon_P (kW)');
    title(AxHandle,'\epsilon_P Distribution by Hour');
    
    for i = 1:12
        set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])
    end

%% Hourly average

end