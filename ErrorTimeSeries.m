function [Ep] = ErrorTimeSeries(D,WindBins,Mean,T)
%ErrorTimeSeries Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n------Error Calcs-------')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

TotalOps   = 28880;
OpsCounter = 0;

[~,~,~,H,MN,~] = datevec(D.Time);

%% Calculations 

% ----------------- Error: Ep -----------------

    for i = 1:size(D.Shear,2)
    
        BinNum = find(D.Shear(T.HubRow,i) <= WindBins,1,'first');
    
        Ep.abs(i)  = abs(D.Power(i) - Mean.All(BinNum))/Mean.All(BinNum);
        Ep.diff(i) = (D.Power(i) - Mean.All(BinNum))/Mean.All(BinNum);
    
        OpsCounter = OpsCounter + 1;
        if mod(OpsCounter,10)==0
            p = OpsCounter/TotalOps*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end

    end

% ----------------- Hourly Ep Distributions -----------------

    for i = 1:24
    
        Indices = (H == i-1);
    
        Ep.dist(i,:) = Indices .* Ep.abs;
        Ep.mean(i)   = mean(nonzeros(Ep.dist(i,:)));

        OpsCounter = OpsCounter + 1;
        if mod(OpsCounter,10)==0
            p = OpsCounter/TotalOps*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end
    
    end

% ----------------- 5 - Min Average of Ep -----------------

    for i = 1:5:size(Ep.abs,2)
        if i > floor(size(Ep.abs,2)/5)*5
            Ep.Avg(i) = mean(Ep.abs(i:end));
        else
            Ep.Avg(i) = mean(Ep.abs(i:i+4));
        end

        OpsCounter = OpsCounter + 1;
        if mod(OpsCounter,10)==0
            p = OpsCounter/TotalOps*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end

    end
    
    Ep.Avg = nonzeros(Ep.Avg);

% ----------------- Time of Day 1 - Min Ep Average -----------------
    
    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);
        
            Ep.TOD(index,:) = Indices .* Ep.abs;
    
            j     = j + 1;
    
            index = index + 1;

            OpsCounter = OpsCounter + 1;
            if mod(OpsCounter,10)==0
                p = OpsCounter/TotalOps*100;
                fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
            end
    
        end
    
        i = i + 1;
        j = 0;

        OpsCounter = OpsCounter + 1;
        if mod(OpsCounter,10)==0
            p = OpsCounter/TotalOps*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end
    
    end
    
    for i = 1:size(Ep.TOD,1)
    
        Ep.TODmean(i) = mean(nonzeros(Ep.TOD(i,:)));

        OpsCounter = OpsCounter + 1;
        if mod(OpsCounter,10)==0
            p = OpsCounter/TotalOps*100;
            fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                              % print progress to screen
        end
    
    end

%% Plots

dateaxis = datetime(D.Time,'ConvertFrom','datenum');
timeaxis = duration(minutes(0:1:1439),'Format','hh:mm');
h = zeros(1,12);

Ep.abs1 = mean(reshape(Ep.abs(1:end-3),10,[]));
dateaxis1 = dateaxis(5:10:end);

% ----------------- Error: Ep -----------------

    figure
        plot(dateaxis1,Ep.abs1,'LineStyle','none','Marker','.','Color','k')
            title(['$\epsilon_P = \frac{|\overline{P} - P(t)|}' ...
                '{\overline{P}}$'],'interpreter','latex','FontSize',18)
            ylabel('$\epsilon_P\;(-)$','interpreter','latex','FontSize',18)
            ylim([0 4])
    
    figure
        plot(dateaxis,Ep.diff,'LineStyle','none','Marker','.','Color','k')
            title('$\epsilon_P = \overline{P} - P(t)$','interpreter',...
                  'latex','FontSize',18)
            ylabel('$\epsilon_P\;(-)$','interpreter','latex','FontSize',18)

% ----------------- Time-of-Day Ep -----------------

    figure
        plot(timeaxis,Ep.TODmean,'LineStyle','none','Marker','.','Color','k')
            xtickformat('hh:mm')
            xlim([timeaxis(1) timeaxis(end)+minutes(1)])
            ylabel('\epsilon_P (kW)')
            title('Average Time-of-Day Error in Power')
%%
% -------------- Hourly Histogram: Day Time --------------

    DayFig   = figure;
    BinWidth = 0.02;

        HourIndices = [5 6 7 8 9  10 11 12 13 14 15 16;
                       6 7 8 9 10 11 12 13 14 15 16 17];
    
        for i = 1:12

            h(i) = subplot(1,12,i);

            [counts,edges,~] = histcounts(nonzeros(Ep.dist(HourIndices(2,i),:)),'BinWidth',BinWidth);
            barh(edges(2:end),counts)

            str = sprintf('N = %3.0f',length(nonzeros(Ep.dist(HourIndices(2,i),:))));
            text(212.5,0.8,str,'HorizontalAlignment','center','FontSize',12)
            xlim([0 425])
            ylim([0 0.9])
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
%%
% -------------- Hourly Histogram: Night Time --------------

    NightFig   = figure;
    BinWidth = 0.02;

%                   [1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]
%                   [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1]

    HourIndices = [17 18 19 20 21 22 23 0 1 2 3 4;
                   18 19 20 21 22 23 24 1 2 3 4 5];

    for i = 1:12

        h(i) = subplot(1,12,i);

        [counts,edges,~] = histcounts(nonzeros(Ep.dist(HourIndices(2,i),:)),'BinWidth',BinWidth);
        barh(edges(2:end),counts)

        str = sprintf('N = %3.0f',length(nonzeros(Ep.dist(HourIndices(2,i),:))));
        text(250,0.8,str,'HorizontalAlignment','center','FontSize',12)
        xlim([0 500])
        ylim([0 0.9])
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
%%
fprintf('\n')





















% figure
% bar(dateaxis(1:5:end),m)


% HourBins = 0:23;
% 
% for i = 1:length(HourBins)
% 
%     Indices = (H == i-1);
% 
%     Ep.dist(i,:) = Indices .* Ep.abs;
%     Ep.mean(i)   = mean(nonzeros(Ep.dist(i,:)));
% 
% end
% 
% HourBins = 0:23;
% MinBins  = 0:59;

end