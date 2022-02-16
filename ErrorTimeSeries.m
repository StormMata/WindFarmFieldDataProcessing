function [Ep] = ErrorTimeSeries(Shear,Power,WindBins,Mean,Time,T)
%ErrorTimeSeries Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n------Error Calcs-------')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

TotalOps   = 28880;
OpsCounter = 0;

%% Calculations 

% ----------------- Error: Ep -----------------

    for i = 1:size(Shear,2)
    
        BinNum = find(Shear(T.HubRow,i) <= WindBins,1,'first');
    
        Ep.abs(i)  = abs(Power(i) - Mean.All(BinNum));
        Ep.diff(i) = Power(i) - Mean.All(BinNum);
    
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

dateaxis = datetime(Time,'ConvertFrom','datenum');
timeaxis = duration(minutes(0:1:1439),'Format','hh:mm');
h = zeros(1,12);

% ----------------- Error: Ep -----------------

    figure
        plot(dateaxis,Ep.abs,'LineStyle','none','Marker','.','Color','k')
            title('$\epsilon_P = |\overline{P} - P(t)|$','interpreter',...
                  'latex','FontSize',18)
            ylabel('$\epsilon_P\;(kW)$','interpreter','latex','FontSize',18)
    
    figure
        plot(dateaxis,Ep.diff,'LineStyle','none','Marker','.','Color','k')
            title('$\epsilon_P = \overline{P} - P(t)$','interpreter',...
                  'latex','FontSize',18)
            ylabel('$\epsilon_P\;(kW)$','interpreter','latex','FontSize',18)

% ----------------- Time-of-Day Ep -----------------

    figure
        plot(timeaxis,Ep.TODmean,'LineStyle','none','Marker','.','Color','k')
            xtickformat('hh:mm')
            xlim([timeaxis(1) timeaxis(end)+minutes(1)])
            ylabel('\epsilon_P (kW)')
            title('Average Time-of-Day Error in Power')

% -------------- Hourly Histogram: Day Time --------------

    DayFig   = figure;
    BinWidth = 20;
    
        for i = 13:24

            h(i-12) = subplot(1,12,i-12);

            [counts,edges,~] = histcounts(nonzeros(Ep.dist(i,:)),'BinWidth',BinWidth);
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
        
        AxHandle=axes(DayFig,'visible','off'); 
            AxHandle.Title.Visible   ='on';
            AxHandle.XLabel.Visible  ='on';
            AxHandle.XLabel.Position = [0.5 0];
            xlabel(AxHandle,'Counts');
            AxHandle.YLabel.Visible  ='on';
            AxHandle.YLabel.Position = [0 0.5];
            ylabel(AxHandle,'\epsilon_P (kW)');
            title(AxHandle,'\epsilon_P Distribution by Hour');
        
        for i = 1:12
    
            set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])
    
        end

% -------------- Hourly Histogram: Night Time --------------

    NightFig = figure;
    BinWidth = 20;

    for i = 1:12

        h(i) = subplot(1,12,i);

        [counts,edges,~] = histcounts(nonzeros(Ep.dist(i,:)),'BinWidth',BinWidth);
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
    
    AxHandle=axes(NightFig,'visible','off'); 
        AxHandle.Title.Visible   ='on';
        AxHandle.XLabel.Visible  ='on';
        AxHandle.XLabel.Position = [0.5 0];
        xlabel(AxHandle,'Counts');
        AxHandle.YLabel.Visible  ='on';
        AxHandle.YLabel.Position = [0 0.5];
        ylabel(AxHandle,'\epsilon_P (kW)');
        title(AxHandle,'\epsilon_P Distribution by Hour');
    
    for i = 1:12

        set(h(i), 'Position', [(1/13*(i))-0.02 .09 1/15 .825])

    end

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