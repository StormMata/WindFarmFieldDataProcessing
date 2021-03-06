z = 43:1:200;

i = 161;

Uold = D.Shear(12,i) .* (z./43).^PLFullold.alpha(i);

Unew = PLFull.Uref(i) .* (z./43).^PLFull.alpha(i);

Uinflec = PLInflec.Uref(i) .* (z./43).^PLInflec.alpha(i);

plot(D.Shear(:,i),flip(T.Heights),'LineWidth',1.5)
hold on
plot(Uold,z,'LineWidth',1.5)
plot(Unew,z,'LineWidth',1.5)
plot(Uinflec,z,'LineWidth',1.5)
yline(T.Hub-T.R,'LineStyle','--')
yline(T.Hub+T.R,'LineStyle','--')

str1 = sprintf('Full Profile Toolbox; R^2 = %2.2f',PLFullold.R(i))
str2 = sprintf('Full Profile Optimization; R^2 = %2.2f',PLFull.R(i))
str3 = sprintf('Partial Profile Optimization; R^2 = %2.2f',PLInflec.R(i))

legend('Shear Profile',str1,str2,str3,'','')

% Full - R2

    histogram(PLFull.R,'binedges',[0:0.05:1],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    hold on
    histogram(PLFullold.R,'binedges',[-1.5:0.05:1],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    legend('Full Profile Optimization','Full Profile Toolbox')
    xlim([-1.5 1])
    xlabel('Coefficient of Determination')
    ylabel('Probability of Occurrence')

% Full - RMSE

    Normer = mean(D.Shear(:,:));

    histogram(PLFull.NRMSE,'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    hold on
    histogram(PLFullold.RMSE./Normer,'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    legend('Full Profile Optimization','Full Profile Toolbox')
%     xlim([-1.5 1])
    xlabel('Coefficient of Determination')
    ylabel('Probability of Occurrence')

% Inflec

    histogram(PLInflec.R,'binedges',[0:0.05:1],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    hold on
    histogram(PLInflecold.R,'binedges',[0:0.05:1],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    legend('Partial Profile Optimization','Partial Profile Toolbox')
    xlim([0 1])
    xlabel('Coefficient of Determination')
    ylabel('Probability of Occurrence')

% Ekman - R2

    histogram(Ekman.NRMSE,'displaystyle','stairs','Normalization','probability',LineWidth=1.5)


    histogram(Ekmanold.R,'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    set(gca,'YScale','log')
    legend('Ekman Toolbox')
    xlabel('Coefficient of Determination')
    ylabel('Probability of Occurrence')
    xlim([-55 1])
    
    y = histogram(Ekman.R,'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
%     legend('Ekman Optimization')
%     xlabel('Coefficient of Determination')
%     ylabel('Probability of Occurrence')
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);

    xlim([-0.05 1])
    hold on

    %left side 
    xs = [-0.05 0  0 0.05 0.05 -0.05 ];
    ys = [y.Values(1) y.Values(1) y.Values(2) y.Values(2) 0 0];
    fill(xs,ys,[0.6350 0.0780 0.1840])

    %right side 
    xs = [0.9 0.95 0.95 1 1 0.9];
    ys = [y.Values(20) y.Values(20) y.Values(21) y.Values(21) 0 0];
    fill(xs,ys,[0.6350 0.0780 0.1840])

% Ekman - RMSE

    histogram(Ekman.NRMSE,'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    legend('Ekman Optimization')
    xlabel('NRMSE')
    ylabel('Probability of Occurrence')
    xlim([-55 1])
    
%     EkmanR = Ekman.R;
    [~,~,~,EkmanTime,~,~] = datevec(D.Time);
    
    filter = (Ekman.NRMSE>=0.0394);
    
    filter = double(filter);
    
    filter(filter == 0) = NaN;
    
    EkmanTime = EkmanTime.*filter;
    
    histogram(EkmanTime,'binedges',[0:1:24],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    title('Hours of the Day for Highest Ekman Fit Normalized RMSE Values')
    xlabel('Hour of the Day')
    ylabel('Probability of Occurrence')
    xlim([0 24])

% Nightime
    
%     EkmanR = Ekman.R;
    [~,~,~,EkmanTime,~,~] = datevec(D.Time);
    
    filter = (Ekman.R<=0.05);
    
    filter = double(filter);
    
    filter(filter == 0) = NaN;
    
%     EkmanR = EkmanR.*filter;
    
    EkmanTime = EkmanTime.*filter;
    
    histogram(EkmanTime,'binedges',[0:1:24],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    title('Hours of the Day for Lowest Ekman Fit R^2 Values')
    xlabel('Hour of the Day')
    ylabel('Probability of Occurrence')
    xlim([0 24])

% Daytime

    EkmanR = Ekman.R;
    [~,~,~,EkmanTime,~,~] = datevec(D.Time);
    
    filter = (Ekman.R>=0.90);
    
    filter = double(filter);
    
    filter(filter == 0) = NaN;
    
    EkmanR = EkmanR.*filter;
    
    EkmanTime = EkmanTime.*filter;
    
    histogram(EkmanTime,'binedges',[0:1:24],'displaystyle','stairs','Normalization','probability',LineWidth=1.5)
    title('Hours of the Day for Highest Ekman Fit R^2 Values')
    xlabel('Hour of the Day')
    ylabel('Probability of Occurrence')
    xlim([0 24])

% Ekman PL Comparison All

    histogram(Ekman.NRMSE,'binwidth',0.01,'Normalization','probability')
    hold on
    histogram(PLFull.NRMSE,'binwidth',0.01,'Normalization','probability')
    str1 = sprintf('Ekman; Median = %2.2f',median(Ekman.NRMSE));
    str2 = sprintf('Power Law; Median = %2.2f',median(PLFull.NRMSE));
    legend(str1,str2)
    xlabel('$RMSE/\overline{y}$','Interpreter','latex')
    ylabel('Probability of Occurrence')
    title('Distribution for All Times')
    xlim([0 0.5])

% Ekman PL Comparison 2 - 8 AM

    [~,~,~,EkmanTime,~,~] = datevec(D.Time);

    filter = (EkmanTime >= 2 & EkmanTime <= 8);

    filter = double(filter);

    filter(filter == 0) = NaN;

    histogram(Ekman.NRMSE.*filter,'binwidth',0.01,'Normalization','probability')
    hold on
    histogram(PLFull.NRMSE.*filter,'binwidth',0.01,'Normalization','probability')
    str1 = sprintf('Ekman; Median = %2.4f',median(Ekman.NRMSE.*filter,'omitnan'));
    str2 = sprintf('Power Law; Median = %2.4f',median(PLFull.NRMSE.*filter,'omitnan'));
    legend(str1,str2)
    xlabel('$RMSE/\overline{y}$','Interpreter','latex')
    ylabel('Probability of Occurrence')
    title('02:00 to 08:00')
    xlim([0 0.5])


% Ekman PL Comparison Hourly

    [~,~,~,EkmanTime,~,~] = datevec(D.Time);

    for i = 0:23

        index = EkmanTime == i;

        ERMSEmeans(i+1) = mean(nonzeros(index .* Ekman.NRMSE));
        PLRMSEmeans(i+1) = mean(nonzeros(index .* PLFull.NRMSE));

    end

    plot([0:23],ERMSEmeans,'LineWidth',1.5)
    hold on
    plot([0:23],PLRMSEmeans,'LineWidth',1.5)
    legend('Ekman','Power Law')
    xlabel('Hour of the Day')
    ylabel('NRMSE')
    xlim([0 23])



