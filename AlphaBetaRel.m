function [AvgPower] = AlphaBetaRel(Shear,PLFull,PDFs,Power,WindBins,T)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

step    = 0.1;

DSrange = 0.6:-step:-0.1;
SSrange = -0.1:step:0.7;

AvgPower = zeros(length(DSrange),length(SSrange));
Counts   = AvgPower;


for i = 1:length(Power)

    BinIndex(i) = find(Shear(T.HubRow,i) <= WindBins,1,'first');
%     BinIndex(i) = find(WindBins <= Shear(T.HubRow,i),1,'first');

end

for i = 1:length(BinIndex)

    PowerAvg(i) = mean(nonzeros((BinIndex == (BinIndex(i))) .* Power));

end


for j = 1:length(DSrange)

    for i = 1:length(SSrange)

        AvgPower(j,i) = mean(nonzeros((PLFull.alpha >= SSrange(i) & PLFull.alpha < (SSrange(i) + step)) .* (PDFs.DSrate >= DSrange(j) & PDFs.DSrate < (DSrange(j) + step)) .* Power ./ PowerAvg));


        
        Counts(j,i)   = sum((PLFull.alpha >= SSrange(i) & PLFull.alpha < (SSrange(i) + step)) .* (PDFs.DSrate >= DSrange(j) & PDFs.DSrate < (DSrange(j) + step)));

    end

end

AvgPower(Counts < 15) = NaN;

% Verify proper orientation
%      AvgPower(1,1)   = 5;
%     AvgPower(1,end) = 5000;
%     AvgPower(8,1)   = 5000;
%     AvgPower(6,1)   = 2500;

    Counts(1,1)   = 2500;
    Counts(1,end) = 5000;
    Counts(8,1)   = 5000;
    Counts(6,1)   = 2500;

%% Plots

figure;
%     imagesc(SSrange,flip(DSrange),AvgPower)
    imagesc(SSrange,DSrange,AvgPower)
     axis xy
    colormap([1 1 1; parula(256)])
    colorbar
    xlabel('Speed Shear (\alpha)')
    ylabel('Direction Shear (\circ m^{-1})')
    hold on
    for i = 1:length(SSrange)
        xline(SSrange(i)-step/2)
    end
    for i = 1:length(DSrange)
        yline(DSrange(i)-step/2)
    end
    x = [SSrange 0.8];
    y = 2/3*x - 0.1;
    plot(x,y,'color','r','LineWidth',2)
    hold off

% figure;
%     HM = heatmap(Counts);
%     HM.NodeChildren(3).YDir='normal';


% % Matrix(isnan(Matrix)) = 0;
% 
% % HM = heatmap(AvgPower)
% % HM.NodeChildren(3).YDir='normal';
% % h.XDisplayLabels = SSrange;
% % h.YDisplayLabels = DSrange;
% HM = imagesc(SSrange,flip(DSrange),AvgPower);
% colormap([1 1 1; parula(256)])
% % HM = imagesc(Matrix);
% % set(gca,'YDir','reverse');
% axis xy
% % YDisplayLabels = flip(DSrange);
% % ax = gca;
% % ax.YDir = 'normal'
% % set(gca,'YDir','normal') 
% % % Initialize a color map array of 256 colors.
% % colorMap = jet(256);
% % % Apply the colormap and show the colorbar
% % colormap(colorMap);
% % colorbar;

end