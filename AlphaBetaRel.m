function [AlphaBeta,AB] = AlphaBetaRel(D,PLFull,PDFs,WindBins,T,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

AB.DSrange = N.DS(2):-N.s:N.DS(1);
AB.SSrange = N.SS(1):N.s:N.SS(2);

[~,~,AlphaBeta.BinIndex] = histcounts(D.Shear(T.HubRow,:),WindBins);
AlphaBeta.BinIndex = AlphaBeta.BinIndex+1;

PowerBinAvg = zeros(1,length(AlphaBeta.BinIndex));

for i = 1:length(AlphaBeta.BinIndex)

    PowerBinAvg(i) = mean(nonzeros((AlphaBeta.BinIndex == (AlphaBeta.BinIndex(i))) .* D.Power));

end

for j = 1:length(AB.DSrange)

    for i = 1:length(AB.SSrange)

        AlphaBeta.Power(j,i) = mean(nonzeros((PLFull.alpha >= AB.SSrange(i) & PLFull.alpha < (AB.SSrange(i) + N.s)) .* (PDFs.DSrate >= AB.DSrange(j) & PDFs.DSrate < (AB.DSrange(j) + N.s)) .* D.Power ./ PowerBinAvg));

        AlphaBeta.Counts(j,i)   = sum((PLFull.alpha >= AB.SSrange(i) & PLFull.alpha < (AB.SSrange(i) + N.s)) .* (PDFs.DSrate >= AB.DSrange(j) & PDFs.DSrate < (AB.DSrange(j) + N.s)));

    end

end

AlphaBeta.Power(AlphaBeta.Counts < 15) = NaN;

AlphaBeta.PowerLH = AlphaBeta.Power;

AB.DS = N.DS;
AB.SS = N.SS;
AB.s  = N.s;

AB.xlow  = AB.SSrange(1) - AB.s/2;                                                % Set x axis minimum
AB.xhigh = AB.SSrange(end) - AB.s/2;                                              % Set x axis maximum

%% Plots - Full Colormap

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

    for i = 1:length(AB.SSrange)                                               % Vertical brid lines
        xline(AB.SSrange(i)-AB.s/2)
    end

    for i = 1:length(AB.DSrange)                                               % Horizontal brid lines
        yline(AB.DSrange(i)-AB.s/2)
    end

    x = [AB.SSrange 0.8];                                                      % Add threshold line
    y = 2/3*x - 0.1;
    plot(x,y,'color','r','LineWidth',2)

    hold off

%% Plots - Simple Above/Below Case

low  = min(AlphaBeta.Power,[],'all');                                              % Find minimum
high = max(AlphaBeta.Power,[],'all');                                              % Find maximum

AlphaBeta.PowerLH(AlphaBeta.Power<1 & ~isnan(AlphaBeta.Power)) = low;                              % Convert all values >1 to maximum
AlphaBeta.PowerLH(AlphaBeta.Power>=1) = high;                                               % Convert all values <1 to minimum

if abs(high-1) > abs(1-low)                                                 % Calibrate colormap values
    AlphaBeta.PowerLH(end,end) = 1-abs(high-1);
else
    AlphaBeta.PowerLH(end,end) = 1+abs(low-1);
end

figure;
    h = imagesc(AB.SSrange,AB.DSrange,AlphaBeta.PowerLH);
    axis xy                                                                 % Flip image
    colormap([0.66 0.66 0.66; .33 .33 .33])                                 % Bi-tone colomap
    colorbar
    set(h,'alphadata',~isnan(AlphaBeta.PowerLH))                                     % Set NaNs to white

    xlabel('Speed Shear (\alpha)')
    ylabel('Direction Shear (\circ m^{-1})')
    ylabel(colorbar('eastoutside'),'Normalized Power (-)')
    xlim([AB.xlow AB.xhigh])
    hold on

    for i = 1:length(AB.SSrange)                                               % Vertical brid lines
        xline(AB.SSrange(i)-AB.s/2)
    end

    for i = 1:length(AB.DSrange)                                               % Horizontal brid lines
        yline(AB.DSrange(i)-AB.s/2)
    end

    x = [AB.SSrange 0.8];                                                      % Add threshold line
    y = 2/3*x - 0.1;
    plot(x,y,'color','r','LineWidth',2)

    hold off

end