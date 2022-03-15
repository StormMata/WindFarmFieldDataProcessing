function [AlphaBeta,AB] = AlphaBetaRel(D,PLFull,PDFs,WindBins,T,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n----Alpha Beta Calcs----')
fprintf('\n------------------------\n')

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

        AlphaBeta.Power(j,i)  = mean(nonzeros((PLFull.alpha >= AB.SSrange(i) & PLFull.alpha < (AB.SSrange(i) + N.s)) .* (PDFs.DSrate >= AB.DSrange(j) & PDFs.DSrate < (AB.DSrange(j) + N.s)) .* D.Power ./ PowerBinAvg));

        AlphaBeta.Counts(j,i) = sum((PLFull.alpha >= AB.SSrange(i) & PLFull.alpha < (AB.SSrange(i) + N.s)) .* (PDFs.DSrate >= AB.DSrange(j) & PDFs.DSrate < (AB.DSrange(j) + N.s)));

    end

end

%% Time Distributions of Three Bins of Interest

UpLeft = nonzeros((PLFull.alpha >= -0.1 & PLFull.alpha < (-0.1 + 0.1)) .* (PDFs.DSrate >= 0.4 & PDFs.DSrate < (0.4 + 0.1)) .* D.Time);

LowLeft = nonzeros((PLFull.alpha >= -0.1 & PLFull.alpha < (-0.1 + 0.1)) .* (PDFs.DSrate >= -0.1 & PDFs.DSrate < (-0.1 + 0.1)) .* D.Time);

LowRight = nonzeros((PLFull.alpha >= 0.7 & PLFull.alpha < (0.7 + 0.1)) .* (PDFs.DSrate >= -0.1 & PDFs.DSrate < (-0.1 + 0.1)) .* D.Time);

[~,~,~,AlphaBeta.UpLeft,~,~] = datevec(UpLeft);

[~,~,~,AlphaBeta.LowLeft,~,~] = datevec(LowLeft);

[~,~,~,AlphaBeta.LowRight,~,~] = datevec(LowRight);

figure;
    histogram(AlphaBeta.UpLeft)
    xlabel('Hour of the Day')
    ylabel('Counts')
    title('Low Speed Shear - High Veer')

figure;
    histogram(AlphaBeta.LowLeft)
    xlabel('Hour of the Day')
    ylabel('Counts')
    title('Low Speed Shear - Low Backing/Low Veer')

figure;
    histogram(AlphaBeta.LowRight)
    xlabel('Hour of the Day')
    ylabel('Counts')
    title('High Speed Shear - Low Backing/Low Veer')

%% Processing

AlphaBeta.Power(AlphaBeta.Counts < 15) = NaN;

AlphaBeta.PowerLH = AlphaBeta.Power;

AB.DS = N.DS;
AB.SS = N.SS;
AB.s  = N.s;

AB.xlow  = AB.SSrange(1) - AB.s/2;                                                % Set x axis minimum
AB.xhigh = AB.SSrange(end) - AB.s/2;                                              % Set x axis maximum

fprintf('\n\nComplete.\n\n')

end