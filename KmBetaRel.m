function [KmBeta,KB] = KmBetaRel(D,PLFull,PDFs,WindBins,T,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n-----Km Beta Calcs------')
fprintf('\n------------------------\n')

KB.DSrange = N.DS(2):-N.s:N.DS(1);
KB.Kmrange = N.Km(1):N.s:N.Km(2);

[~,~,KmBeta.BinIndex] = histcounts(D.Shear(T.HubRow,:),WindBins);
KmBeta.BinIndex = KmBeta.BinIndex+1;

PowerBinAvg = zeros(1,length(KmBeta.BinIndex));

for i = 1:length(KmBeta.BinIndex)

    PowerBinAvg(i) = mean(nonzeros((KmBeta.BinIndex == (KmBeta.BinIndex(i))) .* D.Power));

end

for j = 1:length(KB.DSrange)

    for i = 1:length(KB.Kmrange)

        KmBeta.Power(j,i)  = mean(nonzeros((PLFull.alpha >= KB.Kmrange(i) & PLFull.alpha < (KB.Kmrange(i) + N.s)) .* (PDFs.DSrate >= KB.DSrange(j) & PDFs.DSrate < (KB.DSrange(j) + N.s)) .* D.Power ./ PowerBinAvg));

        KmBeta.Counts(j,i) = sum((PLFull.alpha >= KB.Kmrange(i) & PLFull.alpha < (KB.Kmrange(i) + N.s)) .* (PDFs.DSrate >= KB.DSrange(j) & PDFs.DSrate < (KB.DSrange(j) + N.s)));

    end

end

%% Time Distributions of Three Bins of Interest

% UpLeft = nonzeros((PLFull.alpha >= -0.1 & PLFull.alpha < (-0.1 + 0.1)) .* (PDFs.DSrate >= 0.4 & PDFs.DSrate < (0.4 + 0.1)) .* D.Time);
% 
% LowLeft = nonzeros((PLFull.alpha >= -0.1 & PLFull.alpha < (-0.1 + 0.1)) .* (PDFs.DSrate >= -0.1 & PDFs.DSrate < (-0.1 + 0.1)) .* D.Time);
% 
% LowRight = nonzeros((PLFull.alpha >= 0.7 & PLFull.alpha < (0.7 + 0.1)) .* (PDFs.DSrate >= -0.1 & PDFs.DSrate < (-0.1 + 0.1)) .* D.Time);
% 
% [~,~,~,KmBeta.UpLeft,~,~] = datevec(UpLeft);
% 
% [~,~,~,KmBeta.LowLeft,~,~] = datevec(LowLeft);
% 
% [~,~,~,KmBeta.LowRight,~,~] = datevec(LowRight);
% 
% figure;
%     histogram(KmBeta.UpLeft)
%     xlabel('Hour of the Day')
%     ylabel('Counts')
%     title('Low Speed Shear - High Veer')
% 
% figure;
%     histogram(KmBeta.LowLeft)
%     xlabel('Hour of the Day')
%     ylabel('Counts')
%     title('Low Speed Shear - Low Backing/Low Veer')
% 
% figure;
%     histogram(KmBeta.LowRight)
%     xlabel('Hour of the Day')
%     ylabel('Counts')
%     title('High Speed Shear - Low Backing/Low Veer')

%% Processing

KmBeta.Power(KmBeta.Counts < 15) = NaN;

KmBeta.PowerLH = KmBeta.Power;

KB.DS = N.DS;
KB.Km = N.Km;
KB.s  = N.s;

KB.xlow  = KB.Kmrange(1) - KB.s/2;                                                % Set x axis minimum
KB.xhigh = KB.Kmrange(end) - KB.s/2;                                              % Set x axis maximum

fprintf('\n\nComplete.\n\n')

end