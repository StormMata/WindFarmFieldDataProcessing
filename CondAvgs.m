function [WindBins,Mean,Num,Dist,STD,Sig] = CondAvgs(T,D,AvgComp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n--Conditional Averages--')
fprintf('\n------------------------\n')

WindBins   = linspace(0,20,T.numbins);                                          % x-axis wind speed bins
BinWidth   = (WindBins(2) - WindBins(1))/2;                                     % Find width of bins

for i = 1:length(WindBins)
    % Find indices of values in each bin
        TempIndices  = (D.Shear(T.HubRow,:) > WindBins(i) - BinWidth & D.Shear(T.HubRow,:) <= WindBins(i) + BinWidth);
%         TempIndices2  = (I > 2.5 & I < 7.5);
%         TempIndices   = TempIndices1 .* TempIndices2;

    % Compute averages for each category in each bin
        Mean.All(i)    = median(nonzeros(TempIndices .* D.Power))/2100;
        Mean.Larger(i) = median(nonzeros(AvgComp.Greater .* TempIndices .* D.Power))/2100;
        Mean.Less(i)   = median(nonzeros(AvgComp.Less .* TempIndices .* D.Power))/2100;

    % Compute descriptive stats
        Num.All(i)     = length(nonzeros(TempIndices .* D.Power));
        Num.Larger(i)  = length(nonzeros(AvgComp.Greater .* TempIndices .* D.Power));
        Num.Less(i)    = length(nonzeros(AvgComp.Less .* TempIndices .* D.Power));

        Dist.Larger(i,:)= AvgComp.Greater .* TempIndices .* D.Power;
        Dist.Less(i,:)  = AvgComp.Less .* TempIndices .* D.Power;
        Dist.All(i,:)   = TempIndices .* D.Power;

        STD.All(i)     = std(nonzeros(TempIndices .* D.Power/2100));
        STD.Larger(i)  = std(nonzeros(AvgComp.Greater .* TempIndices .* D.Power));
        STD.Less(i)    = std(nonzeros(AvgComp.Less .* TempIndices .* D.Power));
        
    % Compute inferential stats
        Sig.Large(i)   = ttest2(nonzeros(AvgComp.Greater .* TempIndices .* D.Power),nonzeros(TempIndices .* D.Power));
        Sig.Less(i)    = ttest2(nonzeros(AvgComp.Less .* TempIndices .* D.Power),nonzeros(TempIndices .* D.Power));
        Sig.Both(i)    = ttest2(nonzeros(AvgComp.Less .* TempIndices .* D.Power),nonzeros(AvgComp.Greater .* TempIndices .* D.Power));
end

% for i = 1:length(AlphaBins)
%         TempIndicesA = (alpha > AlphaBins(i) - AlphaBinWidth & alpha <= AlphaBins(i) + AlphaBinWidth);
%         Mean.Alpha(i)  = mean(nonzeros(TempIndicesA .* D.Power));
% end

fprintf('\n\nComplete.\n\n')

end