function [WindBins,Mean,Num,Dist,STD,Sig] = CondAvgs(T,Shear,Power,AvgComp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

WindBins   = linspace(0,20,T.numbins);                                          % x-axis wind speed bins
BinWidth   = (WindBins(2) - WindBins(1))/2;                                     % Find width of bins

for i = 1:length(WindBins)
    % Find indices of values in each bin
        TempIndices  = (Shear(T.HubRow,:) > WindBins(i) - BinWidth & Shear(T.HubRow,:) <= WindBins(i) + BinWidth);
%         TempIndices2  = (I > 2.5 & I < 7.5);
%         TempIndices   = TempIndices1 .* TempIndices2;
    % Compute averages for each category in each bin
        Mean.All(i)    = mean(nonzeros(TempIndices .* Power));
        Mean.Larger(i) = mean(nonzeros(AvgComp.Greater .* TempIndices .* Power));
        Mean.Less(i)   = mean(nonzeros(AvgComp.Less .* TempIndices .* Power));

    % Compute descriptive stats
        Num.All(i)     = length(nonzeros(TempIndices .* Power));
        Num.Larger(i)  = length(nonzeros(AvgComp.Greater .* TempIndices .* Power));
        Num.Less(i)    = length(nonzeros(AvgComp.Less .* TempIndices .* Power));

        Dist.Larger(i,:)= AvgComp.Greater .* TempIndices .* Power;
        Dist.Less(i,:)  = AvgComp.Less .* TempIndices .* Power;
        Dist.All(i,:)   = TempIndices .* Power;

        STD.All(i)     = std(nonzeros(TempIndices .* Power));
        STD.Larger(i)  = std(nonzeros(AvgComp.Greater .* TempIndices .* Power));
        STD.Less(i)    = std(nonzeros(AvgComp.Less .* TempIndices .* Power));
        
    % Compute inferential stats
        Sig.Large(i)   = ttest2(nonzeros(AvgComp.Greater .* TempIndices .* Power),nonzeros(TempIndices .* Power));
        Sig.Less(i)    = ttest2(nonzeros(AvgComp.Less .* TempIndices .* Power),nonzeros(TempIndices .* Power));
        Sig.Both(i)    = ttest2(nonzeros(AvgComp.Less .* TempIndices .* Power),nonzeros(AvgComp.Greater .* TempIndices .* Power));
end

% for i = 1:length(AlphaBins)
%         TempIndicesA = (alpha > AlphaBins(i) - AlphaBinWidth & alpha <= AlphaBins(i) + AlphaBinWidth);
%         Mean.Alpha(i)  = mean(nonzeros(TempIndicesA .* Power));
% end