function [Ep] = ErrorTimeSeries(D,WindBins,Mean,T)
%ErrorTimeSeries Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n------Error Calcs-------')
fprintf('\n------------------------\n')

[~,~,~,H,MN,~] = datevec(D.Time);

%% Calculations 

% ----------------- Error: Ep -----------------

    for i = 1:size(D.Shear,2)
    
        BinNum = find(D.Shear(T.HubRow,i) <= WindBins,1,'first');
%         BinNum = find(WindBins <= D.Shear(T.HubRow,1),1,'last');
    
        Ep.abs(i)  = abs(D.Power(i)/2100 - Mean.All(BinNum))/Mean.All(BinNum);
        Ep.diff(i) = (D.Power(i)/2100 - Mean.All(BinNum))/Mean.All(BinNum);

    end

% ----------------- Hourly Ep Distributions -----------------

    for i = 1:24
    
        Indices = (H == i-1);
    
        Ep.dist(i,:) = Indices .* Ep.abs;
        Ep.mean(i)   = mean(nonzeros(Ep.dist(i,:)));
    
    end

% ----------------- Time of Day 1 - Min Ep Average -----------------
    
    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);
        
%             Ep.TOD(index,:) = Indices .* Ep.abs;
            Ep.TODMagMean(index)  = mean(nonzeros(Indices .* Ep.abs));
            Ep.TODDiffMean(index) = mean(nonzeros(Indices .* Ep.diff));
    
            j     = j + 1;
    
            index = index + 1;
    
        end
    
        i = i + 1;
        j = 0;
    
    end
    
Ep.TODAxis = duration(minutes(linspace(0,1439,length(Ep.TODMagMean))),'Format','hh:mm');

fprintf('\n\nComplete.\n\n')

end