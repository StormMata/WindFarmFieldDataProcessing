function [PDFs] = ShearCharacterization(Veer,T,PLFull,PLInflec,PlotCommand)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Find Direction Shear PDF

for i = 1:size(Veer,2)

    Veer(Veer > 180) = Veer(Veer > 180) - 360;                              % Convert angles >180 to negative angles    [deg]
    
    VeerI = interp1(flip(T.Heights)',Veer(:,i),43:1:169);                   % Interpolate direction shear profile       [deg]
    
    PDFs.DSrate(i) = (VeerI(T.Hub+T.R-43)-VeerI(T.Hub-T.R-43))/(T.R*2);     % Rate of direction shear                   [deg/m]

end

%% Find Direction and Speed Shear Evolution



%% Plots

if strcmp(PlotCommand,'PlotOn') == 1

% Direction Shear
    
    figure;                                                                 % 100 bins
        [PDFs.DSProb1,PDFs.DSEdges1] = histcounts(PDFs.DSrate,110,...
            'Normalization','probability','BinLimits',[-0.6 0.8]);
        plot(PDFs.DSEdges1(1:end-1),PDFs.DSProb1,'Color','k','LineWidth',1.5);
        title('Direction Shear Probability Distribution')
        xlabel('Direction Shear (deg m^{-1})')
        ylabel('Probability of Occurrence')
        xlim([-0.6 0.8])
    
    figure;                                                                 % 350 bins
        [PDFs.DSProb2,PDFs.DSEdges2] = histcounts(PDFs.DSrate,350,...
            'Normalization','probability','BinLimits',[-0.6 0.8]);
        plot(PDFs.DSEdges2(1:end-1),PDFs.DSProb2,'Color','k','LineWidth',1.5);
        title('Direction Shear Probability Distribution')
        xlabel('Direction Shear (deg m^{-1})')
        ylabel('Probability of Occurrence')
        xlim([-0.6 0.8])

% Speed Shear
    
    figure;                                                                 % Full profile fit
        [PDFs.SSProb1,PDFs.SSEdges1] = histcounts(PLFull.alpha,...
            'Normalization','probability','BinLimits',[-0.75 1.5]);
        plot(PDFs.SSEdges1(1:end-1),PDFs.SSProb1,'Color','k','LineWidth',1.5);
        title('Speed Shear Probability Distribution - Full Fit')
        xlabel('Direction Shear (\alpha)')
        ylabel('Probability of Occurrence')
        xlim([-0.75 1.5])

    figure;                                                                 % Inflection profile fit
        [PDFs.SSProb2,PDFs.SSEdges2] = histcounts(PLInflec.alpha,...
            'Normalization','probability','BinLimits',[-0.5 1.5]);
        plot(PDFs.SSEdges2(1:end-1),PDFs.SSProb2,'Color','k','LineWidth',1.5);
        title('Speed Shear Probability Distribution - Inflection Fit')
        xlabel('Direction Shear (\alpha)')
        ylabel('Probability of Occurrence')
        xlim([-0.35 1.4])

end

end