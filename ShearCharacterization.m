function [PDFs] = ShearCharacterization(D,T,PLFull,PLInflec,Ekman)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fprintf('\n------------------------')
fprintf('\n------Shear Calcs-------')
fprintf('\n------------------------\n')

%% Find Direction Shear PDF

for i = 1:size(D.Veer,2)

    D.Veer(D.Veer > 180) = D.Veer(D.Veer > 180) - 360;                          % Convert angles >180 to negative angles    [deg]
    
    VeerI = interp1(flip(T.Heights)',D.Veer(:,i),43:1:169);                     % Interpolate direction shear profile       [deg]
    
    PDFs.DSrate(i) = (VeerI(T.Hub+T.R-43)-VeerI(T.Hub-T.R-43))/(T.R*2);         % Rate of direction shear                   [deg/m]

end

%% Remove NaNs from Fit Data

PLInflec.alpha(isnan(PLInflec.alpha)) = 0;
Ekman.K(isnan(Ekman.K))               = 0;


%% Find Direction Shear Evolution by TOD

    [~,~,~,H,MN,~] = datevec(D.Time);

    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);

            PDFs.DSTOD(index) = mean(nonzeros(Indices .* PDFs.DSrate));
    
            j     = j + 1;
    
            index = index + 1;
    
        end
    
        i = i + 1;
        j = 0;
    
    end

%% Find Speed Shear Evolution by TOD

    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);

            PDFs.SSTODFull(index)   = mean(nonzeros(Indices .* PLFull.alpha));
            PDFs.SSTODInflec(index) = mean(nonzeros(Indices .* PLInflec.alpha));
    
            j     = j + 1;
    
            index = index + 1;
    
        end
    
        i = i + 1;
        j = 0;
    
    end

%% Find Ekman Parameter Evolution by TOD

    index = 1;
    i     = 0;
    j     = 0;
    
    while i <= 23
    
        while j <= 59
    
            Indices = (H == i & MN == j);
        
            PDFs.EkTODmedian(index) = median(nonzeros(Indices .* Ekman.K));
            PDFs.EkTODmean(index)   = mean(nonzeros(Indices .* Ekman.K));
    
            j     = j + 1;
    
            index = index + 1;
    
        end
    
        i = i + 1;
        j = 0;
    
    end

%% Computer Time Averages

PDFs.SSTODAvgFull   = mean(reshape(PDFs.SSTODFull,10,[]));
PDFs.SSTODAvgInflec = mean(reshape(PDFs.SSTODInflec,10,[]));
PDFs.EkTODMed       = mean(reshape(PDFs.EkTODmedian,10,[]));
PDFs.EkTODAvg       = mean(reshape(PDFs.EkTODmean,10,[]));
PDFs.DSTODAvg       = mean(reshape(PDFs.DSTOD,10,[]));

PDFs.TODAxis        = duration(minutes(linspace(0,1439,length(PDFs.DSTOD))),'Format','hh:mm');
PDFs.TODAvgAxis     = duration(minutes(linspace(0,1439,length(PDFs.DSTODAvg))),'Format','hh:mm');

%% Find Direction and Speed Shear Evolution by Height


% Heights = flip(T.Heights)';
% dudz  = gradient(D.Shear(:,160))./gradient(Heights);                % Calculate sign of shear profile
%             index = find(dudz < 0, 1, 'last')

fprintf('\n\nComplete.\n\n')


end