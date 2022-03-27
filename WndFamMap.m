function [WndFamMap] = WndFamMap()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

TimeStart = tic;

fprintf('\n------------------------')
fprintf('\n--Wind Farm Elevations--')
fprintf('\n------------------------\n')
fprintf('\nComplete:             0')

lat  = linspace(23.0420,22.9115,350);
long = linspace(69.8665,70.0040,350);

Elev = zeros(length(lat),length(long));

for i = 1:length(lat)

    for j = 1:length(long)

        Elev(i,j) = elevation(txsite('Latitude',lat(i),'Longitude',long(j)));

    end

    if mod(i,10)==0
        p = i/size(Elev,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % print progress to screen
    end

end

WndFamMap.lat  = lat;                                                       % Store latitudes in structure
WndFamMap.long = long;                                                      % Store longitudes in structure
WndFamMap.El   = Elev;                                                      % Store elevations in structure

TimeEnd = toc(TimeStart);

fprintf('\n\nTime: %10d minutes \n%16.0f seconds\n', floor(TimeEnd/60), rem(TimeEnd,60));

end