lat  = linspace(23.01736,23.00000,400);
long = linspace(69.91132,69.93274,400);

Elev = zeros(length(lat),length(long));

fprintf('\nComplete:             0')

for i = 1:length(lat)

    for j = 1:length(long)

        Elev(i,j) = elevation(txsite('Latitude',lat(i),'Longitude',long(j)));

    end

    if mod(i,10)==0
        p = i/size(Elev,2)*100;
        fprintf('\b\b\b\b\b\b\b%6.1f%%',p)                                  % print progress to screen
    end

end

figure; imagesc(long,lat,Elev);axis xy
colormap summer;
hold on; contour(long,lat,Elev,'showtext','on','linecolor','#525252')
hold on; plot(69.915731,23.014843,'marker','o','markersize',12,'color','k','markerfacecolor','#fc4903')
hold on; plot(69.921821,23.012679,'marker','o','markersize',12,'color','k','markerfacecolor','#fc4903')
hold on; plot(69.927336,23.009665,'marker','o','markersize',12,'color','k','markerfacecolor','#fc4903')
hold on; plot(69.930022,23.006146,'marker','o','markersize',12,'color','k','markerfacecolor','#fc4903')
hold on; plot(69.918603,23.003806,'marker','o','markersize',12,'color','k','markerfacecolor','#fc4903')
hold on; plot(69.931198,23.001170,'marker','o','markersize',12,'color','k','markerfacecolor','#fc4903')
rectangle('Position',[69.925,23.003,0.00245,0.0002],'FaceColor','#a7a7a7','EdgeColor','#a7a7a7','LineWidth',3)
rectangle('Position',[69.92745,23.003,0.00245,0.0002],'FaceColor','#c5c5c5','EdgeColor','#c5c5c5','LineWidth',3)
text(69.925,23.0035,'0','HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
text(69.92745,23.0035,'0.25','HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
text(69.93,23.0035,'0.5 km','HorizontalAlignment','center','FontSize',10,'FontWeight','bold')

centerx = 69.93;
centery = 23.015;

x = [centerx - 0.0003, centerx, centerx];
y = [centery - 0.0003, centery, centery + 0.0006];
c = [0 0 0];
fill(x, y, c)

x = [centerx, centerx, centerx + 0.0003];
y = [centery, centery + 0.0006, centery - 0.0003];
c = [1 1 1];
fill(x, y, c)

text(centerx,centery + 0.00085,'N','HorizontalAlignment','center','FontSize',13,'FontWeight','bold')

centerx = 69.91715;
centery = 23.0155;
side = 0.0001;

x = [centerx + side, centerx + side, centerx - side, centerx - side];
y = [centery - side, centery + side, centery + side, centery - side];
c = [0 0 0];
fill(x, y, c)

centerx = 69.91715;
centery = 23.0155;
side = 0.00015;

x = [centerx, centerx - side, centerx + side];
y = [centery + side, centery + 3*side, centery + 3*side];
c = [0 0 0];
fill(x, y, c)

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
axis equal
% title('Wind Turbine Array')
ylabel(colorbar('eastoutside'),'Elevation (m)')
