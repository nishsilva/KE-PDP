%% Kuroshio indices using Copernicus data and KOCR region
load /projectnb/pdpanalysis/data/Cop_means.mat;
Cop_SST = reshape(Cop_means, 2400*1800*7, 38);
sst_anom = detrend(Cop_SST')'; 						%calculate the anomalies
data_anom = reshape(sst_anom, 2400*1800, 7, 38);
data_anom = data_anom(:,1:6,:);

data_anom = reshape(data_anom, 2400, 1800, 6*38);
data_anom_nan = data_anom;
data_anom_nan(isnan(data_anom_nan)) = 0;     		%remove NaNs

large = imboxfilt(data_anom_nan, 101);              %largescale SSTAs

large_2 = imboxfilt(data_anom_nan, 51);             %filter to capture mesoscale <100km

meso = minus(data_anom_nan,large_2);					%mesoscale SSTAs

%%
lats = ncread('/projectnb/pdpanalysis/data/Copernicus_SST_Jan2021/East/1981_E.nc', 'lat');
lons = ncread('/project/pdpanalysis/nish/data/Copernicus_SST/2000.nc', 'lon');


%%
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 8 12]);
%tiledlayout(1,3,'TileSpacing','Compact', 'Padding', 'loose');

subplot('position', [0.1 0.7 0.8 0.295]);
%nexttile
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 15,...
          'MLabelLocation',25);

limit = 2.5;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-limit limit]); % set the limits of the colorbar

surfm(double(lats), double(lons), data_anom(:,:,134)','Facecolor', 'interp');
%colorbar('westoutside')
load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth',1.5)

textm(60, 122, '(a)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'Unfiltered SSTAs', 'FontSize',20, 'FontWeight', 'bold',...
%    'HorizontalAlignment', 'center');
   
title('Unfiltered SSTAs', 'FontSize',17, 'FontWeight', 'bold')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
tightmap
%hold on;

%la = [34 34 44 44 34]; 
%lo = [142 182 182 142 142];

%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','r', 'Linewidth', 2)


subplot('position', [0.1 0.4 0.8 0.295]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 15,...
          'MLabelLocation',25);

%limit = 3;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-limit limit]); % set the limits of the colorbar

surfm(double(lats), double(lons), large(:,:,134)','Facecolor', 'interp');
%colorbar('southoutside');
load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth',1.5)

textm(60, 122, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'Large-scale SSTAs', 'FontSize',20, 'FontWeight', 'bold',...
%    'HorizontalAlignment', 'center');
   
%title('Correlation b/w SSTAs')    
title('Large-scale SSTAs', 'FontSize',17, 'FontWeight', 'bold')   
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
tightmap
%hold on;

%la = [30 42 42 42 30]; 
%lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','r', 'Linewidth', 2)



subplot('position', [0.1 0.1 0.8 0.295]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 15,...
          'MLabelLocation',25);

%limit = 3;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
cmap = colormap(cmap_mod);
caxis([-limit limit]); % set the limits of the colorbar

surfm(double(lats), double(lons), meso(:,:,134)','Facecolor', 'interp');

load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth',1.5)

textm(60, 122, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'Mesoscale SSTAs', 'FontSize',20, 'FontWeight', 'bold',...
%    'HorizontalAlignment', 'center');
   
%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')   
title('Mesoscale SSTAs', 'FontSize',17, 'FontWeight', 'bold') 
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
tightmap
%subplot('position',[0.1 0 0.8 0.08])

%limit = 3;
caxis([-limit limit]);
h = colorbar('Location', 'Southoutside');
set(h, 'position',[0.1 0.05 0.8 0.012] );
cbfreeze(h);
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')
h.FontWeight = 'bold';
h.FontSize = 16;
%xlabel(h, 'SSTA (K)')
%cb.Layout.Tile = 'south';
%cb.Label.String = 'SSTA (K)';

%print('Nish_092921_KE_PDP_MS_Fig_1_PS','-dpng', '-r500');
%print('Nish_072522_KE_PDP_MS_Fig_1_PS_revised','-dpng', '-r500');