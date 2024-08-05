%% Produce Regression maps for KI Indices

load /projectnb/pdpanalysis/data/Cop_means.mat;
Cop_SST = reshape(Cop_means, 2400*1800*7, 38);
sst_anom = detrend(Cop_SST')';
data_anom = reshape(sst_anom, 2400*1800, 7, 38);
data_anom = data_anom(:,1:6,:);
Seasonal_SSTA = mean(data_anom, 2);
Seasonal_SSTA = reshape(Seasonal_SSTA, 2400*1800*1, 38);

clear('Cop_means', 'Cop_SST', 'sst_anom', 'data_anom');

lons = ncread('/project/pdpanalysis/nish/data/Copernicus_SST/2000.nc', 'lon');
lats = ncread('/projectnb/pdpanalysis/data/Copernicus_SST_Jan2021/East/2000_E.nc', 'lat');

load HVD_Cop_MS_EOF_1;
load HVD_Cop_MS_KI;
load HVD_Cop_LS_EOF_1;
load HVD_Cop_LS_EOF_2;

%%
%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 16 8]);
%figure(1);
%tiledlayout(2,2,'TileSpacing','tight');

KI = mean(reshape(HVD_Cop_MS_KI, 6, 38),1);
Var = Seasonal_SSTA;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 2400, 1800);
SST_reg = reshape(A_regress, 2400, 1800);


[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/38; % convert to power spectral density

rng(1);
rand_t = randn(38, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Seasonal_SSTA*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 2400, 1800);
mask = double(mask);
mask(mask == 0) = NaN;

lats_sparse = lats(1:20:1800);
lons_sparse = lons(1:20:2400);

LA = reshape(double(lats_sparse*ones([1 120])), 90, 120);
LO = reshape(double(lons_sparse*ones([1 90])), 120, 90);

mask = mask(1:20:2400, 1:20:1800);

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 15 9]);

subplot('position', [0.05 0.55 0.4 0.4]);
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 25,...
          'MLabelLocation',25,...
          'FlineWidth', 1);

limit = 0.6;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-0.6 0.6]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 0.5);

%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(60, 122, '(a)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'MS KI', 'FontSize',20, 'FontWeight', 'bold');
title('MS KI', 'FontSize',22, 'FontWeight', 'bold');
%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap

KI = mean(reshape(HVD_Cop_MS_EOF_1, 6, 38),1);
Var = Seasonal_SSTA;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 2400, 1800);
SST_reg = reshape(A_regress, 2400, 1800);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/38; % convert to power spectral density

rng(1);
rand_t = randn(38, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Seasonal_SSTA*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);

mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 2400, 1800);
mask = double(mask);
mask(mask == 0) = NaN;


lats_sparse = lats(1:20:1800);
lons_sparse = lons(1:20:2400);

LA = reshape(double(lats_sparse*ones([1 120])), 90, 120);
LO = reshape(double(lons_sparse*ones([1 90])), 120, 90);

mask = mask(1:20:2400, 1:20:1800);

subplot('position', [0.49 0.55 0.4 0.4]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 25,...
          'MLabelLocation',25,...
          'FlineWidth', 1);

limit = 0.6;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-0.6 0.6]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 0.5);

% [c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
% 				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');
				
%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(60, 122, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'MS EOF 1', 'FontSize',20, 'FontWeight', 'bold');
title('MS EOF 1', 'FontSize',22, 'FontWeight', 'bold');
   
%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap

KI = mean(reshape(HVD_Cop_LS_EOF_1, 6, 38),1);
Var = Seasonal_SSTA;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 2400, 1800);
SST_reg = reshape(A_regress, 2400, 1800);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/38; % convert to power spectral density

rng(1);
rand_t = randn(38, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Seasonal_SSTA*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);

mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 2400, 1800);
mask = double(mask);
mask(mask == 0) = NaN;


lats_sparse = lats(1:20:1800);
lons_sparse = lons(1:20:2400);

LA = reshape(double(lats_sparse*ones([1 120])), 90, 120);
LO = reshape(double(lons_sparse*ones([1 90])), 120, 90);

mask = mask(1:20:2400, 1:20:1800);

subplot('position', [0.05 0.1 0.4 0.4]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 25,...
          'MLabelLocation',25,...
          'FlineWidth', 1);

limit = 0.6;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-0.6 0.6]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 0.05);

%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');
%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(60, 122, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'LS EOF 1', 'FontSize',20, 'FontWeight', 'bold');

title('LS EOF 1', 'FontSize',22, 'FontWeight', 'bold');

%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap

KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 38),1);
Var = Seasonal_SSTA;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 2400, 1800);
SST_reg = reshape(A_regress, 2400, 1800);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/38; % convert to power spectral density

rng(1);
rand_t = randn(38, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Seasonal_SSTA*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);

mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 2400, 1800);
mask = double(mask);
mask(mask == 0) = NaN;

lats_sparse = lats(1:20:1800);
lons_sparse = lons(1:20:2400);

LA = reshape(double(lats_sparse*ones([1 120])), 90, 120);
LO = reshape(double(lons_sparse*ones([1 90])), 120, 90);

mask = mask(1:20:2400, 1:20:1800);

subplot('position', [0.49 0.1 0.4 0.4]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 25,...
          'MLabelLocation',25,...
          'FlineWidth', 1);

limit = 0.6;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-0.6 0.6]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 0.5);

% [c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
% 				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(60, 122, '(d)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'LS EOF 2', 'FontSize',20, 'FontWeight', 'bold');
   
title('LS EOF 2', 'FontSize',22, 'FontWeight', 'bold');
%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap
%subplot('position',[0.1 0 0.8 0.08])

limit = 0.6;
caxis([-limit limit]);
h = colorbar('Location', 'eastoutside', 'XTick', -0.6: 0.1 :0.6);
set(h, 'position',[0.91 0.15 0.012 0.75] );
cbfreeze(h);
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')
h.FontWeight = 'bold';
h.FontSize = 14;
xlabel(h, '')
%cb.Layout.Tile = 'south';
%cb.Label.String = 'SSTA (K)';

%cb = colorbar('FontSize',12,'FontWeight','bold' );
%cb.Layout.Tile = 'south';
%cb.Label.String = 'SSTA (K)';

%print('Nish_092921_KE_PDP_MS_Fig_3_PS','-dpng', '-r500');

print('Nish_030323_KE_PDP_MS_Fig_3_PS_revised','-dpng', '-r500');
