load GPH_850_NP_NA.mat;
GPH_850_NP_NA = GPH_850_NP_NA./9.8;
GPH_850_NP_NA = reshape(GPH_850_NP_NA, 701* 361 *7, 38);
GPH_850_NP_NA = detrend(GPH_850_NP_NA')';
GPH_850_NP_NA = reshape(GPH_850_NP_NA, 701* 361, 7, 38);
GPH_850_NP_NA = GPH_850_NP_NA(:, 3:4, :);
Seasonal_GPH_850 = reshape(mean(GPH_850_NP_NA, 2), 701*361*1, 38);

lons = ncread('/project/pdpanalysis/nish/data/ERA_5_single_level_energy/ERA5_850hPa_GP.nc', 'longitude');
lats = ncread('/project/pdpanalysis/nish/data/ERA_5_single_level_energy/ERA5_850hPa_GP.nc', 'latitude');

load HVD_Cop_MS_EOF_1;
load HVD_Cop_MS_KI;
load HVD_Cop_LS_EOF_1;
load HVD_Cop_LS_EOF_2;

%%
%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 18 6]);
%figure(1);
%tiledlayout(2,2,'TileSpacing','tight');

HVD_Cop_MS_KI = reshape(HVD_Cop_MS_KI, 6, 38);
KI = mean(HVD_Cop_MS_KI(1:2,:),1);
Var = Seasonal_GPH_850;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 701, 361);
SST_reg = reshape(A_regress, 701, 361);

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

T_rand_cov = Var*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 701, 361);
mask = double(mask);
mask(mask == 0) = NaN;

lats_sparse = lats(1:10:361);
lons_sparse = lons(1:10:701);

LA = reshape(double(lats_sparse*ones([1 71])), 37, 71);
LO = reshape(double(lons_sparse*ones([1 37])), 71, 37);

mask = mask(1:10:701, 1:10:361);

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 14 7]);

subplot('position', [0.05 0.55 0.4 0.4]);
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 290];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 30,...
          'FlineWidth', 1);

limit = 14;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-14 14]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

            
% [c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
% 				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');
%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(20, 122, '(a)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'MS KI', 'FontSize',20, 'FontWeight', 'bold');
title('MS KI', 'FontSize',17, 'FontWeight', 'bold');
%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap

%KI = mean(reshape(HVD_Cop_MS_EOF_1, 6, 38),1);
HVD_Cop_MS_EOF_1 = reshape(HVD_Cop_MS_EOF_1, 6, 38);
KI = mean(HVD_Cop_MS_EOF_1(1:2,:),1);
Var = Seasonal_GPH_850;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 701, 361);
SST_reg = reshape(A_regress, 701, 361);

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

T_rand_cov = Var*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);

mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 701, 361);
mask = double(mask);
mask(mask == 0) = NaN;

lats_sparse = lats(1:10:361);
lons_sparse = lons(1:10:701);

LA = reshape(double(lats_sparse*ones([1 71])), 37, 71);
LO = reshape(double(lons_sparse*ones([1 37])), 71, 37);

mask = mask(1:10:701, 1:10:361);

subplot('position', [0.49 0.55 0.4 0.4]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 290];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 14;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-14 14]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

% [c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
% 				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');
				
%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(20, 122, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'MS EOF 1', 'FontSize',20, 'FontWeight', 'bold');
title('MS EOF 1', 'FontSize',17, 'FontWeight', 'bold');

%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap

%KI = mean(reshape(HVD_Cop_LS_EOF_1, 6, 38),1);
HVD_Cop_LS_EOF_1 = reshape(HVD_Cop_LS_EOF_1, 6, 38);
KI = mean(HVD_Cop_LS_EOF_1(1:2,:),1);
Var = Seasonal_GPH_850;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 701, 361);
SST_reg = reshape(A_regress, 701, 361);

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

T_rand_cov = Var*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);

mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 701, 361);
mask = double(mask);
mask(mask == 0) = NaN;

lats_sparse = lats(1:10:361);
lons_sparse = lons(1:10:701);

LA = reshape(double(lats_sparse*ones([1 71])), 37, 71);
LO = reshape(double(lons_sparse*ones([1 37])), 71, 37);

mask = mask(1:10:701, 1:10:361);

subplot('position', [0.05 0.1 0.4 0.4]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 290];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 14;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-14 14]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

% [c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
% 				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');
%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(20, 122, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'LS EOF 1', 'FontSize',20, 'FontWeight', 'bold');
title('LS EOF 1', 'FontSize',17, 'FontWeight', 'bold');
   
%title('Correlation b/w SSTAs')    
%title('Regression b/w SSTAs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 1.5)
tightmap

%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 38),1);
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 38);
KI = mean(HVD_Cop_LS_EOF_2(1:2,:),1);
Var = Seasonal_GPH_850;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 701, 361);
SST_reg = reshape(A_regress, 701, 361);

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
T_rand_cov = Var*rand_t/38;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);

mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 701, 361);
mask = double(mask);
mask(mask == 0) = NaN;

lats_sparse = lats(1:10:361);
lons_sparse = lons(1:10:701);

LA = reshape(double(lats_sparse*ones([1 71])), 37, 71);
LO = reshape(double(lons_sparse*ones([1 37])), 71, 37);

mask = mask(1:10:701, 1:10:361);

subplot('position', [0.49 0.1 0.4 0.4]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [15 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 290];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 14;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-limit limit]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

% [c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
% 				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

%colorbar
load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 1)

textm(20, 122, '(d)', 'FontSize',20, 'FontWeight', 'bold');
%textm(70, 180, 'LS EOF 2', 'FontSize',20, 'FontWeight', 'bold');
title('LS EOF 2', 'FontSize',17, 'FontWeight', 'bold');
   
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

limit = 14;
caxis([-limit limit]);
h = colorbar('Location', 'eastoutside', 'XTick', -14: 2 :14);
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

%print('Nish_092921_KE_PDP_MS_Fig_4_PS','-dpng', '-r500');

print('Nish_030223_KE_PDP_MS_Fig_5_new_lag','-dpng', '-r500');
