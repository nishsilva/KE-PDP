clear all;
close all;
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 15.5 7.5]);
%% GPH Correlation
%% Read 850 GPH data
ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc');
gph = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'hgt');
gph = gph;
lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

%% Etract winter months
gph_winter = reshape(gph, 73 , 37 , 12 , 73);
gph_winter(:,:,5:10,:) = [];
gph_win = reshape(gph_winter, 73 , 37 , 6 * 73);
gph_w =  gph_win(:,:,5:(438-2));
gph_w = reshape(gph_w, 73*37*6, 72);
gph_w = gph_w(:,1:71);
gph_w = detrend(gph_w')';
gph_w = reshape(gph_w, 73*37,6, 71);
gph_w = gph_w(:, 3:4, :);
winter_mean = mean(gph_w, 2);
Seasonal_GPH = reshape(winter_mean, 73*37, 1*71);
Seasonal_GPH = detrend(Seasonal_GPH')';


load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:2, :), 1);

%KI = mean(reshape(HVD_ERSST_EOF2, 6, 71),1);
Var = Seasonal_GPH;

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 73, 37);
SST_reg = reshape(A_regress, 73, 37);

GPH_no_lag_cor = SST_cor;
GPH_no_lag_reg = SST_reg;

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 73, 37);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 73])), 37, 73);
LO = reshape(double(lons*ones([1 37])), 73, 37);

%nexttile
subplot('position', [0.05 0.51 0.42 0.4])
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 12;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
%caxis([-1.2 1.2]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);

%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

q = colorbar('FontSize',12,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -12: 04 :12);
%q.Ruler.TickLabelFormat='%0.1f';
set(q, 'position',[0.48 0.52 0.01 0.385] );
%ylabel(q, 'GPH (m)', 'FontSize',12,'FontWeight','bold')

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
	
%textm(20, 123, '(a) 850hPA GPH', 'FontSize',15, 'FontWeight', 'bold');
textm(20, 123, '(a)', 'FontSize',15, 'FontWeight', 'bold');
%textm(80, 210, '850 GPH', 'FontSize',20, 'FontWeight', 'bold');
   
title('850hPa GPH', 'FontSize',15, 'FontWeight', 'bold')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
tightmap


%% GPH Correlation
%% Read 300 GPH data
ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc');
gph = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc', 'hgt');
gph = gph;
lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc', 'lat');

%% Etract winter months
gph_winter = reshape(gph, 73 , 37 , 12 , 73);
gph_winter(:,:,5:10,:) = [];
gph_win = reshape(gph_winter, 73 , 37 , 6 * 73);
gph_w =  gph_win(:,:,5:(438-2));
gph_w = reshape(gph_w, 73*37*6, 72);
gph_w = gph_w(:,1:71);
gph_w = detrend(gph_w')';
gph_w = reshape(gph_w, 73*37,6, 71);
gph_w = gph_w(:, 3:4, :);
winter_mean = mean(gph_w, 2);
Seasonal_GPH = reshape(winter_mean, 73*37, 1*71);
Seasonal_GPH = detrend(Seasonal_GPH')';

load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:2, :), 1);

%load HVD_ERSST_EOF2.mat;
%LS_EOF_1 = reshape(HVD_ERSST_EOF2, 6, 71);
%Seasonal_LS_EOF_1 = mean(LS_EOF_1, 1);
Seasonal_LS_EOF_1 = detrend(KI')';

KI = Seasonal_LS_EOF_1;
Var = Seasonal_GPH;

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 73, 37);
SST_reg = reshape(A_regress, 73, 37);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 73, 37);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 73])), 37, 73);
LO = reshape(double(lons*ones([1 37])), 73, 37);

GPH_no_lag_cor = SST_cor;
GPH_no_lag_reg = SST_reg;

subplot('position', [0.525 0.51 0.42 0.4])
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'off', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 24;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
%caxis([-2.4 2.4]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);


%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

q = colorbar('FontSize',12,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -24: 08 :24);
%q.Ruler.TickLabelFormat='%0.1f';
set(q, 'position',[0.955 0.52 0.01 0.385] );
%ylabel(q, 'GPH (m)', 'FontSize',12,'FontWeight','bold')

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
	
%textm(80, 123, '(b) 300hPa GPH', 'FontSize',16, 'FontWeight', 'bold');
textm(20, 123, '(b)', 'FontSize',15, 'FontWeight', 'bold');
%textm(80, 210, '850 GPH', 'FontSize',20, 'FontWeight', 'bold');
   
title('300hPa GPH', 'FontSize',15, 'FontWeight', 'bold')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)

la = [37.5 37.5 62.5 62.5 37.5]; 
lo = [250 300 300 250 250];

[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)

tightmap
%% contours for heat flux

gph_winter = reshape(gph, 73 , 37 , 12 , 73);
gph_winter(:,:,5:10,:) = [];
gph_win = reshape(gph_winter, 73 , 37 , 6 * 73);
gph_w =  gph_win(:,:,5:(438-2));
gph_w = reshape(gph_w, 73*37*6, 72);
gph_w = gph_w(:,1:71);
gph_w = detrend(gph_w')';
gph_w = reshape(gph_w, 73*37,6, 71);
gph_w = gph_w(:, 1:4, :);
winter_mean = mean(gph_w, 2);
Seasonal_GPH = reshape(winter_mean, 73*37, 1*71);
Seasonal_GPH = detrend(Seasonal_GPH')';

load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:4, :), 1);

%load HVD_ERSST_EOF2.mat;
%LS_EOF_1 = reshape(HVD_ERSST_EOF2, 6, 71);
%Seasonal_LS_EOF_1 = mean(LS_EOF_1, 1);
Seasonal_LS_EOF_1 = detrend(KI')';

KI = Seasonal_LS_EOF_1;
Var = Seasonal_GPH;

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 73, 37);
SST_reg = reshape(A_regress, 73, 37);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 73, 37);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 73])), 37, 73);
LO = reshape(double(lons*ones([1 37])), 73, 37);

GPH_no_lag_cor = SST_cor;
GPH_no_lag_reg = SST_reg;
%%
ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_LHF.nc');
lhtfl = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_LHF.nc', 'lhtfl');

ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_SHF.nc');
shtfl = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_SHF.nc', 'shtfl');


lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_SHF.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_SHF.nc', 'lat');

lons_1 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats_1 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

thflx = lhtfl + shtfl;

thflx = thflx(:,:, 1:876);
thflx_winter = reshape(thflx, 192 , 94 , 12 , 73);
thflx_winter(:,:,5:10,:) = [];
thflx_win = reshape(thflx_winter, 192 , 94 , 6 * 73);
thflx_w =  thflx_win(:,:,5:(438-2));
thflx_w = reshape(thflx_w, 192*94*6, 72);
thflx_w = thflx_w(:,1:71);
thflx_w = detrend(thflx_w')';
thflx_w = reshape(thflx_w, 192*94,6, 71);
thflx_w = thflx_w(:, 1:4, :);
winter_mean = mean(thflx_w, 2);
Seasonal_thflx = reshape(winter_mean, 192*94, 1*71);
Seasonal_thflx = detrend(Seasonal_thflx')';

load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:4, :), 1);
Seasonal_LS_EOF_2 = detrend(KI')';

KI = Seasonal_LS_EOF_2;
Var = Seasonal_thflx;

%KI = KI(:, 1:70);
%Var = Var(:, 2:end);


y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 192, 94);
SST_reg = reshape(A_regress, 192, 94);

%GPH_no_lag_cor = SST_cor;
%GPH_no_lag_reg = SST_reg;

%load GPH_no_lag_reg.mat;

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 192, 94);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 192])), 94, 192);
LO = reshape(double(lons*ones([1 94])), 192, 94);

subplot('position', [0.05 0.05 0.42 0.4])
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 10;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-10 10]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');
%colorbar('FontSize',12,'FontWeight','bold');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);
%plotm(double(LA), double(LO'), mask', 'k.', 'MarkerSize', 10);

temp4 = double(GPH_no_lag_reg);   %Here you assign the lag-correlation field
plot_data_1 = ...
 reshape(temp4, 73,  37);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

for i=1:6
        [c h1] = contourm(double(lats_1), double(lons_1), double(plot_data_1'), [(i)*6 (i)*6], ...
                'LineWidth',2,'color', [0 1 0]/2	, 'LineStyle', '-');
            %if(~isempty(c));clabel(c,h1);end
        [c h1] = contourm(double(lats_1), double(lons_1), double(plot_data_1'), [(0-i)*6 (0-i)*6], ...
                'LineWidth',2,'color', [0 1 0]/2	, 'LineStyle', '-.');  
            %if(~isempty(c));clabel(c,h1);end
	end
	
q = colorbar('FontSize',12,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -10: 4 :10);
set(q, 'position',[0.48 0.06 0.01 0.385] );

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
	
%textm(80, 123, '(c) v^{\prime}v^\prime', 'FontSize',16, 'FontWeight', 'bold');
textm(20, 123, '(c)', 'FontSize',15, 'FontWeight', 'bold');

%textm(80, 118, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%textm(80, 210, 'Storm Track', 'FontSize',20, 'FontWeight', 'bold');
   
title('Total Heat Flux', 'FontSize',15, 'FontWeight', 'bold') 

hold on;

la = [31 31 44 44 31]; 
lo = [140 182 182 140 140];

%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','g', 'Linewidth', 2)   
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
tightmap
%% contours for v'T'

gph_winter = reshape(gph, 73 , 37 , 12 , 73);
gph_winter(:,:,5:10,:) = [];
gph_win = reshape(gph_winter, 73 , 37 , 6 * 73);
gph_w =  gph_win(:,:,5:(438-2));
gph_w = reshape(gph_w, 73*37*6, 72);
gph_w = gph_w(:,1:71);
gph_w = detrend(gph_w')';
gph_w = reshape(gph_w, 73*37,6, 71);
gph_w = gph_w(:, 1:4, :);
winter_mean = mean(gph_w, 2);
Seasonal_GPH = reshape(winter_mean, 73*37, 1*71);
Seasonal_GPH = detrend(Seasonal_GPH')';

load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:4, :), 1);

%load HVD_ERSST_EOF2.mat;
%LS_EOF_1 = reshape(HVD_ERSST_EOF2, 6, 71);
%Seasonal_LS_EOF_1 = mean(LS_EOF_1, 1);
Seasonal_LS_EOF_1 = detrend(KI')';

KI = Seasonal_LS_EOF_1;
Var = Seasonal_GPH;

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 73, 37);
SST_reg = reshape(A_regress, 73, 37);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 73, 37);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 73])), 37, 73);
LO = reshape(double(lons*ones([1 37])), 73, 37);

GPH_no_lag_cor = SST_cor;
GPH_no_lag_reg = SST_reg;

%%
air_temp_d_1 = ncread('/project/pdpanalysis/nish/data/VT/NCEP_850_T_1948_1988_daily.nc', 'air');
air_temp_d_2 = ncread('/project/pdpanalysis/nish/data/VT/NCEP_850_T_1989_2020_daily.nc', 'air');

air_temp_d = cat(4, air_temp_d_1, air_temp_d_2);

air_temp_d(:,:,:,[60, (1461+60), ((1461*2)+60), ((1461*3)+60), ((1461*4)+60), ((1461*5)+60),...
 ((1461*6)+60),...
 ((1461*7)+60), ((1461*8)+60),...
((1461*9)+60), ((1461*10)+60), ((1461*11)+60), ((1461*12)+60), ((1461*13)+60),...
 ((1461*14)+60), ((1461*15)+60), ((1461*16)+60), ((1461*17)+60), ((1461*18)+60)]) = []; %remove feb 29th of all leap years

% air_temp_d = reshape(air_temp_d, 73, 37, 365, 73);
% air_temp_d = air_temp_d(:,:,1:59, 2:72);
% %air_temp_d(:, :,121:304, :) = [];
% air_temp_d = reshape(air_temp_d, 73, 37, 59 *71);
% %air_temp_d = air_temp_d(:, :, 121:12971);
% air_temp_d = reshape(air_temp_d, 73*37*59, 71);
% air_temp_d_anom = detrend(air_temp_d')';

% air_temp_d = reshape(air_temp_d, 73, 37, 365, 73);
% air_temp_d = air_temp_d(:,:,305:365, 1:71);
% %air_temp_d(:, :,121:304, :) = [];
% air_temp_d = reshape(air_temp_d, 73, 37, 61 *71);
% %air_temp_d = air_temp_d(:, :, 121:12971);
% air_temp_d = reshape(air_temp_d, 73*37*61, 71);
% air_temp_d_anom = detrend(air_temp_d')';



vwnd_d_1 = ncread('/project/pdpanalysis/nish/data/VT/NCEP_850_v_wind_1948_1988_daily.nc', 'vwnd');
vwnd_d_2 = ncread('/project/pdpanalysis/nish/data/VT/NCEP_850_v_wind_1989_2020_daily.nc', 'vwnd');

vwnd_d = cat(4, vwnd_d_1, vwnd_d_2);

vwnd_d(:,:,:,[60, (1461+60), ((1461*2)+60), ((1461*3)+60), ((1461*4)+60), ((1461*5)+60),...
 ((1461*6)+60),...
 ((1461*7)+60), ((1461*8)+60),...
((1461*9)+60), ((1461*10)+60), ((1461*11)+60), ((1461*12)+60), ((1461*13)+60),...
 ((1461*14)+60), ((1461*15)+60), ((1461*16)+60), ((1461*17)+60), ((1461*18)+60)]) = []; %remove feb 29th of all leap years
%%
% vwnd_d = reshape(vwnd_d, 73, 37, 365, 73);
% vwnd_d = vwnd_d(:,:,1:59, 2:72);
% %vwnd_d(:, :,121:304, :) = [];
% vwnd_d = reshape(vwnd_d, 73, 37, 59 *71);
% %vwnd_d = vwnd_d(:, :, 121:12971);
% vwnd_d = reshape(vwnd_d, 73*37*59, 71);
% vwnd_d_anom = detrend(vwnd_d')';
% 
% 
% vT = vwnd_d_anom.*air_temp_d_anom;
% vT = reshape(vT, 73*37, 59, 71);
% 
% vT_Jan = vT(:,1:31,:);
% vT_Jan_mean = mean(vT_Jan, 2);
% vT_Feb = vT(:,32:59,:);
% vT_Feb_mean = mean(vT_Feb, 2);
% %vT_Mar = vT(:,60:90,:);
% %vT_Mar_mean = mean(vT_Mar, 2);
% 
% vT_seas = [vT_Jan_mean vT_Feb_mean];
% 
% vT_seas = reshape(vT_seas, 73*37*2, 71);
% 
% %v_prm_T_prm = vT_seas;
% 
% v_prm_T_prm = detrend(vT_seas')';


% vwnd_d = reshape(vwnd_d, 73, 37, 365, 73);
% vwnd_d = vwnd_d(:,:,305:365, 1:71);
% %vwnd_d(:, :,121:304, :) = [];
% vwnd_d = reshape(vwnd_d, 73, 37, 61 *71);
% %vwnd_d = vwnd_d(:, :, 121:12971);
% vwnd_d = reshape(vwnd_d, 73*37*61, 71);
% vwnd_d_anom = detrend(vwnd_d')';
% 
% 
% vT = vwnd_d_anom.*air_temp_d_anom;
% vT = reshape(vT, 73*37, 61, 71);
% 
% vT_Nov = vT(:,1:30,:);
% vT_Nov_mean = mean(vT_Nov, 2);
% vT_Dec = vT(:,31:61,:);
% vT_Dec_mean = mean(vT_Dec, 2);
% %vT_Mar = vT(:,60:90,:);
% %vT_Mar_mean = mean(vT_Mar, 2);
% 
% vT_seas = [vT_Nov_mean vT_Dec_mean];
% %vT_seas = vT_Dec_mean;
% 
% vT_seas = reshape(vT_seas, 73*37*2, 71);
% 
% %v_prm_T_prm = vT_seas;
% 
% v_prm_T_prm = detrend(vT_seas')';

%%
vwnd_d = reshape(vwnd_d, 73, 37, 365, 73);
vwnd_d = vwnd_d(:,:,:,1:72);
vwnd_d(:, :,60:304, :) = [];
vwnd_d = reshape(vwnd_d, 73,37, 120*72);
vwnd_d = vwnd_d(:,:, 60:8579);

vwnd_d = reshape(vwnd_d, 73*37*120, 71);
vwnd_d_anom = detrend(vwnd_d')';

air_temp_d = reshape(air_temp_d, 73, 37, 365, 73);
air_temp_d = air_temp_d(:,:,:,1:72);
air_temp_d(:, :,60:304, :) = [];
air_temp_d = reshape(air_temp_d, 73,37, 120*72);
air_temp_d = air_temp_d(:,:, 60:8579);

air_temp_d = reshape(air_temp_d, 73*37*120, 71);
air_temp_d_anom = detrend(air_temp_d')';

vT = vwnd_d_anom.*air_temp_d_anom;
vT = reshape(vT, 73*37, 120, 71);

vT_Nov = vT(:,1:30,:);
vT_Nov_mean = mean(vT_Nov, 2);
vT_Dec = vT(:,31:61,:);
vT_Dec_mean = mean(vT_Dec, 2);
vT_Jan = vT(:,62:92,:);
vT_Jan_mean = mean(vT_Jan, 2);
vT_Feb = vT(:,93:120,:);
vT_Feb_mean = mean(vT_Feb, 2);

vT_seas = [vT_Nov_mean vT_Dec_mean vT_Jan_mean vT_Feb_mean];
vT_seas = reshape(vT_seas, 73*37*4, 71);

v_prm_T_prm = detrend(vT_seas')';

%%
%load HVD_ERSST_EOF2.mat;
%LS_EOF_1 = reshape(HVD_ERSST_EOF2, 6, 71);
%Seasonal_LS_EOF_1 = mean(LS_EOF_1, 1);
%Seasonal_LS_EOF_1 = detrend(Seasonal_LS_EOF_1')';

v_prm_T_prm_seas = reshape(mean(reshape(v_prm_T_prm, 73*37, 4, 71),2), 73*37, 1*71);


load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:4, :), 1);

%load HVD_ERSST_EOF2.mat;
%LS_EOF_1 = reshape(HVD_ERSST_EOF2, 6, 71);
%Seasonal_LS_EOF_1 = mean(LS_EOF_1, 1);
Seasonal_LS_EOF_1 = detrend(KI')';

KI = Seasonal_LS_EOF_1;
Var = v_prm_T_prm_seas;

lons = ncread('/project/pdpanalysis/nish/data/VT/NCEP_850_T_1948_1988_daily.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/VT/NCEP_850_T_1948_1988_daily.nc', 'lat');

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 73, 37);
SST_reg = reshape(A_regress, 73, 37);

%GPH_no_lag_cor = SST_cor;
%GPH_no_lag_reg = SST_reg;

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 73, 37);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 73])), 37, 73);
LO = reshape(double(lons*ones([1 37])), 73, 37);


subplot('position', [0.525 0.05 0.42 0.4])
%subplot(2, 2, 1);
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'FontColor', 'k',...
          'LabelFormat', 'compass',...
          'PLabelMeridian', 'west',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);

limit = 3;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
%caxis([-1 1]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);

%colorbar('FontSize',12,'FontWeight','bold');

q = colorbar('FontSize',12,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -4: 1 :4);
set(q, 'position',[0.955 0.06 0.01 0.385] );    

temp4 = double(GPH_no_lag_reg);   %Here you assign the lag-correlation field
plot_data_1 = ...
 reshape(temp4, 73,  37);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

for i=1:6
        [c h1] = contourm(double(lats_1), double(lons_1), double(plot_data_1'), [(i)*6 (i)*6], ...
                'LineWidth',2,'color', [0 1 0]/2	, 'LineStyle', '-');
            %if(~isempty(c));clabel(c,h1);end
        [c h1] = contourm(double(lats_1), double(lons_1), double(plot_data_1'), [(0-i)*6 (0-i)*6], ...
                'LineWidth',2,'color', [0 1 0]/2	, 'LineStyle', '-.');  
            %if(~isempty(c));clabel(c,h1);end
	end

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
	
%textm(20, 120, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%textm(80, 210, '850 GPH', 'FontSize',20, 'FontWeight', 'bold');

textm(20, 123, '(d)', 'FontSize',15, 'FontWeight', 'bold');

   
title('v^{\prime}T^\prime', 'FontSize',15, 'FontWeight', 'bold')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
tightmap




%print('Nish_102523_KE_PDP_MS_Fig6_revised','-dpng', '-r500');