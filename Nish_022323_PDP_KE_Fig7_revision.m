clear all;
close all;
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 17 4.5]);
uwnd = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_uwnd_mon_global.nc', 'uwnd');
uwnd = reshape(uwnd, 144*73, 12, 73);
uwnd = reshape(uwnd, 144*73, 12, 73);
uwnd = uwnd(:, 1:2, 2:72);
uwnd_JF = reshape(uwnd, 144*73, 2*71);
uwnd_JF = reshape(uwnd_JF, 144*73,2, 71);
uwnd_seas = mean(uwnd_JF, 2);
uwnd_seas = reshape(uwnd_seas, 144*73, 1*71);


vwnd = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_vwnd_mon_global.nc', 'vwnd');
vwnd = reshape(vwnd, 144*73, 12, 73);
vwnd = reshape(vwnd, 144*73, 12, 73);
vwnd = vwnd(:, 1:2, 2:72);
vwnd_JF = reshape(vwnd, 144*73, 2*71);
vwnd_JF = reshape(vwnd_JF, 144*73,2, 71);
vwnd_seas = mean(vwnd_JF, 2);
vwnd_seas = reshape(vwnd_seas, 144*73, 1*71);

%uwnd_clim = reshape(squeeze(mean(uwnd_seas,2)), 144, 73);
%vwnd_clim = reshape(squeeze(mean(vwnd_seas,2)), 73, 144);


gph_300 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'hgt');
gph_300 = gph_300;
gph_300 = reshape(gph_300, 144*73, 12, 73);
gph_300 = reshape(gph_300, 144*73, 12, 73);
gph_300 = gph_300(:, 1:2, 2:72);
gph_300_JF = reshape(gph_300, 144*73, 2*71);
gph_300_JF = reshape(gph_300_JF, 144*73,2, 71);
gph_300_seas = mean(gph_300_JF, 2);
gph_300_seas = reshape(gph_300_seas, 144*73, 1*71);

gph_300_seas_anom = detrend(gph_300_seas')';

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lat');

a = 6.37*1E6;    % radius of earth
R = 7.292*1E-5; % rotation rate of earth (omega)

temp = cos(lats*pi/180); % calculate cos(lat)
cos_lat = (temp*ones(1,144));
cos_lat_reshape = reshape((cos_lat),73*144,1);

temp = sin(2*lats*pi/180); % calculate sin(2lat)
temp(1) = temp(2);temp(73) = temp(72); temp(37) = temp(36);
sin2_lat = (temp*ones(1,144));
sin2_lat_reshape = reshape((sin2_lat),73*144,1);

temp = cos(lats*pi/180);  % calculate cos(lat)^-1 for spherical geometry
temp2=temp.\ones(73,1);
inv_cos_lat = (temp2*ones(1,144));

for i=1:71  %loop through each year

%i
  uwnd_t = reshape((reshape(uwnd_seas(:,i), 144, 73*1))', 144*73,1);
  vwnd_t = reshape((reshape(vwnd_seas(:,i), 144, 73*1))', 144*73,1);
  gph_t = reshape((reshape(gph_300_seas(:,i), 144, 73*1))', 144*73,1);

  gph_anom_t = reshape((reshape(gph_300_seas_anom(:,i), 144, 73*1))', 144*73,1);
    
  temp = flipud(reshape(uwnd_t,73,144)); % 'flipud' is used because lat runs from 90 -> -90; if your lat runs S -> N, remove all flipud
  uwnd_zonal_anom = detrend(temp','constant')';	 % remove zonal mean at each lat
  temp = flipud(reshape(vwnd_t,73,144));
  vwnd_zonal_anom = detrend(temp','constant')';	% remove zonal mean at each lat
  temp = ...
	9.8*(flipud(reshape(squeeze(gph_t),73,144)));
  hgt_zonal_anom = flipud(detrend(temp','constant')');% remove zonal mean at each lat	

  temp = 9.8*flipud(reshape(gph_anom_t,73,144));
  gph_zonal_time_anom = flipud(detrend(temp','constant')'); 
  
  clear divu_v divu_u
  temp = [vwnd_zonal_anom vwnd_zonal_anom]; % concatenates along lat circles so can take zonal derivative from lon 360 to 1 
  temp2 = [uwnd_zonal_anom uwnd_zonal_anom];
  temp3 = [hgt_zonal_anom hgt_zonal_anom];
  temp4 = [sin2_lat sin2_lat];

for j = 1:73 % loop through each lat
  temp5 = temp(j,:).*temp3(j,:);
  fx=gradient(temp5,2.5*pi/180);  % take lon derivative, where 2.5 is degree spacing in lon
  fphi = fx./(2*a*R*temp4(j,:));  % account for spherical geometry 
  divx_v(j,:) = circshift(fphi(2:2+144-1),1);  % calculate derivative for zonal component

  temp5 = temp2(j,:).*temp3(j,:);
  fx=gradient(temp5,2.5*pi/180);
  fphi = fx./(2*a*R*temp4(j,:));
  divx_u(j,:) = circshift(fphi(2:2+144-1),1);  % calculate derivative for meridional component 
end

  Fs_uwnd_300_seas(:,i) = (300./1000.)*...  %calculate zonal component of Fs; 300/1000 is 'p'
    ((cos_lat_reshape).*(reshape(flipud(vwnd_zonal_anom.^2),73*144,1)) - ...  % again don't use flipud if your lat runs from S -> N
	(reshape(flipud(divx_v),73*144,1)));
  Fs_vwnd_300_seas(:,i) = (300./1000.)* ... %calculate meridional component of Fs
    ((cos_lat_reshape).* ...
	 (-reshape(flipud(vwnd_zonal_anom.*uwnd_zonal_anom),73*144,1)) + ...
	(reshape(flipud(divx_u),73*144,1)));

  zonal_GPH(:,:,i) = hgt_zonal_anom;
  
  time_zonal_GPH_anom(:,:,i) = gph_zonal_time_anom;

end

%%

Fs_uwnd_300_seas_clim = squeeze(mean(Fs_uwnd_300_seas,2));
Fs_vwnd_300_seas_clim = squeeze(mean(Fs_vwnd_300_seas,2));

% Fs_uwnd_300_seas_plumb = Fs_uwnd_300_seas(:,:);
% Fs_vwnd_300_seas_plumb = Fs_vwnd_300_seas(:,:);
% 
% Fs_u_clim = mean(Fs_uwnd_300_seas,2);
% Fs_v_clim = mean(Fs_vwnd_300_seas,2);

Fs_u_clim = reshape(Fs_uwnd_300_seas_clim, 73, 144);
Fs_v_clim = reshape(Fs_vwnd_300_seas_clim, 73, 144);

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lat');


lats_2 = (lats.*ones(1,length(lons)));
lons_2 = (lons.*ones(1,length(lats)))';
 
zonal_GPH_anom = (reshape(zonal_GPH, 73*144, 71))./9.8;
zonal_GPH_anom_clim = reshape(squeeze(mean(zonal_GPH_anom,2)), 73, 144);

time_zonal_GPH_anom = (reshape(time_zonal_GPH_anom, 73*144, 71))./9.8;
time_zonal_GPH_anom_clim = reshape(squeeze(mean(time_zonal_GPH_anom,2)), 73, 144);

%time_zonal_GPH_anom = detrend(reshape(zonal_GPH_anom, 73*144, 71)')';
%time_zonal_GPH_anom_clim = reshape(squeeze(mean(time_zonal_GPH_anom,2)), 73, 144);
 
 %% plot climatology 

subplot('position', [0.04 0.05 0.42 0.9]) 
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
%h = worldmap;
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth',1);
limit = 300;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
surfm(double(lats), double(lons), zonal_GPH_anom_clim,'Facecolor', 'interp');
colorbar('FontSize',12,'FontWeight','bold');
q = colorbar('FontSize',12,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -300: 100 :300);
set(q, 'position',[0.475 0.152 0.01 0.7] );    
    load coast
   geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
 
x_vec = Fs_u_clim;
plot_data1 = reshape(x_vec, 73, 144);

y_vec = Fs_v_clim;
plot_data2 = reshape(y_vec, 73, 144);

plot_data1_sparse = plot_data1(1:1:73,1:1:144);
plot_data2_sparse = plot_data2(1:1:73,1:1:144);

lat_sparse = double(lats_2(1:1:73,1:1:144));
lon_sparse = double(lons_2(1:1:73,1:1:144));

scale = 1;

lat_temp = lat_sparse; 
lon_temp = lon_sparse;
v_temp = double(plot_data2_sparse)*scale;
u_temp = double(plot_data1_sparse)*scale;

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);
hold on
quiverm(lat_temp, lon_temp, v_temp, u_temp, 'k', 1, 'filled' );

textm(20, 123, '(a)', 'FontSize',16, 'FontWeight', 'bold');
%title('Fs-vectors', 'FontSize',15, 'FontWeight', 'bold');

tightmap
%% KE influence

%% GPH Correlation
%% Read 300 GPH data
ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc');
gph = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc', 'hgt');
gph = gph;
lons_1 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats_1 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc', 'lat');

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

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lat');

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


load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:2, :), 1);
Seasonal_LS_EOF_1 = detrend(KI')';

KI = Seasonal_LS_EOF_1';
Var = Fs_uwnd_300_seas(:,1:71);

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS/length(KI);

E_x_cor = reshape(A_corr, 73, 144);
E_x_reg = reshape(A_regress, 73, 144);

KI = Seasonal_LS_EOF_1';
Var = Fs_vwnd_300_seas(:,1:71);

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS/length(KI);

E_y_cor = reshape(A_corr, 73, 144);
E_y_reg = reshape(A_regress, 73, 144);

%% KE influence on Zonal height anomalies

load HVD_ERSST_EOF2.mat;
HVD_ERSST_EOF2 = reshape(HVD_ERSST_EOF2, 6, 71);
KI = mean(HVD_ERSST_EOF2(1:2, :), 1);
Seasonal_LS_EOF_1 = detrend(KI')';

zonal_GPH_anom = detrend(reshape(zonal_GPH, 73*144, 71)')';

KI = Seasonal_LS_EOF_1';
Var =time_zonal_GPH_anom;


y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS/length(KI);

Zonal_hgt_cor = reshape(A_corr', 73, 144);
Zonal_hgt_reg = reshape(A_regress', 73, 144);

[maxlags,~,~] = size(KI);
[r_KE,lags] = autocorr(KI',maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/71; % convert to power spectral density

rng(1);
rand_t = randn(71, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp99(:,i) = norm;
end

rand_t = temp99;

T_rand_cov = Var*rand_t/71;
T_rand_cov_sort = sort(abs(T_rand_cov), 2);
T_rand_cov_sig  = T_rand_cov_sort(:, 950);
mask = (abs(A_regress) > T_rand_cov_sig);
mask = reshape(mask, 73, 144);
mask = double(mask);
mask(mask == 0) = NaN;

lons = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/300hpa_gph_mon_global.nc', 'lat');

LA = reshape(double(lats*ones([1 144])), 73, 144);
LO = reshape(double(lons*ones([1 73])), 144, 73)';

hgt_anom_std = std(Var', 1);
%%
subplot('position', [0.53 0.05 0.42 0.9])
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 1);
limit = 30;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
surfm(double(lats), double(lons), Zonal_hgt_reg,'Facecolor', 'interp');
colorbar('FontSize',12,'FontWeight','bold');
q = colorbar('FontSize',12,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -30: 10 :30);
set(q, 'position',[0.965 0.152 0.01 0.7] );     

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
	
%textm(80, 120, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%textm(80, 210, '850 GPH', 'FontSize',20, 'FontWeight', 'bold');
textm(20, 123, '(b)', 'FontSize',16, 'FontWeight', 'bold');

x_vec = E_x_reg;
plot_data1 = reshape(x_vec, 73, 144);

y_vec = E_y_reg;
plot_data2 = reshape(y_vec, 73, 144);

plot_data1_sparse = plot_data1(1:1:73,1:1:144);
plot_data2_sparse = plot_data2(1:1:73,1:1:144);

lat_sparse = double(lats_2(1:1:73,1:1:144));
lon_sparse = double(lons_2(1:1:73,1:1:144));


scale = 1;

lat_temp = lat_sparse; lon_temp = lon_sparse;
v_temp = double(plot_data2_sparse)*scale;
u_temp = double(plot_data1_sparse)*scale;

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);
grayColor = [.7 .7 .7];

quiverm(lat_temp, lon_temp, v_temp, u_temp, 'k', 1, 'filled');

temp4 = double(GPH_no_lag_reg);   %Here you assign the lag-correlation field
plot_data_1 = ...
 reshape(temp4, 73,  37);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

lons_1 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lon');
lats_1 = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_850.nc', 'lat');

for i=1:6
        [c h1] = contourm(double(lats_1), double(lons_1), double(plot_data_1'), [(i)*6 (i)*6], ...
                'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');
            %if(~isempty(c));clabel(c,h1);end
        [c h1] = contourm(double(lats_1), double(lons_1), double(plot_data_1'), [(0-i)*6 (0-i)*6], ...
                'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-.');  
            %if(~isempty(c));clabel(c,h1);end
end

plotm(double(LA), double(LO), mask, 'r--', 'LineWidth', 2);
%plotm(double(LA), double(LO), mask, 'color', '#1b7837', 'Marker', '-', 'LineStyle', 'none',...
%    'MarkerSize', 8);

%plotm(double(LA), double(LO), mask, 'color', '#2ca25f',...
%    'Marker', '|', 'MarkerSize', 8, 'LineStyle', 'none');

%title('KE influence Fs-Vectors & Time/Zonal GPH Anomalies', 'FontSize',15, 'FontWeight', 'bold');
tightmap

%print('Nish_082622_KE_PDP_MS_Fig7_revised_new_fig','-dpng', '-r500');

print('Nish_030323_KE_PDP_MS_Fig7_revised_new_fig','-dpng', '-r500');