
load HVD_ERSST_EOF1;
load HVD_ERSST_EOF2;
load HVD_Cop_LS_EOF_1;
load HVD_Cop_LS_EOF_2;

x = NaN(1,198);
y = [x,HVD_Cop_LS_EOF_2];

NOAA_Cop_LS_Indices = [y; HVD_ERSST_EOF2];

save NOAA_Cop_LS_Indices.mat NOAA_Cop_LS_Indices;

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 18 9]);
%tiledlayout(2,2,'TileSpacing','Compact', 'Padding', 'compact');

%nexttile([1 2])
subplot('position', [0.082 0.48 0.85 0.45])

b = plot(NOAA_Cop_LS_Indices(2,:), 'Color', '#2b8cbe', 'LineWidth',2);
hold on
a = plot(NOAA_Cop_LS_Indices(1,:), '-.', 'Color', '#de2d26', 'LineWidth',2);
load dates_vec_6.mat;
dates_vec = extractAfter(dates_vec,["01-"]);
yline(0, '--', 'LineWidth',1.5);
xlim([0 428]);
set(gca,'xtick',(1:70:428),'xticklabel',dates_vec(1:70:428), 'FontSize',20, 'FontWeight', 'bold');
legend([b, a],{'ERSST Index', 'OSTIA Index'},...
    'Location','northeast');
legend boxoff;
text(2.5, 3.5, '(a)', 'FontSize',20, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 16; 
ax.LineWidth = 1;

%% Reading data from monthly files

%tiledlayout(2,2,'TileSpacing','compact', 'Padding', 'compact');

ncdisp('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc');

ncvars =  {'time', 'lon', 'lat', 'sst'};
projectdir = '/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/';
dinfo = dir( fullfile(projectdir, '*.nc') );
num_files = length(dinfo);
filenames = fullfile( projectdir, {dinfo.name} );
startday = '1854-01-01 00:00';
start = datetime(startday,'InputFormat','yyyy-MM-dd HH:mm');

for K = 1 : num_files
  this_file = filenames{K};
  %lats(:,:,K) = ncread(this_file, ncvars{1},441, 240);
  %lons(:,:,K) = ncread(this_file, ncvars{2},481, 480);
  times(K)= ncread(this_file, ncvars{1});
  sst(:,:,K) = ncread(this_file, ncvars{4});
  %sst(:,:,K) = ncread(this_file, ncvars{4}, startLoc, count);
  
  %data_day(:,K) = start + seconds(times);
  %sst(:,:,:,K) = downsample_ts(anom, data_day, 'month', 'mean'); 
  %mean_sst = repmat(means, 1, 1, 2);
  %mean_sst(:,:,[i, i+1]) = means;

end

lons = ncread('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc', 'lat');


%% Extract the winter months and redo the EOF Index
SST_v3 = sst(:,:,1129:1992);                    %selects data from 1948 Jan - 2019 Dec
SST_v3 = reshape(SST_v3, 180*89, 12, 72);
SST_v3 = reshape(SST_v3, 180*89*12, 72);
sst_anom = detrend(SST_v3')';

data_anom = reshape(sst_anom, 180*89, 12, 72);
data_anom = reshape(data_anom, 180, 89, 12 , 72);
data_anom(:,:,5:10,:) = [];
SST_W = reshape(data_anom, 180, 89, 6 * 72);
SST_W = SST_W(:,:,5:(432-2));		%selects 1948 Nov - 2019 May
winter_sst = reshape(SST_W, 180*89,6, 71);
winter_mean = mean(winter_sst, 2);
Seasonal_SSTA = reshape(winter_mean, 180*89, 1*71);
Seasonal_SSTA = detrend(Seasonal_SSTA')';

load HVD_ERSST_EOF2.mat;

KI = mean(reshape(HVD_ERSST_EOF2, 6, 71),1);
Var = Seasonal_SSTA;

lons = ncread('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc', 'lat');

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end));
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';
    
A_corr = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
A_regress = Var*T_norm_LS'/length(KI);

SST_cor = reshape(A_corr, 180, 89);
SST_reg = reshape(A_regress, 180, 89);

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
mask = reshape(mask, 180, 89);
mask = double(mask);
mask(mask == 0) = NaN;

LA = reshape(double(lats*ones([1 180])), 89, 180);
LO = reshape(double(lons*ones([1 89])), 180, 89);

%%
%nexttile
subplot('position', [0.08 0.035 0.32 0.42])

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [19.8 68];  %this is where you set the lat/lon limits you want to plot
LonLim = [120.5 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize',12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
      
      framem('FLineWidth',0.5);

limit = 0.5;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
caxis([-0.5 0.5]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');

%plotm(double(LA), double(LO'), mask', 'k.', 'MarkerSize', 3.3);


%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

q = colorbar('FontSize',14,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -0.5: 0.2 :0.5);
set(q, 'position',[0.405 0.06 0.01 0.32] );    
%ylabel(q, 'SSTA (K)', 'FontSize',12,'FontWeight','bold')

load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth',1)

	
textm(23, 121, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%textm(72.7, 210, 'SSTA', 'FontSize',20, 'FontWeight', 'bold');
   
%title('Regression b/w SSTAs')    
%colorbar('location', 'South');    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
%tightmap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
winter_mean = mean(gph_w, 2);
Seasonal_GPH = reshape(winter_mean, 73*37, 1*71);
Seasonal_GPH = detrend(Seasonal_GPH')';


load HVD_ERSST_EOF2.mat;

KI = mean(reshape(HVD_ERSST_EOF2, 6, 71),1);
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
subplot('position', [0.47 0.02 0.48 0.44])
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [22 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
      
      framem('FLineWidth',0.5);

limit = 12;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
%caxis([-1.2 1.2]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

%plotm(double(LA'), double(LO), mask, 'k.', 'MarkerSize', 3.3);

plotm(double(LA'), double(LO), mask, 'color', '#2ca25f',...
    'Marker', '.', 'MarkerSize', 8, 'LineStyle', 'none');


%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

q = colorbar('FontSize',14,'FontWeight','bold', 'Location', 'eastoutside',...
		'XTick', -12: 4 :12);
set(q, 'position',[0.927 0.062 0.01 0.33] );     
%ylabel(q, 'GPH (m)', 'FontSize',12,'FontWeight','bold')

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1)
	
textm(25, 123, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%textm(80, 210, '850 GPH', 'FontSize',20, 'FontWeight', 'bold');
   
%title('Regression b/w GPHs')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)
%tightmap

%print('Nish_092921_KE_PDP_MS_Fig_5_PS','-dpng', '-r500');
%print('Nish_030323_KE_PDP_MS_Fig_5_PS_revised','-dpng', '-r500');
print('Nish_080524_KE_PDP_MS_Fig_5_PS_revised','-dpng', '-r500');