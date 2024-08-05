%%High Pressure Pattern Index
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 6 8]);
figure(1)
hFig=gcf;
set(hFig);
%% GPH Correlation
%% Read 300 GPH data
ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc');
gph = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc', 'hgt');
gph = gph./(9.80665);
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
gph_win = gph_w(:,3:4,:);
winter_mean = mean(gph_win, 2);
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

rng(1);
rand_t = randn(71, 1000);

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

ax(1) = subplot('position', [0.01 0.5 0.98 00.49]);
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [14.5 75];  %this is where you set the lat/lon limits you want to plot
LonLim = [122 300];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'off', ...
          'FontSize', 8, ...
          'FontWeight', 'normal', ...
          'LabelFormat', 'none');

limit = 2.4;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray')
cmap_mod = flipud(precip_cmap);
colormap(gca, cmap_mod)
%caxis([-2.4 2.4]); % set the limits of the colorbar

surfm(double(lats), double(lons), SST_reg','Facecolor', 'interp');

plotm(double(LA'), double(LO), mask, 'w.', 'MarkerSize', 3.3);


%[c h] = contourm(double(lats), double(lons), double(mask'), 1, ...
%				'LineWidth',2,'color', [0 1 0]/2, 'LineStyle', '-');

cb1 = colorbar('FontSize',15,'FontWeight','bold', 'Location', 'southoutside',...
		'XTick', -2: 0.8 :2);
%ylabel(q, '300hPa GPH (m)', 'FontSize',12,'FontWeight','bold')

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'Linewidth', 1.5)
	
%textm(80, 123, '(b) 300hPa GPH', 'FontSize',16, 'FontWeight', 'bold');
textm(20, 123, '(a)', 'FontSize',15, 'FontWeight', 'bold');
%textm(80, 210, '850 GPH', 'FontSize',20, 'FontWeight', 'bold');
   
%title('300hPa GPH', 'FontSize',15, 'FontWeight', 'bold')    
    
hold on;

%la = [30 30 42 42 30]; 
%lo = [140 172 172 140 140];



%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','k', 'Linewidth', 2)

la = [37.5 37.5 62.5 62.5 37.5]; 
lo = [250 300 300 250 250];

[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)

tightmap


%print('Nish_012222_Monopole_Domain','-dpng', '-r500');

%% Monopole Index Generation

%% GPH Correlation
%% Read 300 GPH data
ncdisp('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc');
gph = ncread('/project/pdpanalysis/nish/data/NCEP_NCAR/NCEP_NCAR_R1_GPH_300.nc', 'hgt');

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
gph_win = gph_w(:,3:4,:);
winter_mean = mean(gph_win, 2);
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

SST_reg_monopole = reshape(SST_reg(53:73, 11:21), 21*11, 1);

gph_w = reshape(gph_w, 73, 37,6, 71);
gph_w_monopole_domain = reshape(gph_w(53:73, 11:21, :,:), 21*11, 6, 71); 

monopole_ew = reshape(mean(gph_w_monopole_domain(:, 1:2, :), 2), 21*11, 1*71);
monopole_lw = reshape(mean(gph_w_monopole_domain(:, 3:4, :), 2), 21*11, 1*71);

for n = 1:71

	rho_ew(:, n) = corr(monopole_ew(:, n), SST_reg_monopole);
	rho_lw(:, n) = corr(monopole_lw(:, n), SST_reg_monopole);
	
end

rho_ew = detrend(rho_ew')';
rho_lw = detrend(rho_lw')';

monopole_index_ew = rho_ew;
monopole_index_lw = rho_lw;

%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 21 6]);
%figure(2)

%b = plot(rho, 'Color', '#2b8cbe', 'LineWidth',1.5);
%load dates_vec_6.mat;
%hold on
%yline(0, '--', 'LineWidth',1.5);
%xlim([0 428]);
%set(gca,'xtick',(1:70:428),'xticklabel',dates_vec(1:70:428), 'FontSize',20, 'FontWeight', 'bold');
%legend([b],{'Monopole Index'},...
%    'Location','northeast');
%legend boxoff;
%text(2.5, 3.5, '(a)', 'FontSize',20, 'FontWeight', 'bold');
%ax = gca;
%ax.FontSize = 20; 
%ax.LineWidth = 2;

	
%print('Nish_012222_Monopole_Index','-dpng', '-r500');

	
%% Monopole Index Regression on SSTs

%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 6 4]);


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
  %times(K)= ncread(this_file, ncvars{1});
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
winter_sst(isnan(winter_sst))=0;

winter_ew = mean(winter_sst(:, 1:2, :), 2);
winter_lw = mean(winter_sst(:, 3:4, :), 2);

winter_mean = mean(winter_sst, 2);
Seasonal_SSTA = reshape(winter_mean, 180*89, 1*71);
Seasonal_SSTA = detrend(Seasonal_SSTA')';


%% Granger Causality

monopole_index_ew = monopole_index_ew;
monopole_index_lw = monopole_index_lw;

%winter_sst(isnan(winter_sst))=0;
SST_ew = winter_ew;
SST_lw = winter_lw;

early_winter_mean = SST_ew;
late_winter_mean = SST_lw;

Seasonal_ew_SST = reshape(early_winter_mean, 180*89, 1*71);
Seasonal_lw_SST = reshape(late_winter_mean, 180*89, 1*71);

Seasonal_ew_SST = detrend(Seasonal_ew_SST')'; %independent variable
Seasonal_lw_SST = detrend(Seasonal_lw_SST')'; 

TI_ew_mean = mean(monopole_index_ew, 1);
TI_lw_mean = mean(monopole_index_lw, 1);  %dependent variable

clear x
clear y

%%
%unrestricted

for i = 1:16020

	x(:,1) = Seasonal_ew_SST(i,:);
	x(:,2) = TI_ew_mean;

    beta = zeros(2,1);
    
	fun = inline('b(1)*x(:,1)+b(2)*x(:,2)', 'b', 'x');
	
	[unres_cov_coeff(:,i), r_unres] = (nlinfit(x,TI_lw_mean', fun, beta));
    
    rss_unres(:,i) = sum(r_unres.^2);
	
% restricted
	
	beta = zeros(1,1);
	fun = inline('b(1)*x(:,2)', 'b', 'x(:,2)');
	
	[res_cov_coeff(:,i), r_res] = (nlinfit(x(:,2),TI_lw_mean', fun, beta));
	
    rss_res(:,i) = sum(r_res.^2);
	
	ru_rr(:,i) = (sum(r_unres.^2)) - sum((r_res.^2));


end

%%
for j = 1:16020
    
    gc_KI_GPH(:,j) = (rss_res(:,j) - rss_unres(:,j))/(rss_unres(:,j)/(71-2));
    
end

%% significance
% the test statistic is evaluated against an F distribution with s (1) and t-k
% (71-2)
% degrees of freedom in the numerator and denominator
% p = fcdf(x,v1 (numerator),v2 (denominator))

for k = 1:16020;
    
    p_val(:,k) = fcdf(gc_KI_GPH(:,k),1, (71-2), 'upper'); 
    
end

mask = p_val < 0.05;
mask = double(mask);
mask(mask == 0) = NaN;

%%

%%
lons = ncread('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc', 'lon');
lats = ncread('/project/pdpanalysis/nish/data/NOAA_Reconstructed_SST_V3/Monthly_Files/ersst.185401.nc', 'lat');


LA = reshape(double(lats*ones([1 180])), 89, 180);
LO = reshape(double(lons*ones([1 89])), 180, 89);

GC_KE = reshape(gc_KI_GPH, 180, 89);
shad = reshape(ru_rr, 180, 89);
mask_sig = reshape(mask, 180, 89);

%%
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 6 4]);
figure(4)

%ax(2) = subplot('position', [0.01 0 0.98 00.5]);

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')
LatLim = [19.8 66];  %this is where you set the lat/lon limits you want to plot
LonLim = [119.5 240];
h = worldmap(LatLim,LonLim);
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'compass',...
           'MLineLocation', 30,...
          'MLabelLocation',30,...
          'FlineWidth', 0.5);

limit = 6.1;  %here you set the max/min value you want shaded
contourcmap([0:limit*2/255:limit],'gray')
cmap_mod_1 = flipud(precip_cmap);
colormap(gca, cmap_mod_1(128:end,:))
%caxis([-10 10]); % set the limits of the colorbar

surfm(double(lats), double(lons), GC_KE','Facecolor', 'interp');
cb2 = colorbar('FontSize',10,'FontWeight','bold', 'Location', 'eastoutside',...
    'XTick', 0: 1 :6.1);

cb2.Ruler.TickLabelFormat='%0.1f';
cb2.Box = 'off';
cb2.TickLength = 0.01;
set(cb2, 'position',[0.905 0.25 0.015 0.52] );
%ylabel(q, 'SST GC Monopole', 'FontSize',12,'FontWeight','bold')

plotm(double(LA), double(LO'), mask_sig', 'g.', 'MarkerSize', 6);


load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'Linewidth', 0.5)
	
%textm(80, 120, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%textm(80, 210, 'KE-ew GC 850hPa-lw', 'FontSize',20, 'FontWeight', 'bold');
   
%title('SST ew GC 300 Monopole lw')    
    
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)

%la = [28 28 35 35 28]; 
%lo = [140 170 170 140 140];

%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','k', 'Linewidth', 2)

%print('Nish_071422_PDP_KE_Fig7_PS_revised','-dpng', '-r1000');

print('Nish_030323_PDP_KE_Fig7_PS_revised','-dpng', '-r1000');