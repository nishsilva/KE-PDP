set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 20 5.75]);
tiledlayout(1, 4);

load PDP_R1_full.txt;

load NCEP_LS_Index.mat;
LS_EOF_1 = reshape(NCEP_LS_Index, 6, 71);
Seasonal_LS_EOF_1 = mean(LS_EOF_1, 1);
Seasonal_LS_EOF_1 = detrend(Seasonal_LS_EOF_1')';

load HVD_ERSST_EOF2.mat;
LS_EOF_2 = reshape(HVD_ERSST_EOF2, 6, 71);
Seasonal_LS_EOF_2 = mean(LS_EOF_2, 1);
Seasonal_LS_EOF_2 = detrend(Seasonal_LS_EOF_2')';

fs = 1; 
N=2;
N20yr = 1/20;
N7yr = 1/7;
[b,a] = butter(N,[N20yr*2./fs N7yr*2./fs]);
t_7_20yr = filter(b,a,Seasonal_LS_EOF_2',[],1);  % here t_anom is the time-series; t_7_20yr is the bandpass filtered time-series
t_7_20yr = normalize(detrend(t_7_20yr));
%Seasonal_LS_EOF_2 = smooth(Seasonal_LS_EOF_2', 3)';

load CP_LS_Index.mat;
CP_EOF_2 = reshape(CP_LS_Index, 6, 71);
Seasonal_CP_EOF_2 = mean(CP_EOF_2, 1);
Seasonal_CP_EOF_2 = detrend(Seasonal_CP_EOF_2')';

Years = 1948:2018;
Years_vec = string(Years);

nexttile([1,3])
x = plot(PDP_R1_full(:,2), '--', 'Color', '#0072BD', 'LineWidth',2.5);
yline(0,'--', 'LineWidth',1.5);
xlim([0 71]);
ylim([-3 3]);
hold on
a = plot(t_7_20yr, '--', 'Color', '#D95319', 'LineWidth',2.5);
%b = plot(Seasonal_KEE_EOF_3, 'Color', '#33a02c', 'LineWidth',1.5);
%c = plot(Seasonal_CP_EOF_2, 'Color', '#3182bd', 'LineWidth',1.5);
set(gca,'xtick',(1:10:72),'xticklabel',Years_vec(1:10:72));
ax = gca;
ax.FontSize = 22; 
ax.LineWidth = 2;
ax.FontWeight = 'bold';
%grid on
%legend([x, a, b, c],{'PDP','KE Index', 'CP EOF2', 'KE EOF 3'},...
%    'Location','NorthEast');
%legend([b, c],{'KE EOF 3', 'CP EOF2'},...
%    'Location','NorthEast');
legend([x, a],{'PDP', 'LS EOF 2'},...
    'Location','NorthWest');
legend boxoff;
hold on;
xlabel('Year', 'FontSize', 22, 'FontWeight', 'bold');
title('(a)', 'FontSize', 22, 'FontWeight', 'bold')	


nexttile

[c,lags] = xcorr(PDP_R1_full(4:end, 2), t_7_20yr(4:65, :), 5, 'normalized');
stem(lags,c, 'filled', 'LineWidth', 2);
ylim([-0.65 0.65]);
hold on
vline(0, '--', 'Color', '#bdbdbd');
text(-4, -0.55, 'PDP leads', 'FontSize',15, 'FontWeight', 'bold')
text(1.5, 0.55, 'KI leads', 'FontSize',15, 'FontWeight', 'bold')
set(gca,'YTickLabel',[]);
yyaxis right
ylim([-0.65 0.65]);
ax = gca;
ax.FontSize = 20; 
ax.LineWidth = 2;
ax.FontWeight = 'bold';
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
xlabel('Lag Year', 'FontSize', 22, 'FontWeight', 'bold');
title('(b)', 'FontSize', 22, 'FontWeight', 'bold')	
%yyaxis right
%print('Nish_093021_KE_PDP_MS_Fig_7_PS','-dpng', '-r300');