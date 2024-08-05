%% Produce the plots
load HVD_Cop_MS_EOF_1;
load HVD_Cop_MS_KI;
load HVD_Cop_LS_EOF_1;
load HVD_Cop_LS_EOF_2;

x = NaN(1,38);

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 22 12]);
load dates_vec_5.mat;
dates_vec = reshape(dates_vec, 6, 38);
dates_vec = vertcat(dates_vec, x);
dates_vec = reshape(dates_vec, 7*38, 1);

HVD_Cop_MS_EOF_1 = reshape(HVD_Cop_MS_EOF_1, 6,38);
HVD_Cop_MS_EOF_1 = vertcat(HVD_Cop_MS_EOF_1, x);
HVD_Cop_MS_EOF_1 = reshape(HVD_Cop_MS_EOF_1, 1, 7*38);

HVD_Cop_MS_KI = reshape(HVD_Cop_MS_KI, 6,38);
HVD_Cop_MS_KI = vertcat(HVD_Cop_MS_KI, x);
HVD_Cop_MS_KI = reshape(HVD_Cop_MS_KI, 1, 7*38);

HVD_Cop_LS_EOF_1 = reshape(HVD_Cop_LS_EOF_1, 6,38);
HVD_Cop_LS_EOF_1 = vertcat(HVD_Cop_LS_EOF_1, x);
HVD_Cop_LS_EOF_1 = reshape(HVD_Cop_LS_EOF_1, 1, 7*38);

HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6,38);
HVD_Cop_LS_EOF_2 = vertcat(HVD_Cop_LS_EOF_2, x);
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 1, 7*38);

%%
tiledlayout(2,2,'TileSpacing','Compact', 'Padding', 'compact');

nexttile
a = plot(HVD_Cop_MS_KI, '-', 'Color', '#0072BD', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',dates_vec(1:63:266));
ax = gca;
ax.FontSize = 20; 
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([a],{'MS KI'},...
    'Location','northeast');
%title('KE Cop MS KI');

nexttile
b = plot(HVD_Cop_MS_EOF_1, '-', 'Color', '#e41a1c', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',dates_vec(1:63:266));
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([b],{'MS EOF 1'},...
    'Location','northeast');
%title('KE Cop MS EOF 1');

nexttile
c = plot(HVD_Cop_LS_EOF_1, '-', 'Color', '#4daf4a', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',dates_vec(1:63:266));
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([c],{'LS EOF 1'},...
    'Location','northeast');
%title('KE Cop LS EOF 1');

nexttile
d = plot(HVD_Cop_LS_EOF_2, '-', 'Color', '#984ea3', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',dates_vec(1:63:266));
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([d],{'LS EOF 2'},...
    'Location','northeast');

print('Nish_072821_KE_PDP_MS_Fig_2','-dpng', '-r500');%% Produce the plots
load HVD_Cop_MS_EOF_1;
load HVD_Cop_MS_KI;
load HVD_Cop_LS_EOF_1;
load HVD_Cop_LS_EOF_2;

x = NaN(1,38);

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 22 12]);
load dates_vec_5.mat;
dates_vec = reshape(dates_vec, 6, 38);
dates_vec = vertcat(dates_vec, x);
dates_vec = reshape(dates_vec, 7*38, 1);

HVD_Cop_MS_EOF_1 = reshape(HVD_Cop_MS_EOF_1, 6,38);
HVD_Cop_MS_EOF_1 = vertcat(HVD_Cop_MS_EOF_1, x);
HVD_Cop_MS_EOF_1 = reshape(HVD_Cop_MS_EOF_1, 1, 7*38);

HVD_Cop_MS_KI = reshape(HVD_Cop_MS_KI, 6,38);
HVD_Cop_MS_KI = vertcat(HVD_Cop_MS_KI, x);
HVD_Cop_MS_KI = reshape(HVD_Cop_MS_KI, 1, 7*38);

HVD_Cop_LS_EOF_1 = reshape(HVD_Cop_LS_EOF_1, 6,38);
HVD_Cop_LS_EOF_1 = vertcat(HVD_Cop_LS_EOF_1, x);
HVD_Cop_LS_EOF_1 = reshape(HVD_Cop_LS_EOF_1, 1, 7*38);

HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6,38);
HVD_Cop_LS_EOF_2 = vertcat(HVD_Cop_LS_EOF_2, x);
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 1, 7*38);

%%
tiledlayout(2,2,'TileSpacing','Compact', 'Padding', 'compact');

nexttile
a = plot(HVD_Cop_MS_KI, '-', 'Color', '#0072BD', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',{[]});
ax = gca;
ax.FontSize = 18; 
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([a],{'MS KI'},'FontSize', 18,...
    'Location','northeast');
legend boxoff;
text(5, 3.25, '(a)', 'FontSize',20, 'FontWeight', 'bold');
%title('KE Cop MS KI');

nexttile
b = plot(HVD_Cop_MS_EOF_1, '-', 'Color', '#e41a1c', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',{[]});
ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([b],{'MS EOF 1'},'FontSize', 18,...
    'Location','northeast');
legend boxoff;
text(5, 3.25, '(b)', 'FontSize',20, 'FontWeight', 'bold');
%title('KE Cop MS EOF 1');

nexttile
c = plot(HVD_Cop_LS_EOF_1, '-', 'Color', '#4daf4a', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',dates_vec(1:63:266));
ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([c],{'LS EOF 1'},'FontSize', 18,...
    'Location','northeast');
legend boxoff;
text(5, 3.25, '(c)', 'FontSize',20, 'FontWeight', 'bold');
%title('KE Cop LS EOF 1');

nexttile
d = plot(HVD_Cop_LS_EOF_2, '-', 'Color', '#984ea3', 'LineWidth',1.5);
xlim([0 228]);
ylim([-4 4]);
set(gca,'xtick',(1:63:266),'xticklabel',dates_vec(1:63:266));
ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold'; 
ax.LineWidth = 2;
%set(gca,'xtick',(1:5:185));
%set(gca,'XTickLabel',[]);
%grid on
hold on
yline(0,'--', 'LineWidth',1.5);
legend([d],{'LS EOF 2'},'FontSize', 18,...
    'Location','northeast');
legend boxoff;
text(5, 3.25, '(d)', 'FontSize',20, 'FontWeight', 'bold');

%print('Nish_080521_KE_PDP_MS_Fig_2','-dpng', '-r500');

print('Nish_071922_KE_PDP_MS_Fig2_PS_revised','-dpng', '-r500');