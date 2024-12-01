%% 用钢厂数据画
load("..\data_set\dataset_cement.mat");
% 创建新图
figure;

linewidth = 1.5;

for idx = 1:5

load_profile = E_primal_days_cv(:, idx); % 第一个时间序列

% 绘制第一个时间序列
plot3(1:24, idx * ones(size(load_profile)), load_profile, 'linewidth', linewidth); 

hold on;

end

fontsize = 25;
% Set axis labels
x1 = xlabel('Hour','FontSize',fontsize,'FontName', 'Times New Roman','FontWeight','bold');
y1 = ylabel('Day','FontSize',fontsize,'FontName', 'Times New Roman','FontWeight','bold');
z1 = zlabel('Load (MW)','FontSize',fontsize,'FontName', 'Times New Roman','FontWeight','bold');
% Figure size
figureUnits = 'centimeters';
figureWidth = 20;
figureHeight = figureWidth * 3 / 4;
set(gcf, 'Units', figureUnits, 'Position', [10 10 figureWidth figureHeight]);

ax.FontName = 'Times New Roman';
% ax.FontSize = fontsize;
% set(gcf, 'PaperSize', [18, 8]);

% 设置坐标轴刻度字体大小
set(gca, 'FontSize', fontsize);

% y1 刻度化为整数
set(gca, 'YTick', 1:5);

% 设置坐标轴
axis tight;

% 关闭绘图保持功能
hold off;

saveas(gcf,'typical_load_data.pdf');
