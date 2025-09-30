% 创建适合IEEE双栏论文的图形（3.5英寸宽）
figure('Units', 'inches', 'Position', [1, 1, 4.5, 2.6], 'Renderer', 'painters');  % 增加高度以适应完整显示

% 定义共同的参数
data_set_names = ["cement", "steelpowder", "steelmaking"];
NOFMODELS_range = 1:3;
BATCH_SIZE_range = 1:4;
colors = {'r', 'g', 'b'};
[X, Y] = meshgrid(1:4, 1:3);

% 创建左侧子图 - 训练时间
subplot(1, 2, 1);
hold on;

% 绘制训练时间曲面
for data_idx = 1:length(data_set_names)
    avg_times = zeros(length(NOFMODELS_range), length(BATCH_SIZE_range));
    for m = NOFMODELS_range
        for b = BATCH_SIZE_range
            filename = "../results/data_" + data_set_names(data_idx) + m + "ALs_" + b + "batch.mat";
            load(filename);
            avg_times(m,b) = sum(result.time);
        end
    end
    s = surf(X, Y, avg_times);
    s.FaceAlpha = 0.7;
    s.FaceColor = colors{data_idx};
    s.EdgeColor = colors{data_idx};
end

% 设置左侧子图属性
xlabel('Batch Size', 'FontSize', 8);
ylabel('NAL', 'FontSize', 8);
zlabel('Time (s)', 'FontSize', 8);
set(gca, 'XTick', 1:4, 'XTickLabel', {'1', '2', '3', '4'}, 'FontSize', 8);
set(gca, 'YTick', 1:3, 'YTickLabel', {'1', '2', '3'}, 'FontSize', 8);
zlim([0 4800]);  % 限制训练时间显示范围
view(45, 30);
grid on;

% 创建右侧子图 - 模型误差
subplot(1, 2, 2);
hold on;

% 绘制误差曲面
for data_idx = 1:length(data_set_names)
    errors = zeros(length(NOFMODELS_range), length(BATCH_SIZE_range));
    load("../data_set/dataset_" + data_set_names(data_idx) + ".mat");
    
    for m = NOFMODELS_range
        for b = BATCH_SIZE_range
            filename = "../results/data_rc_" + data_set_names(data_idx) + m + "ALs.mat";
            load(filename, "E_reduced_constraints");
            mse = mean(mean(abs(E_reduced_constraints - E_primal_days_cv).^2));
            nrmse = sqrt(mse) / max(E_primal_days_cv(:, 1));
            errors(m,b) = nrmse;
        end
    end
    
    s = surf(X, Y, errors * 100);
    s.FaceAlpha = 0.7;
    s.FaceColor = colors{data_idx};
    s.EdgeColor = colors{data_idx};
end

% 设置右侧子图属性
xlabel('Batch Size', 'FontSize', 8);
ylabel('NAL', 'FontSize', 8);
zlabel('NRMSE (%)', 'FontSize', 8);
set(gca, 'XTick', 1:4, 'XTickLabel', {'1', '2', '3', '4'}, 'FontSize', 8);
set(gca, 'YTick', 1:3, 'YTickLabel', {'1', '2', '3'}, 'FontSize', 8);
view(45, 30);
grid on;

% 添加共同的图例，放在图形上
h = legend(data_set_names, 'Orientation', 'horizontal', ...
    'Position', [0.25 0.95 0.5 0.05], ...
    'FontSize', 8, 'Interpreter', 'none');

% 调整图形尺寸和纸张大小
figureUnits = 'inches';  % 改为英寸单位
figureWidth = 4.5;       % IEEE单栏宽度3.5英寸
figureHeight = 2.5;      % 按比例调整高度
set(gcf, 'Units', figureUnits, 'Position', [2 2 figureWidth figureHeight]);
set(gcf, 'PaperUnits', figureUnits);
set(gcf, 'PaperSize', [figureWidth figureHeight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight]);

% 调整子图间距
tight_spacing = 0.12;
pos1 = get(subplot(1,2,1), 'Position');
pos2 = get(subplot(1,2,2), 'Position');
pos1(1) = tight_spacing;
pos2(1) = 1 - tight_spacing - pos2(3);
set(subplot(1,2,1), 'Position', pos1);
set(subplot(1,2,2), 'Position', pos2);

% 使用print函数导出高质量矢量PDF
print(gcf, 'visualization/combined_results.pdf', '-dpdf', '-painters', '-r600');

