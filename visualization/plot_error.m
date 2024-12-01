% 创建图形比较不同数据集、不同设置下的模型误差
data_set_names = ["cement", "steelpowder", "steelmaking"];
NOFMODELS_range = 1:3;
BATCH_SIZE_range = 1:4;

figure('Position', [100, 100, 800, 600]);  % 调整图形大小

% 创建网格数据
[X, Y] = meshgrid(1:4, 1:3);

% 定义三个不同的颜色
colors = {'r', 'g', 'b'};
hold on;

% 为每个数据集创建一个曲面
for data_idx = 1:length(data_set_names)
    % 准备存储每种配置的误差
    errors = zeros(length(NOFMODELS_range), length(BATCH_SIZE_range));
    
    % 加载数据集
    load("data_set/dataset_" + data_set_names(data_idx) + ".mat");
    
    % 计算不同配置下的误差
    for m = NOFMODELS_range
        for b = BATCH_SIZE_range
            % 加载训练结果
            filename = "results/data_rc_" + data_set_names(data_idx) + m + "ALs.mat";
            load(filename, "E_reduced_constraints");
            
            % 计算NRMSE
            mse = mean(mean(abs(E_reduced_constraints - E_primal_days_cv).^2));
            nrmse = sqrt(mse) / max(E_primal_days_cv(:, 1));
            errors(m,b) = nrmse;
        end
    end
    
    % 创建三维曲面图
    s = surf(X, Y, errors * 100);
    s.FaceAlpha = 0.7;  % 设置透明度
    s.FaceColor = colors{data_idx};
    s.EdgeColor = colors{data_idx};
end

% 设置图形属性
title('Model Error Comparison Across Datasets');
xlabel('Batch Size');
ylabel('Number of ALs');
zlabel('NRMSE (%)');
set(gca, 'XTick', 1:4, 'XTickLabel', {'1', '2', '3', '4'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'1', '2', '3'});

% 调整视角
view(45, 30);

% 添加网格线
grid on;

% 添加图例
legend(data_set_names, 'Location', 'northeast', 'Interpreter', 'none');

hold off;
