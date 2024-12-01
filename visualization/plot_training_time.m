% 创建图形比较不同数据集、不同设置下的训练时间
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
    % 准备存储每种配置的平均训练时间
    avg_times = zeros(length(NOFMODELS_range), length(BATCH_SIZE_range));
    
    % 加载并计算平均训练时间
    for m = NOFMODELS_range
        for b = BATCH_SIZE_range
            % 加载结果文件
            filename = "results/data_" + data_set_names(data_idx) + m + "ALs_" + b + "batch.mat";
            load(filename);
            % 计算总训练时间
            avg_times(m,b) = sum(result.time);
        end
    end
    
    % 创建三维曲面图
    s = surf(X, Y, avg_times);
    s.FaceAlpha = 0.7;  % 设置透明度
    s.FaceColor = colors{data_idx};
    s.EdgeColor = colors{data_idx};
end

% 设置图形属性
title('Training Time Comparison Across Datasets');
xlabel('Batch Size');
ylabel('Number of ALs');
zlabel('Training Time (s)');
set(gca, 'XTick', 1:4, 'XTickLabel', {'1', '2', '3', '4'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'1', '2', '3'});

% 最大显示到10000
zlim([0 10000]);

% 调整视角
view(45, 30);

% 添加网格线
grid on;

% 添加图例
legend(data_set_names, 'Location', 'northeast', 'Interpreter', 'none');

hold off;