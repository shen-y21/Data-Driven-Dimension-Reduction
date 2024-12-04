% 创建适合IEEE单栏格式的图形
figureUnits = 'inches';
figureWidth = 3.5;
figureHeight = 2.5;
figure('Units', figureUnits, 'Position', [2 2 figureWidth figureHeight]);

% 定义颜色和线型组合，为不同的NAL和batch组合
colors = {'r', 'g', 'b'};
line_styles = {'-', '--', ':', '-.'};
window_size = 5;  % 移动平均窗口大小
max_iterations = 40;  % 只显示前40次迭代

hold on;
legend_entries = {};
data_set_name = "steelpowder";

% 遍历所有NAL和batch组合
for NOFMODELS = 1:3
    for BATCH_SIZE = 1:4
        % 加载结果
        filename = "results/data_" + data_set_name + ...
                  NOFMODELS + "ALs_" + BATCH_SIZE + "batch.mat";
        load(filename);
        
        % 获取数据并限制迭代次数
        data = result.J_theta(1:min(max_iterations, length(result.J_theta)));
        
        % 计算移动平均
        smooth_data = movmean(data, window_size);
        
        % 绘制曲线
        plot(smooth_data, 'Color', colors{NOFMODELS}, ...
             'LineStyle', line_styles{BATCH_SIZE}, 'LineWidth', 1);
        
        % 添加图例条目
        legend_entries{end+1} = sprintf('NAL=%d,B=%d', NOFMODELS, BATCH_SIZE);
    end
end

% 设置图形属性
xlabel('Iteration', 'FontSize', 8, 'Rotation', 0);
ylabel('J_{\theta}', 'FontSize', 8, 'Rotation', 90);
grid on;
set(gca, 'FontSize', 8);

% 添加图例
legend(legend_entries, 'FontSize', 7, 'Location', 'northeast', ...
       'Interpreter', 'none', 'NumColumns', 2);

% 调整图形尺寸和保存
set(gcf, 'PaperUnits', figureUnits);
set(gcf, 'PaperSize', [figureWidth figureHeight+0.3]);  % 稍微增加高度以容纳图例
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figureWidth figureHeight+0.3]);

% 保存图形
print('visualization/convergence_analysis_steelpowder', '-dpdf', '-r300');