% 创建适合IEEE单栏格式的图形
close;
figureUnits = 'inches';
figureWidth = 3.5;  % 压缩宽度
figureHeight = 1.5;  % 高度保持不变
figure('Units', figureUnits, 'Position', [2 2 figureWidth figureHeight]);

% 定义颜色和线型组合
colors = {'r', 'g'};  % ALs=1,2
line_styles = {'-', '--', ':'};  % batch=1,2,3
window_size = 1;  % 移动平均窗口大小
max_iterations = 40;  % 只显示前40次迭代
data_set_name = "steelpowder";

% 创建两个子图
for subplot_idx = 1:2
    subplot(1, 2, subplot_idx);  % 修改为1行2列
    hold on;
    legend_entries = {};
    
    % 遍历选定的NAL和batch组合
    for NOFMODELS = 1:2
        for BATCH_SIZE = 1:3
            % 加载结果
            filename = "results/data_" + data_set_name + ...
                      NOFMODELS + "ALs_" + BATCH_SIZE + "batch.mat";
            load(filename);
            
            % 选择数据并限制迭代次数
            if subplot_idx == 1
                data = result.P_max_i(1:min(max_iterations, length(result.P_max_i)), 1);
            else
                data = result.E_max_i(1:min(max_iterations, length(result.E_max_i)), 1);
            end
            
            % 计算移动平均
            smooth_data = movmean(data, window_size);
            
            % 绘制曲线
            plot(smooth_data, 'Color', colors{NOFMODELS}, ...
                 'LineStyle', line_styles{BATCH_SIZE}, 'LineWidth', 1);
            
            % 添加图例条目
            legend_entries{end+1} = sprintf('NAL=%d,B=%d', NOFMODELS, BATCH_SIZE);
        end
    end
    
    % 设置子图属性
    xlabel('Iteration', 'FontSize', 5, 'Rotation', 0);
    if subplot_idx == 1
        ylabel('P^{max}', 'FontSize', 5, 'Rotation', 90);
        ylim([0 250]);  % 设置P_max的y轴范围
    else
        ylabel('E^{max}', 'FontSize', 5, 'Rotation', 90);
        ylim([0 4000]);  % 设置E_max的y轴范围 
    end
    grid on;
    set(gca, 'FontSize', 8);
end

% 在图形上方添加图例
h = legend(legend_entries, 'FontSize', 5.5, ...
    'Position', [0.25 0.95 0.5 0.05], ...
    'Interpreter', 'none', 'NumColumns', 3);

% 调整图形尺寸和保存
set(gcf, 'PaperUnits', figureUnits);
set(gcf, 'PaperSize', [figureWidth+5 figureHeight+5]);  % 稍微增加高度以容纳图例
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figureWidth+2 figureHeight+2]);

% 调整子图间距
tight_spacing = 0.09;
pos1 = get(subplot(1,2,1), 'Position');
pos2 = get(subplot(1,2,2), 'Position');
pos1(1) = tight_spacing;
pos2(1) = 1 - tight_spacing - pos2(3);
set(subplot(1,2,1), 'Position', pos1);
set(subplot(1,2,2), 'Position', pos2);

% 保存图形
print('visualization/convergence_analysis_steelpowder', '-dpdf', '-r300');
