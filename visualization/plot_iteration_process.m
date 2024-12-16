% 创建适合IEEE双栏论文的图形
figure('Units', 'inches', 'Position', [1, 1, 4.5, 2.6], 'Renderer', 'painters');  % 增加高度以适应完整显示

% 定义颜色和线型组合
colors = {'r', 'g'};  % ALs=1,2
line_styles = {'-', '--', ':'};  % batch=1,2,3
window_size = 1;  % 移动平均窗口大小
max_iterations = 40;  % 只显示前40次迭代
data_set_name = "steelpowder";

% 创建两个子图
for subplot_idx = 1:2
    subplot(1, 2, subplot_idx);  
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
            
            % 计算移动平���
            smooth_data = movmean(data, window_size);
            
            % 绘制曲线
            plot(smooth_data, 'Color', colors{NOFMODELS}, ...
                 'LineStyle', line_styles{BATCH_SIZE}, 'LineWidth', 1);
            
            % 添加图例条目
            legend_entries{end+1} = sprintf('NAL=%d,B=%d', NOFMODELS, BATCH_SIZE);
        end
    end
    
    % 统一设置子图属性
    xlabel('Iteration', 'FontSize', 8);
    if subplot_idx == 1
        ylabel('P^{max}', 'FontSize', 8);
        ylim([0 250]);  
    else
        ylabel('E^{max}', 'FontSize', 8);
        ylim([0 4000]);  
    end
    grid on;
    set(gca, 'FontSize', 8);
    set(gca, 'XTick', 0:10:40, 'XTickLabel', {'0','10','20','30','40'});
end

% 添加共同的图例
h = legend(legend_entries, 'Orientation', 'horizontal', ...
    'Position', [0.25 0.95 0.5 0.05], ...
    'NumColumns', 3, ...
    'FontSize', 8, 'Interpreter', 'none');

% 调整图形尺寸和纸张大小
figureUnits = 'inches';  
figureWidth = 4.5;       % IEEE双栏宽度
figureHeight = 2.6;      
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
print(gcf, 'visualization/convergence_analysis_steelpowder', '-dpdf', '-painters', '-r600');
