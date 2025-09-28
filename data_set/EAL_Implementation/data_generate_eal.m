%% 使用电解铝负荷优化模型生成多日数据用于逆优化
% 假设内部参数每天相同，但外部电价变化，导致不同的生产调度方案
% 参考STN和RTN的数据生成方法
clc; clear; yalmip("clear"); close all

%% 读取电价数据
% 使用与STN和RTN相同的PJM电价数据
data_price;

%% 生成用电数据
E_primal_days = [];

fprintf('开始生成电解铝负荷历史数据...\n');
fprintf('总天数: %d\n', length(Price_days));

for idx_day = 1 : length(Price_days)
    fprintf('正在处理第 %d/%d 天...\n', idx_day, length(Price_days));
    
    % 读取当天电价数据
    Price = Price_days(:, idx_day);
    
    % 清除之前的变量
    yalmip('clear');
    
    % 运行电解铝优化模型
    cell_linear_optimization_tightening;
    
    % 检查求解是否成功
    if feasible
        % 提取最优功率方案 (kW)
        % P_opt = value(P_input_opt);
        
        % 转换为MW并记录用电数据
        E_primal = P_opt' / 1000; % 转换为MW
        
        % 累积历史数据
        E_primal_days = [E_primal_days, E_primal];
        
        fprintf('  第%d天求解成功，总功率: %.2f MW\n', idx_day, sum(E_primal));
    else
        fprintf('  第%d天求解失败\n', idx_day);
        % 如果求解失败，使用零功率或前一天的功率
        if idx_day == 1
            E_primal = zeros(24, 1); % 第一天失败使用零功率
        else
            E_primal = E_primal_days(:, end); % 使用前一天的功率
        end
        E_primal_days = [E_primal_days, E_primal];
    end
end

%% 区分训练集和交叉验证集
% 前21天用于训练，后10天用于交叉验证
E_primal_days_train = E_primal_days(:, 1 : 21);
E_primal_days_cv = E_primal_days(:, 22 : end);
Price_days_train = Price_days(:, 1 : 21);
Price_days_cv = Price_days(:, 22 : end);

%% 保存数据
save("dataset_eal.mat", "E_primal_days_train", "Price_days_train", ...
    "E_primal_days_cv", "Price_days_cv");

fprintf('\n数据生成完成！\n');
fprintf('训练集: %d天\n', size(E_primal_days_train, 2));
fprintf('验证集: %d天\n', size(E_primal_days_cv, 2));
fprintf('数据已保存到: dataset_eal.mat\n');

%% 数据统计
fprintf('\n=== 数据统计 ===\n');
fprintf('训练集平均日用电量: %.2f MWh\n', mean(sum(E_primal_days_train)));
fprintf('训练集最大日用电量: %.2f MWh\n', max(sum(E_primal_days_train)));
fprintf('训练集最小日用电量: %.2f MWh\n', min(sum(E_primal_days_train)));
fprintf('验证集平均日用电量: %.2f MWh\n', mean(sum(E_primal_days_cv)));
fprintf('验证集最大日用电量: %.2f MWh\n', max(sum(E_primal_days_cv)));
fprintf('验证集最小日用电量: %.2f MWh\n', min(sum(E_primal_days_cv)));

%% 可视化部分数据
figure('Name', '电解铝负荷历史数据示例', 'Position', [100, 100, 1200, 800]);

% 显示前5天的数据
subplot(2,2,1);
plot(1:24, E_primal_days_train(:, 1:5), 'LineWidth', 1.5);
xlabel('时间 (h)');
ylabel('功率 (MW)');
title('训练集前5天用电曲线');
legend('第1天', '第2天', '第3天', '第4天', '第5天', 'Location', 'best');
grid on;

% 显示对应的电价
subplot(2,2,2);
plot(1:24, Price_days_train(:, 1:5), 'LineWidth', 1.5);
xlabel('时间 (h)');
ylabel('电价 (元/kWh)');
title('训练集前5天电价曲线');
legend('第1天', '第2天', '第3天', '第4天', '第5天', 'Location', 'best');
grid on;

% 日总用电量分布
subplot(2,2,3);
histogram(sum(E_primal_days_train), 10, 'FaceAlpha', 0.7);
xlabel('日总用电量 (MWh)');
ylabel('频次');
title('训练集日总用电量分布');
grid on;

% 电价与用电量关系
subplot(2,2,4);
scatter(mean(Price_days_train), sum(E_primal_days_train), 50, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('平均电价 (元/kWh)');
ylabel('日总用电量 (MWh)');
title('电价与用电量关系');
grid on;

fprintf('\n数据生成脚本执行完成！\n');
