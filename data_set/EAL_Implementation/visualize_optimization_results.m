function visualize_optimization_results(solution, params)
% 可视化函数：显示边界收缩算法优化结果
% 输入：
%   solution - 优化结果结构体
%   params - 参数结构体

fprintf('\n=== 开始可视化优化结果 ===\n');

% 时间轴设置
time_hr = 1:params.T;
time_state = (1:params.T) * params.delta_t; % 1小时间隔的时间点

% 第一张图：控制变量和输出变量
figure('Name', '边界收缩算法优化结果 - 控制变量和输出变量', 'Position', [100, 100, 1400, 900]);

% 电价曲线
subplot(2,3,1);
stairs(time_hr, params.price_TOU, 'r-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('电价 (元/kWh)');
title('分时电价曲线');
grid on;

% 电流曲线
subplot(2,3,2);
stairs(time_hr, solution.I, 'b-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('电流 (kA)');
title('优化电流曲线');
grid on;

% 功率曲线
subplot(2,3,3);
stairs(time_hr, solution.P, 'g-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('功率 (kW)');
title('优化功率曲线');
grid on;

% 氧化铝进料曲线
subplot(2,3,4);
stairs(time_hr, solution.m_feed_Al2O3, 'm-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('氧化铝进料速率 (kg/min)');
title('氧化铝进料优化曲线');
grid on;

% 电流效率曲线
subplot(2,3,5);
plot(time_state, solution.g, 'c-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('电流效率');
title('电流效率变化');
grid on;
ylim([0.85, 0.98]);

% 铝产量曲线
subplot(2,3,6);
plot(time_state, solution.Y_Al, 'k-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('铝产量 (kg/h)');
title('铝产量曲线');
grid on;

% 第二张图：状态变量
figure('Name', '边界收缩算法优化结果 - 状态变量', 'Position', [200, 200, 1400, 1000]);

% 电解质温度曲线
subplot(3,3,1);
plot(time_state, solution.T_B(1:end-1), 'r-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('电解质温度 (°C)');
title('电解质温度变化');
grid on;

% 液相温度曲线
subplot(3,3,2);
plot(time_state, solution.T_liq, 'b-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('液相温度 (°C)');
title('液相温度变化');
grid on;

% 过热度曲线
subplot(3,3,3);
plot(time_state, solution.Superheat, 'g-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('过热度 (°C)');
title('过热度变化');
grid on;

% 各组分浓度变化曲线
subplot(3,3,4);
Al2O3_conc = solution.C_Al2O3(1:end-1) * 100;
AlF3_conc = solution.C_AlF3(1:end-1) * 100;
CaF2_conc = solution.C_CaF2(1:end-1) * 100;
LiF_conc = solution.C_LiF(1:end-1) * 100;

plot(time_state, Al2O3_conc, 'm-', 'LineWidth', 2, 'DisplayName', 'Al2O3');
hold on;
plot(time_state, AlF3_conc, 'r-', 'LineWidth', 2, 'DisplayName', 'AlF3');
plot(time_state, CaF2_conc, 'g-', 'LineWidth', 2, 'DisplayName', 'CaF2');
plot(time_state, LiF_conc, 'b-', 'LineWidth', 2, 'DisplayName', 'LiF');
xlabel('时间 (h)');
ylabel('浓度 (wt%)');
title('各组分浓度变化');
legend('show', 'Location', 'best');
grid on;
hold off;

% 浴液质量曲线
subplot(3,3,5);
plot(time_state, solution.m_B(1:end-1), 'c-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('浴液质量 (kg)');
title('浴液质量变化');
grid on;

% Ledge厚度曲线
subplot(3,3,6);
plot(time_state, solution.L_ledge(1:end-1), 'k-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('Ledge厚度 (m)');
title('Ledge厚度变化');
grid on;
    
% 温度差分曲线
subplot(3,3,7);
plot(time_state, solution.dT_values, 'r-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('温度差分 (°C)');
title('温度差分变化');
grid on;

% 质量差分曲线
subplot(3,3,8);
plot(time_state, solution.dm_values, 'b-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('质量差分 (kg)');
title('质量差分变化');
grid on;

% McCormick包络误差
subplot(3,3,9);
bar([solution.mccormick_errors.mB_dT_mean_error, solution.mccormick_errors.Lledge_dm_mean_error], 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'XTickLabel', {'温度差分×浴液质量', 'Ledge厚度×质量差分'});
ylabel('平均误差 (%)');
title('McCormick包络误差');
grid on;

% 第三张图：McCormick包络误差详细分析
figure('Name', 'McCormick包络误差详细分析', 'Position', [300, 300, 1400, 800]);

% 温度差分×浴液质量误差
subplot(2,3,1);
plot(time_state, solution.mccormick_errors.mB_dT_error, 'r-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('相对误差 (%)');
title('温度差分×浴液质量McCormick误差');
grid on;

% Ledge厚度×质量差分误差
subplot(2,3,2);
plot(time_state, solution.mccormick_errors.Lledge_dm_error, 'b-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('相对误差 (%)');
title('Ledge厚度×质量差分McCormick误差');
grid on;

% 电流×电流效率误差
subplot(2,3,3);
plot(time_state, solution.mccormick_errors.I_g_error, 'g-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('相对误差 (%)');
title('电流×电流效率McCormick误差');
grid on;

% Ledge厚度×过热度误差
subplot(2,3,4);
plot(time_state, solution.mccormick_errors.Lledge_superheat_error, 'k-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('相对误差 (%)');
title('Ledge厚度×过热度McCormick误差');
grid on;

% 各组分浓度×浴液质量误差（合并在一张图）
subplot(2,3,5);
plot(time_state, solution.mccormick_errors.C_Al2O3_m_error(1:end-1), 'm-', 'LineWidth', 2, 'DisplayName', 'Al2O3');
hold on;
plot(time_state, solution.mccormick_errors.C_AlF3_m_error(1:end-1), 'r-', 'LineWidth', 2, 'DisplayName', 'AlF3');
plot(time_state, solution.mccormick_errors.C_CaF2_m_error(1:end-1), 'g-', 'LineWidth', 2, 'DisplayName', 'CaF2');
plot(time_state, solution.mccormick_errors.C_LiF_m_error(1:end-1), 'b-', 'LineWidth', 2, 'DisplayName', 'LiF');
xlabel('时间 (h)');
ylabel('相对误差 (%)');
title('各组分浓度×浴液质量McCormick误差');
legend('show', 'Location', 'best');
grid on;
hold off;

% 第四张图：McCormick包络误差统计对比
figure('Name', 'McCormick包络误差统计对比', 'Position', [400, 400, 1200, 600]);

% 平均误差统计
subplot(1,2,1);
error_stats = [solution.mccormick_errors.mB_dT_mean_error, solution.mccormick_errors.Lledge_dm_mean_error, ...
               solution.mccormick_errors.I_g_mean_error, solution.mccormick_errors.Lledge_superheat_mean_error, ...
               solution.mccormick_errors.C_Al2O3_m_mean_error, solution.mccormick_errors.C_AlF3_m_mean_error, ...
               solution.mccormick_errors.C_CaF2_m_mean_error, solution.mccormick_errors.C_LiF_m_mean_error];
error_labels = {'温度差分×浴液质量', 'Ledge厚度×质量差分', '电流×电流效率', 'Ledge厚度×过热度', ...
                'Al2O3浓度×浴液质量', 'AlF3浓度×浴液质量', 'CaF2浓度×浴液质量', 'LiF浓度×浴液质量'};
bar(error_stats, 'FaceColor', [0.6, 0.8, 0.4]);
set(gca, 'XTickLabel', error_labels);
ylabel('平均相对误差 (%)');
title('McCormick包络平均误差统计');
grid on;
xtickangle(45);

% 最大误差统计
subplot(1,2,2);
max_error_stats = [solution.mccormick_errors.mB_dT_max_error, solution.mccormick_errors.Lledge_dm_max_error, ...
                   solution.mccormick_errors.I_g_max_error, solution.mccormick_errors.Lledge_superheat_max_error, ...
                   solution.mccormick_errors.C_Al2O3_m_max_error, solution.mccormick_errors.C_AlF3_m_max_error, ...
                   solution.mccormick_errors.C_CaF2_m_max_error, solution.mccormick_errors.C_LiF_m_max_error];
bar(max_error_stats, 'FaceColor', [0.8, 0.4, 0.4]);
set(gca, 'XTickLabel', error_labels);
ylabel('最大相对误差 (%)');
title('McCormick包络最大误差统计');
grid on;
xtickangle(45);

fprintf('可视化完成！\n');
end
