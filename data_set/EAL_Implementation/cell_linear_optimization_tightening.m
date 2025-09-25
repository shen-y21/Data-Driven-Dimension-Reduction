% 电解铝负荷用能优化模型 - 边界收缩算法简化版本
% 基于分段McCormick包络的边界收缩算法实现
clc; clear; yalmip("clear"); close all

%% 基本参数设置
% 使用参数初始化函数（包含拟合结果加载）
params = f_initialize_parameters();

%% 边界收缩算法实现
fprintf('=== 开始边界收缩算法 ===\n');

% 初始化变量边界
dT_T_min_init = params.dT_T_min_init; 
dT_T_max_init = params.dT_T_max_init;    % 温度差分初始范围
dm_m_min_init = params.dm_m_min_init; 
dm_m_max_init = params.dm_m_max_init;    % 质量差分初始范围
L_ledge_min_init = params.L_ledge_min;
L_ledge_max_init = params.L_ledge_max;   % Ledge厚度初始范围
m_B_min_init = params.m_B_min_init;  % 浴液质量初始下界
m_B_max_init = params.m_B_max_init;  % 浴液质量初始上界

% 当前边界（时变边界）
dT_T_min_current = repmat(dT_T_min_init, 1, params.T);
dT_T_max_current = repmat(dT_T_max_init, 1, params.T);
dm_m_min_current = repmat(dm_m_min_init, 1, params.T);
dm_m_max_current = repmat(dm_m_max_init, 1, params.T);
L_ledge_min_current = repmat(L_ledge_min_init, 1, params.T+1);
L_ledge_max_current = repmat(L_ledge_max_init, 1, params.T+1);
m_B_min_current = repmat(m_B_min_init, 1, params.T+1);
m_B_max_current = repmat(m_B_max_init, 1, params.T+1);

% 存储迭代历史
iteration_history = [];
error_history = [];  % 存储收敛判据误差历史
boundary_history = [];  % 存储边界收缩历史（第21个时间段）
differential_history = [];  % 存储差分值历史（第21个时间段）
ledge_history = [];  % 存储Ledge厚度历史（第21个时间段）
mass_history = [];  % 存储浴液质量历史（第21个时间段）
shrink_factor_history = [];  % 存储收缩因子历史

for iter = 1:params.max_iterations
    fprintf('\n--- 迭代 %d ---\n', iter);
    fprintf('当前温度差分范围: [%.2f, %.2f]°C (平均范围: %.2f)\n', min(dT_T_min_current), max(dT_T_max_current), mean(dT_T_max_current - dT_T_min_current));
    fprintf('当前质量差分范围: [%.0f, %.0f] kg (平均范围: %.0f)\n', min(dm_m_min_current), max(dm_m_max_current), mean(dm_m_max_current - dm_m_min_current));
    fprintf('当前Ledge厚度范围: [%.3f, %.3f] m (平均范围: %.3f)\n', min(L_ledge_min_current), max(L_ledge_max_current), mean(L_ledge_max_current - L_ledge_min_current));
    fprintf('当前浴液质量范围: [%.0f, %.0f] kg (平均范围: %.0f)\n', min(m_B_min_current), max(m_B_max_current), mean(m_B_max_current - m_B_min_current));
    
    % 计算自适应收缩因子
    [shrink_factors, shrink_factor_info] = calculate_adaptive_shrink_factors(iter, error_history, params);
    fprintf('自适应收缩因子: 温度差分%.4f, 质量差分%.4f, Ledge厚度%.4f, 浴液质量%.4f\n', ...
        shrink_factors.dT, shrink_factors.dm, shrink_factors.ledge, shrink_factors.mass);
    
    % 显示第21个时间段的边界信息
    fprintf('第21个时间段边界 - 温度差分: [%.4f, %.4f]°C, 质量差分: [%.4f, %.4f] kg\n', ...
        dT_T_min_current(21), dT_T_max_current(21), dm_m_min_current(21), dm_m_max_current(21));
    fprintf('第21个时间段边界 - Ledge厚度: [%.4f, %.4f] m, 浴液质量: [%.0f, %.0f] kg\n', ...
        L_ledge_min_current(21), L_ledge_max_current(21), m_B_min_current(21), m_B_max_current(21));
    
    % 使用结构体中的分段数
    N_T_segments = params.N_T_segments;
    N_m_segments = params.N_m_segments;
    
    % 根据参考文献的迭代规则：首次迭代使用分段McCormick，后续迭代使用标准McCormick
    if iter == 1
        % 首次迭代：求解分段McCormick模型
        fprintf('使用分段McCormick包络模型\n');
        [obj_value, solution, feasible] = solve_piecewise_mccormick_optimization(...
            dT_T_min_current, dT_T_max_current, dm_m_min_current, dm_m_max_current, ...
            L_ledge_min_current, L_ledge_max_current, m_B_min_current, m_B_max_current, ...
            N_T_segments, N_m_segments, params);
    else
        % 后续迭代：求解标准McCormick松弛模型
        fprintf('使用标准McCormick松弛模型\n');
        [obj_value, solution, feasible] = solve_standard_mccormick_optimization(...
            dT_T_min_current, dT_T_max_current, dm_m_min_current, dm_m_max_current, ...
            L_ledge_min_current, L_ledge_max_current, m_B_min_current, m_B_max_current, params);
    end
    
    if ~feasible
        fprintf('当前迭代无可行解，停止算法\n');
        break;
    end
    
    % 记录当前解
    iteration_history = [iteration_history; iter, obj_value, dT_T_max_current - dT_T_min_current, dm_m_max_current - dm_m_min_current, ...
        L_ledge_max_current - L_ledge_min_current, m_B_max_current - m_B_min_current];
    
    % 记录收敛判据误差
    mB_dT_error = solution.mccormick_errors.mB_dT_mean_error;
    Lledge_dm_error = solution.mccormick_errors.Lledge_dm_mean_error;
    Lledge_superheat_error = solution.mccormick_errors.Lledge_superheat_mean_error;
    C_Al2O3_m_error = solution.mccormick_errors.C_Al2O3_m_mean_error;
    C_AlF3_m_error = solution.mccormick_errors.C_AlF3_m_mean_error;
    C_CaF2_m_error = solution.mccormick_errors.C_CaF2_m_mean_error;
    C_LiF_m_error = solution.mccormick_errors.C_LiF_m_mean_error;
    error_history = [error_history; iter, mB_dT_error, Lledge_dm_error, Lledge_superheat_error, ...
        C_Al2O3_m_error, C_AlF3_m_error, C_CaF2_m_error, C_LiF_m_error];
    
    % 记录第21个时间段的边界和变量值（用于分析）
    boundary_history = [boundary_history; iter, dT_T_min_current(21), dT_T_max_current(21), dm_m_min_current(21), dm_m_max_current(21)];
    differential_history = [differential_history; iter, solution.dT_values(21), solution.dm_values(21)];
    ledge_history = [ledge_history; iter, L_ledge_min_current(21), L_ledge_max_current(21), solution.L_ledge(21)];
    mass_history = [mass_history; iter, m_B_min_current(21), m_B_max_current(21), solution.m_B(21)];
    
    % 记录收缩因子历史
    shrink_factor_history = [shrink_factor_history; iter, shrink_factors.dT, shrink_factors.dm, shrink_factors.ledge, shrink_factors.mass, shrink_factor_info.adaptation_reason];
    
    % 输出第21个时间段的详细信息
    fprintf('第21个时间段 - 温度差分: %.4f°C, 边界: [%.4f, %.4f]°C\n', solution.dT_values(21), dT_T_min_current(21), dT_T_max_current(21));
    fprintf('第21个时间段 - 质量差分: %.4f kg, 边界: [%.4f, %.4f] kg\n', solution.dm_values(21), dm_m_min_current(21), dm_m_max_current(21));
    
    % 根据迭代类型处理边界更新
    if iter == 1
        % 首次迭代：定义分段点并从分段McCormick结果中提取每个时段的边界
        dT_T_segments = create_segments_with_zero(dT_T_min_current, dT_T_max_current, N_T_segments);
        dm_m_segments = create_segments_with_zero(dm_m_min_current, dm_m_max_current, N_m_segments);
        
        [dT_T_min_new, dT_T_max_new, dm_m_min_new, dm_m_max_new, L_ledge_min_new, L_ledge_max_new, m_B_min_new, m_B_max_new] = extract_time_varying_bounds_from_piecewise(...
            solution, dT_T_segments, dm_m_segments, N_T_segments, N_m_segments, shrink_factors, params);
    else
        % 后续迭代：基于当前解进行边界收缩
        [dT_T_min_new, dT_T_max_new, dm_m_min_new, dm_m_max_new, L_ledge_min_new, L_ledge_max_new, m_B_min_new, m_B_max_new] = shrink_time_varying_bounds(...
            solution, dT_T_min_current, dT_T_max_current, dm_m_min_current, dm_m_max_current, ...
            L_ledge_min_current, L_ledge_max_current, m_B_min_current, m_B_max_current, shrink_factors, params);
    end
    
    % 计算边界变化
    dT_range_change = mean(abs((dT_T_max_new - dT_T_min_new) - (dT_T_max_current - dT_T_min_current)));
    dm_range_change = mean(abs((dm_m_max_new - dm_m_min_new) - (dm_m_max_current - dm_m_min_current)));
    L_ledge_range_change = mean(abs((L_ledge_max_new - L_ledge_min_new) - (L_ledge_max_current - L_ledge_min_current)));
    m_B_range_change = mean(abs((m_B_max_new - m_B_min_new) - (m_B_max_current - m_B_min_current)));
    
    % 检查收敛性（基于McCormick包络误差）
    fprintf('收敛判据误差: 浴液质量×温度差分 %.4f%%, ledge厚度×质量差分 %.4f%%, ledge厚度×过热度 %.4f%%\n', ...
        solution.mccormick_errors.mB_dT_mean_error, solution.mccormick_errors.Lledge_dm_mean_error, solution.mccormick_errors.Lledge_superheat_mean_error);
    fprintf('收敛判据误差: Al2O3×浴液质量 %.4f%%, AlF3×浴液质量 %.4f%%, CaF2×浴液质量 %.4f%%, LiF×浴液质量 %.4f%%\n', ...
        solution.mccormick_errors.C_Al2O3_m_mean_error, solution.mccormick_errors.C_AlF3_m_mean_error, ...
        solution.mccormick_errors.C_CaF2_m_mean_error, solution.mccormick_errors.C_LiF_m_mean_error);
    fprintf('边界变化: 温度差分 %.4f, 质量差分 %.4f, Ledge厚度 %.4f, 浴液质量 %.4f\n', dT_range_change, dm_range_change, L_ledge_range_change, m_B_range_change);
    
    % 主要收敛判据：所有双线性项的包络误差都小于阈值
    mB_dT_error = solution.mccormick_errors.mB_dT_mean_error;
    Lledge_dm_error = solution.mccormick_errors.Lledge_dm_mean_error;
    Lledge_superheat_error = solution.mccormick_errors.Lledge_superheat_mean_error;
    C_Al2O3_m_error = solution.mccormick_errors.C_Al2O3_m_mean_error;
    C_AlF3_m_error = solution.mccormick_errors.C_AlF3_m_mean_error;
    C_CaF2_m_error = solution.mccormick_errors.C_CaF2_m_mean_error;
    C_LiF_m_error = solution.mccormick_errors.C_LiF_m_mean_error;
    
    all_errors = [mB_dT_error, Lledge_dm_error, Lledge_superheat_error, C_Al2O3_m_error, C_AlF3_m_error, C_CaF2_m_error, C_LiF_m_error];
    max_error = max(all_errors);
    
    if max_error < params.tolerance
        fprintf('算法收敛！所有双线性项包络误差都小于阈值 %.4f%%\n', params.tolerance);
        fprintf('最大误差: %.4f%% (浴液质量×温度差分: %.4f%%, ledge厚度×质量差分: %.4f%%, ledge厚度×过热度: %.4f%%)\n', ...
            max_error, mB_dT_error, Lledge_dm_error, Lledge_superheat_error);
        fprintf('浓度×浴液质量误差: Al2O3 %.4f%%, AlF3 %.4f%%, CaF2 %.4f%%, LiF %.4f%%\n', ...
            C_Al2O3_m_error, C_AlF3_m_error, C_CaF2_m_error, C_LiF_m_error);
        break;
    end
    
    % 更新时变边界
    dT_T_min_current = dT_T_min_new;
    dT_T_max_current = dT_T_max_new;
    dm_m_min_current = dm_m_min_new;
    dm_m_max_current = dm_m_max_new;
    L_ledge_min_current = L_ledge_min_new;
    L_ledge_max_current = L_ledge_max_new;
    m_B_min_current = m_B_min_new;
    m_B_max_current = m_B_max_new;
    
    % 输出更新后的第21个时间段边界信息
    fprintf('更新后第21个时间段边界 - 温度差分: [%.4f, %.4f]°C, 质量差分: [%.4f, %.4f] kg\n', ...
        dT_T_min_current(21), dT_T_max_current(21), dm_m_min_current(21), dm_m_max_current(21));
    fprintf('更新后第21个时间段边界 - Ledge厚度: [%.4f, %.4f] m, 浴液质量: [%.0f, %.0f] kg\n', ...
        L_ledge_min_current(21), L_ledge_max_current(21), m_B_min_current(21), m_B_max_current(21));
    
end

%% 输出最终结果
fprintf('\n=== 边界收缩算法完成 ===\n');
fprintf('总迭代次数: %d\n', iter);
fprintf('最终目标函数值: %.2f\n', obj_value);
fprintf('最终温度差分范围: [%.2f, %.2f]°C (平均范围: %.2f)\n', min(dT_T_min_current), max(dT_T_max_current), mean(dT_T_max_current - dT_T_min_current));
fprintf('最终质量差分范围: [%.0f, %.0f] kg (平均范围: %.0f)\n', min(dm_m_min_current), max(dm_m_max_current), mean(dm_m_max_current - dm_m_min_current));
fprintf('最终Ledge厚度范围: [%.3f, %.3f] m (平均范围: %.3f)\n', min(L_ledge_min_current), max(L_ledge_max_current), mean(L_ledge_max_current - L_ledge_min_current));
fprintf('最终浴液质量范围: [%.0f, %.0f] kg (平均范围: %.0f)\n', min(m_B_min_current), max(m_B_max_current), mean(m_B_max_current - m_B_min_current));

% 输出每个时间段的详细边界信息
fprintf('\n=== 各时间段详细边界信息 ===\n');
fprintf('时间段\t温度差分范围(°C)\t\t质量差分范围(kg)\t\tLedge厚度范围(m)\t\t浴液质量范围(kg)\n');
fprintf('------\t----------------\t\t----------------\t\t----------------\t\t----------------\n');
for t = 1:params.T
    fprintf('%2d\t[%.3f, %.3f]\t\t[%.1f, %.1f]\t\t[%.4f, %.4f]\t\t[%.0f, %.0f]\n', t, ...
        dT_T_min_current(t), dT_T_max_current(t), dm_m_min_current(t), dm_m_max_current(t), ...
        L_ledge_min_current(t), L_ledge_max_current(t), m_B_min_current(t), m_B_max_current(t));
end

% 显示最终McCormick包络误差
if ~isempty(solution) && isfield(solution, 'mccormick_errors')
    fprintf('\n=== 最终McCormick包络误差 ===\n');
    
    % 分段McCormick包络误差
    fprintf('分段McCormick包络平均误差: \n');
    fprintf('  - 温度差分×浴液质量平均误差: %.4f%%\n', solution.mccormick_errors.mB_dT_mean_error);
    fprintf('  - Ledge厚度×质量差分平均误差: %.4f%%\n', solution.mccormick_errors.Lledge_dm_mean_error);
    
    % 浓度与浴液质量相乘的误差
    fprintf('浓度×浴液质量平均误差: \n');
    fprintf('  - Al2O3浓度×浴液质量平均误差: %.4f%%\n', solution.mccormick_errors.C_Al2O3_m_mean_error);
    fprintf('  - AlF3浓度×浴液质量平均误差: %.4f%%\n', solution.mccormick_errors.C_AlF3_m_mean_error);
    fprintf('  - CaF2浓度×浴液质量平均误差: %.4f%%\n', solution.mccormick_errors.C_CaF2_m_mean_error);
    fprintf('  - LiF浓度×浴液质量平均误差: %.4f%%\n', solution.mccormick_errors.C_LiF_m_mean_error);
    
    % 其他全局McCormick包络误差
    fprintf('电流×电流效率平均误差: %.4f%%\n', solution.mccormick_errors.I_g_mean_error);
    fprintf('Ledge厚度×过热度平均误差: %.4f%%\n', solution.mccormick_errors.Lledge_superheat_mean_error);
    
end

% 绘制收敛历史
if ~isempty(iteration_history)
    figure('Name', '边界收缩算法收敛历史', 'Position', [100, 100, 1800, 900]);
    
    subplot(3,3,1);
    plot(iteration_history(:,1), iteration_history(:,2), 'b-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('目标函数值');
    title('目标函数收敛历史');
    grid on;
    
    subplot(3,3,2);
    plot(iteration_history(:,1), iteration_history(:,3), 'r-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('温度差分范围');
    title('温度差分范围收缩历史');
    grid on;
    
    subplot(3,3,3);
    plot(iteration_history(:,1), iteration_history(:,4), 'g-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('质量差分范围');
    title('质量差分范围收缩历史');
    grid on;
    
    % 收敛判据误差历史
    subplot(3,3,4);
    plot(error_history(:,1), error_history(:,2), 'm-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('浴液质量×温度差分误差 (%)');
    title('浴液质量×温度差分误差收敛历史');
    grid on;
    yline(params.tolerance, 'r--', 'LineWidth', 2, 'DisplayName', '收敛阈值');
    legend('误差', '阈值', 'Location', 'best');
    
    subplot(3,3,5);
    plot(error_history(:,1), error_history(:,3), 'c-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('Ledge厚度×质量差分误差 (%)');
    title('Ledge厚度×质量差分误差收敛历史');
    grid on;
    yline(params.tolerance, 'r--', 'LineWidth', 2, 'DisplayName', '收敛阈值');
    legend('误差', '阈值', 'Location', 'best');
    
    % 新增：Ledge厚度×过热度误差
    subplot(3,3,6);
    plot(error_history(:,1), error_history(:,4), 'k-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('Ledge厚度×过热度误差 (%)');
    title('Ledge厚度×过热度误差收敛历史');
    grid on;
    yline(params.tolerance, 'r--', 'LineWidth', 2, 'DisplayName', '收敛阈值');
    legend('误差', '阈值', 'Location', 'best');
    
    % 新增：浴液质量×组分浓度误差（Al2O3）
    subplot(3,3,7);
    plot(error_history(:,1), error_history(:,5), 'y-o', 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('Al2O3浓度×浴液质量误差 (%)');
    title('Al2O3浓度×浴液质量误差收敛历史');
    grid on;
    yline(params.tolerance, 'r--', 'LineWidth', 2, 'DisplayName', '收敛阈值');
    legend('误差', '阈值', 'Location', 'best');
    
    % 综合收敛判据图（包含所有主要误差）
    subplot(3,3,8);
    plot(error_history(:,1), error_history(:,2), 'm-o', 'LineWidth', 2, 'DisplayName', '浴液质量×温度差分');
    hold on;
    plot(error_history(:,1), error_history(:,3), 'c-o', 'LineWidth', 2, 'DisplayName', 'Ledge厚度×质量差分');
    plot(error_history(:,1), error_history(:,4), 'k-o', 'LineWidth', 2, 'DisplayName', 'Ledge厚度×过热度');
    plot(error_history(:,1), error_history(:,5), 'y-o', 'LineWidth', 2, 'DisplayName', 'Al2O3浓度×浴液质量');
    xlabel('迭代次数');
    ylabel('McCormick包络误差 (%)');
    title('主要收敛判据误差综合历史');
    grid on;
    yline(params.tolerance, 'r--', 'LineWidth', 2, 'DisplayName', '收敛阈值');
    legend('Location', 'best');
    hold off;
    
    % 所有组分浓度×浴液质量误差综合图
    subplot(3,3,9);
    plot(error_history(:,1), error_history(:,5), 'y-o', 'LineWidth', 2, 'DisplayName', 'Al2O3');
    hold on;
    plot(error_history(:,1), error_history(:,6), 'r-o', 'LineWidth', 2, 'DisplayName', 'AlF3');
    plot(error_history(:,1), error_history(:,7), 'g-o', 'LineWidth', 2, 'DisplayName', 'CaF2');
    plot(error_history(:,1), error_history(:,8), 'b-o', 'LineWidth', 2, 'DisplayName', 'LiF');
    xlabel('迭代次数');
    ylabel('组分浓度×浴液质量误差 (%)');
    title('各组分浓度×浴液质量误差收敛历史');
    grid on;
    yline(params.tolerance, 'r--', 'LineWidth', 2, 'DisplayName', '收敛阈值');
    legend('Location', 'best');
    hold off;
end

% 绘制收缩因子历史
if ~isempty(shrink_factor_history)
    figure('Name', '自适应收缩因子历史', 'Position', [300, 300, 1200, 600]);
    
    subplot(2,2,1);
    plot(shrink_factor_history(:,1), shrink_factor_history(:,2), 'r-o', 'LineWidth', 2, 'DisplayName', '温度差分');
    hold on;
    plot(shrink_factor_history(:,1), shrink_factor_history(:,3), 'g-o', 'LineWidth', 2, 'DisplayName', '质量差分');
    xlabel('迭代次数');
    ylabel('收缩因子');
    title('差分变量收缩因子历史');
    legend('Location', 'best');
    grid on;
    hold off;
    
    subplot(2,2,2);
    plot(shrink_factor_history(:,1), shrink_factor_history(:,4), 'b-o', 'LineWidth', 2, 'DisplayName', 'Ledge厚度');
    hold on;
    plot(shrink_factor_history(:,1), shrink_factor_history(:,5), 'm-o', 'LineWidth', 2, 'DisplayName', '浴液质量');
    xlabel('迭代次数');
    ylabel('收缩因子');
    title('状态变量收缩因子历史');
    legend('Location', 'best');
    grid on;
    hold off;
    
    subplot(2,2,3);
    plot(shrink_factor_history(:,1), shrink_factor_history(:,2), 'r-o', 'LineWidth', 2, 'DisplayName', '温度差分');
    hold on;
    plot(shrink_factor_history(:,1), shrink_factor_history(:,3), 'g-o', 'LineWidth', 2, 'DisplayName', '质量差分');
    plot(shrink_factor_history(:,1), shrink_factor_history(:,4), 'b-o', 'LineWidth', 2, 'DisplayName', 'Ledge厚度');
    plot(shrink_factor_history(:,1), shrink_factor_history(:,5), 'm-o', 'LineWidth', 2, 'DisplayName', '浴液质量');
    xlabel('迭代次数');
    ylabel('收缩因子');
    title('所有收缩因子综合历史');
    legend('Location', 'best');
    grid on;
    hold off;
    
    subplot(2,2,4);
    % 显示收缩因子变化原因
    adaptation_reasons = shrink_factor_history(:,6);
    unique_reasons = unique(adaptation_reasons);
    reason_counts = zeros(length(unique_reasons), 1);
    for i = 1:length(unique_reasons)
        reason_counts(i) = sum(adaptation_reasons == unique_reasons(i));
    end
    pie(reason_counts, unique_reasons);
    title('收缩因子自适应原因分布');
end

% 绘制第21个时间段的边界和差分值收敛历史
if ~isempty(boundary_history) && ~isempty(differential_history)
    figure('Name', '第21个时间段边界和差分值收敛历史', 'Position', [200, 200, 1200, 800]);
    
    % 温度差分边界和值的历史
    subplot(2,2,1);
    plot(boundary_history(:,1), boundary_history(:,2), 'r-', 'LineWidth', 2, 'DisplayName', '温度差分下界');
    hold on;
    plot(boundary_history(:,1), boundary_history(:,3), 'r--', 'LineWidth', 2, 'DisplayName', '温度差分上界');
    plot(differential_history(:,1), differential_history(:,2), 'bo-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '温度差分解');
    xlabel('迭代次数');
    ylabel('温度差分 (°C)');
    title('第21个时间段温度差分边界和解的收敛历史');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % 质量差分边界和值的历史
    subplot(2,2,2);
    plot(boundary_history(:,1), boundary_history(:,4), 'g-', 'LineWidth', 2, 'DisplayName', '质量差分下界');
    hold on;
    plot(boundary_history(:,1), boundary_history(:,5), 'g--', 'LineWidth', 2, 'DisplayName', '质量差分上界');
    plot(differential_history(:,1), differential_history(:,3), 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '质量差分解');
    xlabel('迭代次数');
    ylabel('质量差分 (kg)');
    title('第21个时间段质量差分边界和解的收敛历史');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % 温度差分边界范围的历史
    subplot(2,2,3);
    dT_range = boundary_history(:,3) - boundary_history(:,2);
    plot(boundary_history(:,1), dT_range, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('迭代次数');
    ylabel('温度差分边界范围 (°C)');
    title('第21个时间段温度差分边界范围收缩历史');
    grid on;
    
    % 质量差分边界范围的历史
    subplot(2,2,4);
    dm_range = boundary_history(:,5) - boundary_history(:,4);
    plot(boundary_history(:,1), dm_range, 'm-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('迭代次数');
    ylabel('质量差分边界范围 (kg)');
    title('第21个时间段质量差分边界范围收缩历史');
    grid on;
end

% 保存结果
save('boundary_tightening_results.mat', 'solution', 'iteration_history', 'error_history', ...
     'boundary_history', 'differential_history', 'ledge_history', 'mass_history', ...
     'shrink_factor_history', 'dT_T_min_current', 'dT_T_max_current', 'dm_m_min_current', 'dm_m_max_current', ...
     'L_ledge_min_current', 'L_ledge_max_current', 'm_B_min_current', 'm_B_max_current', ...
     'obj_value', 'iter');

fprintf('结果已保存到 boundary_tightening_results.mat\n');

% 可视化最后一次求解结果
if ~isempty(solution) && feasible
    visualize_optimization_results(solution, params);
end

%% 辅助函数：求解分段McCormick包络优化问题
function [obj_value, solution, feasible] = solve_piecewise_mccormick_optimization(...
    dT_T_min, dT_T_max, dm_m_min, dm_m_max, L_ledge_min, L_ledge_max, m_B_min, m_B_max, ...
    N_T_segments, N_m_segments, params)
    
    % 定义分段点（包含0）
    dT_T_segments = create_segments_with_zero(dT_T_min, dT_T_max, N_T_segments);
    dm_m_segments = create_segments_with_zero(dm_m_min, dm_m_max, N_m_segments);
    
    % 构建优化问题
    [constraints, variables, obj] = build_optimization_problem(...
        dT_T_min, dT_T_max, dm_m_min, dm_m_max, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params, ...
        'piecewise', N_T_segments, N_m_segments, dT_T_segments, dm_m_segments);
    
    % 求解优化问题
    [obj_value, solution, feasible] = solve_optimization_problem(constraints, obj, variables, params);
    
    % 添加分段相关的结果
    if feasible
        solution.lambda_T = value(variables.lambda_T);
        solution.lambda_m = value(variables.lambda_m);
    end
end

%% 辅助函数：创建包含0的分段点
function segments = create_segments_with_zero(min_val, max_val, N_segments)
    % 确保输入为标量
    min_val = min(min_val);
    max_val = max(max_val);
    
    if min_val < 0 && max_val > 0
        % 包含0的情况
        negative_segments = linspace(min_val, 0, ceil(N_segments/2) + 1);
        positive_segments = linspace(0, max_val, floor(N_segments/2) + 1);
        segments = unique([negative_segments, positive_segments]);
    else
        % 不包含0的情况
        segments = linspace(min_val, max_val, N_segments + 1);
    end
end

%% 辅助函数：添加全局McCormick包络约束
function constraints = add_global_mccormick_constraints(constraints, I, g, C_Al2O3, C_AlF3, C_CaF2, C_LiF, m_B, L_ledge, T_B, T_liq, ...
    Z_I_g, Z_C_Al2O3_m, Z_C_AlF3_m, Z_C_CaF2_m, Z_C_LiF_m, Z_Lledge_superheat, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params)
    
    % 1. 电流与电流效率相乘的McCormick包络
    constraints = [constraints, Z_I_g >= params.I_min * g + params.g_min * I - params.I_min * params.g_min];
    constraints = [constraints, Z_I_g >= params.I_max * g + params.g_max * I - params.I_max * params.g_max];
    constraints = [constraints, Z_I_g <= params.I_min * g + params.g_max * I - params.I_min * params.g_max];
    constraints = [constraints, Z_I_g <= params.I_max * g + params.g_min * I - params.I_max * params.g_min];
    
    % 2. 浓度与浴液质量相乘的McCormick包络

    for t = 1:params.T+1
        % Al2O3浓度
        constraints = [constraints, Z_C_Al2O3_m(t) >= params.C_Al2O3_min * m_B(t) + m_B_min(t) * C_Al2O3(t) - params.C_Al2O3_min * m_B_min(t)];
        constraints = [constraints, Z_C_Al2O3_m(t) >= params.C_Al2O3_max * m_B(t) + m_B_max(t) * C_Al2O3(t) - params.C_Al2O3_max * m_B_max(t)];
        constraints = [constraints, Z_C_Al2O3_m(t) <= params.C_Al2O3_min * m_B(t) + m_B_max(t) * C_Al2O3(t) - params.C_Al2O3_min * m_B_max(t)];
        constraints = [constraints, Z_C_Al2O3_m(t) <= params.C_Al2O3_max * m_B(t) + m_B_min(t) * C_Al2O3(t) - params.C_Al2O3_max * m_B_min(t)];
        
        % AlF3浓度
        constraints = [constraints, Z_C_AlF3_m(t) >= params.C_AlF3_min * m_B(t) + m_B_min(t) * C_AlF3(t) - params.C_AlF3_min * m_B_min(t)];
        constraints = [constraints, Z_C_AlF3_m(t) >= params.C_AlF3_max * m_B(t) + m_B_max(t) * C_AlF3(t) - params.C_AlF3_max * m_B_max(t)];
        constraints = [constraints, Z_C_AlF3_m(t) <= params.C_AlF3_min * m_B(t) + m_B_max(t) * C_AlF3(t) - params.C_AlF3_min * m_B_max(t)];
        constraints = [constraints, Z_C_AlF3_m(t) <= params.C_AlF3_max * m_B(t) + m_B_min(t) * C_AlF3(t) - params.C_AlF3_max * m_B_min(t)];
        
        % CaF2浓度
        constraints = [constraints, Z_C_CaF2_m(t) >= params.C_CaF2_min * m_B(t) + m_B_min(t) * C_CaF2(t) - params.C_CaF2_min * m_B_min(t)];
        constraints = [constraints, Z_C_CaF2_m(t) >= params.C_CaF2_max * m_B(t) + m_B_max(t) * C_CaF2(t) - params.C_CaF2_max * m_B_max(t)];
        constraints = [constraints, Z_C_CaF2_m(t) <= params.C_CaF2_min * m_B(t) + m_B_max(t) * C_CaF2(t) - params.C_CaF2_min * m_B_max(t)];
        constraints = [constraints, Z_C_CaF2_m(t) <= params.C_CaF2_max * m_B(t) + m_B_min(t) * C_CaF2(t) - params.C_CaF2_max * m_B_min(t)];
        
        % LiF浓度
        constraints = [constraints, Z_C_LiF_m(t) >= params.C_LiF_min * m_B(t) + m_B_min(t) * C_LiF(t) - params.C_LiF_min * m_B_min(t)];
        constraints = [constraints, Z_C_LiF_m(t) >= params.C_LiF_max * m_B(t) + m_B_max(t) * C_LiF(t) - params.C_LiF_max * m_B_max(t)];
        constraints = [constraints, Z_C_LiF_m(t) <= params.C_LiF_min * m_B(t) + m_B_max(t) * C_LiF(t) - params.C_LiF_min * m_B_max(t)];
        constraints = [constraints, Z_C_LiF_m(t) <= params.C_LiF_max * m_B(t) + m_B_min(t) * C_LiF(t) - params.C_LiF_max * m_B_min(t)];
    end
    
    % 3. Ledge厚度与过热度相乘的McCormick包络
    for t = 1:params.T
        superheat_t = T_B(t) - T_liq(t);  % 过热度 = T_B(t) - T_liq(t)
        constraints = [constraints, Z_Lledge_superheat(t) >= L_ledge_min(t) * superheat_t + params.Superheat_min * L_ledge(t) - L_ledge_min(t) * params.Superheat_min];
        constraints = [constraints, Z_Lledge_superheat(t) >= L_ledge_max(t) * superheat_t + params.Superheat_max * L_ledge(t) - L_ledge_max(t) * params.Superheat_max];
        constraints = [constraints, Z_Lledge_superheat(t) <= L_ledge_min(t) * superheat_t + params.Superheat_max * L_ledge(t) - L_ledge_min(t) * params.Superheat_max];
        constraints = [constraints, Z_Lledge_superheat(t) <= L_ledge_max(t) * superheat_t + params.Superheat_min * L_ledge(t) - L_ledge_max(t) * params.Superheat_min];
    end
end

%% 辅助函数：添加分段McCormick包络约束
function constraints = add_piecewise_mccormick_constraints(constraints, T_B, m_B, L_ledge, ...
    Z_mB_dT, Z_Lledge_dm, lambda_T, lambda_m, ...
    dT_T_segments, dm_m_segments, params, N_T_segments, N_m_segments)
    
    % 定义变量范围
    m_B_min = params.m_B_min_init; m_B_max = params.m_B_max_init;

    % 定义子变量
    dT_sub = sdpvar(params.T, N_T_segments, 'full');
    m_B_sub = sdpvar(params.T, N_T_segments, 'full');
    dm_sub = sdpvar(params.T, N_m_segments, 'full');
    L_ledge_sub = sdpvar(params.T, N_m_segments, 'full');
    
    % 分段McCormick包络
    for t = 1:params.T
        % 变量分解约束
        constraints = [constraints, T_B(t+1) - T_B(t) == sum(dT_sub(t,:))];
        constraints = [constraints, m_B(t+1) - m_B(t) == sum(dm_sub(t,:))];
        constraints = [constraints, m_B(t) == sum(m_B_sub(t,:))];
        constraints = [constraints, L_ledge(t) == sum(L_ledge_sub(t,:))];

        % 分段选择约束
        constraints = [constraints, sum(lambda_T(t, :)) == 1];
        constraints = [constraints, sum(lambda_m(t, :)) == 1];

        for s = 1:N_T_segments
            dT_s_min = dT_T_segments(s);
            dT_s_max = dT_T_segments(s+1);
            
            % 分段区间约束
            constraints = [constraints, dT_s_min * lambda_T(t,s) <= dT_sub(t,s) <= dT_s_max * lambda_T(t,s)];
            constraints = [constraints, m_B_min * lambda_T(t,s) <= m_B_sub(t,s) <= m_B_max * lambda_T(t,s)];

            % 分段内的McCormick包络约束
            est_T_1(s) = m_B_min * dT_sub(t,s) + dT_s_min * m_B_sub(t,s) - m_B_min * dT_s_min * lambda_T(t,s);
            est_T_2(s) = m_B_max * dT_sub(t,s) + dT_s_max * m_B_sub(t,s) - m_B_max * dT_s_max * lambda_T(t,s);
            est_T_3(s) = m_B_min * dT_sub(t,s) + dT_s_max * m_B_sub(t,s) - m_B_min * dT_s_max * lambda_T(t,s);
            est_T_4(s) = m_B_max * dT_sub(t,s) + dT_s_min * m_B_sub(t,s) - m_B_max * dT_s_min * lambda_T(t,s);
        end
        constraints = [constraints, Z_mB_dT(t) >= sum(est_T_1)];
        constraints = [constraints, Z_mB_dT(t) >= sum(est_T_2)];
        constraints = [constraints, Z_mB_dT(t) <= sum(est_T_3)];
        constraints = [constraints, Z_mB_dT(t) <= sum(est_T_4)];

        for s = 1:N_m_segments
            dm_s_min = dm_m_segments(s);
            dm_s_max = dm_m_segments(s+1);
            
            % 分段区间约束
            constraints = [constraints, dm_s_min * lambda_m(t,s) <= dm_sub(t,s) <= dm_s_max * lambda_m(t,s)];
            constraints = [constraints, params.L_ledge_min * lambda_m(t,s) <= L_ledge_sub(t,s) <= params.L_ledge_max * lambda_m(t,s)];
           
            % 分段内的McCormick包络约束
            est_m_1(s) = params.L_ledge_min * dm_sub(t,s) + dm_s_min * L_ledge_sub(t,s) - params.L_ledge_min * dm_s_min * lambda_m(t,s);
            est_m_2(s) = params.L_ledge_max * dm_sub(t,s) + dm_s_max * L_ledge_sub(t,s) - params.L_ledge_max * dm_s_max * lambda_m(t,s);
            est_m_3(s) = params.L_ledge_min * dm_sub(t,s) + dm_s_max * L_ledge_sub(t,s) - params.L_ledge_min * dm_s_max * lambda_m(t,s);
            est_m_4(s) = params.L_ledge_max * dm_sub(t,s) + dm_s_min * L_ledge_sub(t,s) - params.L_ledge_max * dm_s_min * lambda_m(t,s);
        end
        constraints = [constraints, Z_Lledge_dm(t) >= sum(est_m_1)];
        constraints = [constraints, Z_Lledge_dm(t) >= sum(est_m_2)];
        constraints = [constraints, Z_Lledge_dm(t) <= sum(est_m_3)];
        constraints = [constraints, Z_Lledge_dm(t) <= sum(est_m_4)];
    end

end

%% 辅助函数：添加动态约束
function constraints = add_dynamic_constraints(constraints, T_B, L_ledge, m_B, C_Al2O3, C_AlF3, C_CaF2, C_LiF, ...
    T_liq, Superheat, g, I, m_feed_Al2O3, Z_I_g, Z_C_Al2O3_m, Z_C_AlF3_m, Z_C_CaF2_m, Z_C_LiF_m, ...
    Z_mB_dT, Z_Lledge_dm, Z_Lledge_superheat, params)

    % 液相温度计算
    constraints = [constraints, ...
        T_liq == params.beta_linear(1) * C_Al2O3(1:params.T)*100 + ...
                 params.beta_linear(2) * C_AlF3(1:params.T)*100 + ...
                 params.beta_linear(3) * C_CaF2(1:params.T)*100 + ...
                 params.beta_linear(4)];
    
    % 过热度计算
    constraints = [constraints, Superheat == T_B(1:params.T) - T_liq];
    
    % 电流效率计算
    constraints = [constraints, g == params.coeff_global(1) + ...
                                params.coeff_global(2)*Superheat + ...
                                params.coeff_global(3)*C_AlF3(1:params.T)*100 + ...
                                params.coeff_global(4)*C_LiF(1:params.T)*100];
    
    % 功率计算
    P_N = (params.EMF + params.R_cell*params.I_N*1000)*params.I_N;
    dP_dI_N = params.EMF + 2*params.R_cell*params.I_N*1000;
    P = P_N + dP_dI_N*(I - repmat(params.I_N,1,params.T));
    
    % 浴液温度更新约束
    constraints = [constraints, ...
        Z_mB_dT * params.c_B + params.m_M * params.c_M * (T_B(2:params.T+1) - T_B(1:params.T)) == ...
        3600*params.delta_t * (P(1:params.T)*1000 - ...
        (T_B(1:params.T) - params.T_gas)/params.R_B_gas - ...
        (T_B(1:params.T) - T_liq)/params.R_B_ledge - ...
        (T_B(1:params.T) - params.T_amb)/params.R_bottom_total - ...
        0.596*60*1e3*m_feed_Al2O3(1:params.T) - ...
        5.474*3600*1e3*(Z_I_g*1000*params.M_Al)/(params.F*params.z))];
    
    % 浴液质量更新约束
    R_shell_side = 0.25/(47 * params.A_ledge);
    R_shell_side_amb = 1/(22 * params.A_ledge);
    coeff_m_B = params.delta_t*3600 / params.lf;
    constraints = [constraints, ...
        Z_Lledge_dm / (1.5 * params.A_ledge) + (R_shell_side + R_shell_side_amb) * (m_B(2:params.T+1) - m_B(1:params.T)) == ...
        coeff_m_B * (Z_Lledge_superheat / (1.5 * params.A_ledge * params.R_B_ledge) + ...
        (R_shell_side + R_shell_side_amb) / params.R_B_ledge * (T_B(1:params.T) - T_liq) - (T_liq - params.T_amb))];
    
    % ledge厚度更新约束
    coeff_L_ledge = 1 / (params.rho_ledge * params.A_ledge);
    constraints = [constraints, ...
        L_ledge(2:params.T+1) == L_ledge(1:params.T) - ...
        coeff_L_ledge * (m_B(2:params.T+1) - m_B(1:params.T))];
    
    % 氧化铝浓度更新约束
    coeff_m_Al2O3_feed = 60*params.delta_t;
    coeff_m_Al2O3_reaction = -3600*params.delta_t * params.M_Al2O3 / (params.F*2*params.z);
    constraints = [constraints, ...
        Z_C_Al2O3_m(2:params.T+1) == Z_C_Al2O3_m(1:params.T) + ...
        coeff_m_Al2O3_feed * m_feed_Al2O3(1:params.T) + ...
        coeff_m_Al2O3_reaction * (Z_I_g*1000)];
    
    % 其他组分浓度更新约束
    constraints = [constraints, Z_C_AlF3_m(2:params.T+1) == params.C_AlF3_init*params.m_B_init];
    constraints = [constraints, Z_C_CaF2_m(2:params.T+1) == params.C_CaF2_init*params.m_B_init];
    constraints = [constraints, Z_C_LiF_m(2:params.T+1) == params.C_LiF_init*params.m_B_init];
end

%% 辅助函数：添加调节约束
function constraints = add_regulation_constraints(constraints, I, lamda_keep, lamda_up, lamda_down, params)
    % 调节状态唯一性约束
    constraints = [constraints, lamda_keep + lamda_up + lamda_down == ones(1, params.T-1)];
    
    % 功率爬坡约束
    P_N = (params.EMF + params.R_cell*params.I_N*1000)*params.I_N;
    P = P_N + (params.EMF + 2*params.R_cell*params.I_N*1000)*(I - repmat(params.I_N,1,params.T));
    
    constraints = [constraints, P(2:end) - P(1:end-1) >= -params.R_down * lamda_down * P_N];
    constraints = [constraints, P(2:end) - P(1:end-1) <= params.R_up * lamda_up * P_N];
    
    % 最大调节次数约束
    constraints = [constraints, sum(lamda_up + lamda_down) <= params.N_rel];
    
    % 连续调节约束
    if params.T > 2
        constraints = [constraints, lamda_up(2:end) + lamda_down(2:end) + ...
            lamda_up(1:end-1) + lamda_down(1:end-1) <= ones(1, params.T-2)];
    end
end

%% 辅助函数：计算McCormick包络误差
function mccormick_errors = calculate_mccormick_errors(I, g, T_B, m_B, L_ledge, T_liq, ...
    Z_I_g, Z_C_Al2O3_m, Z_C_AlF3_m, Z_C_CaF2_m, Z_C_LiF_m, Z_mB_dT, Z_Lledge_dm, Z_Lledge_superheat, variables, params)
    
    % 提取优化结果
    I_opt = value(I);
    g_opt = value(g);
    T_B_opt = value(T_B);
    m_B_opt = value(m_B);
    L_ledge_opt = value(L_ledge);
    T_liq_opt = value(T_liq);
    Z_I_g_opt = value(Z_I_g);
    Z_C_Al2O3_m_opt = value(Z_C_Al2O3_m);
    Z_C_AlF3_m_opt = value(Z_C_AlF3_m);
    Z_C_CaF2_m_opt = value(Z_C_CaF2_m);
    Z_C_LiF_m_opt = value(Z_C_LiF_m);
    Z_mB_dT_opt = value(Z_mB_dT);
    Z_Lledge_dm_opt = value(Z_Lledge_dm);
    Z_Lledge_superheat_opt = value(Z_Lledge_superheat);
    
    % 计算真实乘积值
    I_g_product = I_opt .* g_opt;
    
    % 各组分浓度与浴液质量相乘的真实值
    % 注意：需要从优化结果中提取实际的浓度值，而不是使用固定初始值
    C_Al2O3_opt = value(variables.C_Al2O3);
    C_AlF3_opt = value(variables.C_AlF3);
    C_CaF2_opt = value(variables.C_CaF2);
    C_LiF_opt = value(variables.C_LiF);
    
    C_Al2O3_m_product = C_Al2O3_opt .* m_B_opt;
    C_AlF3_m_product = C_AlF3_opt .* m_B_opt;
    C_CaF2_m_product = C_CaF2_opt .* m_B_opt;
    C_LiF_m_product = C_LiF_opt .* m_B_opt;
    
    mB_dT_product = m_B_opt(1:end-1) .* (T_B_opt(2:end) - T_B_opt(1:end-1));
    Lledge_dm_product = L_ledge_opt(1:end-1) .* (m_B_opt(2:end) - m_B_opt(1:end-1));
    Lledge_superheat_product = L_ledge_opt(1:end-1) .* (T_B_opt(1:end-1) - T_liq_opt);
    
    % 计算相对误差
    eps_val = 1e-10;  % 避免除零
    
    % 电流与电流效率相乘的误差
    I_g_error = abs(Z_I_g_opt - I_g_product) ./ (abs(I_g_product) + eps_val) * 100;
    
    % 各组分浓度与浴液质量相乘的误差
    C_Al2O3_m_error = abs(Z_C_Al2O3_m_opt - C_Al2O3_m_product) ./ (abs(C_Al2O3_m_product) + eps_val) * 100;
    C_AlF3_m_error = abs(value(Z_C_AlF3_m_opt) - C_AlF3_m_product) ./ (abs(C_AlF3_m_product) + eps_val) * 100;
    C_CaF2_m_error = abs(value(Z_C_CaF2_m_opt) - C_CaF2_m_product) ./ (abs(C_CaF2_m_product) + eps_val) * 100;
    C_LiF_m_error = abs(value(Z_C_LiF_m_opt) - C_LiF_m_product) ./ (abs(C_LiF_m_product) + eps_val) * 100;
    
    % 温度差分与浴液质量相乘的误差（分段McCormick）
    mB_dT_error = abs(Z_mB_dT_opt - mB_dT_product) ./ (abs(mB_dT_product) + eps_val) * 100;
    
    % Ledge厚度与质量差分相乘的误差（分段McCormick）
    Lledge_dm_error = abs(Z_Lledge_dm_opt - Lledge_dm_product) ./ (abs(Lledge_dm_product) + eps_val) * 100;
    
    % Ledge厚度与过热度相乘的误差
    Lledge_superheat_error = abs(Z_Lledge_superheat_opt - Lledge_superheat_product) ./ (abs(Lledge_superheat_product) + eps_val) * 100;
    
    % 计算平均相对误差
    mccormick_errors = struct();
    
    % 记录每个时段的误差（用于可视化）
    mccormick_errors.I_g_error = I_g_error;
    mccormick_errors.C_Al2O3_m_error = C_Al2O3_m_error;
    mccormick_errors.C_AlF3_m_error = C_AlF3_m_error;
    mccormick_errors.C_CaF2_m_error = C_CaF2_m_error;
    mccormick_errors.C_LiF_m_error = C_LiF_m_error;
    mccormick_errors.mB_dT_error = mB_dT_error;
    mccormick_errors.Lledge_dm_error = Lledge_dm_error;
    mccormick_errors.Lledge_superheat_error = Lledge_superheat_error;
    
    % 电流与电流效率相乘的误差
    mccormick_errors.I_g_mean_error = mean(I_g_error);
    mccormick_errors.I_g_max_error = max(I_g_error);
    
    % 各组分浓度与浴液质量相乘的误差
    mccormick_errors.C_Al2O3_m_mean_error = mean(C_Al2O3_m_error);
    mccormick_errors.C_Al2O3_m_max_error = max(C_Al2O3_m_error);
    mccormick_errors.C_AlF3_m_mean_error = mean(C_AlF3_m_error);
    mccormick_errors.C_AlF3_m_max_error = max(C_AlF3_m_error);
    mccormick_errors.C_CaF2_m_mean_error = mean(C_CaF2_m_error);
    mccormick_errors.C_CaF2_m_max_error = max(C_CaF2_m_error);
    mccormick_errors.C_LiF_m_mean_error = mean(C_LiF_m_error);
    mccormick_errors.C_LiF_m_max_error = max(C_LiF_m_error);
    
    % 分段McCormick包络的误差
    mccormick_errors.mB_dT_mean_error = mean(mB_dT_error);
    mccormick_errors.mB_dT_max_error = max(mB_dT_error);
    mccormick_errors.Lledge_dm_mean_error = mean(Lledge_dm_error);
    mccormick_errors.Lledge_dm_max_error = max(Lledge_dm_error);
    
    % Ledge厚度与过热度相乘的误差
    mccormick_errors.Lledge_superheat_mean_error = mean(Lledge_superheat_error);
    mccormick_errors.Lledge_superheat_max_error = max(Lledge_superheat_error);
    
    end

%% 辅助函数：求解标准McCormick松弛优化问题
function [obj_value, solution, feasible] = solve_standard_mccormick_optimization(...
    dT_T_min, dT_T_max, dm_m_min, dm_m_max, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params)
    
    % 构建优化问题
    [constraints, variables, obj] = build_optimization_problem(...
        dT_T_min, dT_T_max, dm_m_min, dm_m_max, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params, ...
        'standard');
    
    % 求解优化问题
    [obj_value, solution, feasible] = solve_optimization_problem(constraints, obj, variables, params);
end

%% 辅助函数：添加标准McCormick包络约束
function constraints = add_standard_mccormick_constraints(constraints, T_B, m_B, L_ledge, ...
    Z_mB_dT, Z_Lledge_dm, dT_T_min, dT_T_max, dm_m_min, dm_m_max, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params)
    
    % 1. 温度差分与浴液质量相乘的标准McCormick包络
    for t = 1:params.T
        dT_t = T_B(t+1) - T_B(t);
        
        % 标准McCormick包络约束（使用时变边界）
        constraints = [constraints, Z_mB_dT(t) >= m_B_min(t) * dT_t + dT_T_min(t) * m_B(t) - m_B_min(t) * dT_T_min(t)];
        constraints = [constraints, Z_mB_dT(t) >= m_B_max(t) * dT_t + dT_T_max(t) * m_B(t) - m_B_max(t) * dT_T_max(t)];
        constraints = [constraints, Z_mB_dT(t) <= m_B_min(t) * dT_t + dT_T_max(t) * m_B(t) - m_B_min(t) * dT_T_max(t)];
        constraints = [constraints, Z_mB_dT(t) <= m_B_max(t) * dT_t + dT_T_min(t) * m_B(t) - m_B_max(t) * dT_T_min(t)];
    end
    
    % 2. Ledge厚度与质量差分相乘的标准McCormick包络
    for t = 1:params.T
        dm_t = m_B(t+1) - m_B(t);
        
        % 标准McCormick包络约束（使用时变边界）
        constraints = [constraints, Z_Lledge_dm(t) >= L_ledge_min(t) * dm_t + dm_m_min(t) * L_ledge(t) - L_ledge_min(t) * dm_m_min(t)];
        constraints = [constraints, Z_Lledge_dm(t) >= L_ledge_max(t) * dm_t + dm_m_max(t) * L_ledge(t) - L_ledge_max(t) * dm_m_max(t)];
        constraints = [constraints, Z_Lledge_dm(t) <= L_ledge_min(t) * dm_t + dm_m_max(t) * L_ledge(t) - L_ledge_min(t) * dm_m_max(t)];
        constraints = [constraints, Z_Lledge_dm(t) <= L_ledge_max(t) * dm_t + dm_m_min(t) * L_ledge(t) - L_ledge_max(t) * dm_m_min(t)];
    end
end

%% 辅助函数：从分段McCormick结果中提取时变边界
function [dT_T_min_new, dT_T_max_new, dm_m_min_new, dm_m_max_new, L_ledge_min_new, L_ledge_max_new, m_B_min_new, m_B_max_new] = extract_time_varying_bounds_from_piecewise(...
    solution, dT_T_segments, dm_m_segments, N_T_segments, N_m_segments, shrink_factors, params)
    
    % 初始化时变边界
    dT_T_min_new = zeros(1, params.T);
    dT_T_max_new = zeros(1, params.T);
    dm_m_min_new = zeros(1, params.T);
    dm_m_max_new = zeros(1, params.T);
    L_ledge_min_new = zeros(1, params.T+1);
    L_ledge_max_new = zeros(1, params.T+1);
    m_B_min_new = zeros(1, params.T+1);
    m_B_max_new = zeros(1, params.T+1);
    
    % 提取每个时段的边界
    for t = 1:params.T
        % 找到温度差分对应的分段
        lambda_T_t = solution.lambda_T(t, :);
        [~, active_segment_T] = max(lambda_T_t);
        
        % 提取温度差分边界
        dT_T_min_new(t) = dT_T_segments(active_segment_T);
        dT_T_max_new(t) = dT_T_segments(active_segment_T + 1);
        
        % 检查温度差分解是否刚好落在分段边界处
        dT_current = solution.dT_values(t);
        if abs(dT_current - dT_T_min_new(t)) < 1e-10 || abs(dT_current - dT_T_max_new(t)) < 1e-10
            fprintf('警告：时段%d温度差分解 %.4f 刚好落在分段边界 [%.4f, %.4f] 处\n', ...
                t, dT_current, dT_T_min_new(t), dT_T_max_new(t));
        end
        
        % 找到质量差分对应的分段
        lambda_m_t = solution.lambda_m(t, :);
        [~, active_segment_m] = max(lambda_m_t);
        
        % 提取质量差分边界
        dm_m_min_new(t) = dm_m_segments(active_segment_m);
        dm_m_max_new(t) = dm_m_segments(active_segment_m + 1);
        
        % 检查质量差分解是否刚好落在分段边界处
        dm_current = solution.dm_values(t);
        if abs(dm_current - dm_m_min_new(t)) < 1e-10 || abs(dm_current - dm_m_max_new(t)) < 1e-10
            fprintf('警告：时段%d质量差分解 %.4f 刚好落在分段边界 [%.4f, %.4f] 处\n', ...
                t, dm_current, dm_m_min_new(t), dm_m_max_new(t));
        end
    end
    
    % 对于Ledge厚度和浴液质量，使用当前解进行边界收缩
    for t = 1:params.T+1
        % Ledge厚度边界收缩
        L_ledge_current = solution.L_ledge(t);
        L_ledge_range = params.L_ledge_max - params.L_ledge_min;
        L_ledge_range_target = (1 - shrink_factors.ledge) * L_ledge_range;
        L_ledge_range_target = max(L_ledge_range_target, 1e-6);
        
        L_ledge_center = max(params.L_ledge_min, min(params.L_ledge_max, L_ledge_current));
        L_ledge_min_new(t) = max(params.L_ledge_min, L_ledge_center - 0.5 * L_ledge_range_target);
        L_ledge_max_new(t) = min(params.L_ledge_max, L_ledge_center + 0.5 * L_ledge_range_target);
        
        % 浴液质量边界收缩
        m_B_current = solution.m_B(t);
        m_B_range = params.m_B_max_init - params.m_B_min_init;  % 固定范围
        m_B_range_target = (1 - shrink_factors.mass) * m_B_range;
        m_B_range_target = max(m_B_range_target, 1e-6);
        
        m_B_center = max(params.m_B_min_init, min(params.m_B_max_init, m_B_current));
        m_B_min_new(t) = max(params.m_B_min_init, m_B_center - 0.5 * m_B_range_target);
        m_B_max_new(t) = min(params.m_B_max_init, m_B_center + 0.5 * m_B_range_target);
    end
    
    fprintf('从分段McCormick结果中提取时变边界完成\n');
    fprintf('温度差分边界范围: [%.2f, %.2f]°C\n', min(dT_T_min_new), max(dT_T_max_new));
    fprintf('质量差分边界范围: [%.0f, %.0f] kg\n', min(dm_m_min_new), max(dm_m_max_new));
    fprintf('Ledge厚度边界范围: [%.3f, %.3f] m\n', min(L_ledge_min_new), max(L_ledge_max_new));
    fprintf('浴液质量边界范围: [%.0f, %.0f] kg\n', min(m_B_min_new), max(m_B_max_new));
end

%% 辅助函数：计算自适应收缩因子
function [shrink_factors, shrink_factor_info] = calculate_adaptive_shrink_factors(iter, error_history, params)
    % 自适应收缩因子计算策略
    % 基于迭代次数和收敛误差历史动态调整收缩因子
    
    % 基础收缩因子（从参数中获取）
    base_dT = params.shrink_factor_dT;
    base_dm = params.shrink_factor_dm;
    base_ledge = params.shrink_factor_ledge;
    base_mass = params.shrink_factor_mass;
    
    % 迭代次数衰减因子（随迭代次数增加而减小）
    max_iter = params.max_iterations;
    iteration_decay = max(0.1, 1 - (iter - 1) / max_iter);  % 从1衰减到0.1
    
    % 误差自适应因子（基于最近几轮的误差变化）
    error_adaptation = 1.0;  % 默认值
    if iter > 3 && ~isempty(error_history)
        % 计算最近3轮的误差变化趋势
        recent_errors = error_history(max(1, end-2):end, 2:end);  % 排除迭代次数列
        if size(recent_errors, 1) >= 2
            % 计算误差变化率
            error_change = mean(recent_errors(end, :)) - mean(recent_errors(1, :));
            if error_change > 0
                % 误差增加，减小收缩因子（更保守）
                error_adaptation = max(0.3, 1 - 0.5 * min(1, error_change / 10));
            else
                % 误差减小，保持收缩因子不变（不增加）
                error_adaptation = 1.0;
            end
        end
    end
    
    % 边界范围自适应因子（如果边界已经很窄，进一步减小收缩因子）
    boundary_adaptation = 1.0;  % 默认值
    if iter > 1
        % 这里可以根据当前边界范围进一步调整
        % 暂时使用固定值，后续可以根据实际边界范围调整
    end
    
    % 计算最终的自适应收缩因子
    shrink_factors = struct();
    shrink_factors.dT = base_dT * iteration_decay * error_adaptation * boundary_adaptation;
    shrink_factors.dm = base_dm * iteration_decay * error_adaptation * boundary_adaptation;
    shrink_factors.ledge = base_ledge * iteration_decay * error_adaptation * boundary_adaptation;
    shrink_factors.mass = base_mass * iteration_decay * error_adaptation * boundary_adaptation;
    
    % 确保收缩因子在合理范围内
    min_shrink = 0.001;  % 最小收缩因子
    max_shrink = 0.2;    % 最大收缩因子
    
    shrink_factors.dT = max(min_shrink, min(max_shrink, shrink_factors.dT));
    shrink_factors.dm = max(min_shrink, min(max_shrink, shrink_factors.dm));
    shrink_factors.ledge = max(min_shrink, min(max_shrink, shrink_factors.ledge));
    shrink_factors.mass = max(min_shrink, min(max_shrink, shrink_factors.mass));
    
    % 记录自适应原因
    shrink_factor_info = struct();
    shrink_factor_info.iteration_decay = iteration_decay;
    shrink_factor_info.error_adaptation = error_adaptation;
    shrink_factor_info.boundary_adaptation = boundary_adaptation;
    
    % 确定主要自适应原因
    if iteration_decay < 0.5
        shrink_factor_info.adaptation_reason = 1;  % 迭代次数衰减
    elseif error_adaptation < 0.8
        shrink_factor_info.adaptation_reason = 2;  % 误差增加，减小收缩因子
    else
        shrink_factor_info.adaptation_reason = 3;  % 正常，保持收缩因子
    end
end

%% 辅助函数：收缩时变边界
function [dT_T_min_new, dT_T_max_new, dm_m_min_new, dm_m_max_new, L_ledge_min_new, L_ledge_max_new, m_B_min_new, m_B_max_new] = shrink_time_varying_bounds(...
    solution, dT_T_min_current, dT_T_max_current, dm_m_min_current, dm_m_max_current, ...
    L_ledge_min_current, L_ledge_max_current, m_B_min_current, m_B_max_current, shrink_factors, params)
    
    % 提取当前解中的差分值
    dT_values = solution.dT_values;
    dm_values = solution.dm_values;
    L_ledge_values = solution.L_ledge;
    m_B_values = solution.m_B;
    
    % 初始化新边界
    dT_T_min_new = zeros(1, params.T);
    dT_T_max_new = zeros(1, params.T);
    dm_m_min_new = zeros(1, params.T);
    dm_m_max_new = zeros(1, params.T);
    L_ledge_min_new = zeros(1, params.T+1);
    L_ledge_max_new = zeros(1, params.T+1);
    m_B_min_new = zeros(1, params.T+1);
    m_B_max_new = zeros(1, params.T+1);
    
    % 对每个时段进行边界收缩
    for t = 1:params.T
        % 当前时段的边界范围
        dT_range_current = dT_T_max_current(t) - dT_T_min_current(t);
        dm_range_current = dm_m_max_current(t) - dm_m_min_current(t);

        % 目标收缩范围（确保单调收缩）
        dT_range_target = (1 - shrink_factors.dT) * dT_range_current;
        dm_range_target = (1 - shrink_factors.dm) * dm_range_current;

        % 确保目标范围为正（最小范围阈值）
        dT_range_target = max(dT_range_target, 1e-6);
        dm_range_target = max(dm_range_target, 1e-6);

        % 改进的边界收缩策略：确保当前解在新边界内
        dT_current = dT_values(t);
        dm_current = dm_values(t);
        
        % 计算当前解在旧边界中的相对位置
        dT_relative_pos = (dT_current - dT_T_min_current(t)) / dT_range_current;
        dm_relative_pos = (dm_current - dm_m_min_current(t)) / dm_range_current;
        
        % 确保相对位置在[0,1]范围内
        dT_relative_pos = max(0, min(1, dT_relative_pos));
        dm_relative_pos = max(0, min(1, dm_relative_pos));
        
        % 基于相对位置计算新边界，确保当前解在新边界内
        dT_T_min_new(t) = dT_current - dT_relative_pos * dT_range_target;
        dT_T_max_new(t) = dT_current + (1 - dT_relative_pos) * dT_range_target;
        dm_m_min_new(t) = dm_current - dm_relative_pos * dm_range_target;
        dm_m_max_new(t) = dm_current + (1 - dm_relative_pos) * dm_range_target;

        % 确保单调收缩：新边界必须在旧边界内
        dT_T_min_new(t) = max(dT_T_min_current(t), dT_T_min_new(t));
        dT_T_max_new(t) = min(dT_T_max_current(t), dT_T_max_new(t));
        dm_m_min_new(t) = max(dm_m_min_current(t), dm_m_min_new(t));
        dm_m_max_new(t) = min(dm_m_max_current(t), dm_m_max_new(t));
        
        % 最终验证：确保当前解在新边界内（安全保护）
        if dT_current < dT_T_min_new(t) || dT_current > dT_T_max_new(t)
            fprintf('警告：第%d时段温度差分解 %.4f 不在新边界 [%.4f, %.4f] 内，强制调整\n', ...
                t, dT_current, dT_T_min_new(t), dT_T_max_new(t));
            dT_T_min_new(t) = min(dT_T_min_new(t), dT_current);
            dT_T_max_new(t) = max(dT_T_max_new(t), dT_current);
        end
        
        if dm_current < dm_m_min_new(t) || dm_current > dm_m_max_new(t)
            fprintf('警告：第%d时段质量差分解 %.4f 不在新边界 [%.4f, %.4f] 内，强制调整\n', ...
                t, dm_current, dm_m_min_new(t), dm_m_max_new(t));
            dm_m_min_new(t) = min(dm_m_min_new(t), dm_current);
            dm_m_max_new(t) = max(dm_m_max_new(t), dm_current);
        end

        % 确保不突破初始物理边界
        dT_T_min_new(t) = max(params.dT_T_min_init, dT_T_min_new(t));
        dT_T_max_new(t) = min(params.dT_T_max_init, dT_T_max_new(t));
        dm_m_min_new(t) = max(params.dm_m_min_init, dm_m_min_new(t));
        dm_m_max_new(t) = min(params.dm_m_max_init, dm_m_max_new(t));

        % 确保下界严格小于上界
        if dT_T_min_new(t) >= dT_T_max_new(t)
            dT_T_min_new(t) = dT_T_max_new(t) - 1e-6;
        end
        if dm_m_min_new(t) >= dm_m_max_new(t)
            dm_m_min_new(t) = dm_m_max_new(t) - 1e-6;
        end
    end
    
    % 对Ledge厚度和浴液质量进行边界收缩
    for t = 1:params.T+1
        % Ledge厚度边界收缩（使用更保守的收缩因子）
        L_ledge_range_current = L_ledge_max_current(t) - L_ledge_min_current(t);
        L_ledge_range_target = (1 - shrink_factors.ledge) * L_ledge_range_current;
        % 设置合理的最小范围阈值（1cm）
        L_ledge_range_target = max(L_ledge_range_target, 1e-6);
        
        L_ledge_current = L_ledge_values(t);
        
        % 计算当前解在旧边界中的相对位置
        L_ledge_relative_pos = (L_ledge_current - L_ledge_min_current(t)) / L_ledge_range_current;
        L_ledge_relative_pos = max(0, min(1, L_ledge_relative_pos));
        
        % 基于相对位置计算新边界，确保当前解在新边界内
        L_ledge_min_new(t) = L_ledge_current - L_ledge_relative_pos * L_ledge_range_target;
        L_ledge_max_new(t) = L_ledge_current + (1 - L_ledge_relative_pos) * L_ledge_range_target;
        
        % 确保单调收缩
        L_ledge_min_new(t) = max(L_ledge_min_current(t), L_ledge_min_new(t));
        L_ledge_max_new(t) = min(L_ledge_max_current(t), L_ledge_max_new(t));
        
        % 最终验证：确保当前解在新边界内（安全保护）
        if L_ledge_current < L_ledge_min_new(t) || L_ledge_current > L_ledge_max_new(t)
            fprintf('警告：第%d时段Ledge厚度解 %.4f 不在新边界 [%.4f, %.4f] 内，强制调整\n', ...
                t, L_ledge_current, L_ledge_min_new(t), L_ledge_max_new(t));
            L_ledge_min_new(t) = min(L_ledge_min_new(t), L_ledge_current);
            L_ledge_max_new(t) = max(L_ledge_max_new(t), L_ledge_current);
        end
        
        % 确保不突破初始物理边界
        L_ledge_min_new(t) = max(params.L_ledge_min, L_ledge_min_new(t));
        L_ledge_max_new(t) = min(params.L_ledge_max, L_ledge_max_new(t));
        
        % 确保下界严格小于上界
        if L_ledge_min_new(t) >= L_ledge_max_new(t)
            L_ledge_min_new(t) = L_ledge_max_new(t) - 1e-6;
        end
        
        % 浴液质量边界收缩（使用更保守的收缩因子）
        m_B_range_current = m_B_max_current(t) - m_B_min_current(t);
        m_B_range_target = (1 - shrink_factors.mass) * m_B_range_current;
        % 设置合理的最小范围阈值（500kg）
        m_B_range_target = max(m_B_range_target, 1e-6);
        
        m_B_current = m_B_values(t);
        
        % 计算当前解在旧边界中的相对位置
        m_B_relative_pos = (m_B_current - m_B_min_current(t)) / m_B_range_current;
        m_B_relative_pos = max(0, min(1, m_B_relative_pos));
        
        % 基于相对位置计算新边界，确保当前解在新边界内
        m_B_min_new(t) = m_B_current - m_B_relative_pos * m_B_range_target;
        m_B_max_new(t) = m_B_current + (1 - m_B_relative_pos) * m_B_range_target;
        
        % 确保单调收缩
        m_B_min_new(t) = max(m_B_min_current(t), m_B_min_new(t));
        m_B_max_new(t) = min(m_B_max_current(t), m_B_max_new(t));
        
        % 最终验证：确保当前解在新边界内（安全保护）
        if m_B_current < m_B_min_new(t) || m_B_current > m_B_max_new(t)
            fprintf('警告：第%d时段浴液质量解 %.0f 不在新边界 [%.0f, %.0f] 内，强制调整\n', ...
                t, m_B_current, m_B_min_new(t), m_B_max_new(t));
            m_B_min_new(t) = min(m_B_min_new(t), m_B_current);
            m_B_max_new(t) = max(m_B_max_new(t), m_B_current);
        end
        
        % 确保不突破初始物理边界
        m_B_min_new(t) = max(params.m_B_min_init, m_B_min_new(t));
        m_B_max_new(t) = min(params.m_B_max_init, m_B_max_new(t));
        
        % 确保下界严格小于上界
        if m_B_min_new(t) >= m_B_max_new(t)
            m_B_min_new(t) = m_B_max_new(t) - 1e-6;
        end
    end
    
    fprintf('时变边界收缩完成\n');
    fprintf('温度差分边界范围: [%.2f, %.2f]°C (平均范围: %.2f)\n', min(dT_T_min_new), max(dT_T_max_new), mean(dT_T_max_new - dT_T_min_new));
    fprintf('质量差分边界范围: [%.0f, %.0f] kg (平均范围: %.0f)\n', min(dm_m_min_new), max(dm_m_max_new), mean(dm_m_max_new - dm_m_min_new));
    fprintf('Ledge厚度边界范围: [%.3f, %.3f] m (平均范围: %.3f)\n', min(L_ledge_min_new), max(L_ledge_max_new), mean(L_ledge_max_new - L_ledge_min_new));
    fprintf('浴液质量边界范围: [%.0f, %.0f] kg (平均范围: %.0f)\n', min(m_B_min_new), max(m_B_max_new), mean(m_B_max_new - m_B_min_new));
end

%% 通用函数：构建优化问题
function [constraints, variables, obj] = build_optimization_problem(...
    dT_T_min, dT_T_max, dm_m_min, dm_m_max, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params, model_type, varargin)
    
    % 解析输入参数
    if strcmp(model_type, 'piecewise')
        N_T_segments = varargin{1};
        N_m_segments = varargin{2};
        dT_T_segments = varargin{3};
        dm_m_segments = varargin{4};
    end
    
    % 定义决策变量
    variables = struct();
    variables.I = sdpvar(1, params.T, 'full');
    variables.m_feed_Al2O3 = sdpvar(1, params.T, 'full');
    variables.T_B = sdpvar(1, params.T+1, 'full');
    variables.L_ledge = sdpvar(1, params.T+1, 'full');
    variables.m_B = sdpvar(1, params.T+1, 'full');
    variables.C_Al2O3 = sdpvar(1, params.T+1, 'full');
    variables.C_AlF3 = sdpvar(1, params.T+1, 'full');
    variables.C_CaF2 = sdpvar(1, params.T+1, 'full');
    variables.C_LiF = sdpvar(1, params.T+1, 'full');
    variables.T_liq = sdpvar(1, params.T, 'full');
    variables.Superheat = sdpvar(1, params.T, 'full');
    variables.g = sdpvar(1, params.T, 'full');
    
    % McCormick包络变量
    variables.Z_I_g = sdpvar(1, params.T, 'full');
    variables.Z_C_Al2O3_m = sdpvar(1, params.T+1, 'full');
    variables.Z_C_AlF3_m = sdpvar(1, params.T+1, 'full');
    variables.Z_C_CaF2_m = sdpvar(1, params.T+1, 'full');
    variables.Z_C_LiF_m = sdpvar(1, params.T+1, 'full');
    variables.Z_mB_dT = sdpvar(1, params.T, 'full');
    variables.Z_Lledge_dm = sdpvar(1, params.T, 'full');
    variables.Z_Lledge_superheat = sdpvar(1, params.T, 'full');
    
    % 调节动作变量
    variables.lamda_keep = binvar(1, params.T-1, 'full');
    variables.lamda_up = binvar(1, params.T-1, 'full');
    variables.lamda_down = binvar(1, params.T-1, 'full');
    
    % 分段变量（仅用于分段模型）
    if strcmp(model_type, 'piecewise')
        variables.lambda_T = binvar(params.T, N_T_segments, 'full');
        variables.lambda_m = binvar(params.T, N_m_segments, 'full');
    end
    
    % 定义约束
    constraints = [];
    
    % 初始条件
    constraints = [constraints, variables.I(1) == params.I_N];
    constraints = [constraints, variables.T_B(1) == params.T_B_init];
    constraints = [constraints, variables.L_ledge(1) == params.L_ledge_init];
    constraints = [constraints, variables.m_B(1) == params.m_B_init];
    constraints = [constraints, variables.C_Al2O3(1) == params.C_Al2O3_init];
    constraints = [constraints, variables.C_AlF3(1) == params.C_AlF3_init];
    constraints = [constraints, variables.C_CaF2(1) == params.C_CaF2_init];
    constraints = [constraints, variables.C_LiF(1) == params.C_LiF_init];
    
    % 基本约束
    constraints = [constraints, params.I_min <= variables.I <= params.I_max];
    constraints = [constraints, params.m_feed_Al2O3_min <= variables.m_feed_Al2O3 <= params.m_feed_Al2O3_max];
    constraints = [constraints, params.g_min <= variables.g <= params.g_max];
    constraints = [constraints, params.T_B_min <= variables.T_B <= params.T_B_max];
    constraints = [constraints, params.Superheat_min <= variables.Superheat <= params.Superheat_max];
    constraints = [constraints, params.C_Al2O3_min <= variables.C_Al2O3 <= params.C_Al2O3_max];
    constraints = [constraints, params.L_ledge_min <= variables.L_ledge <= params.L_ledge_max];
    
    % % 差分变量边界约束
    % for t = 1:params.T
    %     constraints = [constraints, dT_T_min(t) <= variables.T_B(t+1) - variables.T_B(t) <= dT_T_max(t)];
    %     constraints = [constraints, dm_m_min(t) <= variables.m_B(t+1) - variables.m_B(t) <= dm_m_max(t)];
    % end
    
    % 全局McCormick包络约束
    constraints = add_global_mccormick_constraints(constraints, variables.I, variables.g, variables.C_Al2O3, variables.C_AlF3, variables.C_CaF2, variables.C_LiF, variables.m_B, variables.L_ledge, variables.T_B, variables.T_liq, ...
        variables.Z_I_g, variables.Z_C_Al2O3_m, variables.Z_C_AlF3_m, variables.Z_C_CaF2_m, variables.Z_C_LiF_m, variables.Z_Lledge_superheat, L_ledge_min, L_ledge_max, m_B_min, m_B_max, params);
    
    % McCormick包络约束（根据模型类型选择）
    if strcmp(model_type, 'piecewise')
        % 分段McCormick包络约束
        constraints = add_piecewise_mccormick_constraints(constraints, variables.T_B, variables.m_B, variables.L_ledge, ...
            variables.Z_mB_dT, variables.Z_Lledge_dm, variables.lambda_T, variables.lambda_m, ...
            dT_T_segments, dm_m_segments, params, N_T_segments, N_m_segments);
    else
        % 标准McCormick包络约束
        constraints = add_standard_mccormick_constraints(constraints, variables.T_B, variables.m_B, variables.L_ledge, ...
            variables.Z_mB_dT, variables.Z_Lledge_dm, dT_T_min, dT_T_max, dm_m_min, dm_m_max, ...
            L_ledge_min, L_ledge_max, m_B_min, m_B_max, params);
    end
    
    % 动态约束
    constraints = add_dynamic_constraints(constraints, variables.T_B, variables.L_ledge, variables.m_B, variables.C_Al2O3, variables.C_AlF3, variables.C_CaF2, variables.C_LiF, ...
        variables.T_liq, variables.Superheat, variables.g, variables.I, variables.m_feed_Al2O3, variables.Z_I_g, variables.Z_C_Al2O3_m, variables.Z_C_AlF3_m, variables.Z_C_CaF2_m, variables.Z_C_LiF_m, ...
        variables.Z_mB_dT, variables.Z_Lledge_dm, variables.Z_Lledge_superheat, params);
    
    % 调节约束
    constraints = add_regulation_constraints(constraints, variables.I, variables.lamda_keep, variables.lamda_up, variables.lamda_down, params);
    
    % 目标函数
    P_N = (params.EMF + params.R_cell*params.I_N*1000)*params.I_N;
    dP_dI_N = params.EMF + 2*params.R_cell*params.I_N*1000;
    P = P_N + dP_dI_N*(variables.I - repmat(params.I_N,1,params.T));
    
    Cost_elec = sum(P .* params.price_TOU') * params.delta_t;
    Cost_Al2O3 = sum(variables.m_feed_Al2O3 .* params.price_Al2O3) * params.delta_t;
    Y_Al = 3600*1e3 * (variables.Z_I_g * params.M_Al) / (params.F*params.z);
    Revenue_Al = sum(Y_Al .* params.price_AL) * params.delta_t;
    obj = Cost_elec + Cost_Al2O3 - Revenue_Al;
end

%% 通用函数：求解优化问题
function [obj_value, solution, feasible] = solve_optimization_problem(constraints, obj, variables, params)
    
    % 求解
    ops = sdpsettings('solver', 'gurobi', 'verbose', 1);
    result = optimize(constraints, obj, ops);
    
    if result.problem == 0
        feasible = true;
        obj_value = value(obj);
        
        % 提取差分值
        dT_values = value(variables.T_B(2:end)) - value(variables.T_B(1:end-1));
        dm_values = value(variables.m_B(2:end)) - value(variables.m_B(1:end-1));
        
        % 计算McCormick包络误差
        mccormick_errors = calculate_mccormick_errors(variables.I, variables.g, variables.T_B, variables.m_B, variables.L_ledge, variables.T_liq, ...
            variables.Z_I_g, variables.Z_C_Al2O3_m, variables.Z_C_AlF3_m, variables.Z_C_CaF2_m, variables.Z_C_LiF_m, variables.Z_mB_dT, variables.Z_Lledge_dm, variables.Z_Lledge_superheat, variables, params);
        
        % 构建解结构
        solution = struct();
        solution.I = value(variables.I);
        solution.m_feed_Al2O3 = value(variables.m_feed_Al2O3);
        solution.T_B = value(variables.T_B);
        solution.L_ledge = value(variables.L_ledge);
        solution.m_B = value(variables.m_B);
        solution.C_Al2O3 = value(variables.C_Al2O3);
        solution.C_AlF3 = value(variables.C_AlF3);
        solution.C_CaF2 = value(variables.C_CaF2);
        solution.C_LiF = value(variables.C_LiF);
        solution.T_liq = value(variables.T_liq);
        solution.Superheat = value(variables.Superheat);
        solution.g = value(variables.g);
        solution.mccormick_errors = mccormick_errors;
        
        % McCormick包络变量
        solution.Z_I_g = value(variables.Z_I_g);
        solution.Z_C_Al2O3_m = value(variables.Z_C_Al2O3_m);
        solution.Z_C_AlF3_m = value(variables.Z_C_AlF3_m);
        solution.Z_C_CaF2_m = value(variables.Z_C_CaF2_m);
        solution.Z_C_LiF_m = value(variables.Z_C_LiF_m);
        solution.Z_mB_dT = value(variables.Z_mB_dT);
        solution.Z_Lledge_dm = value(variables.Z_Lledge_dm);
        solution.Z_Lledge_superheat = value(variables.Z_Lledge_superheat);
        
        % 计算功率和产量
        P_N = (params.EMF + params.R_cell*params.I_N*1000)*params.I_N;
        dP_dI_N = params.EMF + 2*params.R_cell*params.I_N*1000;
        P = P_N + dP_dI_N*(solution.I - repmat(params.I_N,1,params.T));
        Y_Al = 3600*1e3 * (value(variables.Z_I_g) * params.M_Al) / (params.F*params.z);
        
        solution.P = P;
        solution.Y_Al = Y_Al;
        solution.dT_values = dT_values;
        solution.dm_values = dm_values;
    else
        feasible = false;
        obj_value = inf;
        solution = [];
    end
end