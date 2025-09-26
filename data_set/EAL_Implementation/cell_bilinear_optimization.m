% 电解铝负荷用能优化模型（双线性版本）
% 基于简化电解槽线性仿真模型的负荷优化，保留原始双线性约束
% clc; clear; yalmip("clear"); 
close all

%% 初始化参数
% 运行参数初始化脚本
run('initialize_parameters.m');

% 计算初始质量
m_Al2O3_init = C_Al2O3_init * m_B_init;
m_AlF3_init = C_AlF3_init * m_B_init;
m_CaF2_init = C_CaF2_init * m_B_init;
m_LiF_init = C_LiF_init * m_B_init;

%% 定义决策变量
I = sdpvar(1, T, 'full');              % 电流变量 (kA)
m_feed_Al2O3 = sdpvar(1, T, 'full');   % 氧化铝进料速率 (kg/min)

% 状态变量（每小时更新一次）
T_B = sdpvar(1, T+1, 'full');          % 浴液温度 (°C)
L_ledge = sdpvar(1, T+1, 'full');      % ledge厚度 (m)
m_B = sdpvar(1, T+1, 'full');          % 浴液质量 (kg)
C_Al2O3 = sdpvar(1, T+1, 'full');      % 氧化铝浓度
C_AlF3 = sdpvar(1, T+1, 'full');       % 氟化铝浓度
C_CaF2 = sdpvar(1, T+1, 'full');       % 氟化钙浓度
C_LiF = sdpvar(1, T+1, 'full');        % 氟化锂浓度

% 辅助变量
T_liq = sdpvar(1, T, 'full');          % 液相温度 (°C)（每小时更新一次）
Superheat = sdpvar(1, T, 'full');      % 过热度 (°C)（每小时更新一次）

% 电流效率变量（温度相关）
g = sdpvar(1, T, 'full');              % 电流效率（每小时更新一次）

% 调节动作变量
lamda_keep = binvar(1, T-1, 'full');   % 保持状态
lamda_up = binvar(1, T-1, 'full');     % 上调状态
lamda_down = binvar(1, T-1, 'full');   % 下调状态

%% 约束条件
constraints = [];

% 初始条件约束
constraints = [constraints, I(1) == I_N];
constraints = [constraints, T_B(1) == T_B_init];
constraints = [constraints, L_ledge(1) == L_ledge_init];
constraints = [constraints, m_B(1) == m_B_init];
constraints = [constraints, C_Al2O3(1) == C_Al2O3_init];
constraints = [constraints, C_AlF3(1) == C_AlF3_init];
constraints = [constraints, C_CaF2(1) == C_CaF2_init];
constraints = [constraints, C_LiF(1) == C_LiF_init];

% 电流上下限约束
constraints = [constraints, I_min <= I <= I_max];

% 氧化铝进料约束
constraints = [constraints, m_feed_Al2O3_min <= m_feed_Al2O3 <= m_feed_Al2O3_max];

% 电流效率范围约束（从参数文件读取）
constraints = [constraints, g_min <= g <= g_max];

% 调节状态唯一性约束
constraints = [constraints, lamda_keep + lamda_up + lamda_down == ones(1, T-1)];

% 功率计算（线性化处理，在额定电流附近泰勒展开）
P_N = (EMF + R_cell*I_N*1000)*I_N;  % 额定功率 (kW)
dP_dI_N = EMF + 2*R_cell*I_N*1000;  % 功率对电流的导数在额定点的值
P = P_N + dP_dI_N*(I - repmat(I_N,1,T)); % 线性化功率表达式

% 验证线性化精度（可选，用于调试）
fprintf('=== 功率线性化验证 ===\n');
fprintf('额定电流: %.1f kA\n', I_N);
fprintf('额定功率: %.2f kW\n', P_N);
fprintf('功率对电流导数: %.2f kW/kA\n', dP_dI_N);

% 计算线性化误差示例
I_test = [466, 527];
P_exact = (EMF + R_cell*I_test*1000).*I_test;
P_linear = P_N + dP_dI_N*(I_test - I_N);
error_pct = abs(P_exact - P_linear) ./ P_exact * 100;

fprintf('电流范围: %.1f - %.1f kA\n', I_min, I_max);
fprintf('最大线性化误差: %.2f%%\n', max(error_pct));
fprintf('========================\n\n');

% 功率爬坡约束（基于每小时功率变化）
constraints = [constraints, P(2:end) - P(1:end-1) >= -R_down * lamda_down * P_N];
constraints = [constraints, P(2:end) - P(1:end-1) <= R_up * lamda_up * P_N];

% 最大调节次数约束
constraints = [constraints, sum(lamda_up + lamda_down) <= N_rel];

% 连续调节约束
constraints = [constraints, lamda_up(2:end) + lamda_down(2:end) + ...
    lamda_up(1:end-1) + lamda_down(1:end-1) <= ones(1, T-2)];

% 动态约束（矩阵形式）
% 定义系数矩阵和常数向量
coeff_T_liq = beta_linear; % [Al2O3, AlF3, CaF2, 常数项]

% 液相温度计算（矩阵形式）
constraints = [constraints, ...
    T_liq == coeff_T_liq(1) * C_Al2O3(1:T)*100 + ...
         coeff_T_liq(2) * C_AlF3(1:T)*100 + ...
         coeff_T_liq(3) * C_CaF2(1:T)*100 + ...
         coeff_T_liq(4)];

% 过热度计算（矩阵形式）
constraints = [constraints, Superheat == T_B(1:T) - T_liq];

% 电流效率计算（温度相关，矩阵形式）
constraints = [constraints, g == coeff_global(1) + ...
                                coeff_global(2)*Superheat + ...
                                coeff_global(3)*C_AlF3(1:T)*100 + ...
                                coeff_global(4)*C_LiF(1:T)*100];

% 状态变量更新约束（使用原始双线性表达式）
% 浴液温度更新约束（保留双线性项）
constraints = [constraints, ...
    m_B(1:T) .* c_B .* (T_B(2:T+1) - T_B(1:T)) + m_M * c_M * (T_B(2:T+1) - T_B(1:T)) == ...
    3600*delta_t * (P(1:T)*1000 - ...
    (T_B(1:T) - T_gas)/R_B_gas - ...
    (T_B(1:T) - T_liq)/R_B_ledge - ...
    (T_B(1:T) - T_amb)/R_bottom_total - ...
    0.596*60*1e3*m_feed_Al2O3(1:T) - ...
    5.474*3600*1e3*(I(1:T) .* g(1:T)*1000*M_Al)/(F*z))];

% 浴液质量更新约束（保留双线性项）
% 计算固定的热阻分量
R_shell_side = 0.25/(47 * A_ledge); % 槽壳侧部导热 (K/W)
R_shell_side_amb = 1/(22 * A_ledge); % 槽壳侧部到环境对流 (K/W)

% 计算每个时刻的ledge导热热阻
R_ledge_t = L_ledge(1:T) / (1.5 * A_ledge); % 每个时刻的ledge导热热阻 (K/W)
R_ledge_amb_t = R_ledge_t + R_shell_side + R_shell_side_amb; % 每个时刻的总ledge到环境热阻 (K/W)

% 使用原始双线性表达式重新整理浴液质量更新约束
coeff_m_B = delta_t*3600 / lf;
constraints = [constraints, ...
    R_ledge_amb_t .* (m_B(2:T+1) - m_B(1:T)) == ...
    coeff_m_B * (R_ledge_amb_t / R_B_ledge .* (T_B(1:T) - T_liq) - (T_liq - T_amb))];

% ledge厚度更新约束
coeff_L_ledge = 1 / (rho_ledge * A_ledge);
constraints = [constraints, ...
    L_ledge(2:T+1) == L_ledge(1:T) - ...
    coeff_L_ledge * (m_B(2:T+1) - m_B(1:T))];

% 各组分浓度更新约束（保留双线性项）
coeff_m_Al2O3_feed = 60*delta_t;
coeff_m_Al2O3_reaction = -3600*delta_t * M_Al2O3 / (F*2*z);

% 氧化铝浓度更新约束（保留双线性项）
constraints = [constraints, ...
    C_Al2O3(2:T+1) .* m_B(2:T+1) == C_Al2O3(1:T) .* m_B(1:T) + ...
    coeff_m_Al2O3_feed * m_feed_Al2O3(1:T) + ...
    coeff_m_Al2O3_reaction * (I(1:T) .* g(1:T)*1000)];

% 其他组分浓度更新约束
constraints = [constraints, C_AlF3(2:T+1) .* m_B(2:T+1) == C_AlF3_init*m_B_init];
constraints = [constraints, C_CaF2(2:T+1) .* m_B(2:T+1) == C_CaF2_init*m_B_init];
constraints = [constraints, C_LiF(2:T+1) .* m_B(2:T+1) == C_LiF_init*m_B_init];

% 温度约束
constraints = [constraints, T_B_min <= T_B <= T_B_max];

% 过热度约束
constraints = [constraints, Superheat_min <= Superheat <= Superheat_max];

% 氧化铝浓度约束
constraints = [constraints, C_Al2O3_min <= C_Al2O3 <= C_Al2O3_max];

% ledge厚度约束
constraints = [constraints, L_ledge_min <= L_ledge <= L_ledge_max];

%% 目标函数
% 电费成本（功率每小时更新一次）
Cost_elec = sum(P .* price_TOU') * delta_t;

% 氧化铝用料成本
Cost_Al2O3 = sum(m_feed_Al2O3 .* price_Al2O3) * delta_t;

% 产品收益（使用原始双线性表达式）
Y_Al = 3600*1e3 * (I .* g * M_Al) / (F*z); % kg/h
Revenue_Al = sum(Y_Al .* price_AL) * delta_t;

% 总目标函数（最小化成本）
obj = Cost_elec + Cost_Al2O3 - Revenue_Al;

%% 求解
ops = sdpsettings('solver', 'gurobi', 'verbose', 1);
result = optimize(constraints, obj, ops);

%% 结果分析
if result.problem == 0
    fprintf('优化求解成功！\n');
    
    % 提取优化结果
    I_opt = value(I);
    m_feed_Al2O3_opt = value(m_feed_Al2O3);
    T_B_opt = value(T_B);
    L_ledge_opt = value(L_ledge);
    m_B_opt = value(m_B);
    C_Al2O3_opt = value(C_Al2O3);
    C_AlF3_opt = value(C_AlF3);
    C_CaF2_opt = value(C_CaF2);
    C_LiF_opt = value(C_LiF);
    P_input_opt = value(P);
    Y_Al_opt = value(Y_Al);
    T_liq_opt = value(T_liq);
    Superheat_opt = value(Superheat);
    g_opt = value(g);
    lamda_up_opt = value(lamda_up);
    lamda_down_opt = value(lamda_down);
    lamda_keep_opt = value(lamda_keep);
    
    % 计算各项成本
    elec_cost = value(Cost_elec);
    Al2O3_cost = value(Cost_Al2O3);
    Al_revenue = value(Revenue_Al);
    total_cost = elec_cost + Al2O3_cost - Al_revenue;
    
    fprintf('\n=== 优化结果 ===\n');
    fprintf('电费成本：%.2f元\n', elec_cost);
    fprintf('氧化铝用料成本：%.2f元\n', Al2O3_cost);
    fprintf('铝产品收益：%.2f元\n', Al_revenue);
    fprintf('总运行成本：%.2f元\n', total_cost);
    fprintf('总铝产量：%.2f kg\n', sum(Y_Al_opt));
    fprintf('平均电流效率：%.4f\n', mean(g_opt));
    fprintf('电流效率范围：%.4f - %.4f\n', min(g_opt), max(g_opt));
    
    % 保存优化结果
    save('electrolysis_bilinear_optimization_results.mat', 'I_opt', 'm_feed_Al2O3_opt', 'T_B_opt', ...
         'L_ledge_opt', 'm_B_opt', 'C_Al2O3_opt', 'C_AlF3_opt', 'C_CaF2_opt', 'C_LiF_opt', ...
         'P_input_opt', 'Y_Al_opt', 'T_liq_opt', 'Superheat_opt', 'g_opt', ...
         'lamda_up_opt', 'lamda_down_opt', 'lamda_keep_opt', ...
         'elec_cost', 'Al2O3_cost', 'Al_revenue', 'total_cost');
    
    fprintf('优化结果已保存到 electrolysis_bilinear_optimization_results.mat\n');
    
    %% 可视化结果
    time_hr = 1:T;
    time_state = (1:T) * delta_t; % 1小时间隔的时间点
    
    % 第一张图：电价、电流、功率、氧化铝进料、电流效率、产量
    figure('Name', '双线性优化结果 - 控制变量和输出变量', 'Position', [100, 100, 1200, 800]);
    
    % 电价曲线
    subplot(3,2,1);
    stairs(time_hr, price_TOU, 'r-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('电价 (元/kWh)');
    title('分时电价曲线');
    grid on;
    
    % 电流曲线
    subplot(3,2,2);
    stairs(time_hr, I_opt, 'b-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('电流 (kA)');
    title('优化电流曲线（每小时更新）');
    grid on;
    
    % 功率曲线
    subplot(3,2,3);
    stairs(time_hr, P_input_opt, 'g-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('功率 (kW)');
    title('优化功率曲线（每小时更新）');
    grid on;
    
    % 氧化铝进料曲线
    subplot(3,2,4);
    stairs(time_hr, m_feed_Al2O3_opt, 'm-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('氧化铝进料速率 (kg/min)');
    title('氧化铝进料优化曲线（每小时更新）');
    grid on;
    
    % 电流效率曲线（简化模型 vs 复杂模型）
    subplot(3,2,5);
    
    % 计算真实电流效率（复杂模型）
    % 使用优化结果中的过热度、AlF3浓度和LiF浓度
    g_real = 1 - 10.^(0.0095*Superheat_opt - ...
        0.019*C_AlF3_opt(1:end-1)*100 - ...
        0.06*C_LiF_opt(1:end-1)*100 + ...
        0.982151 - 2);
    
    % 绘制两条曲线进行对比
    plot(time_state, g_opt, 'c-', 'LineWidth', 2, 'DisplayName', '简化模型');
    hold on;
    plot(time_state, g_real, 'r--', 'LineWidth', 2, 'DisplayName', '复杂模型');
    xlabel('时间 (h)');
    ylabel('电流效率');
    title('电流效率变化对比（每小时更新）');
    legend('show', 'Location', 'best');
    grid on;
    ylim([0.85, 0.98]);
    hold off;
    
    % 铝产量曲线
    subplot(3,2,6);
    plot(time_state, Y_Al_opt, 'k-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('铝产量 (kg/h)');
    title('铝产量曲线（每小时更新）');
    grid on;
    
    % 第二张图：电解质温度、液相温度、过热度、氧化铝浓度、氧化铝质量、浴液质量、ledge厚度
    figure('Name', '双线性优化结果 - 状态变量', 'Position', [100, 100, 1400, 1000]);
    
    % 电解质温度曲线
    subplot(3,2,1);
    plot(time_state, T_B_opt(1:end-1), 'r-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('电解质温度 (°C)');
    title('电解质温度变化（每小时更新）');
    grid on;
    
    % 液相温度曲线（简化模型 vs 复杂模型）
    subplot(3,2,2);
    
    % 计算真实液相温度（复杂模型）
    % 各组分的化学式代表其浓度的100倍
    Al2O3_conc = C_Al2O3_opt(1:end-1) * 100;  % 转换为100倍浓度
    AlF3_conc = C_AlF3_opt(1:end-1) * 100;    % 转换为100倍浓度
    CaF2_conc = C_CaF2_opt(1:end-1) * 100;    % 转换为100倍浓度
    LiF_conc = C_LiF_opt(1:end-1) * 100;      % 转换为100倍浓度
    MgF2_conc = 0.01 * ones(size(Al2O3_conc)) * 100; % 假设MgF2浓度固定为0.01%
    
    % 复杂液相温度模型
    T_liq_real = 1011 + 0.5*AlF3_conc - 0.13*AlF3_conc.^2.2 - ...
                  (3.45*CaF2_conc)./(1+0.0173*CaF2_conc) + ...
                  0.124*CaF2_conc.*AlF3_conc - ...
                  0.00542*(CaF2_conc.*AlF3_conc).^1.5 - ...
                  (7.93*Al2O3_conc)./(1+0.0936*Al2O3_conc-0.0017*Al2O3_conc.^2-0.0023*AlF3_conc.*Al2O3_conc) - ...
                  3.95*MgF2_conc;
    
    % 绘制两条曲线进行对比
    plot(time_state, T_liq_opt, 'b-', 'LineWidth', 2, 'DisplayName', '简化模型');
    hold on;
    plot(time_state, T_liq_real, 'r--', 'LineWidth', 2, 'DisplayName', '复杂模型');
    xlabel('时间 (h)');
    ylabel('液相温度 (°C)');
    title('液相温度变化对比（每小时更新）');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
    
    % 过热度曲线
    subplot(3,2,3);
    plot(time_state, Superheat_opt, 'g-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('过热度 (°C)');
    title('过热度变化（每小时更新）');
    grid on;
    
    % 各组分浓度变化曲线
    subplot(3,2,4);
    plot(time_state, Al2O3_conc, 'm-', 'LineWidth', 2, 'DisplayName', 'Al2O3');
    hold on;
    plot(time_state, AlF3_conc, 'r-', 'LineWidth', 2, 'DisplayName', 'AlF3');
    plot(time_state, CaF2_conc, 'g-', 'LineWidth', 2, 'DisplayName', 'CaF2');
    plot(time_state, LiF_conc, 'b-', 'LineWidth', 2, 'DisplayName', 'LiF');
    xlabel('时间 (h)');
    ylabel('浓度 (wt%)');
    title('各组分浓度变化（每小时更新）');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
    
    % 浴液质量曲线
    subplot(3,2,5);
    plot(time_state, m_B_opt(1:end-1), 'c-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('浴液质量 (kg)');
    title('浴液质量变化（每小时更新）');
    grid on;
    
    % Ledge厚度曲线
    subplot(3,2,6);
    plot(time_state, L_ledge_opt(1:end-1), 'k-', 'LineWidth', 2);
    xlabel('时间 (h)');
    ylabel('Ledge厚度 (m)');
    title('Ledge厚度变化（每小时更新）');
    grid on;
    
    fprintf('双线性优化模型求解完成！\n');
    
else
    fprintf('优化求解失败，错误代码：%d\n', result.problem);
    fprintf('错误信息：%s\n', result.info);
end