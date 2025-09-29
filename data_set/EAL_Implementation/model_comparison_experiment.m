% 电解铝模型对比实验脚本
% 比较基准模型和当前模型的效果
% 基准模型：不考虑ledge厚度、氧化铝浓度、过热度等状态变量的约束
% 当前模型：考虑所有状态变量的完整约束，与cell_bilinear_optimization保持一致
clc; clear; yalmip("clear"); close all

%% 初始化参数
% 运行参数初始化脚本
run('initialize_parameters.m');

% 计算初始质量
m_Al2O3_init = C_Al2O3_init * m_B_init;
m_AlF3_init = C_AlF3_init * m_B_init;
m_CaF2_init = C_CaF2_init * m_B_init;
m_LiF_init = C_LiF_init * m_B_init;

%% 第一步：求解基准模型（简化模型）
fprintf('=== 求解基准模型（简化模型）===\n');

% 定义决策变量
I_benchmark = sdpvar(1, T, 'full');              % 电流变量 (kA)
m_feed_Al2O3_benchmark = sdpvar(1, T, 'full');   % 氧化铝进料速率 (kg/min)

% 状态变量（只考虑浴液温度）
T_B_benchmark = sdpvar(1, T+1, 'full');          % 浴液温度 (°C)

% 调节动作变量
lamda_keep = binvar(1, T-1, 'full');   % 保持状态
lamda_up = binvar(1, T-1, 'full');     % 上调状态
lamda_down = binvar(1, T-1, 'full');   % 下调状态

% 约束条件
constraints_benchmark = [];

% 初始条件约束
constraints_benchmark = [constraints_benchmark, I_benchmark(1) == I_N];
constraints_benchmark = [constraints_benchmark, T_B_benchmark(1) == T_B_init];

% 电流上下限约束
constraints_benchmark = [constraints_benchmark, I_min <= I_benchmark <= I_max];

% 氧化铝进料约束
constraints_benchmark = [constraints_benchmark, m_feed_Al2O3_min <= m_feed_Al2O3_benchmark <= m_feed_Al2O3_max];

% 功率计算（线性化处理）
P_N = (EMF + R_cell*I_N*1000)*I_N;  % 额定功率 (kW)
dP_dI_N = EMF + 2*R_cell*I_N*1000;  % 功率对电流的导数在额定点的值
P_benchmark = P_N + dP_dI_N*(I_benchmark - repmat(I_N,1,T)); % 线性化功率表达式

% 浴液温度更新约束（基准模型，使用简化的热阻）
% 计算熔融区总质量和平均比热容
m_melt = m_B_init + m_M;  % 熔融区总质量 (kg)
c_melt = (m_B_init*c_B + m_M*c_M)/m_melt; % 熔融区平均比热容 (J/kg·K)

% 简化的热阻（不考虑ledge厚度变化）
R_ledge_amb_simple = 0.005879;  % 简化的ledge到环境热阻 (K/W)

coeff_T_B = 3600*delta_t / (m_melt * c_melt);
constraints_benchmark = [constraints_benchmark, ...
    T_B_benchmark(2:T+1) == T_B_benchmark(1:T) + ...
    coeff_T_B * (P_benchmark(1:T)*1000 - ...
    (T_B_benchmark(1:T) - T_gas)/R_B_gas - ...
    (T_B_benchmark(1:T) - T_amb)/(R_B_ledge+R_ledge_amb_simple) - ...
    (T_B_benchmark(1:T) - T_amb)/R_bottom_total - ...
    0.596*60*1e3*m_feed_Al2O3_benchmark(1:T) - ...
    5.474*3600*1e3*(I_benchmark(1:T)*1000.*g_N*M_Al)/(F*z))];

% 温度约束
constraints_benchmark = [constraints_benchmark, T_B_min <= T_B_benchmark <= T_B_max];

% 功率爬坡约束（基于每小时功率变化）
constraints_benchmark = [constraints_benchmark, P_benchmark(2:end) - P_benchmark(1:end-1) >= -R_down * lamda_down * P_N];
constraints_benchmark = [constraints_benchmark, P_benchmark(2:end) - P_benchmark(1:end-1) <= R_up * lamda_up * P_N];

% 最大调节次数约束
constraints_benchmark = [constraints_benchmark, sum(lamda_up + lamda_down) <= N_rel];

% 连续调节约束
constraints_benchmark = [constraints_benchmark, lamda_up(2:end) + lamda_down(2:end) + ...
    lamda_up(1:end-1) + lamda_down(1:end-1) <= ones(1, T-2)];

% 目标函数（最小化成本）
Cost_elec_benchmark = sum(P_benchmark .* price_TOU') * delta_t;
Cost_Al2O3_benchmark = sum(m_feed_Al2O3_benchmark .* price_Al2O3) * delta_t;
Y_Al_benchmark = 3600*1e3 * (I_benchmark(1:T) .* g_N * M_Al) / (F*z); % kg/h
Revenue_Al_benchmark = sum(Y_Al_benchmark .* price_AL) * delta_t;
obj_benchmark = Cost_elec_benchmark + Cost_Al2O3_benchmark - Revenue_Al_benchmark;

% 求解基准模型
ops = sdpsettings('solver', 'gurobi', 'verbose', 1);
result_benchmark = optimize(constraints_benchmark, obj_benchmark, ops);

if result_benchmark.problem == 0
    fprintf('基准模型求解成功！\n');
    
    % 提取基准模型结果
    I_benchmark_opt = value(I_benchmark);
    m_feed_Al2O3_benchmark_opt = value(m_feed_Al2O3_benchmark);
    T_B_benchmark_opt = value(T_B_benchmark);
    
    % 计算基准模型成本
    elec_cost_benchmark = value(Cost_elec_benchmark);
    Al2O3_cost_benchmark = value(Cost_Al2O3_benchmark);
    Al_revenue_benchmark = value(Revenue_Al_benchmark);
    total_cost_benchmark = value(obj_benchmark);
    
    fprintf('基准模型总运行成本：%.2f元\n', total_cost_benchmark);
    
else
    fprintf('基准模型求解失败，错误代码：%d\n', result_benchmark.problem);
    return;
end

%% 第二步：直接求解当前模型
% 这一步的目的是获得当前模型的全局最优解
% 电流和氧化铝进料都可以自由优化，不受基准模型约束
fprintf('\n=== 直接求解当前模型 ===\n');

% 定义当前模型的决策变量（电流和氧化铝进料都可优化）
I_current = sdpvar(1, T, 'full');        % 电流 (kA)
m_feed_Al2O3_current = sdpvar(1, T, 'full');   % 氧化铝进料速率 (kg/min)

% 状态变量（每小时更新一次）
T_B_current = sdpvar(1, T+1, 'full');          % 浴液温度 (°C)
L_ledge_current = sdpvar(1, T+1, 'full');      % ledge厚度 (m)
m_B_current = sdpvar(1, T+1, 'full');          % 浴液质量 (kg)
C_Al2O3_current = sdpvar(1, T+1, 'full');      % 氧化铝浓度
C_AlF3_current = sdpvar(1, T+1, 'full');       % 氟化铝浓度
C_CaF2_current = sdpvar(1, T+1, 'full');       % 氟化钙浓度
C_LiF_current = sdpvar(1, T+1, 'full');        % 氟化锂浓度

% 辅助变量
T_liq_current = sdpvar(1, T, 'full');          % 液相温度 (°C)（每小时更新一次）
Superheat_current = sdpvar(1, T, 'full');      % 过热度 (°C)（每小时更新一次）

% 电流效率变量（温度相关）
g_current = sdpvar(1, T, 'full');              % 电流效率（每小时更新一次）

% 引入McCormick包络的辅助变量
Z_I_g_current = sdpvar(1, T, 'full');          % Z_I_g(t) = I(t) * g(t)

Z_C_Al2O3_m_current = sdpvar(1, T+1, 'full');  % Z_C_Al2O3_m(t) = C_Al2O3(t) * m_B(t)
Z_C_AlF3_m_current = sdpvar(1, T+1, 'full');   % Z_C_AlF3_m(t) = C_AlF3(t) * m_B(t)
Z_C_CaF2_m_current = sdpvar(1, T+1, 'full');   % Z_C_CaF2_m(t) = C_CaF2(t) * m_B(t)
Z_C_LiF_m_current = sdpvar(1, T+1, 'full');    % Z_C_LiF_m(t) = C_LiF(t) * m_B(t)

% 调节动作变量
lamda_keep_current = binvar(1, T-1, 'full');   % 保持状态
lamda_up_current = binvar(1, T-1, 'full');     % 上调状态
lamda_down_current = binvar(1, T-1, 'full');   % 下调状态

% 约束条件
constraints_current = [];

% 初始条件约束
constraints_current = [constraints_current, I_current(1) == I_N];
constraints_current = [constraints_current, T_B_current(1) == T_B_init];
constraints_current = [constraints_current, L_ledge_current(1) == L_ledge_init];
constraints_current = [constraints_current, m_B_current(1) == m_B_init];
constraints_current = [constraints_current, C_Al2O3_current(1) == C_Al2O3_init];
constraints_current = [constraints_current, C_AlF3_current(1) == C_AlF3_init];
constraints_current = [constraints_current, C_CaF2_current(1) == C_CaF2_init];
constraints_current = [constraints_current, C_LiF_current(1) == C_LiF_init];

% 初始McCormick包络变量值
constraints_current = [constraints_current, Z_C_Al2O3_m_current(1) == C_Al2O3_init * m_B_init];
constraints_current = [constraints_current, Z_C_AlF3_m_current(1) == C_AlF3_init * m_B_init];
constraints_current = [constraints_current, Z_C_CaF2_m_current(1) == C_CaF2_init * m_B_init];
constraints_current = [constraints_current, Z_C_LiF_m_current(1) == C_LiF_init * m_B_init];

% 电流上下限约束
constraints_current = [constraints_current, I_min <= I_current <= I_max];

% 氧化铝进料约束
constraints_current = [constraints_current, m_feed_Al2O3_min <= m_feed_Al2O3_current <= m_feed_Al2O3_max];

% 电流效率范围约束
constraints_current = [constraints_current, g_min <= g_current <= g_max];

% McCormick包络约束
% 电流效率与电流相乘
constraints_current = [constraints_current, Z_I_g_current >= I_min * g_current + g_min * I_current - I_min * g_min];
constraints_current = [constraints_current, Z_I_g_current >= I_max * g_current + g_max * I_current - I_max * g_max];
constraints_current = [constraints_current, Z_I_g_current <= I_min * g_current + g_max * I_current - I_min * g_max];
constraints_current = [constraints_current, Z_I_g_current <= I_max * g_current + g_min * I_current - I_max * g_min];

% 浓度与浴液质量相乘（使用从参数文件读取的范围）
for t = 1:T+1
    % Al2O3浓度与浴液质量的McCormick包络
    constraints_current = [constraints_current, Z_C_Al2O3_m_current(t) >= C_Al2O3_min * m_B_current(t) + m_B_min * C_Al2O3_current(t) - C_Al2O3_min * m_B_min];
    constraints_current = [constraints_current, Z_C_Al2O3_m_current(t) >= C_Al2O3_max * m_B_current(t) + m_B_max * C_Al2O3_current(t) - C_Al2O3_max * m_B_max];
    constraints_current = [constraints_current, Z_C_Al2O3_m_current(t) <= C_Al2O3_min * m_B_current(t) + m_B_max * C_Al2O3_current(t) - C_Al2O3_min * m_B_max];
    constraints_current = [constraints_current, Z_C_Al2O3_m_current(t) <= C_Al2O3_max * m_B_current(t) + m_B_min * C_Al2O3_current(t) - C_Al2O3_max * m_B_min];
    
    % AlF3浓度与浴液质量的McCormick包络
    constraints_current = [constraints_current, Z_C_AlF3_m_current(t) >= C_AlF3_min * m_B_current(t) + m_B_min * C_AlF3_current(t) - C_AlF3_min * m_B_min];
    constraints_current = [constraints_current, Z_C_AlF3_m_current(t) >= C_AlF3_max * m_B_current(t) + m_B_max * C_AlF3_current(t) - C_AlF3_max * m_B_max];
    constraints_current = [constraints_current, Z_C_AlF3_m_current(t) <= C_AlF3_min * m_B_current(t) + m_B_max * C_AlF3_current(t) - C_AlF3_min * m_B_max];
    constraints_current = [constraints_current, Z_C_AlF3_m_current(t) <= C_AlF3_max * m_B_current(t) + m_B_min * C_AlF3_current(t) - C_AlF3_max * m_B_min];
    
    % CaF2浓度与浴液质量的McCormick包络
    constraints_current = [constraints_current, Z_C_CaF2_m_current(t) >= C_CaF2_min * m_B_current(t) + m_B_min * C_CaF2_current(t) - C_CaF2_min * m_B_min];
    constraints_current = [constraints_current, Z_C_CaF2_m_current(t) >= C_CaF2_max * m_B_current(t) + m_B_max * C_CaF2_current(t) - C_CaF2_max * m_B_max];
    constraints_current = [constraints_current, Z_C_CaF2_m_current(t) <= C_CaF2_min * m_B_current(t) + m_B_max * C_CaF2_current(t) - C_CaF2_min * m_B_max];
    constraints_current = [constraints_current, Z_C_CaF2_m_current(t) <= C_CaF2_max * m_B_current(t) + m_B_min * C_CaF2_current(t) - C_CaF2_max * m_B_min];
    
    % LiF浓度与浴液质量的McCormick包络
    constraints_current = [constraints_current, Z_C_LiF_m_current(t) >= C_LiF_min * m_B_current(t) + m_B_min * C_LiF_current(t) - C_LiF_min * m_B_min];
    constraints_current = [constraints_current, Z_C_LiF_m_current(t) >= C_LiF_max * m_B_current(t) + m_B_max * C_LiF_current(t) - C_LiF_max * m_B_max];
    constraints_current = [constraints_current, Z_C_LiF_m_current(t) <= C_LiF_min * m_B_current(t) + m_B_max * C_LiF_current(t) - C_LiF_min * m_B_max];
    constraints_current = [constraints_current, Z_C_LiF_m_current(t) <= C_LiF_max * m_B_current(t) + m_B_min * C_LiF_current(t) - C_LiF_max * m_B_min];
end

% 调节状态唯一性约束
constraints_current = [constraints_current, lamda_keep_current + lamda_up_current + lamda_down_current == ones(1, T-1)];

% 功率计算
P_current = P_N + dP_dI_N*(I_current - repmat(I_N,1,T));

% 功率爬坡约束（基于每小时功率变化）
constraints_current = [constraints_current, P_current(2:end) - P_current(1:end-1) >= -R_down * lamda_down_current * P_N];
constraints_current = [constraints_current, P_current(2:end) - P_current(1:end-1) <= R_up * lamda_up_current * P_N];

% 最大调节次数约束
constraints_current = [constraints_current, sum(lamda_up_current + lamda_down_current) <= N_rel];

% 连续调节约束
constraints_current = [constraints_current, lamda_up_current(2:end) + lamda_down_current(2:end) + ...
    lamda_up_current(1:end-1) + lamda_down_current(1:end-1) <= ones(1, T-2)];

% 动态约束（矩阵形式）
% 液相温度计算（矩阵形式）
constraints_current = [constraints_current, ...
    T_liq_current == beta_linear(1) * C_Al2O3_current(1:T)*100 + ...
         beta_linear(2) * C_AlF3_current(1:T)*100 + ...
         beta_linear(3) * C_CaF2_current(1:T)*100 + ...
         beta_linear(4)];

% 过热度计算（矩阵形式）
constraints_current = [constraints_current, Superheat_current == T_B_current(1:T) - T_liq_current];

% 电流效率计算（温度相关，矩阵形式）
constraints_current = [constraints_current, g_current == coeff_global(1) + ...
                                coeff_global(2)*Superheat_current + ...
                                coeff_global(3)*C_AlF3_current(1:T)*100 + ...
                                coeff_global(4)*C_LiF_current(1:T)*100];

% 状态变量更新约束（使用McCormick包络变量和双线性表达式）
% 浴液温度更新约束（直接使用双线性表达式）
constraints_current = [constraints_current, ...
    m_B_current(1:T) .* c_B .* (T_B_current(2:T+1) - T_B_current(1:T)) + m_M * c_M * (T_B_current(2:T+1) - T_B_current(1:T)) == ...
    3600*delta_t * (P_current(1:T)*1000 - ...
    (T_B_current(1:T) - T_gas)/R_B_gas - ...
    (T_B_current(1:T) - T_liq_current)/R_B_ledge - ...
    (T_B_current(1:T) - T_amb)/R_bottom_total - ...
    0.596*60*1e3*m_feed_Al2O3_current(1:T) - ...
    5.474*3600*1e3*(Z_I_g_current*1000*M_Al)/(F*z))];

% 浴液质量更新约束（使用双线性表达式）
% 计算固定的热阻分量
R_shell_side = 0.25/(47 * A_ledge); % 槽壳侧部导热 (K/W)
R_shell_side_amb = 1/(22 * A_ledge); % 槽壳侧部到环境对流 (K/W)

% 计算每个时刻的ledge导热热阻
R_ledge_t = L_ledge_current(1:T) / (1.5 * A_ledge); % 每个时刻的ledge导热热阻 (K/W)
R_ledge_amb_t = R_ledge_t + R_shell_side + R_shell_side_amb; % 每个时刻的总ledge到环境热阻 (K/W)

% 使用双线性表达式重新整理浴液质量更新约束
coeff_m_B = delta_t*3600 / lf;
constraints_current = [constraints_current, ...
    R_ledge_amb_t .* (m_B_current(2:T+1) - m_B_current(1:T)) == ...
    coeff_m_B * (R_ledge_amb_t / R_B_ledge .* (T_B_current(1:T) - T_liq_current) - (T_liq_current - T_amb))];

% ledge厚度更新约束
coeff_L_ledge = 1 / (rho_ledge * A_ledge);
constraints_current = [constraints_current, ...
    L_ledge_current(2:T+1) == L_ledge_current(1:T) - ...
    coeff_L_ledge * (m_B_current(2:T+1) - m_B_current(1:T))];

% 各组分浓度更新约束（使用McCormick包络变量）
coeff_m_Al2O3_feed = 60*delta_t;
coeff_m_Al2O3_reaction = -3600*delta_t * M_Al2O3 / (F*2*z);

% 氧化铝浓度更新约束（使用McCormick包络变量）
constraints_current = [constraints_current, ...
    Z_C_Al2O3_m_current(2:T+1) == Z_C_Al2O3_m_current(1:T) + ...
    coeff_m_Al2O3_feed * m_feed_Al2O3_current(1:T) + ...
    coeff_m_Al2O3_reaction * (Z_I_g_current*1000)];

% 其他组分浓度更新约束
constraints_current = [constraints_current, Z_C_AlF3_m_current(2:T+1) == C_AlF3_init*m_B_init];
constraints_current = [constraints_current, Z_C_CaF2_m_current(2:T+1) == C_CaF2_init*m_B_init];
constraints_current = [constraints_current, Z_C_LiF_m_current(2:T+1) == C_LiF_init*m_B_init];

% 温度约束
constraints_current = [constraints_current, T_B_min <= T_B_current <= T_B_max];

% 过热度约束
constraints_current = [constraints_current, Superheat_min <= Superheat_current <= Superheat_max];

% 氧化铝浓度约束
constraints_current = [constraints_current, C_Al2O3_min <= C_Al2O3_current <= C_Al2O3_max];

% ledge厚度约束
constraints_current = [constraints_current, L_ledge_min <= L_ledge_current <= L_ledge_max];

% 总产量约束
Y_Al_current = 3600*1e3 * (Z_I_g_current * M_Al) / (F*z); % kg/h
Y_total_current = sum(Y_Al_current); % 24小时总产量 (kg)
constraints_current = [constraints_current, Y_total_min <= Y_total_current <= Y_total_max];

% 目标函数：最小化总成本
% 电费成本
Cost_elec_current = sum(price_TOU' .* P_current .* delta_t);
% 氧化铝成本
Cost_Al2O3_current = sum(price_Al2O3 * m_feed_Al2O3_current * delta_t * 60);
% 铝收入
Revenue_Al_current = sum(Y_Al_current .* price_AL) * delta_t;
% 总成本
obj_current = Cost_elec_current + Cost_Al2O3_current - Revenue_Al_current;

% 求解当前模型
result_current = optimize(constraints_current, obj_current, ops);

if result_current.problem == 0
    fprintf('当前模型直接求解成功！\n');
    
    % 提取当前模型结果
    I_current_opt = value(I_current);
    m_feed_Al2O3_current_opt = value(m_feed_Al2O3_current);
    T_B_current_opt = value(T_B_current);
    L_ledge_current_opt = value(L_ledge_current);
    m_B_current_opt = value(m_B_current);
    C_Al2O3_current_opt = value(C_Al2O3_current);
    C_AlF3_current_opt = value(C_AlF3_current);
    C_CaF2_current_opt = value(C_CaF2_current);
    C_LiF_current_opt = value(C_LiF_current);
    T_liq_current_opt = value(T_liq_current);
    Superheat_current_opt = value(Superheat_current);
    g_current_opt = value(g_current);
    Y_Al_current_opt = value(Y_Al_current);
    
    % 计算成本
    elec_cost_current = value(Cost_elec_current);
    Al2O3_cost_current = value(Cost_Al2O3_current);
    Al_revenue_current = value(Revenue_Al_current);
    total_cost_current = value(obj_current);
    
    fprintf('当前模型总运行成本：%.2f元\n', total_cost_current);
    
else
    fprintf('当前模型直接求解失败，错误代码：%d\n', result_current.problem);
    return;
end

%% 第三步：固定基准模型控制变量在当前模型中计算状态变量
% 这一步的目的是计算当电流和氧化铝进料都固定为基准模型结果时，
% 当前模型（考虑所有状态变量约束）中各状态变量的变化
% 使用仿真模型进行小时级状态变量计算
fprintf('\n=== 固定基准模型控制变量在当前模型中计算状态变量 ===\n');

% 固定控制变量
I_fixed = I_benchmark_opt;  % 固定电流为基准模型结果
m_feed_Al2O3_fixed = m_feed_Al2O3_current_opt;  % 固定氧化铝进料为当前模型结果

% 初始化状态变量（小时级）
T_B_fixed_opt = zeros(1, T+1);
L_ledge_fixed_opt = zeros(1, T+1);
m_B_fixed_opt = zeros(1, T+1);
C_Al2O3_fixed_opt = zeros(1, T+1);
C_AlF3_fixed_opt = zeros(1, T+1);
C_CaF2_fixed_opt = zeros(1, T+1);
C_LiF_fixed_opt = zeros(1, T+1);
T_liq_fixed_opt = zeros(1, T);
Superheat_fixed_opt = zeros(1, T);
g_fixed_opt = zeros(1, T);

% 初始条件
T_B_fixed_opt(1) = T_B_init;
L_ledge_fixed_opt(1) = L_ledge_init;
m_B_fixed_opt(1) = m_B_init;
C_Al2O3_fixed_opt(1) = C_Al2O3_init;
C_AlF3_fixed_opt(1) = C_AlF3_init;
C_CaF2_fixed_opt(1) = C_CaF2_init;
C_LiF_fixed_opt(1) = C_LiF_init;

% 计算固定的热阻分量
R_shell_side = 0.25/(47 * A_ledge); % 槽壳侧部导热 (K/W)
R_shell_side_amb = 1/(22 * A_ledge); % 槽壳侧部到环境对流 (K/W)

% 小时级仿真循环
for t = 1:T
    % 液相温度计算（使用线性模型）
    T_liq_fixed_opt(t) = beta_linear(1) * C_Al2O3_fixed_opt(t)*100 + ...
                         beta_linear(2) * C_AlF3_fixed_opt(t)*100 + ...
                         beta_linear(3) * C_CaF2_fixed_opt(t)*100 + ...
                         beta_linear(4);
    
    % 过热度计算
    Superheat_fixed_opt(t) = T_B_fixed_opt(t) - T_liq_fixed_opt(t);
    
    % 电流效率计算（温度相关）
    g_fixed_opt(t) = coeff_global(1) + ...
                     coeff_global(2)*Superheat_fixed_opt(t) + ...
                     coeff_global(3)*C_AlF3_fixed_opt(t)*100 + ...
                     coeff_global(4)*C_LiF_fixed_opt(t)*100;
    
    % 功率计算（基于固定的电流）
    P_fixed_t = P_N + dP_dI_N*(I_fixed(t) - I_N);
    
    % 浴液温度更新（使用双线性表达式）
    % 计算ledge导热热阻
    R_ledge_t = L_ledge_fixed_opt(t) / (1.5 * A_ledge);
    R_ledge_amb_t = R_ledge_t + R_shell_side + R_shell_side_amb;
    
    % 浴液温度更新
    dT_B = 3600*delta_t / (m_B_fixed_opt(t) * c_B + m_M * c_M) * ...
           (P_fixed_t*1000 - ...
           (T_B_fixed_opt(t) - T_gas)/R_B_gas - ...
           (T_B_fixed_opt(t) - T_liq_fixed_opt(t))/R_B_ledge - ...
           (T_B_fixed_opt(t) - T_amb)/R_bottom_total - ...
           0.596*60*1e3*m_feed_Al2O3_fixed(t) - ...
           5.474*3600*1e3*(I_fixed(t)*1000*g_fixed_opt(t)*M_Al)/(F*z));
    
    T_B_fixed_opt(t+1) = T_B_fixed_opt(t) + dT_B;
    
    % 浴液质量更新（使用双线性表达式）
    coeff_m_B = delta_t*3600 / lf;
    dmB = coeff_m_B * ((T_B_fixed_opt(t) - T_liq_fixed_opt(t)) / R_B_ledge - ...
                       (T_liq_fixed_opt(t) - T_amb) / R_ledge_amb_t);
    m_B_fixed_opt(t+1) = m_B_fixed_opt(t) + dmB;
    
    % ledge厚度更新
    coeff_L_ledge = 1 / (rho_ledge * A_ledge);
    L_ledge_fixed_opt(t+1) = L_ledge_fixed_opt(t) - coeff_L_ledge * dmB;
    
    % 各组分浓度更新
    % 氧化铝质量更新
    dm_Al2O3_feed = 60*delta_t * m_feed_Al2O3_fixed(t);
    dm_Al2O3_reaction = -3600*delta_t * (I_fixed(t)*1000*g_fixed_opt(t)*M_Al2O3) / (F*2*z);
    m_Al2O3_fixed = C_Al2O3_fixed_opt(t) * m_B_fixed_opt(t) + dm_Al2O3_feed + dm_Al2O3_reaction;
    
    % 各组分浓度
    C_Al2O3_fixed_opt(t+1) = m_Al2O3_fixed / m_B_fixed_opt(t+1);
    C_AlF3_fixed_opt(t+1) = C_AlF3_init * m_B_init / m_B_fixed_opt(t+1);
    C_CaF2_fixed_opt(t+1) = C_CaF2_init * m_B_init / m_B_fixed_opt(t+1);
    C_LiF_fixed_opt(t+1) = C_LiF_init * m_B_init / m_B_fixed_opt(t+1);
end

% 计算固定控制变量情况下的成本
P_fixed = P_N + dP_dI_N*(I_fixed - repmat(I_N,1,T));
Cost_elec_fixed = sum(price_TOU' .* P_fixed .* delta_t);
Cost_Al2O3_fixed = sum(price_Al2O3 * m_feed_Al2O3_fixed .* delta_t * 60);
Y_Al_fixed = 3600*1e3 * (I_fixed*1000 .* g_fixed_opt * M_Al) / (F*z); % kg/h
Revenue_Al_fixed = sum(Y_Al_fixed .* price_AL) * delta_t;
total_cost_fixed = Cost_elec_fixed + Cost_Al2O3_fixed - Revenue_Al_fixed;

fprintf('固定控制变量的当前模型总运行成本：%.2f元\n', total_cost_fixed);

%% 第四步：结果对比分析
% 对比两个结果：
% 1. 固定基准模型控制变量在当前模型中计算状态变量：电流和氧化铝进料固定，计算各状态变量变化
% 2. 当前模型（直接求解）：考虑状态变量约束，电流和进料都可优化
fprintf('\n=== 模型对比结果分析 ===\n');

% 计算各状态变量的变化
time_hr = 1:T;
time_state = (1:T) * delta_t;

% 计算氧化铝浓度（转换为百分比）
C_Al2O3_fixed = C_Al2O3_fixed_opt(1:end-1) * 100;  % 固定控制变量的当前模型
C_Al2O3_current = C_Al2O3_current_opt(1:end-1) * 100;  % 当前模型（直接求解）

%% 可视化对比结果
% 第一张图：电流和氧化铝进料对比
figure('Name', '模型对比 - 控制变量', 'Position', [100, 100, 1200, 400]);

% 电流曲线对比
subplot(1,2,1);
stairs(time_hr, I_fixed, 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
stairs(time_hr, I_current_opt, 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('电流 (kA)');
title('电流曲线对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% 氧化铝进料曲线对比
subplot(1,2,2);
stairs(time_hr, m_feed_Al2O3_fixed, 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
stairs(time_hr, m_feed_Al2O3_current_opt, 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('氧化铝进料速率 (kg/min)');
title('氧化铝进料曲线对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% 第二张图：状态变量对比
figure('Name', '模型对比 - 状态变量', 'Position', [100, 100, 1400, 1000]);

% 电解质温度对比
subplot(3,2,1);
plot(time_state, T_B_fixed_opt(1:end-1), 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
plot(time_state, T_B_current_opt(1:end-1), 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('电解质温度 (°C)');
title('电解质温度变化对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% 液相温度对比
subplot(3,2,2);
plot(time_state, T_liq_fixed_opt, 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
plot(time_state, T_liq_current_opt, 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('液相温度 (°C)');
title('液相温度变化对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% 过热度对比
subplot(3,2,3);
plot(time_state, Superheat_fixed_opt, 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
plot(time_state, Superheat_current_opt, 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('过热度 (°C)');
title('过热度变化对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% 氧化铝浓度对比
subplot(3,2,4);
plot(time_state, C_Al2O3_fixed, 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
plot(time_state, C_Al2O3_current, 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('氧化铝浓度 (wt%)');
title('氧化铝浓度变化对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% 浴液质量对比
subplot(3,2,5);
plot(time_state, m_B_fixed_opt(1:end-1), 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
plot(time_state, m_B_current_opt(1:end-1), 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('浴液质量 (kg)');
title('浴液质量变化对比');
legend('show', 'Location', 'best');
grid on;
hold off;

% Ledge厚度对比
subplot(3,2,6);
plot(time_state, L_ledge_fixed_opt(1:end-1), 'r--', 'LineWidth', 2, 'DisplayName', '固定控制变量');
hold on;
plot(time_state, L_ledge_current_opt(1:end-1), 'g:', 'LineWidth', 2, 'DisplayName', '当前模型（直接求解）');
xlabel('时间 (h)');
ylabel('Ledge厚度 (m)');
title('Ledge厚度变化对比');
legend('show', 'Location', 'best');
grid on;
hold off;

%% 数值对比分析
fprintf('\n=== 数值对比分析 ===\n');

% 成本对比
fprintf('成本对比：\n');
fprintf('  固定控制变量总成本：%.2f元\n', total_cost_fixed);
fprintf('  当前模型（直接求解）总成本：%.2f元\n', total_cost_current);
fprintf('  直接求解相比固定控制变量节省：%.2f元 (%.2f%%)\n', ...
    total_cost_fixed - total_cost_current, ...
    (total_cost_fixed - total_cost_current) / total_cost_fixed * 100);

% 温度对比
fprintf('\n电解质温度差异（直接求解 vs 固定控制变量）：\n');
T_B_diff = T_B_current_opt(1:end-1) - T_B_fixed_opt(1:end-1);
fprintf('  最大差异：%.2f°C\n', max(abs(T_B_diff)));
fprintf('  平均差异：%.2f°C\n', mean(abs(T_B_diff)));
fprintf('  标准差：%.2f°C\n', std(T_B_diff));

% 过热度对比
fprintf('\n过热度差异（直接求解 vs 固定控制变量）：\n');
Superheat_diff = Superheat_current_opt - Superheat_fixed_opt;
fprintf('  最大差异：%.2f°C\n', max(abs(Superheat_diff)));
fprintf('  平均差异：%.2f°C\n', mean(abs(Superheat_diff)));
fprintf('  标准差：%.2f°C\n', std(Superheat_diff));

% 氧化铝浓度对比
fprintf('\n氧化铝浓度差异（直接求解 vs 固定控制变量）：\n');
C_Al2O3_diff = C_Al2O3_current - C_Al2O3_fixed;
fprintf('  最大差异：%.4f wt%%\n', max(abs(C_Al2O3_diff)));
fprintf('  平均差异：%.4f wt%%\n', mean(abs(C_Al2O3_diff)));
fprintf('  标准差：%.4f wt%%\n', std(C_Al2O3_diff));

% 浴液质量对比
fprintf('\n浴液质量差异（直接求解 vs 固定控制变量）：\n');
m_B_diff = m_B_current_opt(1:end-1) - m_B_fixed_opt(1:end-1);
fprintf('  最大差异：%.2f kg\n', max(abs(m_B_diff)));
fprintf('  平均差异：%.2f kg\n', mean(abs(m_B_diff)));
fprintf('  标准差：%.2f kg\n', std(m_B_diff));

% Ledge厚度对比
fprintf('\nLedge厚度差异（直接求解 vs 固定控制变量）：\n');
L_ledge_diff = L_ledge_current_opt(1:end-1) - L_ledge_fixed_opt(1:end-1);
fprintf('  最大差异：%.4f m\n', max(abs(L_ledge_diff)));
fprintf('  平均差异：%.4f m\n', mean(abs(L_ledge_diff)));
fprintf('  标准差：%.4f m\n', std(L_ledge_diff));

%% 保存对比结果
save('model_comparison_results.mat', 'I_benchmark_opt', 'm_feed_Al2O3_benchmark_opt', ...
     'T_B_benchmark_opt', 'I_fixed', 'm_feed_Al2O3_fixed', 'T_B_fixed_opt', 'L_ledge_fixed_opt', 'm_B_fixed_opt', ...
     'C_Al2O3_fixed_opt', 'C_AlF3_fixed_opt', 'C_CaF2_fixed_opt', 'C_LiF_fixed_opt', ...
     'T_liq_fixed_opt', 'Superheat_fixed_opt', 'g_fixed_opt', 'Y_Al_fixed', ...
     'C_Al2O3_fixed', 'elec_cost_benchmark', 'Al2O3_cost_benchmark', ...
     'Al_revenue_benchmark', 'total_cost_benchmark', ...
     'I_current_opt', 'm_feed_Al2O3_current_opt', 'T_B_current_opt', 'L_ledge_current_opt', ...
     'm_B_current_opt', 'C_Al2O3_current_opt', 'C_AlF3_current_opt', 'C_CaF2_current_opt', 'C_LiF_current_opt', ...
     'T_liq_current_opt', 'Superheat_current_opt', 'g_current_opt', 'Y_Al_current_opt', ...
     'C_Al2O3_current', 'elec_cost_current', 'Al2O3_cost_current', 'Al_revenue_current', 'total_cost_current', ...
     'total_cost_fixed', 'Cost_elec_fixed', 'Cost_Al2O3_fixed', 'Revenue_Al_fixed');

fprintf('\n对比结果已保存到 model_comparison_results.mat\n');
fprintf('对比实验完成！\n');
