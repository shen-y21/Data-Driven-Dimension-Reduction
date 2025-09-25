% 电解铝负荷用能优化模型
% 基于简化电解槽线性仿真模型的负荷优化
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

% 引入McCormick包络的辅助变量
Z_I_g = sdpvar(1, T, 'full');          % Z_I_g(t) = I(t) * g(t)

Z_C_Al2O3_m = sdpvar(1, T+1, 'full');  % Z_C_Al2O3_m(t) = C_Al2O3(t) * m_B(t)
Z_C_AlF3_m = sdpvar(1, T+1, 'full');   % Z_C_AlF3_m(t) = C_AlF3(t) * m_B(t)
Z_C_CaF2_m = sdpvar(1, T+1, 'full');   % Z_C_CaF2_m(t) = C_CaF2(t) * m_B(t)
Z_C_LiF_m = sdpvar(1, T+1, 'full');    % Z_C_LiF_m(t) = C_LiF(t) * m_B(t)

Z_mB_dT = sdpvar(1, T, 'full');        % Z_mB_dT(t) = m_B(t) * (T_B(t+1) - T_B(t))
Z_Lledge_dm = sdpvar(1, T, 'full');    % Z_Lledge_dm(t) = L_ledge(t) * (m_B(t+1) - m_B(t))

Z_Lledge_superheat = sdpvar(1, T, 'full');  % Z_Lledge_superheat(t) = L_ledge(t) * (T_B(t) - T_liq(t))

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

% 初始McCormick包络变量值
constraints = [constraints, Z_C_Al2O3_m(1) == C_Al2O3_init * m_B_init];
constraints = [constraints, Z_C_AlF3_m(1) == C_AlF3_init * m_B_init];
constraints = [constraints, Z_C_CaF2_m(1) == C_CaF2_init * m_B_init];
constraints = [constraints, Z_C_LiF_m(1) == C_LiF_init * m_B_init];

% 电流上下限约束
constraints = [constraints, I_min <= I <= I_max];

% 氧化铝进料约束
constraints = [constraints, m_feed_Al2O3_min <= m_feed_Al2O3 <= m_feed_Al2O3_max];

% 电流效率范围约束（从参数文件读取）
constraints = [constraints, g_min <= g <= g_max];

% McCormick包络约束
% 电流效率与电流相乘
constraints = [constraints, Z_I_g >= I_min * g + g_min * I - I_min * g_min];
constraints = [constraints, Z_I_g >= I_max * g + g_max * I - I_max * g_max];
constraints = [constraints, Z_I_g <= I_min * g + g_max * I - I_min * g_max];
constraints = [constraints, Z_I_g <= I_max * g + g_min * I - I_max * g_min];

% 浓度与浴液质量相乘（使用从参数文件读取的范围）

for t = 1:T+1
    % Al2O3浓度与浴液质量的McCormick包络
    constraints = [constraints, Z_C_Al2O3_m(t) >= C_Al2O3_min * m_B(t) + m_B_min * C_Al2O3(t) - C_Al2O3_min * m_B_min];
    constraints = [constraints, Z_C_Al2O3_m(t) >= C_Al2O3_max * m_B(t) + m_B_max * C_Al2O3(t) - C_Al2O3_max * m_B_max];
    constraints = [constraints, Z_C_Al2O3_m(t) <= C_Al2O3_min * m_B(t) + m_B_max * C_Al2O3(t) - C_Al2O3_min * m_B_max];
    constraints = [constraints, Z_C_Al2O3_m(t) <= C_Al2O3_max * m_B(t) + m_B_min * C_Al2O3(t) - C_Al2O3_max * m_B_min];
    
    % AlF3浓度与浴液质量的McCormick包络
    constraints = [constraints, Z_C_AlF3_m(t) >= C_AlF3_min * m_B(t) + m_B_min * C_AlF3(t) - C_AlF3_min * m_B_min];
    constraints = [constraints, Z_C_AlF3_m(t) >= C_AlF3_max * m_B(t) + m_B_max * C_AlF3(t) - C_AlF3_max * m_B_max];
    constraints = [constraints, Z_C_AlF3_m(t) <= C_AlF3_min * m_B(t) + m_B_max * C_AlF3(t) - C_AlF3_min * m_B_max];
    constraints = [constraints, Z_C_AlF3_m(t) <= C_AlF3_max * m_B(t) + m_B_min * C_AlF3(t) - C_AlF3_max * m_B_min];
    
    % CaF2浓度与浴液质量的McCormick包络
    constraints = [constraints, Z_C_CaF2_m(t) >= C_CaF2_min * m_B(t) + m_B_min * C_CaF2(t) - C_CaF2_min * m_B_min];
    constraints = [constraints, Z_C_CaF2_m(t) >= C_CaF2_max * m_B(t) + m_B_max * C_CaF2(t) - C_CaF2_max * m_B_max];
    constraints = [constraints, Z_C_CaF2_m(t) <= C_CaF2_min * m_B(t) + m_B_max * C_CaF2(t) - C_CaF2_min * m_B_max];
    constraints = [constraints, Z_C_CaF2_m(t) <= C_CaF2_max * m_B(t) + m_B_min * C_CaF2(t) - C_CaF2_max * m_B_min];
    
    % LiF浓度与浴液质量的McCormick包络
    constraints = [constraints, Z_C_LiF_m(t) >= C_LiF_min * m_B(t) + m_B_min * C_LiF(t) - C_LiF_min * m_B_min];
    constraints = [constraints, Z_C_LiF_m(t) >= C_LiF_max * m_B(t) + m_B_max * C_LiF(t) - C_LiF_max * m_B_max];
    constraints = [constraints, Z_C_LiF_m(t) <= C_LiF_min * m_B(t) + m_B_max * C_LiF(t) - C_LiF_min * m_B_max];
    constraints = [constraints, Z_C_LiF_m(t) <= C_LiF_max * m_B(t) + m_B_min * C_LiF(t) - C_LiF_max * m_B_min];
end

% 分段McCormick包络约束
% 创建包含0的分段点：将范围分为负值段和正值段
dT_negative_segments = linspace(dT_T_min, 0, ceil(N_T_segments/2) + 1);  % 负值段
dT_positive_segments = linspace(0, dT_T_max, floor(N_T_segments/2) + 1);  % 正值段
dT_T_segments = unique([dT_negative_segments, dT_positive_segments]);  % 合并并去重

% 创建包含0的分段点：将范围分为负值段和正值段
dm_negative_segments = linspace(dm_m_min, 0, ceil(N_m_segments/2) + 1);  % 负值段
dm_positive_segments = linspace(0, dm_m_max, floor(N_m_segments/2) + 1);  % 正值段
dm_m_segments = unique([dm_negative_segments, dm_positive_segments]);  % 合并并去重

% 更新实际分段数（可能因去重而改变）
N_T_segments_actual = length(dT_T_segments) - 1;
N_m_segments_actual = length(dm_m_segments) - 1;

fprintf('分段区间定义（包含0作为分段点）:\n');
fprintf('  温度差分分段点: [%.2f', dT_T_segments(1));
for i = 2:length(dT_T_segments)
    fprintf(', %.2f', dT_T_segments(i));
end
fprintf('] (%d段)\n', N_T_segments_actual);

fprintf('  质量差分分段点: [%.0f', dm_m_segments(1));
for i = 2:length(dm_m_segments)
    fprintf(', %.0f', dm_m_segments(i));
end
fprintf('] (%d段)\n', N_m_segments_actual);

% 重新定义二进制变量，使用实际分段数
lambda_T = binvar(T, N_T_segments_actual, 'full');  % 温度差分分段选择变量
lambda_m = binvar(T, N_m_segments_actual, 'full');  % 质量差分分段选择变量

fprintf('二进制变量重新定义: lambda_T(%dx%d), lambda_m(%dx%d)\n', T, N_T_segments_actual, T, N_m_segments_actual);

% 验证分段区间包含0
if ~any(abs(dT_T_segments) < 1e-10)
    warning('温度差分分段区间不包含0，可能影响McCormick包络精度！');
end
if ~any(abs(dm_m_segments) < 1e-10)
    warning('质量差分分段区间不包含0，可能影响McCormick包络精度！');
end

% 定义子变量（用于分段McCormick包络约束）
dT_sub = sdpvar(T, N_T_segments_actual, 'full');
m_B_sub = sdpvar(T, N_T_segments_actual, 'full');
dm_sub = sdpvar(T, N_m_segments_actual, 'full');
L_ledge_sub = sdpvar(T, N_m_segments_actual, 'full');

% 温度差分分段约束（使用变量分解方法）
for t = 1:T
    % 变量分解约束
    constraints = [constraints, T_B(t+1) - T_B(t) == sum(dT_sub(t,:))];
    constraints = [constraints, m_B(t) == sum(m_B_sub(t,:))];
    
    % 分段选择约束
    constraints = [constraints, sum(lambda_T(t, :)) == 1];
    
    for s = 1:N_T_segments_actual
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
end

% 质量差分分段约束（使用变量分解方法）
for t = 1:T
    % 变量分解约束
    constraints = [constraints, m_B(t+1) - m_B(t) == sum(dm_sub(t,:))];
    constraints = [constraints, L_ledge(t) == sum(L_ledge_sub(t,:))];
    
    % 分段选择约束
    constraints = [constraints, sum(lambda_m(t, :)) == 1];
    
    for s = 1:N_m_segments_actual
        dm_s_min = dm_m_segments(s);
        dm_s_max = dm_m_segments(s+1);
        
        % 分段区间约束
        constraints = [constraints, dm_s_min * lambda_m(t,s) <= dm_sub(t,s) <= dm_s_max * lambda_m(t,s)];
        constraints = [constraints, L_ledge_min * lambda_m(t,s) <= L_ledge_sub(t,s) <= L_ledge_max * lambda_m(t,s)];
        
        % 分段内的McCormick包络约束
        est_m_1(s) = L_ledge_min * dm_sub(t,s) + dm_s_min * L_ledge_sub(t,s) - L_ledge_min * dm_s_min * lambda_m(t,s);
        est_m_2(s) = L_ledge_max * dm_sub(t,s) + dm_s_max * L_ledge_sub(t,s) - L_ledge_max * dm_s_max * lambda_m(t,s);
        est_m_3(s) = L_ledge_min * dm_sub(t,s) + dm_s_max * L_ledge_sub(t,s) - L_ledge_min * dm_s_max * lambda_m(t,s);
        est_m_4(s) = L_ledge_max * dm_sub(t,s) + dm_s_min * L_ledge_sub(t,s) - L_ledge_max * dm_s_min * lambda_m(t,s);
    end
    constraints = [constraints, Z_Lledge_dm(t) >= sum(est_m_1)];
    constraints = [constraints, Z_Lledge_dm(t) >= sum(est_m_2)];
    constraints = [constraints, Z_Lledge_dm(t) <= sum(est_m_3)];
    constraints = [constraints, Z_Lledge_dm(t) <= sum(est_m_4)];
end

% Ledge厚度与过热度相乘的McCormick包络约束
for t = 1:T
    superheat_t = T_B(t) - T_liq(t);  % 过热度 = T_B(t) - T_liq(t)
    
    % 全局McCormick包络约束：Z_Lledge_superheat(t) = L_ledge(t) * superheat_t
    constraints = [constraints, Z_Lledge_superheat(t) >= L_ledge_min * superheat_t + Superheat_min * L_ledge(t) - L_ledge_min * Superheat_min];
    constraints = [constraints, Z_Lledge_superheat(t) >= L_ledge_max * superheat_t + Superheat_max * L_ledge(t) - L_ledge_max * Superheat_max];
    constraints = [constraints, Z_Lledge_superheat(t) <= L_ledge_min * superheat_t + Superheat_max * L_ledge(t) - L_ledge_min * Superheat_max];
    constraints = [constraints, Z_Lledge_superheat(t) <= L_ledge_max * superheat_t + Superheat_min * L_ledge(t) - L_ledge_max * Superheat_min];
end

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

% 调节状态唯一性约束
constraints = [constraints, lamda_keep + lamda_up + lamda_down == ones(1, T-1)];

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
% constraints = [constraints, g == g_base + g_temp_coeff * (T_B(1:T) - T_base)];
constraints = [constraints, g == coeff_global(1) + ...
                                coeff_global(2)*Superheat + ...
                                coeff_global(3)*C_AlF3(1:T)*100 + ...
                                coeff_global(4)*C_LiF(1:T)*100];

% 状态变量更新约束（矩阵形式）
% 浴液温度更新约束（使用分段McCormick包络）
constraints = [constraints, ...
    Z_mB_dT * c_B + m_M * c_M * (T_B(2:T+1) - T_B(1:T)) == ...
    3600*delta_t * (P(1:T)*1000 - ...
    (T_B(1:T) - T_gas)/R_B_gas - ...
    (T_B(1:T) - T_liq)/R_B_ledge - ...
    (T_B(1:T) - T_amb)/R_bottom_total - ...
    0.596*60*1e3*m_feed_Al2O3(1:T) - ...
    5.474*3600*1e3*(Z_I_g*1000*M_Al)/(F*z))];

% 浴液质量更新约束（使用分段McCormick包络）
% 计算固定的热阻分量
R_shell_side = 0.25/(47 * A_ledge); % 槽壳侧部导热 (K/W)
R_shell_side_amb = 1/(22 * A_ledge); % 槽壳侧部到环境对流 (K/W)

% 计算每个时刻的ledge导热热阻
R_ledge_t = L_ledge(1:T) / (1.5 * A_ledge); % 每个时刻的ledge导热热阻 (K/W)
R_ledge_amb_t = R_ledge_t + R_shell_side + R_shell_side_amb; % 每个时刻的总ledge到环境热阻 (K/W)

% 使用分段McCormick包络重新整理浴液质量更新约束
coeff_m_B = delta_t*3600 / lf;
constraints = [constraints, ...
    Z_Lledge_dm / (1.5 * A_ledge) + (R_shell_side + R_shell_side_amb) * (m_B(2:T+1) - m_B(1:T)) == ...
    coeff_m_B * (Z_Lledge_superheat / (1.5 * A_ledge * R_B_ledge) + ...
    (R_shell_side + R_shell_side_amb) / R_B_ledge * (T_B(1:T) - T_liq) - (T_liq - T_amb))];

% ledge厚度更新约束
coeff_L_ledge = 1 / (rho_ledge * A_ledge);
constraints = [constraints, ...
    L_ledge(2:T+1) == L_ledge(1:T) - ...
    coeff_L_ledge * (m_B(2:T+1) - m_B(1:T))];

% 各组分浓度更新约束（使用McCormick包络变量）
coeff_m_Al2O3_feed = 60*delta_t;
coeff_m_Al2O3_reaction = -3600*delta_t * M_Al2O3 / (F*2*z);

% 氧化铝浓度更新约束
constraints = [constraints, ...
    Z_C_Al2O3_m(2:T+1) == Z_C_Al2O3_m(1:T) + ...
    coeff_m_Al2O3_feed * m_feed_Al2O3(1:T) + ...
    coeff_m_Al2O3_reaction * (Z_I_g*1000)];

% 其他组分浓度更新约束
constraints = [constraints, Z_C_AlF3_m(2:T+1) == C_AlF3_init*m_B_init];
constraints = [constraints, Z_C_CaF2_m(2:T+1) == C_CaF2_init*m_B_init];
constraints = [constraints, Z_C_LiF_m(2:T+1) == C_LiF_init*m_B_init];

% 温度约束
constraints = [constraints, T_B_min <= T_B <= T_B_max];
% constraints = [constraints, T_B(T) == T_B_init];

% 过热度约束
constraints = [constraints, Superheat_min <= Superheat <= Superheat_max];

% 氧化铝浓度约束
constraints = [constraints, C_Al2O3_min <= C_Al2O3 <= C_Al2O3_max];

% ledge厚度约束
constraints = [constraints, L_ledge_min <= L_ledge <= L_ledge_max];
% constraints = [constraints, L_ledge(T) == L_ledge_init];

%% 目标函数
% 电费成本（功率每小时更新一次）
Cost_elec = sum(P .* price_TOU') * delta_t;

% 氧化铝用料成本
Cost_Al2O3 = sum(m_feed_Al2O3 .* price_Al2O3) * delta_t;

% 产品收益（状态变量和决策变量都是每小时更新一次，直接计算）
Y_Al = 3600*1e3 * (Z_I_g * M_Al) / (F*z); % kg/h
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
    Z_I_g_opt = value(Z_I_g);
    Z_C_Al2O3_m_opt = value(Z_C_Al2O3_m);
    Z_C_AlF3_m_opt = value(Z_C_AlF3_m);
    Z_C_CaF2_m_opt = value(Z_C_CaF2_m);
    Z_C_LiF_m_opt = value(Z_C_LiF_m);
    Z_mB_dT_opt = value(Z_mB_dT);
    Z_Lledge_dm_opt = value(Z_Lledge_dm);
    Z_Lledge_superheat_opt = value(Z_Lledge_superheat);
    lambda_T_opt = value(lambda_T);
    lambda_m_opt = value(lambda_m);
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
    
    % 验证McCormick包络约束
    fprintf('\n=== McCormick包络验证 ===\n');
    
    % 验证电流效率与电流相乘的McCormick包络
    I_g_product = I_opt .* g_opt;
    mccormick_error_I_g = abs(Z_I_g_opt - I_g_product);
    mccormick_rel_error_I_g = mccormick_error_I_g ./ (abs(I_g_product) + eps) * 100; % 相对误差（%）
    fprintf('电流效率与电流相乘McCormick包络最大相对误差: %.4f%%\n', max(mccormick_rel_error_I_g));
    fprintf('电流效率与电流相乘McCormick包络平均相对误差: %.4f%%\n', mean(mccormick_rel_error_I_g));
    
    % 验证浓度与浴液质量相乘的McCormick包络
    C_Al2O3_m_product = C_Al2O3_opt .* m_B_opt;
    C_AlF3_m_product = C_AlF3_opt .* m_B_opt;
    C_CaF2_m_product = C_CaF2_opt .* m_B_opt;
    C_LiF_m_product = C_LiF_opt .* m_B_opt;
    
    mccormick_error_C_Al2O3 = abs(Z_C_Al2O3_m_opt - C_Al2O3_m_product);
    mccormick_error_C_AlF3 = abs(Z_C_AlF3_m_opt - C_AlF3_m_product);
    mccormick_error_C_CaF2 = abs(Z_C_CaF2_m_opt - C_CaF2_m_product);
    mccormick_error_C_LiF = abs(Z_C_LiF_m_opt - C_LiF_m_product);
    
    % 计算相对误差（%）
    mccormick_rel_error_C_Al2O3 = mccormick_error_C_Al2O3 ./ (abs(C_Al2O3_m_product) + eps) * 100;
    mccormick_rel_error_C_AlF3 = mccormick_error_C_AlF3 ./ (abs(C_AlF3_m_product) + eps) * 100;
    mccormick_rel_error_C_CaF2 = mccormick_error_C_CaF2 ./ (abs(C_CaF2_m_product) + eps) * 100;
    mccormick_rel_error_C_LiF = mccormick_error_C_LiF ./ (abs(C_LiF_m_product) + eps) * 100;
    
    fprintf('Al2O3浓度与浴液质量相乘McCormick包络最大相对误差: %.4f%%\n', max(mccormick_rel_error_C_Al2O3));
    fprintf('AlF3浓度与浴液质量相乘McCormick包络最大相对误差: %.4f%%\n', max(mccormick_rel_error_C_AlF3));
    fprintf('CaF2浓度与浴液质量相乘McCormick包络最大相对误差: %.4f%%\n', max(mccormick_rel_error_C_CaF2));
    fprintf('LiF浓度与浴液质量相乘McCormick包络最大相对误差: %.4f%%\n', max(mccormick_rel_error_C_LiF));
    
    % 验证分段McCormick包络
    fprintf('\n=== 分段McCormick包络验证 ===\n');
    
    % 验证浴液质量与温度差分相乘的分段McCormick包络
    mB_dT_product = m_B_opt(1:end-1) .* (T_B_opt(2:end) - T_B_opt(1:end-1));
    segmented_error_mB_dT = abs(Z_mB_dT_opt - mB_dT_product);
    segmented_rel_error_mB_dT = segmented_error_mB_dT ./ (abs(mB_dT_product) + eps) * 100;
    fprintf('浴液质量与温度差分相乘分段McCormick包络最大相对误差: %.4f%%\n', max(segmented_rel_error_mB_dT));
    fprintf('浴液质量与温度差分相乘分段McCormick包络平均相对误差: %.4f%%\n', mean(segmented_rel_error_mB_dT));
    
    % 验证Ledge厚度与浴液质量差分相乘的分段McCormick包络
    Lledge_dm_product = L_ledge_opt(1:end-1) .* (m_B_opt(2:end) - m_B_opt(1:end-1));
    segmented_error_Lledge_dm = abs(Z_Lledge_dm_opt - Lledge_dm_product);
    segmented_rel_error_Lledge_dm = segmented_error_Lledge_dm ./ (abs(Lledge_dm_product) + eps) * 100;
    fprintf('Ledge厚度与浴液质量差分相乘分段McCormick包络最大相对误差: %.4f%%\n', max(segmented_rel_error_Lledge_dm));
    fprintf('Ledge厚度与浴液质量差分相乘分段McCormick包络平均相对误差: %.4f%%\n', mean(segmented_rel_error_Lledge_dm));
    
    % 验证Ledge厚度与过热度相乘的McCormick包络
    Lledge_superheat_product = L_ledge_opt(1:end-1) .* (T_B_opt(1:end-1) - T_liq_opt(1:end));
    mccormick_error_Lledge_superheat = abs(Z_Lledge_superheat_opt - Lledge_superheat_product);
    mccormick_rel_error_Lledge_superheat = mccormick_error_Lledge_superheat ./ (abs(Lledge_superheat_product) + eps) * 100;
    fprintf('Ledge厚度与过热度相乘McCormick包络最大相对误差: %.4f%%\n', max(mccormick_rel_error_Lledge_superheat));
    fprintf('Ledge厚度与过热度相乘McCormick包络平均相对误差: %.4f%%\n', mean(mccormick_rel_error_Lledge_superheat));
    
    % 分析分段选择情况
    fprintf('\n=== 分段选择分析 ===\n');
    fprintf('温度差分分段选择统计:\n');
    for s = 1:N_T_segments_actual
        count = sum(lambda_T_opt(:, s));
        fprintf('  分段%d (%.2f ~ %.2f°C): %d次\n', s, dT_T_segments(s), dT_T_segments(s+1), count);
    end
    
    fprintf('质量差分分段选择统计:\n');
    for s = 1:N_m_segments_actual
        count = sum(lambda_m_opt(:, s));
        fprintf('  分段%d (%.0f ~ %.0f kg): %d次\n', s, dm_m_segments(s), dm_m_segments(s+1), count);
    end

    fprintf('========================\n');
    
    % 保存优化结果
    save('electrolysis_optimization_results.mat', 'I_opt', 'm_feed_Al2O3_opt', 'T_B_opt', ...
         'L_ledge_opt', 'm_B_opt', 'C_Al2O3_opt', 'C_AlF3_opt', 'C_CaF2_opt', 'C_LiF_opt', ...
         'P_input_opt', 'Y_Al_opt', 'T_liq_opt', 'Superheat_opt', 'g_opt', 'Z_I_g_opt', ...
         'Z_C_Al2O3_m_opt', 'Z_C_AlF3_m_opt', 'Z_C_CaF2_m_opt', 'Z_C_LiF_m_opt', ...
         'Z_mB_dT_opt', 'Z_Lledge_dm_opt', 'Z_Lledge_superheat_opt', 'lambda_T_opt', 'lambda_m_opt', ...
         'lamda_up_opt', 'lamda_down_opt', 'lamda_keep_opt', ...
         'elec_cost', 'Al2O3_cost', 'Al_revenue', 'total_cost');
    
    %% 可视化结果
    time_hr = 1:T;
    time_state = (1:T) * delta_t; % 1小时间隔的时间点
    
    % 第一张图：电价、电流、功率、氧化铝进料、电流效率、产量
    figure('Name', '优化结果 - 控制变量和输出变量', 'Position', [100, 100, 1200, 800]);
    
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
    figure('Name', '优化结果 - 状态变量', 'Position', [100, 100, 1400, 1000]);
    
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
    
else
    fprintf('优化求解失败，错误代码：%d\n', result.problem);
end