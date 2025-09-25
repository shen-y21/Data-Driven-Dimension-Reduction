function params = f_initialize_parameters()
% 初始化电解铝负荷用能优化模型的所有参数
% 返回包含所有参数的结构体

params = struct();

%% 加载拟合结果
% 加载电流效率拟合结果
if exist('current_efficiency_fitting_results.mat', 'file')
    load('current_efficiency_fitting_results.mat');
    params.coeff_global = coeff_global;  % 电流效率线性模型系数
    fprintf('已加载电流效率拟合结果\n');
    fprintf('电流效率线性模型: g = %.4f + %.4f*Superheat + %.4f*AlF3 + %.4f*LiF\n', ...
        coeff_global(1), coeff_global(2), coeff_global(3), coeff_global(4));
else
    fprintf('未找到电流效率拟合结果文件，请先运行 current_efficiency_fitting_simple.m\n');
    error('缺少电流效率拟合结果文件');
end

% 加载液相温度计算式系数
if exist('liquidus_temperature_linear.mat', 'file')
    load('liquidus_temperature_linear.mat');
    params.beta_linear = beta_linear;  % 液相温度线性模型系数
    fprintf('已加载液相温度计算式系数\n');
    fprintf('液相温度线性模型: T_liq = %.6f*Al2O3 + %.6f*AlF3 + %.6f*CaF2 + %.3f\n', ...
        beta_linear(1), beta_linear(2), beta_linear(3), beta_linear(4));
else
    fprintf('未找到液相温度系数文件，请先运行 Tliq_data_driven.m\n');
    error('缺少液相温度系数文件');
end

%% 时间参数
params.T = 24;                    % 优化总时长：24小时
params.delta_t = 1;               % 优化时间步长：1小时

%% 边界收缩算法参数
params.max_iterations = 50;       % 最大迭代次数
params.tolerance = 3;             % 收敛容差（McCormick包络误差阈值，%）
params.shrink_factor_dT = 0.05;   % 温度差分收缩因子
params.shrink_factor_dm = 0.05;   % 质量差分收缩因子
params.shrink_factor_ledge = 0.05; % Ledge厚度收缩因子（更保守）
params.shrink_factor_mass = 0.05;  % 浴液质量收缩因子（更保守）

%% 电价和价格参数
Price = [0.6982,0.6982,0.2329,0.2329,0.2329,0.2329,0.6982,0.9312,0.9312,0.6982,0.6982,0.6982, ...
    0.2329,0.2329,0.2329,0.2329,0.6982,0.6982,0.9312,0.9312,0.9312,0.9312,0.9312,0.9312]';
params.price_TOU = Price;
params.price_AL = 19;              % 铝现货价格（元/kg）
params.price_Al2O3 = 2.5;          % 氧化铝价格（元/kg）

%% 电解槽几何参数
params.L = 16;                     % 长度 (m)
params.W = 4;                      % 宽度 (m)
params.H_B = 0.18;                 % 浴液高度 (m)
params.H_M = 0.22;                 % 铝液高度 (m)
params.H_melt = params.H_B + params.H_M;  % 熔融区总高度 (m)
params.ACD = 0.03;                 % 极距 (m)

% 计算各面积
params.A_bottom = params.L * params.W;                     % 底部面积 (m²)
params.A_anode = params.A_bottom * 0.85;                   % 阳极面积 (m²)
params.A_B = 2*(params.L+params.W)*params.H_B;             % 浴液侧面积 (m²)
params.A_M = 2*(params.L+params.W)*params.H_M;             % 铝液侧面积 (m²)
params.A_ledge = 2*(params.L+params.W)*params.H_melt;      % 熔融区侧面积 (m²)

%% 热阻参数
params.R_B_gas = 0.001769;         % 顶部散热热阻 (K/W)
params.R_B_ledge = 0.000076;       % 侧部散热热阻 (K/W)
params.R_bottom_total = 0.003749;  % 底部散热热阻 (K/W)

%% 热容参数
params.rho_B = 2100;               % 浴液密度 (kg/m³)
params.V_B = params.L*params.W*params.H_B-params.A_anode*(params.H_B-params.ACD);        
params.m_B_init = params.rho_B * params.V_B;   % 浴液初始质量 (kg)
params.c_B = 1200;                 % 浴液比热容 (J/kg·K)

params.rho_M = 2300;               % 铝液密度 (kg/m³)
params.V_M = params.L*params.W*params.H_M;        
params.m_M = params.rho_M * params.V_M;        % 铝液质量 (kg)
params.c_M = 880;                  % 铝液比热容 (J/kg·K)

params.rho_ledge = 3000;           % ledge密度 (kg/m³)
params.lf = 310000;                % 潜热 (J/kg)

%% 电化学参数
params.F = 96485;                  % 法拉第常数
params.z = 3;                      % 电子转移数
params.M_Al2O3 = 101.96 / 1000;   % 氧化铝摩尔质量 (kg/mol)
params.M_Al = 26.98 / 1000;       % 铝摩尔质量 (kg/mol)
params.EMF = 2.16;                 % 反电动势 (V)
params.R_cell = 3.58e-6;           % 槽电阻 (Ω)

%% 环境参数
params.T_amb = 30;                 % 环境温度 (°C)
params.T_gas = 50;                 % 顶部气体温度 (°C)

%% 初始条件
params.T_B_init = 965;             % 浴液初始温度 (°C)
params.C_Al2O3_init = 0.03;        % 氧化铝初始浓度
params.C_AlF3_init = 0.11;         % 氟化铝初始浓度
params.C_CaF2_init = 0.06;         % 氟化钙初始浓度
params.C_LiF_init = 0.015;         % 氟化锂初始浓度
params.L_ledge_init = 0.065;       % 初始ledge厚度 (m)

%% 约束参数
params.I_N = 500;                  % 额定电流 (kA)
params.I_min = 0.7 * params.I_N;   % 电流下限 (kA)
params.I_max = 1.1 * params.I_N;   % 电流上限 (kA)
params.g_min = 0.92;
params.g_max = 0.96;

params.T_B_min = 950;              % 浴液温度下限 (°C)
params.T_B_max = 980;              % 浴液温度上限 (°C)
params.Superheat_min = 5;          % 最小过热度 (°C)
params.Superheat_max = 15;         % 最大过热度 (°C)

params.C_Al2O3_min = 0.02;         % 氧化铝浓度下限
params.C_Al2O3_max = 0.04;         % 氧化铝浓度上限
params.L_ledge_min = 0.05;         % ledge厚度下限 (m)
params.L_ledge_max = 0.10;         % ledge厚度上限 (m)

params.R_up = 0.2;                 % 功率上爬坡速率
params.R_down = 0.2;               % 功率下爬坡速率
params.N_rel = 10;                 % 最大调节次数

params.m_feed_Al2O3_min = 3;       % 最小进料速率 (kg/min)
params.m_feed_Al2O3_max = 8;       % 最大进料速率 (kg/min)

%% 分段McCormick包络参数
params.N_T_segments = 4;           % 温度差分分段数
params.N_m_segments = 4;           % 质量差分分段数

% 各组分浓度范围（用于McCormick包络）
params.C_AlF3_min = 0.093;         % AlF3浓度下限
params.C_AlF3_max = 0.122;         % AlF3浓度上限
params.C_CaF2_min = 0.051;         % CaF2浓度下限
params.C_CaF2_max = 0.067;         % CaF2浓度上限
params.C_LiF_min = 0.012;          % LiF浓度下限
params.C_LiF_max = 0.017;          % LiF浓度上限

%% 边界收缩算法边界参数
params.dT_T_min_init = -15;        % 温度差分初始下界
params.dT_T_max_init = 15;         % 温度差分初始上界
params.dm_m_min_init = -500;       % 质量差分初始下界
params.dm_m_max_init = 500;        % 质量差分初始上界
params.m_B_min_init = 5000;        % 浴液质量初始下界
params.m_B_max_init = 9000;        % 浴液质量初始上界

fprintf('参数结构体初始化完成\n');
fprintf('优化时长: %d小时, 时间步长: %d小时\n', params.T, params.delta_t);
fprintf('McCormick包络参数: 温度差分分段数=%d, 质量差分分段数=%d\n', params.N_T_segments, params.N_m_segments);
fprintf('边界收缩参数: 最大迭代%d次, 容差%.0e\n', ...
    params.max_iterations, params.tolerance);
fprintf('收缩因子: 温度差分%.2f, 质量差分%.2f, Ledge厚度%.2f, 浴液质量%.2f\n', ...
    params.shrink_factor_dT, params.shrink_factor_dm, params.shrink_factor_ledge, params.shrink_factor_mass);

end
