% 初始化电解铝负荷用能优化模型的所有参数
% 脚本形式，直接定义所有参数变量

%% 加载拟合结果
% 加载电流效率拟合结果
if exist('current_efficiency_fitting_results.mat', 'file')
    load('current_efficiency_fitting_results.mat');
    % coeff_global 变量已从文件加载
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
    % beta_linear 变量已从文件加载
    fprintf('已加载液相温度计算式系数\n');
    fprintf('液相温度线性模型: T_liq = %.6f*Al2O3 + %.6f*AlF3 + %.6f*CaF2 + %.3f\n', ...
        beta_linear(1), beta_linear(2), beta_linear(3), beta_linear(4));
else
    fprintf('未找到液相温度系数文件，请先运行 Tliq_data_driven.m\n');
    error('缺少液相温度系数文件');
end

%% 时间参数
T = 24;                    % 优化总时长：24小时
delta_t = 1;               % 优化时间步长：1小时

%% 电价和价格参数
% 如果外部没有提供电价数据，使用默认电价
if ~exist('Price', 'var') || isempty(Price)
    Price = [0.6982,0.6982,0.2329,0.2329,0.2329,0.2329,0.6982,0.9312,0.9312,0.6982,0.6982,0.6982, ...
        0.2329,0.2329,0.2329,0.2329,0.6982,0.6982,0.9312,0.9312,0.9312,0.9312,0.9312,0.9312]';
end
price_TOU = Price;
price_AL = 19;              % 铝现货价格（元/kg）
price_Al2O3 = 2.5;          % 氧化铝价格（元/kg）

%% 电解槽几何参数
L = 16;                     % 长度 (m)
W = 4;                      % 宽度 (m)
H_B = 0.18;                 % 浴液高度 (m)
H_M = 0.22;                 % 铝液高度 (m)
H_melt = H_B + H_M;         % 熔融区总高度 (m)
ACD = 0.03;                 % 极距 (m)

% 计算各面积
A_bottom = L * W;                     % 底部面积 (m²)
A_anode = A_bottom * 0.85;            % 阳极面积 (m²)
A_B = 2*(L+W)*H_B;                    % 浴液侧面积 (m²)
A_M = 2*(L+W)*H_M;                    % 铝液侧面积 (m²)
A_ledge = 2*(L+W)*H_melt;             % 熔融区侧面积 (m²)

%% 热阻参数
R_B_gas = 0.001769;         % 顶部散热热阻 (K/W)
R_B_ledge = 0.000076;       % 侧部散热热阻 (K/W)
R_bottom_total = 0.003749;  % 底部散热热阻 (K/W)

%% 热容参数
rho_B = 2100;               % 浴液密度 (kg/m³)
V_B = L*W*H_B-A_anode*(H_B-ACD);        
m_B_init = rho_B * V_B;     % 浴液初始质量 (kg)
c_B = 1200;                 % 浴液比热容 (J/kg·K)

rho_M = 2300;               % 铝液密度 (kg/m³)
V_M = L*W*H_M;        
m_M = rho_M * V_M;          % 铝液质量 (kg)
c_M = 880;                  % 铝液比热容 (J/kg·K)

rho_ledge = 3000;           % ledge密度 (kg/m³)
lf = 310000;                % 潜热 (J/kg)

%% 电化学参数
F = 96485;                  % 法拉第常数
z = 3;                      % 电子转移数
M_Al2O3 = 101.96 / 1000;   % 氧化铝摩尔质量 (kg/mol)
M_Al = 26.98 / 1000;       % 铝摩尔质量 (kg/mol)
EMF = 2.16;                 % 反电动势 (V)
R_cell = 3.58e-6;           % 槽电阻 (Ω)

%% 环境参数
T_amb = 30;                 % 环境温度 (°C)
T_gas = 50;                 % 顶部气体温度 (°C)

%% 初始条件
T_B_init = 965;             % 浴液初始温度 (°C)
C_Al2O3_init = 0.03;        % 氧化铝初始浓度
C_AlF3_init = 0.11;         % 氟化铝初始浓度
C_CaF2_init = 0.06;         % 氟化钙初始浓度
C_LiF_init = 0.015;         % 氟化锂初始浓度
L_ledge_init = 0.065;       % 初始ledge厚度 (m)

%% 约束参数
I_N = 500;                  % 额定电流 (kA)
I_min = 0.7 * I_N;          % 电流下限 (kA)
I_max = 1.1 * I_N;          % 电流上限 (kA)
g_min = 0.92;
g_max = 0.96;

T_B_min = 950;              % 浴液温度下限 (°C)
T_B_max = 980;              % 浴液温度上限 (°C)
Superheat_min = 5;          % 最小过热度 (°C)
Superheat_max = 15;         % 最大过热度 (°C)

C_Al2O3_min = 0.02;        % 氧化铝浓度下限
C_Al2O3_max = 0.04;        % 氧化铝浓度上限
L_ledge_min = 0.05;         % ledge厚度下限 (m)
L_ledge_max = 0.10;         % ledge厚度上限 (m)

R_up = 0.2;                 % 功率上爬坡速率
R_down = 0.2;               % 功率下爬坡速率
N_rel = 10;                 % 最大调节次数

m_feed_Al2O3_min = 4;       % 最小进料速率 (kg/min)
m_feed_Al2O3_max = 6;       % 最大进料速率 (kg/min)

%% 产量约束参数
g_N = 0.94;                 % 额定电流效率
Y_N = 3600*1e3 * (I_N * g_N * M_Al) / (F*z); % 额定产量 (kg/h)
Y_total_min_ratio = 0.9;   % 总产量最小比例（相对于额定产量）
Y_total_max_ratio = 1.1;   % 总产量最大比例（相对于额定产量）
Y_total_min = Y_N * T * Y_total_min_ratio; % 24小时最小总产量 (kg)
Y_total_max = Y_N * T * Y_total_max_ratio; % 24小时最大总产量 (kg)

%% McCormick包络参数
% 分段McCormick包络参数
N_T_segments = 4;           % 温度差分分段数
N_m_segments = 4;           % 质量差分分段数

% 温度差分分段范围
dT_T_min = -10;             % 温度差分最小值 (°C)
dT_T_max = 10;              % 温度差分最大值 (°C)

% 质量差分分段范围
dm_m_min = -300;            % 质量差分最小值 (kg)
dm_m_max = 300;             % 质量差分最大值 (kg)

% 各组分浓度范围（用于McCormick包络）
C_AlF3_min = 0.093;         % AlF3浓度下限
C_AlF3_max = 0.122;         % AlF3浓度上限
C_CaF2_min = 0.051;         % CaF2浓度下限
C_CaF2_max = 0.067;         % CaF2浓度上限
C_LiF_min = 0.012;          % LiF浓度下限
C_LiF_max = 0.017;          % LiF浓度上限

% 浴液质量范围（用于McCormick包络）
m_B_min = 5000;             % 浴液质量下限 (kg)
m_B_max = 9000;             % 浴液质量上限 (kg)

fprintf('参数初始化完成\n');
fprintf('优化时长: %d小时, 时间步长: %d小时\n', T, delta_t);
fprintf('McCormick包络参数: 温度差分分段数=%d, 质量差分分段数=%d\n', N_T_segments, N_m_segments);
