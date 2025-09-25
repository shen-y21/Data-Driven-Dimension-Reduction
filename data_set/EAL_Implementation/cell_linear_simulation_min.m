clc; clear; close all

%% 加载拟合结果
% 加载电流效率拟合结果
if exist('current_efficiency_fitting_results.mat', 'file')
    load('current_efficiency_fitting_results.mat');
    % coeff_global 变量已从文件加载
    fprintf('已加载电流效率拟合结果\n');
    fprintf('电流效率线性模型: g = %.4f + %.4f*Superheat + %.4f*AlF3 + %.4f*LiF\n', ...
        coeff_global(1), coeff_global(2), coeff_global(3), coeff_global(4));
else
    fprintf('未找到电流效率拟合结果文件，使用默认系数\n');
    coeff_global = [0.95, 0, 0, 0]; % 默认值
end

% 加载液相温度计算式系数
if exist('liquidus_temperature_linear.mat', 'file')
    load('liquidus_temperature_linear.mat');
    % beta_linear 变量已从文件加载
    fprintf('已加载液相温度计算式系数\n');
    fprintf('液相温度线性模型: T_liq = %.6f*Al2O3 + %.6f*AlF3 + %.6f*CaF2 + %.3f\n', ...
        beta_linear(1), beta_linear(2), beta_linear(3), beta_linear(4));
else
    fprintf('未找到液相温度系数文件，使用默认系数\n');
    beta_linear = [-5.752324, -4.378122, -2.201580, 1031.732]; % 默认值
end

%% 参数定义
% 注意：本模型所有温度变量均使用摄氏度（°C）
% 热阻单位为K/W，但由于温度差在°C和K中数值相同，计算不受影响
T = 5;                     % 总模拟时长：24小时
delta_t = 1;              % 时间步长，分钟
N = T * 60/ delta_t;     % 步数

% 常量参数
F = 96485; z = 3;
M_Al2O3 = 101.96 / 1000;
M_Al = 26.98 / 1000;
M_AlF3 = 83.98 / 1000;
MgF2 = 1.0;              

%% 几何参数
L = 16;       % 长度
W = 4;        % 宽度
H_B = 0.18;    % 浴液高度
H_M = 0.22;    % 铝液高度
H_melt = H_B + H_M; % 熔融区总高度
ACD = 0.03;

L_ledge_init = 0.065;

% 计算各面积
A_bottom = L * W;                     % 底部面积
A_anode = A_bottom * 0.85;            % 阳极面积
A_B = 2*(L+W)*H_B;                    % 浴液侧面积
A_M = 2*(L+W)*H_M;                    % 铝液侧面积
A_ledge = 2*(L+W)*H_melt;             % 熔融区侧面积

%% 简化热阻计算
% 注意：热阻单位为K/W，但温度使用°C
% 由于温度差在°C和K中数值相同，热阻计算不受影响
% 顶部散热热阻 (保持不变)
R_B_anode = 1/(1100 * A_anode);
R_anode = 0.4/(6.6 * A_anode);
R_crust = 0.04/(2.6 * A_bottom);
R_crust_gas = 1/(40 * A_bottom);
R_crust_total = R_crust + R_crust_gas;
R_rod = 2.6/(210 * 0.03);
R_rod_gas = 1/(40 * 1.7);
R_rod_total = R_rod + R_rod_gas;
R_B_gas = R_B_anode + R_anode + 1/(1/R_crust_total+1/R_rod_total);
% 与优化模型保持一致，使用固定值
R_B_gas = 0.001769;

% 底部散热整体热阻：熔融区 → 底部环境
R_B_lining = 1/(1250 * A_bottom);           % 熔融区到内衬对流
R_lining = 0.025/(0.14 * A_bottom); % 内衬导热
R_shell_bottom = 0.25/(31 * A_bottom);     % 槽壳导热
R_shell_bottom_amb = 1/(20 * A_bottom);  % 槽壳到环境对流
R_bottom_total = R_B_lining + R_lining + R_shell_bottom + R_shell_bottom_amb;
R_bottom_total = 0.003749;

% 侧部散热整体热阻：熔融区 → ledge → 侧部环境
R_B_ledge = 1/(600*A_B+1100*A_M);      % 熔融区到ledge
R_B_ledge = 0.000076;
R_ledge = L_ledge_init / (1.5 * A_ledge); % ledge导热
R_shell_side = 0.25/(47 * A_ledge); % 槽壳侧部导热
R_shell_side_amb = 1/(22 * A_ledge); % 槽壳侧部到环境对流
R_ledge_amb = R_ledge + R_shell_side + R_shell_side_amb;
% R_ledge_amb = 0.005879;

fprintf('=== 简化热阻分析 ===\n');
fprintf('底部散热热阻: %.6f K/W\n', R_bottom_total);
fprintf('侧部散热热阻: %.6f K/W\n', R_B_ledge);
fprintf('顶部散热热阻: %.6f K/W\n', R_B_gas);
fprintf('侧部环境热阻: %.6f K/W\n', R_ledge_amb);

%% 热容参数计算
% 熔融区（浴液 + 铝液）
rho_B = 2100;         
V_B = L*W*H_B-A_anode*(H_B-ACD);
m_B_init = rho_B * V_B;
c_B = 1200;           

rho_M = 2300;         
V_M = L*W*H_M;        
m_M = rho_M * V_M;    
c_M = 880;           

m_melt = m_B_init + m_M;   
c_melt = (m_B_init*c_B + m_M*c_M)/m_melt;

% ledge参数
rho_ledge = 3000;   
lf = 310000;

%% 初始条件（所有温度使用摄氏度）
T_amb = 30;  % 环境温度，°C
T_amb = T_amb * ones(N,1);

T_gas = 50;  % 顶部气体温度，°C
T_gas = T_gas * ones(N,1);

T_B_init = 965;  % 浴液初始温度，°C

C_Al2O3_init = 0.03;
C_AlF3_init = 0.11;
C_CaF2_init = 0.06;
C_LiF_init = 0.015;
C_MgF2_init = 0.01;
m_Al2O3_init = C_Al2O3_init*m_B_init;
m_AlF3_init = C_AlF3_init*m_B_init;
m_CaF2_init = C_CaF2_init*m_B_init;
m_LiF_init = C_LiF_init*m_B_init;

%% 加载优化结果作为输入
% 尝试从优化结果文件加载电流和氧化铝进料速率
if exist('electrolysis_optimization_results.mat', 'file')
    load('electrolysis_optimization_results.mat');
    fprintf('已加载优化结果文件\n');
    
    % 检查优化结果的时间长度
    if length(I_opt) >= T
        % 使用优化结果的前T小时数据
        I_opt_hourly = I_opt(1:T);  % kA
        m_feed_Al2O3_opt_hourly = m_feed_Al2O3_opt(1:T);  % kg/min
        
        % 将小时数据扩展为分钟数据（阶梯函数）
        % 每个小时的值在60分钟内保持不变
        I_line = zeros(N, 1);
        m_feed_al2o3 = zeros(N, 1);
        
        for hour = 1:T
            start_min = (hour-1)*60 + 1;
            end_min = hour*60;
            I_line(start_min:end_min) = I_opt_hourly(hour) * 1000;  % 转换为A
            m_feed_al2o3(start_min:end_min) = m_feed_Al2O3_opt_hourly(hour);  % kg/min
        end
        
        fprintf('使用优化结果作为输入：\n');
        fprintf('  电流范围: %.1f - %.1f kA\n', min(I_line)/1000, max(I_line)/1000);
        fprintf('  氧化铝进料范围: %.1f - %.1f kg/min\n', min(m_feed_al2o3), max(m_feed_al2o3));
    else
        fprintf('优化结果时间长度不足，使用默认输入\n');
        m_feed_al2o3 = 5 * ones(N,1);  % 固定加料
        I_line = 500000 * ones(N,1);
        I_line(120+1:end) = 480000;  % 第121分钟起突变为480 kA
    end
else
    fprintf('未找到优化结果文件，使用默认输入\n');
    m_feed_al2o3 = 5 * ones(N,1);  % 固定加料
    I_line = 500000 * ones(N,1);
    I_line(120+1:end) = 480000;  % 第121分钟起突变为480 kA
end

V = 3.95;     % 槽电压，V    
EMF = 2.16;
R_cell = 3.58e-6;

%% 初始化变量
% 状态变量
C_Al2O3 = zeros(N+1,1);    % 氧化铝浓度
C_AlF3 = zeros(N+1,1);     % 氟化铝浓度
C_CaF2 = zeros(N+1,1);     % 氟化钙浓度
C_LiF = zeros(N+1,1);      % 氟化锂浓度
m_Al2O3 = zeros(N+1,1);    % 氧化铝质量
m_B = zeros(N+1,1);        % 浴液质量
T_B = zeros(N+1,1);        % 浴液温度
L_ledge = zeros(N+1,1);    % Ledge厚度

T_liq = zeros(N,1);        % 液相温度
T_liq_real = zeros(N,1);
Superheat = zeros(N,1);    % 过热度
Superheat_real = zeros(N,1);
g = zeros(N,1);            % 电流效率

% 初始值
C_Al2O3(1) = C_Al2O3_init;
C_AlF3(1) = C_AlF3_init;
C_CaF2(1) = C_CaF2_init;
C_LiF(1) = C_LiF_init;
m_Al2O3(1) = m_Al2O3_init;
m_B(1) = m_B_init;
T_B(1) = T_B_init;
L_ledge(1) = L_ledge_init;

% 记录能量量
E_diss = zeros(N,1);
E_proc = zeros(N,1);
E_loss = zeros(N,1);
E_net = zeros(N,1);

% 记录热流量
Q_top = zeros(N,1);
Q_side = zeros(N,1);
Q_bottom = zeros(N,1);
Q_ledge_amb = zeros(N,1);

% 记录氧化铝质量变化
dm_Al2O3_feed = zeros(N,1);
dm_Al2O3_reaction = zeros(N,1);
dm_Al2O3 = zeros(N,1);

%% 主循环
for k = 1:N

    % 液相温度、过热度（°C）

    % 使用质量与浴液质量的比值计算浓度
    AlF3 = C_AlF3(k) * 100; 
    CaF2 = C_CaF2(k) * 100;
    Al2O3 = C_Al2O3(k) * 100;
    LiF = C_LiF(k) * 100;
    MgF2 = C_MgF2_init * 100;

    T_liq_real(k) = 1011 + 0.5*AlF3 - 0.13*AlF3^2.2 ...
        - (3.45*CaF2)/(1+0.0173*CaF2) + 0.124*CaF2*AlF3 ...
        - 0.00542*(CaF2*AlF3)^1.5 ...
        - (7.93*Al2O3)/(1+0.0936*Al2O3-0.0017*Al2O3^2-0.0023*AlF3*Al2O3)...
        - 3.95*MgF2; % 液相温度，°C
    Superheat_real(k) = T_B(k) - T_liq_real(k);

    % 使用从文件加载的系数计算液相温度
    T_liq(k) = beta_linear(1)*Al2O3 + beta_linear(2)*AlF3 + beta_linear(3)*CaF2 + beta_linear(4);
    Superheat(k) = T_B(k) - T_liq(k);
    
    % 计算电流效率
    g(k) = coeff_global(1) + coeff_global(2)*Superheat(k) + coeff_global(3)*AlF3 + coeff_global(4)*LiF;

    % 输入能量 (W)
    % I0 = I_line(1);
    % E_diss0 = V * I0 - R_ext * I0^2;
    % dE_diss_dI = V - 2 * R_ext * I0;
    % E_diss(k) = E_diss0 + dE_diss_dI * (I_line(k) - I0);
    % E_diss(k) = V * I_line(k);
    E_diss(k) = (EMF+R_cell*I_line(k))*I_line(k);

    % 反应能耗 (W)
    E_proc(k) = 0.596 * 60*1e3 * m_feed_al2o3(k) + 5.474 * 3600*1e3 * (I_line(k)*g(k)*M_Al) / (F*z);
    
    % 简化散热损失 (W) - 使用整体热阻
    Q_top(k) = (T_B(k)-T_gas(k)) / R_B_gas;                    % 顶部散热（温度差为°C，热阻为K/W，结果正确）
    Q_side(k) = (T_B(k)-T_liq(k)) / R_B_ledge;                 % 侧部散热（浴液-ledge）
    Q_bottom(k) = (T_B(k)-T_amb(k)) / R_bottom_total;          % 底部散热（整体热阻）
    Q_ledge_amb(k) = (T_liq(k)-T_amb(k))/R_ledge_amb;
    
    E_loss(k) = Q_top(k) + Q_side(k) + Q_bottom(k);

    % 净能量输入 (W)
    E_net(k) = E_diss(k) - E_proc(k) - E_loss(k);

    % 氧化铝质量变化（与优化模型保持一致）
    % 注意：delta_t是分钟，需要转换为小时用于与优化模型一致
    dm_Al2O3_feed(k) = delta_t * m_feed_al2o3(k);  % 分钟 * kg/min = kg
    dm_Al2O3_reaction(k) = -60*delta_t * (I_line(k)*g(k)*M_Al2O3) / (F*2*z);  % 分钟 * kg/min = kg
    dm_Al2O3(k) = dm_Al2O3_feed(k) + dm_Al2O3_reaction(k);
    m_Al2O3(k+1) = m_Al2O3(k) + dm_Al2O3(k);

    % 浴液质量变化（与优化模型保持一致）
    % 基于ledge厚度与过热度的关系
    % 计算ledge导热热阻
    R_ledge_t = L_ledge(k) / (1.5 * A_ledge);
    R_ledge_amb_t = R_ledge_t + R_shell_side + R_shell_side_amb;
    
    % 浴液质量变化公式（与优化模型一致）
    coeff_m_B = delta_t*60 / lf;  % 注意：delta_t是分钟，需要转换为秒
    dmB = coeff_m_B * ((T_B(k) - T_liq(k)) / R_B_ledge - (T_liq(k) - T_amb(k)) / R_ledge_amb_t);
    m_B(k+1) = m_B(k) + dmB;

    % 各组分浓度（由质量计算）
    C_Al2O3(k+1) = m_Al2O3(k+1) / m_B(k+1);
    C_AlF3(k+1) = m_AlF3_init / m_B(k+1);
    C_CaF2(k+1) = m_CaF2_init / m_B(k+1);
    C_LiF(k+1) = m_LiF_init / m_B(k+1);

    % ledge厚度变化
    dLL = - dmB / (rho_ledge * A_ledge);
    L_ledge(k+1) = L_ledge(k) + dLL;

    % 浴液温度变化（与优化模型保持一致）
    % 使用浴液和铝液分别的热容，而不是总熔融区热容
    dT_B = 60*delta_t / (m_B(k) * c_B + m_M * c_M) * (E_diss(k) - E_loss(k) - E_proc(k));
    T_B(k+1) = T_B(k) + dT_B;

end

%% 可视化
time_hr = (0:N) * (delta_t / 60);

% 主要物理量变化图
figure;
subplot(3,2,1); 
plot(time_hr, T_B, 'r', 'LineWidth', 1.5); 
title('浴液温度 T_B (°C)'); ylabel('°C'); grid on;

subplot(3,2,2); 
plot(time_hr(1:N), T_liq, 'b', 'LineWidth', 1.5); hold on;
plot(time_hr(1:N), T_liq_real, 'r--', 'LineWidth', 1.5);
title('液相温度 T_{liq} (°C)'); ylabel('°C'); 
legend('简化模型', '实际模型', 'Location', 'best'); grid on;

subplot(3,2,3); 
plot(time_hr(1:N), Superheat, 'k', 'LineWidth', 1.5); hold on;
plot(time_hr(1:N), Superheat_real, 'g--', 'LineWidth', 1.5);
title('过热度 Superheat (°C)'); ylabel('°C'); 
legend('简化模型', '实际模型', 'Location', 'best'); grid on;

subplot(3,2,4); 
plot(time_hr, C_Al2O3*100, 'b', 'LineWidth', 1.5); hold on;
plot(time_hr, C_AlF3*100, 'r', 'LineWidth', 1.5);
plot(time_hr, C_CaF2*100, 'g', 'LineWidth', 1.5);
plot(time_hr, C_LiF*100, 'm', 'LineWidth', 1.5);
title('各组分浓度 (wt%)'); ylabel('wt%'); 
legend('Al2O3', 'AlF3', 'CaF2', 'LiF', 'Location', 'best'); grid on;

subplot(3,2,5); 
plot(time_hr, m_B, 'g', 'LineWidth', 1.5);
title('浴液质量 m_B (kg)'); ylabel('kg'); grid on;

subplot(3,2,6); 
plot(time_hr, L_ledge*1000, 'b', 'LineWidth', 1.5);
title('Ledge厚度 L_L (mm)'); ylabel('mm'); xlabel('时间 (h)'); grid on;

% 创建一个新的图形来显示氧化铝质量变化分解
figure;
subplot(2,1,1);
plot(time_hr(1:N), dm_Al2O3_feed, 'g-', 'LineWidth', 1.5); hold on;
plot(time_hr(1:N), dm_Al2O3_reaction, 'r-', 'LineWidth', 1.5);
plot(time_hr(1:N), dm_Al2O3, 'k-', 'LineWidth', 2);
title('氧化铝质量变化分解 (每分钟变化)'); ylabel('kg/min'); 
legend('加料贡献', '反应消耗', '质量变化', '总变化', 'Location', 'best'); 
grid on;

subplot(2,1,2);
plot(time_hr(1:N), cumsum(dm_Al2O3_feed), 'g-', 'LineWidth', 1.5); hold on;
plot(time_hr(1:N), cumsum(dm_Al2O3_reaction), 'r-', 'LineWidth', 1.5);
plot(time_hr(1:N), cumsum(dm_Al2O3), 'k-', 'LineWidth', 2);
title('氧化铝质量变化累积贡献'); ylabel('kg'); xlabel('时间 (h)');
legend('加料累积', '反应累积', '质量变化累积', '总变化累积', 'Location', 'best'); 
grid on;

% 电流和电流效率变化图
figure;
subplot(2,1,1);
plot(time_hr(1:N), I_line/1000, 'k', 'LineWidth', 2);
ylabel('电流 (kA)');
title('电流随时间变化');
grid on;

subplot(2,1,2);
plot(time_hr(1:N), g, 'b', 'LineWidth', 2);
ylabel('电流效率');
xlabel('时间 (h)');
title('电流效率随时间变化');
grid on;

% 如果使用了优化结果，添加对比图
if exist('electrolysis_optimization_results.mat', 'file') && exist('I_opt', 'var')
    figure('Name', '优化结果与仿真结果对比', 'Position', [100, 100, 1200, 800]);
    
     % 电流对比
     subplot(2,2,1);
     stairs(time_hr(1:N), I_line/1000, 'b-', 'LineWidth', 2, 'DisplayName', '仿真结果'); hold on;
     if length(I_opt) >= T
         % 优化结果从第1小时开始，需要转换为从0时刻开始
         time_opt = 0:T-1;  % 0, 1, 2, 3, 4小时
         stairs(time_opt, I_opt(1:T), 'r--', 'LineWidth', 2, 'DisplayName', '优化结果');
     end
     ylabel('电流 (kA)');
     title('电流对比');
     legend('show', 'Location', 'best');
     grid on;
     
     % 氧化铝进料对比
     subplot(2,2,2);
     stairs(time_hr(1:N), m_feed_al2o3, 'b-', 'LineWidth', 2, 'DisplayName', '仿真结果'); hold on;
     if exist('m_feed_Al2O3_opt', 'var') && length(m_feed_Al2O3_opt) >= T
         % 优化结果从第1小时开始，需要转换为从0时刻开始
         time_opt = 0:T-1;  % 0, 1, 2, 3, 4小时
         stairs(time_opt, m_feed_Al2O3_opt(1:T), 'r--', 'LineWidth', 2, 'DisplayName', '优化结果');
     end
     ylabel('氧化铝进料速率 (kg/min)');
     title('氧化铝进料对比');
     legend('show', 'Location', 'best');
     grid on;
     
     % 浴液温度对比
     subplot(2,2,3);
     plot(time_hr, T_B, 'b-', 'LineWidth', 2, 'DisplayName', '仿真结果'); hold on;
     if exist('T_B_opt', 'var') && length(T_B_opt) >= T+1
         % 优化结果从第1小时开始，需要转换为从0时刻开始
         time_opt = 0:T;  % 0, 1, 2, 3, 4, 5小时
         plot(time_opt, T_B_opt(1:T+1), 'r--', 'LineWidth', 2, 'DisplayName', '优化结果');
     end
     ylabel('浴液温度 (°C)');
     title('浴液温度对比');
     legend('show', 'Location', 'best');
     grid on;
     
     % 氧化铝浓度对比
     subplot(2,2,4);
     plot(time_hr, C_Al2O3*100, 'b-', 'LineWidth', 2, 'DisplayName', '仿真结果'); hold on;
     if exist('C_Al2O3_opt', 'var') && length(C_Al2O3_opt) >= T+1
         % 优化结果从第1小时开始，需要转换为从0时刻开始
         time_opt = 0:T;  % 0, 1, 2, 3, 4, 5小时
         plot(time_opt, C_Al2O3_opt(1:T+1)*100, 'r--', 'LineWidth', 2, 'DisplayName', '优化结果');
     end
     ylabel('氧化铝浓度 (wt%)');
     xlabel('时间 (h)');
     title('氧化铝浓度对比');
     legend('show', 'Location', 'best');
     grid on;
end

% 能量量曲线图
figure;
plot(time_hr(1:N), E_diss/1e6, 'r-', 'LineWidth', 1.5); hold on;
plot(time_hr(1:N), E_proc/1e6, 'b-', 'LineWidth', 1.5);
plot(time_hr(1:N), E_loss/1e6, 'k-', 'LineWidth', 1.5);
plot(time_hr(1:N), E_net/1e6, 'g-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('能量 (MW)');
title('能量项随时间变化');
legend('E_{diss}','E_{proc}','E_{loss}','E_{net}');
grid on;

% 净能量输入单独图
figure;
plot(time_hr(1:N), E_net/1e6, 'g-', 'LineWidth', 2); hold on;
plot(time_hr(1:N), zeros(size(time_hr(1:N))), 'k--', 'LineWidth', 1);  % 零线
xlabel('时间 (h)');
ylabel('净能量 (MW)');
title('净能量输入随时间变化');
legend('E_{net}', '零线');
grid on;

% 散热分析
figure;
plot(time_hr(1:N), Q_top/1e3, 'r-', 'LineWidth', 1.5); hold on;
plot(time_hr(1:N), Q_side/1e3, 'b-', 'LineWidth', 1.5);
plot(time_hr(1:N), Q_bottom/1e3, 'g-', 'LineWidth', 1.5);
plot(time_hr(1:N), Q_ledge_amb/1e3, 'm-', 'LineWidth', 1.5);
plot(time_hr(1:N), E_loss/1e3, 'k-', 'LineWidth', 2);
xlabel('时间 (h)');
ylabel('散热功率 (kW)');
title('简化模型散热分解');
legend('顶部散热', '侧部散热', '底部散热', '侧壁-环境散热', '总散热', 'Location', 'best');
grid on;

%% 散热统计分析
fprintf('\n=== 简化模型散热分析 ===\n');
fprintf('平均散热功率 (kW):\n');
fprintf('  顶部散热: %.1f ± %.1f\n', mean(Q_top)/1e3, std(Q_top)/1e3);
fprintf('  侧部散热: %.1f ± %.1f\n', mean(Q_side)/1e3, std(Q_side)/1e3);
fprintf('  底部散热: %.1f ± %.1f\n', mean(Q_bottom)/1e3, std(Q_bottom)/1e3);
fprintf('  总散热: %.1f ± %.1f\n', mean(E_loss)/1e3, std(E_loss)/1e3);

% 计算散热百分比
Q_top_percent = Q_top ./ E_loss * 100;
Q_side_percent = Q_side ./ E_loss * 100;
Q_bottom_percent = Q_bottom ./ E_loss * 100;

fprintf('\n平均散热百分比:\n');
fprintf('  顶部散热: %.1f%%\n', mean(Q_top_percent));
fprintf('  侧部散热: %.1f%%\n', mean(Q_side_percent));
fprintf('  底部散热: %.1f%%\n', mean(Q_bottom_percent));

%% 热阻参数标定
fprintf('\n=== 热阻参数标定 ===\n');

% 1. 根据初始时刻E_net=0，计算初始时刻的E_loss
E_net_initial = 0;
E_loss_initial = E_diss(1) - E_proc(1) - E_net_initial;

% 2. 按照48%、43%、9%的比例分配散热
Q_top_target = E_loss_initial * 0.56;
Q_side_target = E_loss_initial * 0.17;
Q_bottom_target = E_loss_initial * 0.27;

% 3. 反推热阻值
% 初始温度（°C）
T_B_initial = T_B(1);
T_gas_initial = T_gas(1);
T_amb_initial = T_amb(1);
T_liq_initial = T_liq(1);

% 反推热阻（温度差为°C，热阻单位为K/W，数值上等价）
R_B_gas_calibrated = (T_B_initial - T_gas_initial) / Q_top_target;
R_B_ledge_calibrated = (T_B_initial - T_liq_initial) / Q_side_target;
R_bottom_total_calibrated = (T_B_initial - T_amb_initial) / Q_bottom_target;

% 4. 标定R_ledge_amb：利用初始时刻dmB=0的条件
% 根据 dmB = -60*delta_t/lf * (Q_side - (T_liq-T_amb)/R_ledge_amb)
% 初始时刻 dmB = 0，所以 Q_side(1) = (T_liq(1)-T_amb(1))/R_ledge_amb
% 因此 R_ledge_amb = (T_liq(1)-T_amb(1))/Q_side(1)
R_ledge_amb_calibrated = (T_liq_initial - T_amb_initial) / Q_side_target;

fprintf('\n标定后的热阻值:\n');
fprintf('  底部散热热阻: %.6f K/W\n', R_bottom_total_calibrated);
fprintf('  侧部散热热阻: %.6f K/W\n', R_B_ledge_calibrated);
fprintf('  顶部散热热阻: %.6f K/W\n', R_B_gas_calibrated);
fprintf('  侧部环境热阻: %.6f K/W\n', R_ledge_amb_calibrated);

%% 优化结果与仿真结果对比分析
if exist('electrolysis_optimization_results.mat', 'file') && exist('I_opt', 'var')
    fprintf('\n=== 优化结果与仿真结果对比分析 ===\n');
    
    % 电流对比
    if length(I_opt) >= T
        I_sim_avg = mean(I_line(1:T*60)/1000);  % 仿真平均电流 (kA)
        I_opt_avg = mean(I_opt(1:T));  % 优化平均电流 (kA)
        I_error = abs(I_sim_avg - I_opt_avg) / I_opt_avg * 100;
        fprintf('电流对比:\n');
        fprintf('  仿真平均电流: %.2f kA\n', I_sim_avg);
        fprintf('  优化平均电流: %.2f kA\n', I_opt_avg);
        fprintf('  相对误差: %.2f%%\n', I_error);
    end
    
    % 氧化铝进料对比
    if exist('m_feed_Al2O3_opt', 'var') && length(m_feed_Al2O3_opt) >= T
        feed_sim_avg = mean(m_feed_al2o3(1:T*60));  % 仿真平均进料 (kg/min)
        feed_opt_avg = mean(m_feed_Al2O3_opt(1:T));  % 优化平均进料 (kg/min)
        feed_error = abs(feed_sim_avg - feed_opt_avg) / feed_opt_avg * 100;
        fprintf('\n氧化铝进料对比:\n');
        fprintf('  仿真平均进料: %.2f kg/min\n', feed_sim_avg);
        fprintf('  优化平均进料: %.2f kg/min\n', feed_opt_avg);
        fprintf('  相对误差: %.2f%%\n', feed_error);
    end
    
    % 浴液温度对比
    if exist('T_B_opt', 'var') && length(T_B_opt) >= T+1
        T_sim_avg = mean(T_B(1:T*60+1));  % 仿真平均温度 (°C)
        T_opt_avg = mean(T_B_opt(1:T+1));  % 优化平均温度 (°C)
        T_error = abs(T_sim_avg - T_opt_avg) / T_opt_avg * 100;
        fprintf('\n浴液温度对比:\n');
        fprintf('  仿真平均温度: %.2f °C\n', T_sim_avg);
        fprintf('  优化平均温度: %.2f °C\n', T_opt_avg);
        fprintf('  相对误差: %.2f%%\n', T_error);
    end
    
    % 氧化铝浓度对比
    if exist('C_Al2O3_opt', 'var') && length(C_Al2O3_opt) >= T+1
        C_sim_avg = mean(C_Al2O3(1:T*60+1)*100);  % 仿真平均浓度 (wt%)
        C_opt_avg = mean(C_Al2O3_opt(1:T+1)*100);  % 优化平均浓度 (wt%)
        C_error = abs(C_sim_avg - C_opt_avg) / C_opt_avg * 100;
        fprintf('\n氧化铝浓度对比:\n');
        fprintf('  仿真平均浓度: %.3f wt%%\n', C_sim_avg);
        fprintf('  优化平均浓度: %.3f wt%%\n', C_opt_avg);
        fprintf('  相对误差: %.2f%%\n', C_error);
    end
    
    fprintf('\n========================\n');
end