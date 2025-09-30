%% 电解铝负荷降维可行域模型直接求解程序（虚拟储能模型）
% 基于逆优化框架直接求解电解铝负荷的虚拟储能降维模型
% 将280个电解槽视为一个虚拟储能系统，捕捉状态变量时段耦合约束
% 使用所有21天训练数据，直接求解最小化平均最小二乘误差

clc; clear;

%% 选择电解铝数据集
data_set_name = "eal";

% 加载数据：电价和电解槽用电数据
load("data_set/EAL_Implementation/dataset_eal.mat");

% 将单个电解槽数据乘以280倍（总共有280个电解槽）
E_primal_days_train = E_primal_days_train * 280;
E_primal_days_cv = E_primal_days_cv * 280;

% 训练集大小
NOFTRAIN = 21;

% 求解设置
TimeLimit = 300; % 增加求解时间限制，因为使用所有数据
% 时间间隔长度
delta_t = 1;
% 时间间隔数量
NOFINTERVALS = 24;
% 批次大小设为所有训练数据
BATCH_SIZE = NOFTRAIN;

fprintf('开始直接求解电解铝虚拟储能逆优化问题...\n');
fprintf('训练集大小: %d天\n', NOFTRAIN);
fprintf('时间间隔: %d小时\n', NOFINTERVALS);
fprintf('求解时间限制: %d秒\n', TimeLimit);

%% 变量定义和初始化
yalmip('clear');

% 虚拟储能模型参数
p_max = sdpvar(1, 1, 'full'); % 功率上限 (MW)
p_min = sdpvar(1, 1, 'full'); % 功率下限 (MW)
e_max = sdpvar(1, 1, 'full'); % 能量上限 (MWh)
e_min = sdpvar(1, 1, 'full'); % 能量下限 (MWh)
theta = sdpvar(1, 1, 'full'); % 状态转移参数 (自放电系数)
w = sdpvar(1, 1, 'full');     % 基础功率消耗 (MW)

% 原始问题变量 - 使用所有训练数据
p_t = sdpvar(NOFINTERVALS, BATCH_SIZE, 'full'); % 每小时功率消耗 (MW)
e_t = sdpvar(NOFINTERVALS+1, BATCH_SIZE, 'full'); % 能量状态 (MWh)

% 辅助变量 (用于KKT条件的互补松弛)
z_p_max_t = binvar(NOFINTERVALS, BATCH_SIZE, 'full'); % 功率上限整数变量
z_p_min_t = binvar(NOFINTERVALS, BATCH_SIZE, 'full'); % 功率下限整数变量
z_e_max_t = binvar(NOFINTERVALS+1, BATCH_SIZE, 'full'); % 能量上限整数变量
z_e_min_t = binvar(NOFINTERVALS+1, BATCH_SIZE, 'full'); % 能量下限整数变量

% 对偶变量
mu_p_min_t = sdpvar(NOFINTERVALS, BATCH_SIZE, 'full');
mu_p_max_t = sdpvar(NOFINTERVALS, BATCH_SIZE, 'full');
mu_e_min_t = sdpvar(NOFINTERVALS+1, BATCH_SIZE, 'full');
mu_e_max_t = sdpvar(NOFINTERVALS+1, BATCH_SIZE, 'full');
lambda_t = sdpvar(NOFINTERVALS, BATCH_SIZE, 'full'); % 状态转移约束的对偶变量

%% 参数初始化
% 基于历史数据初始化参数
p_max_val = max(max(E_primal_days_train)) * 1.1; % 功率上限
p_min_val = min(min(E_primal_days_train)) * 0.9; % 功率下限
e_max_val = sum(E_primal_days_train(:, 1)) * 1.2; % 能量上限
e_min_val = sum(E_primal_days_train(:, 1)) * 0.8; % 能量下限
theta_val = 0.95; % 自放电系数 (假设每小时5%的自放电)
w_val = mean(mean(E_primal_days_train)) * 0.1; % 基础功率消耗

% 设置变量初始值
assign(p_max, p_max_val);
assign(p_min, p_min_val);
assign(e_max, e_max_val);
assign(e_min, e_min_val);
assign(theta, theta_val);
assign(w, w_val);

fprintf('参数初始化完成:\n');
fprintf('  p_max: %.2f MW\n', p_max_val);
fprintf('  p_min: %.2f MW\n', p_min_val);
fprintf('  e_max: %.2f MWh\n', e_max_val);
fprintf('  e_min: %.2f MWh\n', e_min_val);
fprintf('  theta: %.4f\n', theta_val);
fprintf('  w: %.2f MW\n', w_val);

%% 电解铝相关参数定义
% 电化学参数
F = 96485;                  % 法拉第常数
z = 3;                      % 电子转移数
M_Al = 26.98 / 1000;       % 铝摩尔质量 (kg/mol)
g_N = 0.94;                 % 额定电流效率
EMF = 2.16;                 % 反电动势 (V)
R_cell = 3.58e-6;           % 槽电阻 (Ω)
I_N = 500;                  % 额定电流 (kA)
P_N = (EMF + R_cell*I_N*1000)*I_N;  % 额定功率 (kW)

% 价格参数
price_AL = 19;              % 铝现货价格（元/kg）

% 计算收益系数
dP_dI_N = EMF + 2*R_cell*I_N*1000;  % 功率对电流的导数在额定点的值
Y_Al_coeff = 3600*1e3 * g_N * M_Al / (F*z);  % 产量系数
Y_Al_offset = Y_Al_coeff * (I_N - P_N / dP_dI_N);  % 产量偏移项
Y_Al_slope = Y_Al_coeff / dP_dI_N * 1000/280;  % 产量对功率的斜率
Revenue_coeff = Y_Al_slope * price_AL * delta_t;  % 收益系数
Revenue_offset = Y_Al_offset * price_AL * delta_t;  % 收益偏移项

%% 约束条件
Constraints = [];

%% KKT条件

% 1. 平稳性条件 (Stationarity)
% 对功率变量 p_t 的梯度
% 修正后的电价系数：price_e_corrected = price_e - Revenue_coeff
price_e_corrected = Price_days_train - Revenue_coeff;

for t = 1:NOFINTERVALS
    for n = 1:BATCH_SIZE
        % 目标函数对 p_t 的偏导数: price_e_corrected(t,n) + lambda_t(t,n) - lambda_t(t+1,n)*theta
        if t < NOFINTERVALS
            Constraints = [Constraints, price_e_corrected(t,n) + lambda_t(t,n) - lambda_t(t+1,n)*theta - mu_p_min_t(t,n) + mu_p_max_t(t,n) == 0];
        else
            % 最后一个时间步
            Constraints = [Constraints, price_e_corrected(t,n) + lambda_t(t,n) - mu_p_min_t(t,n) + mu_p_max_t(t,n) == 0];
        end
    end
end

% 对能量状态变量 e_t 的梯度
for t = 1:NOFINTERVALS+1
    for n = 1:BATCH_SIZE
        if t == 1
            % 初始状态
            Constraints = [Constraints, -lambda_t(1,n) - mu_e_min_t(t,n) + mu_e_max_t(t,n) == 0];
        elseif t <= NOFINTERVALS
            % 中间状态: e_{t+1} = theta * e_t + p_t * delta_t + w * delta_t
            Constraints = [Constraints, lambda_t(t-1,n) - lambda_t(t,n)*theta - mu_e_min_t(t,n) + mu_e_max_t(t,n) == 0];
        else
            % 最终状态
            Constraints = [Constraints, lambda_t(t-1,n) - mu_e_min_t(t,n) + mu_e_max_t(t,n) == 0];
        end
    end
end

%% 原始问题约束

% 2. 功率约束: p_t ∈ [p_min, p_max]
Constraints = [Constraints, p_t >= p_min];
Constraints = [Constraints, p_t <= p_max];

% 3. 能量状态约束: e_t ∈ [e_min, e_max]
Constraints = [Constraints, e_t >= e_min];
Constraints = [Constraints, e_t <= e_max];

% 4. 状态转移约束: e_{t+1} = theta * e_t + p_t * delta_t + w * delta_t
for t = 1:NOFINTERVALS
    for n = 1:BATCH_SIZE
        Constraints = [Constraints, e_t(t+1,n) == theta * e_t(t,n) + p_t(t,n) * delta_t + w * delta_t];
    end
end

%% 对偶可行性
Constraints = [Constraints, mu_p_min_t >= 0];
Constraints = [Constraints, mu_p_max_t >= 0];
Constraints = [Constraints, mu_e_min_t >= 0];
Constraints = [Constraints, mu_e_max_t >= 0];

%% 互补松弛条件 (Fortuny-Amat变换)
M = 1e3;

% 功率上限约束的互补松弛
Constraints = [Constraints, -p_t + p_max <= M * z_p_max_t];
Constraints = [Constraints, mu_p_max_t <= M * (1 - z_p_max_t)];

% 功率下限约束的互补松弛
Constraints = [Constraints, p_t - p_min <= M * z_p_min_t];
Constraints = [Constraints, mu_p_min_t <= M * (1 - z_p_min_t)];

% 能量上限约束的互补松弛
Constraints = [Constraints, -e_t + e_max <= M * z_e_max_t];
Constraints = [Constraints, mu_e_max_t <= M * (1 - z_e_max_t)];

% 能量下限约束的互补松弛
Constraints = [Constraints, e_t - e_min <= M * z_e_min_t];
Constraints = [Constraints, mu_e_min_t <= M * (1 - z_e_min_t)];

%% 先验假设和参数约束

% 参数合理性约束
Constraints = [Constraints, p_max <= M];
Constraints = [Constraints, p_min >= 0];
Constraints = [Constraints, e_max <= M];
Constraints = [Constraints, e_min >= 0];
Constraints = [Constraints, theta >= 0.8]; % 自放电系数不能太小
Constraints = [Constraints, theta <= 1.0]; % 自放电系数不能大于1
Constraints = [Constraints, w >= 0];

% 对偶变量约束，避免数值问题
Constraints = [Constraints, mu_p_min_t <= M * 0.5];
Constraints = [Constraints, mu_p_max_t <= M * 0.5];
Constraints = [Constraints, mu_e_min_t <= M * 0.5];
Constraints = [Constraints, mu_e_max_t <= M * 0.5];

%% 逆优化目标函数
% 损失函数：所有训练场景的平均最小二乘误差
% J_theta = (1/(N*T)) * sum_{n=1}^N sum_{t=1}^T (p_t(t,n) - e_true(t,n))^2

% 计算预测误差
error = reshape(p_t - E_primal_days_train, NOFINTERVALS * BATCH_SIZE, 1);
J_theta = sum(error.^2) / (NOFINTERVALS * BATCH_SIZE);

fprintf('\n开始求解逆优化问题...\n');
fprintf('变量数量: %d\n', length(depends(J_theta)));
fprintf('约束数量: %d\n', length(Constraints));

%% 求解
objective = J_theta;
ops = sdpsettings('debug', 1, 'solver', 'gurobi', ...
    'verbose', 1, ...
    'gurobi.NonConvex', 2, ...
    'allownonconvex', 1, ...
    'gurobi.TimeLimit', TimeLimit, ...
    'usex0', 1, ...
    'gurobi.MIPGap', 1e-4, ...
    'gurobi.IntFeasTol', 1e-6);
ops.gurobi.TuneTimeLimit = TimeLimit;

tic;
sol = optimize(Constraints, objective, ops);
solve_time = toc;

%% 结果分析
if sol.problem == 0
    fprintf('\n✓ 逆优化问题求解成功！\n');
    fprintf('求解时间: %.2f 秒\n', solve_time);
    
    % 提取最优参数
    p_max_opt = value(p_max);
    p_min_opt = value(p_min);
    e_max_opt = value(e_max);
    e_min_opt = value(e_min);
    theta_opt = value(theta);
    w_opt = value(w);
    J_theta_opt = value(J_theta);
    
    fprintf('\n=== 最优参数 ===\n');
    fprintf('p_max: %.2f MW\n', p_max_opt);
    fprintf('p_min: %.2f MW\n', p_min_opt);
    fprintf('e_max: %.2f MWh\n', e_max_opt);
    fprintf('e_min: %.2f MWh\n', e_min_opt);
    fprintf('theta: %.4f\n', theta_opt);
    fprintf('w: %.2f MW\n', w_opt);
    fprintf('最终损失函数值: %.6f\n', J_theta_opt);
    
    % 计算拟合误差统计
    p_pred = value(p_t);
    error_matrix = p_pred - E_primal_days_train;
    mse = mean(error_matrix(:).^2);
    rmse = sqrt(mse);
    mae = mean(abs(error_matrix(:)));
    
    fprintf('\n=== 拟合误差统计 ===\n');
    fprintf('均方误差 (MSE): %.6f\n', mse);
    fprintf('均方根误差 (RMSE): %.6f MW\n', rmse);
    fprintf('平均绝对误差 (MAE): %.6f MW\n', mae);
    
    % 保存结果
    result.p_max = p_max_opt;
    result.p_min = p_min_opt;
    result.e_max = e_max_opt;
    result.e_min = e_min_opt;
    result.theta = theta_opt;
    result.w = w_opt;
    result.J_theta = J_theta_opt;
    result.solve_time = solve_time;
    result.mse = mse;
    result.rmse = rmse;
    result.mae = mae;
    result.p_pred = p_pred;
    result.error_matrix = error_matrix;
    
    save("results/direct_inverse_optimization_eal.mat", "result");
    fprintf('\n结果已保存到: results/direct_inverse_optimization_eal.mat\n');
    
else
    fprintf('\n✗ 逆优化问题求解失败！\n');
    fprintf('错误代码: %d\n', sol.problem);
    fprintf('错误信息: %s\n', sol.info);
    
    if sol.problem == 1
        fprintf('问题不可行，请检查约束条件\n');
    elseif sol.problem == 2
        fprintf('目标函数无界，请检查目标函数定义\n');
    elseif sol.problem == 3
        fprintf('求解时间超限\n');
    elseif sol.problem == 4
        fprintf('数值问题，请检查模型参数\n');
    end
end

fprintf('\n直接求解逆优化问题完成！\n');
