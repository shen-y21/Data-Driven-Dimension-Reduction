%% 基于虚拟储能模型的逆优化
% 拟合虚拟储能模型参数，基于电表数据求解最优参数

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
if ~exist('price_AL', 'var')
    price_AL = 19;          % 铝现货价格（元/kg）
end

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
price_e_corrected = price_e - Revenue_coeff;

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
% 损失函数：历史能耗数据的拟合程度
% J_theta = (1/N) * sum_{n=1}^N ||p^(n) - p_n||^2

% 计算预测误差
error = reshape(p_t - e_true, NOFINTERVALS * BATCH_SIZE, 1);
J_theta = sum(error.^2) / (NOFINTERVALS * BATCH_SIZE);

%% 求解
objective = J_theta;
ops = sdpsettings('debug', 1, 'solver', 'gurobi', ...
    'verbose', 0, ...
    'gurobi.NonConvex', 2, ...
    'allownonconvex', 1, ...
    'gurobi.TimeLimit', TimeLimit, 'usex0', 1);
ops.gurobi.TuneTimeLimit = TimeLimit;
sol = optimize(Constraints, objective, ops);

%% 显示训练过程
disp("Number of iterations: " + idx_itr)
disp("Training loss function J_theta: " + value(J_theta))
disp("Virtual battery parameters:")
disp("  p_max: " + value(p_max))
disp("  p_min: " + value(p_min))
disp("  e_max: " + value(e_max))
disp("  e_min: " + value(e_min))
disp("  theta: " + value(theta))
disp("  w: " + value(w))
