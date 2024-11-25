%% Fitting ALF model parameters on the grid side based on electricity meter data
% Solve for optimal parameters using inverse optimization

%% Constraints
Constraints = [];

%% KKT Conditions, for each batch
% Stationarity
Constraints = [Constraints, repmat(price_e, 1, NOFMODELS, 1) - mu_p_min_ti + mu_p_max_ti ...
    - repmat(mu_e_min_i, NOFINTERVALS, 1, 1) + repmat(mu_e_max_i, NOFINTERVALS, 1, 1) == 0];

% Constraints of the original problem
% Power upper and lower limits (MW). NOFMODELS * NOFINTERVALS
Constraints = [Constraints, p_ti <= repmat(P_max_i, NOFINTERVALS, 1, BATCH_SIZE)];
Constraints = [Constraints, repmat(P_min_i, NOFINTERVALS, 1, BATCH_SIZE) <= p_ti];

% Energy in the last time period is between the maximum and minimum values. NOFMODELS
Constraints = [Constraints, sum(p_ti, 1) <= repmat(E_max_i, 1, 1, BATCH_SIZE)];
Constraints = [Constraints, repmat(E_min_i, 1, 1, BATCH_SIZE) <= sum(p_ti, 1)];

% Dual feasibility
Constraints = [Constraints, mu_e_min_i >= 0];
Constraints = [Constraints, mu_e_max_i >= 0];
Constraints = [Constraints, mu_p_min_ti >= 0];
Constraints = [Constraints, mu_p_max_ti >= 0];

% Complementary slackness (Fortuny-Amat transformation)
M = 1e3;
% P_max
Constraints = [Constraints, -p_ti + repmat(P_max_i, NOFINTERVALS, 1, BATCH_SIZE) <= M * z_p_max_ti];
Constraints = [Constraints, mu_p_max_ti <= M * (1 - z_p_max_ti)];
% P_min
Constraints = [Constraints, p_ti - repmat(P_min_i, NOFINTERVALS, 1, BATCH_SIZE) <= M * z_p_min_ti];
Constraints = [Constraints, mu_p_min_ti <= M * (1 - z_p_min_ti)];
% E_max
Constraints = [Constraints, -sum(p_ti, 1) + repmat(E_max_i, 1, 1, BATCH_SIZE) <= M * z_e_max_i];
Constraints = [Constraints, mu_e_max_i <= M * (1 - z_e_max_i)];
% E_min
Constraints = [Constraints, sum(p_ti, 1) - repmat(E_min_i, 1, 1, BATCH_SIZE) <= M * z_e_min_i];
Constraints = [Constraints, mu_e_min_i <= M * (1 - z_e_min_i)];

%% Prior Assumptions

% Limit variable sizes
% Basic parameter constraints (not necessary to limit non-negativity in the presence of original problem constraints)
Constraints = [Constraints, P_max_i <= M];
Constraints = [Constraints, P_min_i >= 0];
Constraints = [Constraints, E_max_i / NOFINTERVALS <= M];
Constraints = [Constraints, E_min_i >= 0];

% Order constraints, for avoiding multiple solutions
for idx = 1:NOFMODELS - 1
    Constraints = [Constraints, P_max_i(1, idx) <= P_max_i(1, idx + 1)];
end

% Dual variable constraints, for avoiding numerical issues
Constraints = [Constraints, mu_p_min_ti <= M * 0.5];
Constraints = [Constraints, mu_p_max_ti <= M * 0.5];
Constraints = [Constraints, mu_e_min_i <= M * 0.5];
Constraints = [Constraints, mu_e_max_i <= M * 0.5];

%% Inverse Problem Objective Function
% Loss function: fitting degree of historical energy consumption data
% J_theta计算过程说明：
% 1. sum(p_ti, 2)：
%    - p_ti是三维矩阵 (NOFINTERVALS × NOFMODELS × BATCH_SIZE)
%    - 沿着第2维求和（即把所有模型的功率加起来）
%    - 结果是 (NOFINTERVALS × 1 × BATCH_SIZE)
%
% 2. reshape(..., NOFINTERVALS, BATCH_SIZE, 1)：
%    - 将上一步的结果重新整形为 (NOFINTERVALS × BATCH_SIZE × 1)
%    - 这代表了每个时间段的总功率预测值
%
% 3. (...) - e_true：
%    - e_true是实际观测到的用电量数据
%    - 做差得到每个时间点的预测误差
%
% 4. norm(...)^2：
%    - norm()计算的是Frobenius范数（默认情况下）
%    - 相当于将所有误差平方后求和再开根号
%    - 即：J_θ = √∑(预测值-真实值)²
%
% 这是一个典型的最小二乘（least squares）目标函数，用于衡量预测值与真实值之间的整体偏差程度
error = reshape(reshape(sum(p_ti, 2), NOFINTERVALS, BATCH_SIZE, 1) - e_true, ...
    NOFINTERVALS * BATCH_SIZE, 1);
J_theta = sum(error.^2);

%% Solve

objective = J_theta;
ops = sdpsettings('debug', 1, 'solver', 'gurobi', ...
    'verbose', 0, ...
    'gurobi.NonConvex', 2, ...
    'allownonconvex', 1, ...
    'gurobi.TimeLimit', TimeLimit, 'usex0', 1); % 'gurobi.NonConvex', 2, ..., 'usex0', 1
ops.gurobi.TuneTimeLimit = TimeLimit;
sol = optimize(Constraints, objective, ops);

%% Display the training process
disp("Number of iterations: " + idx_itr)
disp("Training loss function J_theta: " + value(J_theta))
