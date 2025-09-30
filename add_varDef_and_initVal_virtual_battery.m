%% 虚拟储能模型变量定义和初始化
% 基于逆优化框架的虚拟储能模型变量设置
% 将280个电解槽视为一个虚拟储能系统

%% 变量设置

% 虚拟储能模型参数
p_max = sdpvar(1, 1, 'full'); % 功率上限 (MW)
p_min = sdpvar(1, 1, 'full'); % 功率下限 (MW)
e_max = sdpvar(1, 1, 'full'); % 能量上限 (MWh)
e_min = sdpvar(1, 1, 'full'); % 能量下限 (MWh)
theta = sdpvar(1, 1, 'full'); % 状态转移参数 (自放电系数)
w = sdpvar(1, 1, 'full');     % 基础功率消耗 (MW)

% 原始问题变量
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

% 初始化结果存储
result.p_max = [];
result.p_min = [];
result.e_max = [];
result.e_min = [];
result.theta = [];
result.w = [];
result.J_theta = [];
result.time = [];
