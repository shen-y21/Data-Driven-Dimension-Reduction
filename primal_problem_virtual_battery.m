%% 虚拟储能模型原始问题求解
% 基于虚拟储能模型求解最优功率消耗
% 状态转移方程: e_{t+1} = theta * e_t + p_t * delta_t + w * delta_t

% 输入: 电价 price_e, 虚拟储能模型参数
% 输出: 最优功率消耗 p_val

%% 变量定义
% 功率变量
p_t = sdpvar(NOFINTERVALS, BATCH_SIZE, 'full'); % 每小时功率消耗 (MW)
e_t = sdpvar(NOFINTERVALS+1, BATCH_SIZE, 'full'); % 能量状态 (MWh)

%% 约束条件
Constraints = [];

% 1. 功率约束: p_t ∈ [p_min, p_max]
Constraints = [Constraints, p_t >= p_min];
Constraints = [Constraints, p_t <= p_max];

% 2. 能量状态约束: e_t ∈ [e_min, e_max]
Constraints = [Constraints, e_t >= e_min];
Constraints = [Constraints, e_t <= e_max];

% 3. 状态转移约束: e_{t+1} = theta * e_t + p_t * delta_t + w * delta_t
for t = 1:NOFINTERVALS
    for n = 1:BATCH_SIZE
        Constraints = [Constraints, e_t(t+1,n) == theta * e_t(t,n) + p_t(t,n) * delta_t + w * delta_t];
    end
end

% 4. 初始状态约束: e_1 = e_initial (假设初始能量状态为0)
% for n = 1:BATCH_SIZE
%     Constraints = [Constraints, e_t(1,n) == 0];
% end

%% 目标函数: 最小化总用电成本
objective = sum(sum(price_e .* p_t));

%% 求解设置
TimeLimit = 120;
ops = sdpsettings('debug', 1, 'solver', 'gurobi', ...
    'verbose', 0, ...
    'gurobi.NonConvex', 2, ...
    'allownonconvex', 1, ...
    'gurobi.TimeLimit', TimeLimit, 'usex0', 1);
ops.gurobi.TuneTimeLimit = TimeLimit;

%% 求解
sol = optimize(Constraints, objective, ops);

%% 记录最优解
p_val = value(p_t);
e_val = value(e_t);

% 检查求解状态
if sol.problem == 0
    fprintf('虚拟储能模型求解成功\n');
else
    fprintf('虚拟储能模型求解失败，错误代码: %d\n', sol.problem);
    fprintf('错误信息: %s\n', sol.info);
end
