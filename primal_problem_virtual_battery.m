%% 虚拟储能模型原始问题求解
% 基于虚拟储能模型求解最优功率消耗
% 状态转移方程: e_{t+1} = theta * e_t + p_t * delta_t + w * delta_t

% 输入: 电价 price_e, 虚拟储能模型参数
% 输出: 最优功率消耗 p_val

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
price_AL = 19;          % 铝现货价格（元/kg）

%% 变量定义
% 功率变量
p_t = sdpvar(NOFINTERVALS, BATCH_SIZE, 'full'); % 每小时功率消耗 (MW)
e_t = sdpvar(NOFINTERVALS+1, BATCH_SIZE, 'full'); % 能量状态 (MWh)

% 注意：电流、铝产量和产品收益现在通过功率直接计算，不再作为决策变量

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

% 5. 功率与电流关系约束: P = P_N + (EMF + 2*R_cell*I_N*1000)*(I - I_N)
% 其中 p_t 是 P 的 280 倍 (p_t = P * 280 / 1000, 单位从kW转换为MW)
dP_dI_N = EMF + 2*R_cell*I_N*1000;  % 功率对电流的导数在额定点的值

% 从功率关系推导电流表达式: I = I_N + (p_t * 1000/280 - P_N) / dP_dI_N
% 铝产量表达式: Y_Al = 3600*1e3 * I * g_N * M_Al / (F*z)
% 将电流表达式代入产量公式，得到产量与功率的关系
Y_Al_coeff = 3600*1e3 * g_N * M_Al / (F*z);  % 产量系数
Y_Al_offset = Y_Al_coeff * (I_N - P_N / dP_dI_N);  % 产量偏移项
Y_Al_slope = Y_Al_coeff / dP_dI_N * 1000/280;  % 产量对功率的斜率

% 产品收益系数: Revenue_Al = sum(Y_Al .* price_AL) * delta_t
% 将产量表达式代入收益公式
Revenue_coeff = Y_Al_slope * price_AL * delta_t;  % 收益系数
Revenue_offset = Y_Al_offset * price_AL * delta_t;  % 收益偏移项

%% 目标函数: 最小化用电成本与产品收益的差值
% 修正后的电价系数：考虑产品收益后的净成本系数
% 净成本系数 = 电价 - 产品收益系数
price_e_corrected = price_e - Revenue_coeff;

% 目标函数：最小化修正后的总成本（包含收益偏移项）
objective = sum(sum(price_e_corrected .* p_t)) - sum(Revenue_offset) * BATCH_SIZE;

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

% 通过功率计算电流、产量和收益
I_val = I_N + (p_val * 1000/280 - P_N) / dP_dI_N;
Y_Al_val = Y_Al_offset + Y_Al_slope * p_val;
Revenue_Al_val = Revenue_offset + Revenue_coeff * p_val;

% 计算成本
Cost_elec_val = sum(sum(price_e .* p_val));
Total_Revenue_Al_val = sum(sum(Revenue_Al_val));

% 检查求解状态
if sol.problem == 0
    fprintf('虚拟储能模型求解成功\n');
    fprintf('用电成本: %.2f 元\n', Cost_elec_val);
    fprintf('产品收益: %.2f 元\n', Total_Revenue_Al_val);
    fprintf('净成本: %.2f 元\n', Cost_elec_val - Total_Revenue_Al_val);
    fprintf('总铝产量: %.2f kg\n', sum(sum(Y_Al_val)));
else
    fprintf('虚拟储能模型求解失败，错误代码: %d\n', sol.problem);
    fprintf('错误信息: %s\n', sol.info);
end
