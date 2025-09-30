%% 电解铝负荷降维可行域模型训练程序（虚拟储能模型）
% 基于逆优化框架训练电解铝负荷的虚拟储能降维模型
% 将280个电解槽视为一个虚拟储能系统，捕捉状态变量时段耦合约束
% 数据集包含21天训练集和10天测试集，单个电解槽数据需要乘以280倍

clc; clear;

%% 选择电解铝数据集
data_set_name = "eal";

% 加载数据：电价和电解槽用电数据
load("data_set/EAL_Implementation/dataset_eal.mat");

% 将单个电解槽数据乘以280倍（总共有280个电解槽）
E_primal_days_train = E_primal_days_train * 280;
E_primal_days_cv = E_primal_days_cv * 280;

% 训练集大小
data_set.NOFTRAIN = 21;

% 求解设置
TimeLimit = 120; % 每次迭代的最大求解时间
% 最大迭代次数
max_itr = 40;
% 时间间隔长度
delta_t = 1;
% 时间间隔数量
NOFINTERVALS = 24;

% 虚拟储能模型：NOFMODELS固定为1（单个虚拟储能）
NOFMODELS = 1;

for BATCH_SIZE = 1 : 4 % 批次大小，每次迭代使用的天数

    yalmip('clear');
    
    % 变量命名和初始化（虚拟储能模型）
    add_varDef_and_initVal_virtual_battery;

    % 逆优化迭代
    for idx_itr = 0:max_itr
        % 求解逆问题，随机使用batch_size天的数据
        idx_days = randi([1, data_set.NOFTRAIN], 1, BATCH_SIZE);

        % 读取这些天的电价信息和用电表数据（MW）
        price_e = Price_days_train(:, idx_days);
        e_true = E_primal_days_train(:, idx_days);

        % 求解逆问题更新val
        tic;
        inverse_optimization_virtual_battery;
        result.time = [result.time; toc];

        % 学习率，当迭代次数大于20时自适应更新
        alpha = 0.05;
        if idx_itr > 20
            alpha = 0.05 * 1 / (idx_itr - 20)^0.5;
        end

        % 更新虚拟储能模型参数
        p_max_val = (1 - alpha) * p_max_val + alpha * value(p_max);
        assign(p_max, p_max_val);
        p_min_val = (1 - alpha) * p_min_val + alpha * value(p_min);
        assign(p_min, p_min_val);
        e_max_val = (1 - alpha) * e_max_val + alpha * value(e_max);
        assign(e_max, e_max_val);
        e_min_val = (1 - alpha) * e_min_val + alpha * value(e_min);
        assign(e_min, e_min_val);
        theta_val = (1 - alpha) * theta_val + alpha * value(theta);
        assign(theta, theta_val);
        w_val = (1 - alpha) * w_val + alpha * value(w);
        assign(w, w_val);

        % 记录迭代过程
        if exist('J_theta', 'var')
            result.J_theta = [result.J_theta; value(J_theta)];
        end
        result.p_max = [result.p_max; p_max_val];
        result.p_min = [result.p_min; p_min_val];
        result.e_max = [result.e_max; e_max_val];
        result.e_min = [result.e_min; e_min_val];
        result.theta = [result.theta; theta_val];
        result.w = [result.w; w_val];

        % 收敛准则：过去NOFDAYS期间损失函数的变化
        err = 1e-3;
        if idx_itr > data_set.NOFTRAIN && ...
                mean(abs(result.J_theta(end - data_set.NOFTRAIN + 1:end))) < err
            break;
        end
    end

    % 保存结果
    save("results/data_" + data_set_name + "_virtual_battery_" + BATCH_SIZE + "batch.mat", "result");
    
    % 显示训练完成信息
    fprintf('电解铝虚拟储能模型训练完成 - 批次大小: %d\n', BATCH_SIZE);
    fprintf('最终损失函数值: %.6f\n', result.J_theta(end));
    fprintf('训练时间: %.2f 秒\n', sum(result.time));
    fprintf('迭代次数: %d\n', length(result.J_theta));
    fprintf('最终参数:\n');
    fprintf('  p_max: %.2f MW\n', p_max_val);
    fprintf('  p_min: %.2f MW\n', p_min_val);
    fprintf('  e_max: %.2f MWh\n', e_max_val);
    fprintf('  e_min: %.2f MWh\n', e_min_val);
    fprintf('  theta: %.4f\n', theta_val);
    fprintf('  w: %.2f MW\n', w_val);
    fprintf('----------------------------------------\n');

end
fprintf('所有电解铝虚拟储能降维模型训练完成！\n');
