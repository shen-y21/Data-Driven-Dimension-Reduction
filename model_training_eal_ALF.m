%% 电解铝负荷降维可行域模型训练程序
% 基于逆优化框架训练电解铝负荷的降维可行域模型
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

% 根据代理模型中可转移负荷数量区分不同模型
for NOFMODELS = 1:3

    for BATCH_SIZE = 1 : 4 % 批次大小，每次迭代使用的天数

        yalmip('clear');
        
        % 变量命名和初始化
        add_varDef_and_initVal;

        % 逆优化迭代
        for idx_itr = 0:max_itr
            % 求解逆问题，随机使用batch_size天的数据
            idx_days = randi([1, data_set.NOFTRAIN], 1, BATCH_SIZE);

            % 读取这些天的电价信息和用电表数据（MW）
            price_e = Price_days_train(:, idx_days);
            e_true = E_primal_days_train(:, idx_days);

            % 求解逆问题更新val
            tic;
            inverse_optimization;
            result.time = [result.time; toc];

            % 学习率，当迭代次数大于40时自适应更新
            alpha = 0.05;
            if idx_itr > 20
                alpha = 0.05 * 1 / (idx_itr - 20)^0.5;
            end

            % 更新代理模型参数和初始值
            P_max_i_val = (1 - alpha) * P_max_i_val + alpha * value(P_max_i);
            assign(P_max_i, P_max_i_val);
            P_min_i_val = (1 - alpha) * P_min_i_val + alpha * value(P_min_i);
            assign(P_min_i, P_min_i_val);
            E_max_i_val = (1 - alpha) * E_max_i_val + alpha * value(E_max_i);
            assign(E_max_i, E_max_i_val);
            E_min_i_val = (1 - alpha) * E_min_i_val + alpha * value(E_min_i);
            assign(E_min_i, E_min_i_val);

            % 记录迭代过程
            if exist('J_theta', 'var')
                result.J_theta = [result.J_theta; value(J_theta)];
            end
            result.P_max_i = [result.P_max_i; P_max_i_val];
            result.P_min_i = [result.P_min_i; P_min_i_val];
            result.E_max_i = [result.E_max_i; E_max_i_val];
            result.E_min_i = [result.E_min_i; E_min_i_val];

            % 收敛准则：过去NOFDAYS期间损失函数的变化
            err = 1e-3;
            if idx_itr > data_set.NOFTRAIN && ...
                    mean(abs(result.J_theta(end - data_set.NOFTRAIN + 1:end))) < err
                break;
            end
        end

        % 保存结果
        save("results/data_" + data_set_name + NOFMODELS + "ALs_" + BATCH_SIZE + "batch.mat", "result");
        
        % 显示训练完成信息
        fprintf('电解铝负荷模型训练完成 - 模型数: %d, 批次大小: %d\n', NOFMODELS, BATCH_SIZE);
        fprintf('最终损失函数值: %.6f\n', result.J_theta(end));
        fprintf('训练时间: %.2f 秒\n', sum(result.time));
        fprintf('迭代次数: %d\n', length(result.J_theta));
        fprintf('----------------------------------------\n');

    end
end

fprintf('所有电解铝负荷降维模型训练完成！\n');
