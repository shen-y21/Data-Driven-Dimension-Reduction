%% 电解铝负荷降维可行域（虚拟储能）模型测试
% 使用训练好的虚拟储能模型参数在测试集上计算最优能耗结果

clc; clear;

%% 选择电解铝数据集
data_set_name = "eal";

% 加载数据：电价和电解槽用电数据
load("data_set/EAL_Implementation/dataset_eal.mat");

% 将单个电解槽数据乘以280倍（总共有280个电解槽）
E_primal_days_cv = E_primal_days_cv * 280;
E_primal_days_train = E_primal_days_train * 280;

% 时间间隔设置
NOFINTERVALS = 24;
delta_t = 1; % 时间间隔长度（小时）

% 虚拟储能模型：NOFMODELS固定为1（单个虚拟储能）
NOFMODELS = 1;

%% 测试不同批次大小的模型
for BATCH_SIZE = 1 : 4

    fprintf('测试虚拟储能模型 - 批次大小: %d\n', BATCH_SIZE);
    
    % 加载训练好的虚拟储能模型参数
    load("results/data_" + data_set_name + "_virtual_battery_" + BATCH_SIZE + "batch.mat", "result");

    % 虚拟储能模型参数
    p_max_val = result.p_max(end); % 功率上限
    p_min_val = result.p_min(end); % 功率下限
    e_max_val = result.e_max(end); % 能量上限
    e_min_val = result.e_min(end); % 能量下限
    theta_val = result.theta(end); % 状态转移参数
    w_val = result.w(end);         % 基础功率消耗

    % 存储测试结果
    E_reduced_constraints = [];
    E_actual = [];
    optimization_times = [];
    optimization_status = [];

    %% 在测试集上测试（10天）
    for idx_day = 1:10
        fprintf('  测试第 %d 天...\n', idx_day);
        
        % 读取电价信息
        price_e = Price_days_cv(:, idx_day);
        
        % 设置虚拟储能模型参数
        p_max = p_max_val;
        p_min = p_min_val;
        e_max = e_max_val;
        e_min = e_min_val;
        theta = theta_val;
        w = w_val;
        
        % 求解虚拟储能模型原始问题
        tic;
        primal_problem_virtual_battery;
        solve_time = toc;
        
        % 记录结果
        E_reduced_constraints = [E_reduced_constraints, p_val];
        E_actual = [E_actual, E_primal_days_cv(:, idx_day)];
        optimization_times = [optimization_times, solve_time];
        
        % 检查求解状态
        if exist('sol', 'var')
            optimization_status = [optimization_status, sol.problem];
        else
            optimization_status = [optimization_status, -1];
        end
        
        % 显示当天结果
        fprintf('    求解时间: %.3f 秒\n', solve_time);
        fprintf('    求解状态: %d\n', optimization_status(end));
        if optimization_status(end) == 0
            fprintf('    预测总能耗: %.2f MWh\n', sum(p_val));
            fprintf('    实际总能耗: %.2f MWh\n', sum(E_primal_days_cv(:, idx_day)));
            fprintf('    预测误差: %.2f%%\n', abs(sum(p_val) - sum(E_primal_days_cv(:, idx_day))) / sum(E_primal_days_cv(:, idx_day)) * 100);
        end
    end

    %% 计算整体性能指标
    successful_days = sum(optimization_status == 0);
    avg_optimization_time = mean(optimization_times);
    
    % 计算预测精度（仅对成功求解的天数）
    if successful_days > 0
        successful_idx = optimization_status == 0;
        prediction_error = abs(E_reduced_constraints(:, successful_idx) - E_actual(:, successful_idx));
        mae = mean(prediction_error(:));
        E_actual_successful = E_actual(:, successful_idx);
        % 避免除零错误，添加小的正数
        mape = mean(prediction_error(:) ./ (E_actual_successful(:) + eps)) * 100;
        rmse = sqrt(mean(prediction_error(:).^2));
    else
        mae = inf;
        mape = inf;
        rmse = inf;
    end

    %% 保存测试结果
    test_results = struct();
    test_results.E_reduced_constraints = E_reduced_constraints;
    test_results.E_actual = E_actual;
    test_results.optimization_times = optimization_times;
    test_results.optimization_status = optimization_status;
    test_results.successful_days = successful_days;
    test_results.avg_optimization_time = avg_optimization_time;
    test_results.mae = mae;
    test_results.mape = mape;
    test_results.rmse = rmse;
    test_results.model_params = struct('p_max', p_max_val, 'p_min', p_min_val, ...
                                     'e_max', e_max_val, 'e_min', e_min_val, ...
                                     'theta', theta_val, 'w', w_val);
    
    save("results/test_virtual_battery_" + data_set_name + "_" + BATCH_SIZE + "batch.mat", "test_results");
    
    %% 显示测试结果摘要
    fprintf('\n=== 虚拟储能模型测试结果摘要 (批次大小: %d) ===\n', BATCH_SIZE);
    fprintf('成功求解天数: %d/%d\n', successful_days, 10);
    fprintf('平均求解时间: %.3f 秒\n', avg_optimization_time);
    if successful_days > 0
        fprintf('平均绝对误差 (MAE): %.2f MW\n', mae);
        fprintf('平均绝对百分比误差 (MAPE): %.2f%%\n', mape);
        fprintf('均方根误差 (RMSE): %.2f MW\n', rmse);
    else
        fprintf('所有测试天数求解失败\n');
    end
    fprintf('模型参数:\n');
    fprintf('  p_max: %.2f MW\n', p_max_val);
    fprintf('  p_min: %.2f MW\n', p_min_val);
    fprintf('  e_max: %.2f MWh\n', e_max_val);
    fprintf('  e_min: %.2f MWh\n', e_min_val);
    fprintf('  theta: %.4f\n', theta_val);
    fprintf('  w: %.2f MW\n', w_val);
    fprintf('==========================================\n\n');

end

fprintf('电解铝虚拟储能模型测试完成！\n');
fprintf('测试结果已保存到 results/ 目录下\n');
