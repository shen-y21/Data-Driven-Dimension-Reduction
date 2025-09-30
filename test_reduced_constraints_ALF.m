%% Calculate the optimal energy consumption results on the test set using the identified parameters for the ALF model.

for data_set_name = ["cement", "steelpowder", "steelmaking", "eal"]

    % Load data: prices and electricity meter data for past time periods
    if strcmp(data_set_name, "eal")
        load("data_set/EAL_Implementation/dataset_" + data_set_name + ".mat");
        % 电解铝数据需要乘以280倍（总共有280个电解槽）
        E_primal_days_cv = E_primal_days_cv * 280;
        E_primal_days_train = E_primal_days_train * 280;
    else
        load("data_set/dataset_" + data_set_name + ".mat");
    end

    % Generate bidding results using the fitted model

    % Input: electricity prices, problem parameters (already input in the upper-level function)
    % Output: optimal energy consumption

    % Number of time intervals
    NOFINTERVALS = 24;

    for NOFMODELS = 1:3

        for BATCH_SIZE = 1 : 4

            % Load parameters of the trained model
            load("results/data_" + data_set_name + NOFMODELS + "ALs_" + BATCH_SIZE + "batch.mat", "result");

            % Problem parameters
            P_max_i_val = result.P_max_i(end, :); % Aggregated power upper limit
            P_min_i_val = result.P_min_i(end, :);
            E_min_i_val = result.E_min_i(end, :);
            E_max_i_val = result.E_max_i(end, :);

            E_reduced_constraints = [];

            for idx_day = 1:10
                % Read price information for the day
                price_e = Price_days_cv(:, idx_day);

                % Original problem
                primal_problem;

                % Record results
                E_reduced_constraints = [E_reduced_constraints, p_val * ones(NOFMODELS, 1)];

            end

            save("results/data_rc_" + data_set_name + NOFMODELS + "ALs.mat", "E_reduced_constraints");

        end

    end

end
