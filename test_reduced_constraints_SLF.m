%% Calculate the optimal energy consumption results on the test set using a simple model (SAL).

data_set_names = ["cement", "steelpowder", "steelmaking", "eal"];
for name_idx = 1:4

    data_set_name = data_set_names(name_idx);
    % Load data: prices and electricity meter data for past time periods
    if strcmp(data_set_name, "eal")
        load("data_set/EAL_Implementation/dataset_" + data_set_name + ".mat");
        % 电解铝数据需要乘以280倍（总共有280个电解槽）
        E_primal_days_cv = E_primal_days_cv * 280;
        E_primal_days_train = E_primal_days_train * 280;
    else
        load("data_set/dataset_" + data_set_name + ".mat");
    end

    % Input: electricity prices, problem parameters (already input in the upper-level function)
    % Output: optimal energy consumption

    % Number of time intervals
    NOFINTERVALS = 24;

    NOFMODELS = 1;

    % Parameters of the trained model
    e_true = E_primal_days_train(:, 1);
    % Problem parameters
    P_max_i_val = max(e_true); % Aggregated power upper limit
    P_min_i_val = min(e_true);
    E_min_i_val = sum(e_true);
    E_max_i_val = max(e_true) * NOFINTERVALS;

    E_reduced_constraints = [];

    for idx_day = 1:10
        % Read price information for the day
        price_e = Price_days_cv(:, idx_day);

        % Original problem
        primal_problem;

        % Record results
        E_reduced_constraints = [E_reduced_constraints, p_val * ones(NOFMODELS)];

    end

    save("results/data_rc_" + data_set_name + "_SALs.mat", "E_reduced_constraints");

end
