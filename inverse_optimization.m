%% Fitting ALF model parameters on the grid side based on electricity meter data
% Solve for optimal parameters using inverse optimization

%% Constraints
Constraints = [];

%% KKT Conditions, for each batch
% Stationarity
Constraints = [Constraints, repmat(price_e, 1, NOFMODELS, 1) - mu_p_min_ti + mu_p_max_ti ...
    - repmat(mu_e_min_i, NOFINTERVALS, 1, 1) + repmat(mu_e_max_i, NOFINTERVALS, 1, 1) == 0];

% Constraints of the original problem
% Power upper and lower limits (MW). NOFMODELS * NOFINTERVALS * BATCH_SIZE
Constraints = [Constraints, p_ti <= repmat(P_max_i, NOFINTERVALS, 1, BATCH_SIZE)];
Constraints = [Constraints, repmat(P_min_i, NOFINTERVALS, 1, BATCH_SIZE) <= p_ti];

% Energy in the last time period is between the maximum and minimum values. NOFMODELS * BATCH_SIZE
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
% Explanation of J_theta calculation process:
% 1. sum(p_ti, 2):
%    - p_ti is a 3D matrix (NOFINTERVALS × NOFMODELS × BATCH_SIZE)
%    - Sum along dimension 2 (summing power across all models)
%    - Result is (NOFINTERVALS × 1 × BATCH_SIZE)
%
% 2. reshape(..., NOFINTERVALS, BATCH_SIZE, 1):
%    - Reshapes the previous result into (NOFINTERVALS × BATCH_SIZE × 1)
%    - Represents total power prediction for each time interval
%
% 3. (...) - e_true:
%    - e_true is the actual observed electricity consumption data
%    - Subtraction gives prediction error at each time point
%
% 4. norm(...)^2:
%    - norm() calculates the Frobenius norm (by default)
%    - Equivalent to square root of sum of squared errors
%    - i.e.: J_θ = √∑(predicted - actual)²
%
% This is a typical least squares objective function, measuring 
% the overall deviation between predicted and actual values
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
