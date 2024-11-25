%% Setting up variables for inverse optimization, this program needs to be separated from constraint construction, so that different constraints can be built repeatedly

%% Variable Setup

% Surrogate model parameters
P_max_i = sdpvar(1, NOFMODELS, 'full'); % Aggregate power upper limit
P_min_i = sdpvar(1, NOFMODELS, 'full'); % Aggregate power lower limit
E_min_i = sdpvar(1, NOFMODELS, 'full'); % Aggregate energy upper limit
E_max_i = sdpvar(1, NOFMODELS, 'full'); % Aggregate energy lower limit

% Original problem variables: Energy (kW), NOFINTERVALS * NOFMODELS * BATCH_SIZE
p_ti = sdpvar(NOFINTERVALS, NOFMODELS, BATCH_SIZE, 'full'); % Hourly energy consumption MWh

% Auxiliary variables
z_p_max_ti = binvar(NOFINTERVALS, NOFMODELS, BATCH_SIZE, 'full'); % Power limit integer variables
z_p_min_ti = binvar(NOFINTERVALS, NOFMODELS, BATCH_SIZE, 'full'); % Power limit integer variables
z_e_max_i = binvar(1, NOFMODELS, BATCH_SIZE, 'full'); % Energy upper limit integer variables
z_e_min_i = binvar(1, NOFMODELS, BATCH_SIZE, 'full'); % Energy lower limit integer variables

% Dual variables
mu_p_min_ti = sdpvar(NOFINTERVALS, NOFMODELS, BATCH_SIZE, 'full');
mu_p_max_ti = sdpvar(NOFINTERVALS, NOFMODELS, BATCH_SIZE, 'full');
mu_e_min_i = sdpvar(1, NOFMODELS, BATCH_SIZE, 'full');
mu_e_max_i = sdpvar(1, NOFMODELS, BATCH_SIZE, 'full');

%% Iterative parameter solving

% Initialization
P_max_i_val = ones(1, NOFMODELS) * max(max(E_primal_days_train)) / NOFMODELS;
P_min_i_val = ones(1, NOFMODELS) * min(min(E_primal_days_train)) / NOFMODELS;
E_max_i_val = zeros(1, NOFMODELS);
E_min_i_val = ones(1, NOFMODELS) * sum(E_primal_days_train(:, 1)) / NOFMODELS;

% Setting variable initial values
assign(P_max_i, P_max_i_val);
assign(P_min_i, P_min_i_val);
assign(E_max_i, E_max_i_val);
assign(E_min_i, E_min_i_val);

% Initialization, variable storage
result.P_max_i = [];
result.P_min_i = [];
result.E_min_i = [];
result.E_max_i = [];
result.J_theta = [];
% 记录每次迭代的计算时间
result.time = [];
