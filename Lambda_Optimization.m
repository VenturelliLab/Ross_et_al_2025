clc; clear; close all;
%% Setup

% Define which dataset to use for parameter fitting
fit_data = '20221227_Data.xlsx';

% Import the model equations and define an ODE solver
model  = @model_full;
solver = @ode89;

% Import the experimental data
x0_vec   = table2array(readtable(fit_data, 'Sheet', 'IC'));
x_vec.x1 = table2array(readtable(fit_data, 'Sheet', 'Data_tyrA'));
x_vec.x2 = table2array(readtable(fit_data, 'Sheet', 'Data_pheA'));


% Create a structure storing experimental details
exp_para.time    = 24;                      % Length of a single passage [hr]
exp_para.x0_vec  = x0_vec;                  % Initial conditions
exp_para.num_ini = size(x0_vec, 2);         % Number of initial conditions
exp_para.num_pas = size(x_vec.x1, 1);       % Number of passages


% Define the range and granularity of the regularization parameter set
lambdas = logspace(-6, -1, 100);


% Define the range and size of the LHS parameter space
n_LHS = 100;
bounds = [0.34	1.50	10.40	0.03	19.60	22.23	0.83	7.96	16.87	2.68	23.82	42.94;...
          1.82	39.37	56.12	2.16	64.96	121.16	1.96	53.16	99.19	6.91	40.04	158.90];
LHS_thetas = getLHS(n_LHS, 12, bounds);


% Set some linear inequality constraints
A = [1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0;... % mu1_max < 2
     0, 0, 0, 0, 0, 0,  1, 0, 0, 0, 0, 0];   % mu2_max < 2
b = [2; 2];


% Set some lower bounds for the parameters
LB = zeros(1, 12);


% Allocate memory for the fmincon results
opt_errors = zeros([length(lambdas), n_LHS]);
opt_thetas = zeros([length(lambdas), n_LHS, 12]);


% Set the optimization options for the fmincon algorithm
fmcOpts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'Iter', ...
    'FiniteDifferenceStepSize',1e-5, ...
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations',5000,...
    'MaxIterations', 200);


%% Optimization

% Run fmincon while looping through lambdas and initial parameter sets
for i = 1:length(lambdas)
    obj_fcn = @(theta) objective(theta, x_vec, model, lambdas(i), exp_para, solver);    
    
    for j = 1:n_LHS        
        disp([num2str(i), ': ', num2str(lambdas(i)), ' - ', num2str(j)])

        theta_opt = fmincon(obj_fcn, LHS_thetas(:, j), A, b, [], [], LB, [], [], fmcOpts);
        opt_errors(i, j) = mse(theta_opt, x_vec, model, exp_para, solver, false);
        opt_thetas(i, j, :) = theta_opt;
    end   
end


%% Cleanup

% Find the best performing parameter sets for each lambda
[min_errors_inLH, I] = min(opt_errors, [], 2);
min_thetas_inLH = zeros(12, length(lambdas));

for i = 1:length(lambdas)
    min_thetas_inLH(:, i) = squeeze(opt_thetas(i, I(i), :));
end
