function y = objective(theta, data, model, lambda, exp_para, solver)
% theta: a vector of model parameters
% data: a 1x1 struct with k fields, each field is an mxn matrix
%   - k is the number of species
%   - m is the number of passages
%   - n is the number of experiments (including initial conditions)
% model: a function handle to simulate dynamics in each passage
% lambda: regularization parameter (here uses L2 regularization)
% exp_para: a 1x1 struct with four fields
%   - time is the duration of a single batch
%   - x0_vec is the initial conditions for eaxh experiment
%   - num_ini is the number of initial conditions
%   - num_pass is the number of passages
% solver: a function handle for the ode solver 


% Extract the relevant experimental parameters
num_ini = exp_para.num_ini;
num_pas = exp_para.num_pas;
period = exp_para.time;


% Initialize the error variable
error = 0;


% Simulate the experiment and compute the error accumulation
for i = 1:num_ini
    init = exp_para.x0_vec(:, i);
    targ_OD = sum(init(1:2));
    x_mod = passage(model, theta, init, period, num_pas, targ_OD, solver);

    x1_exp = data.x1(:, i);
    x2_exp = data.x2(:, i);

    x1_mod = x_mod(2:end, 1);
    x2_mod = x_mod(2:end, 2);

    error = error...
        + (norm(x1_mod - x1_exp, 2) + norm(x2_mod - x2_exp, 2));
end


% Define the weights for parameter normalization prior to regularization
weights = repmat([0.9713, 14.5848, 15.6201, 2.7055, 27.2983, 32.8543], [1, 2]);


% Compute the full error term with regularization
y = error/num_ini + lambda*norm(theta ./ weights, 2);