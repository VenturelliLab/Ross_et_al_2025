function x_pass = passage(model, theta, init, period, n, init_OD, solver)
% model: a function handle to simulate dynamics in each passage
% theta: a vector of model parameters
% init: a 1xn vector of the initial conditions dor the system
% period: an integer describing the batch duration
% n: an integer describing the number of batches to simulate
% init_OD: a number describing the initial OD at each batch
% solver: a function handle for the ode solver 

    if not(size(init, 1) == 1)
        init = init';
    end

    x_pass = zeros(n + 1, length(init));
    x_pass(1, :) = init;

    p_model = @(t, x) model(t, x, theta, zeros(1, 4));


    for i = 1:n
        [~, x] = solver(p_model, [0, period], init);

        OD_com = x(end, 1) + x(end, 2);
        d_factor = init_OD / OD_com;

        init([1, 2]) = x(end, [1, 2]) * d_factor;
        init = max([init; [0, 0, 0, 0, 0]]);

        x_pass(i + 1, :) = max(x(end, :), [0, 0, 0, 0, 0]);
    end
    

end
