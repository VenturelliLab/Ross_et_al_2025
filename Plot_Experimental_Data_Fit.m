%%
clc; clear; close all;

%% Setup
% Load model parameters and data
theta = load('theta_full.mat').theta;
model = @model_full;
solver = @ode89;

x0_table = readtable('20221227_Data.xlsx', 'Sheet', 'IC');
x1_table = readtable('20221227_Data.xlsx', 'Sheet', 'Data_tyrA');
x2_table = readtable('20221227_Data.xlsx', 'Sheet', 'Data_pheA');

x0_vec = table2array(x0_table);
x1_data = table2array(x1_table);
x2_data = table2array(x2_table);

initial_aminoAcid_fractions = [0, 0.05, 0.1];
num_conditions = numel(initial_aminoAcid_fractions);
num_passages = 10;
time_points = 0:num_passages;

errors = zeros(2 * num_passages, num_conditions); % 2 strains, 10 passages

colors = struct(...
    'data_tyrA', '#00d5ff', ...
    'data_pheA', '#ff4f00', ...
    'model_tyrA', '#0095b3', ...
    'model_pheA', '#b33700' ...
);

%% Simulation and Plotting
figure;
for idx = 1:num_conditions
    X = initial_aminoAcid_fractions(idx);

    % Initial conditions
    x0 = [0.09, 0.01, X * 400, X * 200, 11.10];
    
    % Simulate
    x_sim = passage(model, theta, x0, 24, num_passages, 0.1, solver);

    % Extract model and experimental data
    x1_sim = x_sim(:, 1);
    x2_sim = x_sim(:, 2);
    x1_true = [x0(1); x1_data(:, idx)];
    x2_true = [x0(2); x2_data(:, idx)];

    % Plot
    subplot(1, 3, idx); hold on; box on;
    plot(time_points, x1_true, '-', 'LineWidth', 2, 'Color', colors.data_tyrA);
    plot(time_points, x2_true, '-', 'LineWidth', 2, 'Color', colors.data_pheA);
    plot(time_points, x1_sim, '--', 'LineWidth', 2, 'Color', colors.model_tyrA);
    plot(time_points, x2_sim, '--', 'LineWidth', 2, 'Color', colors.model_pheA);
    
    xlim([-0.5, 10.5]);
    ylim([-0.1, 3.1]);
    xlabel('Passage');
    ylabel('OD');
    title(sprintf('%.f/%.f [Tyr]/[Phe]', x0(4), x0(3)));
    
    if idx == 3 % Add legend only to the last subplot
        legend({'Data N1', 'Data N2', 'Model N1', 'Model N2'});
    end
end