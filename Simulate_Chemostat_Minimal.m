clc; clear; close all;
%% Main

% Define a simplified symmetrical parameter set and extract some useful
% parameters for later
theta = [0.1, 0.1, 10, 20, 0.1, 0.1, 10, 20, 2.62, 6.14];
chemostat_params = 0.2;

a12 = theta(1); a13 = theta(2);
a21 = theta(5); a23 = theta(6);
N_tot = theta(9);


% Define the initial conditions and the timespan of the simulation
x_init = [1, 1, 0.51];
t_span = linspace(0, 200, 10000);


% Run the model simulation and extract the state variables
model = @(t, x) model_minimal(t, x, theta, chemostat_params);
[t, x_mod] = ode45(model, t_span, x_init);

N1 = (0 + x_mod(:,3))*N_tot; N2 = (1 - x_mod(:,3))*N_tot;
R1 = x_mod(:,1); R2 = x_mod(:,2); R3 = ones(size(t))*theta(10);


% Compute the resource specific growth rates over time
mu12 = a12*R2; mu13 = a13*R3;
mu21 = a21*R1; mu23 = a23*R3;


%% Data plotting
figure;

% Plot the auxotroph absolute abundance
subplot(3, 1, 1)
plot(t, N1, 'linewidth',2, 'linestyle', '-', 'Color','#00d5ff', 'DisplayName','\Delta{\ittyrA}'); hold on;
plot(t, N2, 'linewidth',2, 'linestyle', '-', 'Color','#ff4f00', 'DisplayName','\Delta{\itpheA}')
xlim([159.4368, 188.2088]); ylim([-0.1, 3.1]); yticks([0, 1, 2, 3]);
title('Minimal Model Simulation'); ylabel('Abundance'); legend;

% Plot the resource concentrations
subplot(3, 1, 2)
plot(t, R1, 'linewidth',2, 'linestyle', '-', 'Color','#ff4f00', 'DisplayName','Phe'); hold on;
plot(t, R2, 'linewidth',2, 'linestyle', '-', 'Color','#00d5ff', 'DisplayName','Tyr')
plot(t, R3, 'linewidth',2, 'linestyle', '-', 'Color','#ffd700', 'DisplayName','Glu')
xlim([159.4368, 188.2088]); ylim([-5, 65]); yticks([0, 30, 60]);
ylabel('Concentration'); legend;

% Plot the intervals of amino acid limitation
subplot(3, 1, 3)
plot(t, mu12<mu13, 'linewidth',2, 'linestyle', '-', 'Color','#0095b3', 'DisplayName','Tyr Limited.'); hold on;
plot(t, mu21<mu23, 'linewidth',2, 'linestyle', '-', 'Color','#b33700', 'DisplayName','Phe Limited');
xlim([159.4368, 188.2088]); ylim([-0.1, 1.1]); yticks([0, 1]);
xlabel('Time (hr.)'), ylabel('True/False'); legend;