clc; clear; close all;
%% Main

% Define a simplified symmetrical parameter set and extract some useful
% parameters
theta = [1, 10, 10, 1, 10, 20, 1, 10, 10, 1, 10, 20];
chemostat_params = [0.2, 1, 1, 11.10];

k13 = theta(02); k15 = theta(03);
k23 = theta(08); k24 = theta(09);


% Define the initial conditions and the timespan of the simulation
x_init = [0.1, 0.11, 1, 1, 11.10];
t_span = linspace(0, 200, 100000);


% Run the model simulation and extract the state variables
model = @(t, x) model_full(t, x, theta, chemostat_params);
[t, x_mod] = ode45(model, t_span, x_init);

N1 = x_mod(:,1); N2 = x_mod(:,2);
R1 = x_mod(:,3); R2 = x_mod(:,4); R3 = x_mod(:,5);


% Compute the resource specific growth rates over time
mu13 = R3./(k13 + R3); mu15 = R2./(k15 + R2);
mu23 = R3./(k23 + R3); mu24 = R1./(k24 + R1);



%% Data plotting
figure;

% Plot the auxotroph absolute abundance
subplot(3, 1, 1);
plot(t, N1, 'linewidth',2, 'linestyle', '-', 'Color','#00D5FF', 'DisplayName','\Delta{\ittyrA}'); hold on;
plot(t, N2, 'linewidth',2, 'linestyle', '-', 'Color','#FF4F00', 'DisplayName','\Delta{\itpheA}')
xlim([151.2352, 183.4325]); ylim([-0.1, 3.1]); yticks([0, 1, 2, 3]);
title('Full Model Simulation'); ylabel('Abundance'); legend;

% Plot the resource concentrations
subplot(3, 1, 2)
plot(t, R1, 'linewidth',2, 'linestyle', '-', 'Color','#FF4F00', 'DisplayName','Phe'); hold on;
plot(t, R2, 'linewidth',2, 'linestyle', '-', 'Color','#00D5FF', 'DisplayName','Tyr')
plot(t, R3, 'linewidth',2, 'linestyle', '-', 'Color','#ffd700', 'DisplayName','Glu')
xlim([151.2352, 183.4325]); ylim([-5, 65]); yticks([0, 30, 60]);
ylabel('Concentration'); legend;

% Plot the intervals of amino acid limitation
subplot(3, 1, 3)
plot(t, mu15<mu13, 'linewidth',2, 'linestyle', '-', 'Color','#0095b3', 'DisplayName','Tyr Limited'); hold on;
plot(t, mu24<mu23, 'linewidth',2, 'linestyle', '-', 'Color','#b33700', 'DisplayName','Phe Limited');
xlim([151.2352, 183.4325]); ylim([-0.1, 1.1]); yticks([0, 1]);
xlabel('Time (hr.)'), ylabel('True/False'); legend;
