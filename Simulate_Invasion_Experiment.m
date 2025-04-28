clc; clear; close all;
%%

% Simplified symmetrical parameter set
theta = [1, 10, 10, 1, 10, 20, 1, 10, 10, 1, 10, 20];

% No amino acid inflow
chemo_params = [0.1, 1, 1, 11.1];

% Initial conditions
x_init = [0.1, 0.11, 0, 1, 1, 11.1];

model = @(t, x) model_cheater(t, x, theta, chemo_params);
[t, x] = ode89(model, [0, 200], x_init);


plot(t, x(:,1), 'linewidth',2, 'Color','#00d5ff', 'HandleVisibility','off'); hold on
plot(t, x(:,2), 'linewidth',2, 'Color','#ff4f00', 'HandleVisibility','off');
plot(t, x(:,3), 'linewidth',2, 'Color','#32cd32', 'HandleVisibility','off');


%%
x_init = x(end,:);
x_init(3) = 0.4;
model = @(t, x) model_cheater(t, x, theta, chemo_params);
[t, x] = ode89(model, [200, 1000], x_init);


plot(t, x(:,1), 'linewidth',2, 'Color','#00d5ff', 'DisplayName','\Delta tyrA');
plot(t, x(:,2), 'linewidth',2, 'Color','#ff4f00', 'DisplayName','\Delta pheA');
plot(t, x(:,3), 'linewidth',2, 'Color','#32cd32', 'DisplayName','Cheater');

ylim([-0.1, 3.6])
legend;

xlabel('Time')
ylabel('Abundance')



%% The following simulations will show the cheater oscillating with the resident community
% clc; clear; close all;
% %%
% 
% % Simplified symmetrical parameter set
% theta = [1, 10, 10, 1, 10, 20, 1, 10, 10, 1, 10, 20];
% theta = [0.7045, 3.7551, 20.8669, 1.2770, 43.3019, 80.3646, 0.9517, 10.6978, 8.6906, 3.5704, 39.0734, 46.6017]; % Inferred
% 
% % No amino acid inflow
% chemo_params = [0.1, 40, 100, 11.1];
% 
% % Initial conditions
% x_init = [0.11, 0.11, 0, 1, 1, 11.1];
% 
% model = @(t, x) model_cheater(t, x, theta, chemo_params);
% [t, x] = ode89(model, [0, 200], x_init);
% 
% 
% plot(t, x(:,1), 'linewidth',2, 'Color','#00d5ff', 'HandleVisibility','off'); hold on
% plot(t, x(:,2), 'linewidth',2, 'Color','#ff4f00', 'HandleVisibility','off');
% plot(t, x(:,3), 'linewidth',2, 'Color','#32cd32', 'HandleVisibility','off');
% 
% 
% %%
% x_init = x(end,:);
% x_init(3) = 0.4;
% model = @(t, x) model_cheater(t, x, theta, chemo_params);
% [t, x] = ode89(model, [200, 1000], x_init);
% 
% 
% plot(t, x(:,1), 'linewidth',2, 'Color','#00d5ff', 'DisplayName','\Delta tyrA');
% plot(t, x(:,2), 'linewidth',2, 'Color','#ff4f00', 'DisplayName','\Delta pheA');
% plot(t, x(:,3), 'linewidth',2, 'Color','#32cd32', 'DisplayName','Cheater');
% 
% ylim([-0.1, 3.6])
% legend;
% 
% xlabel('Time')
% ylabel('Abundance')