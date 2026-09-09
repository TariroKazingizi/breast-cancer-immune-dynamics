function breast_cancer_coexistence_short()

clc;
clear;
close all;

% Parameters adjusted for coexistence
r = 0.45;
K = 1e9;

mu_y = 0.0412;
mu_N = 0.0412;
mu_G = 0.0412;

eta = 1e-9;
lambda = 4.16;

s_y = 1.3e4;
s_N = 1.3e4;
s_G = 0.1;

gamma1 = 1e-10;
gamma2 = 1e-9;
gamma3 = 1e-9;

kappa_x = 0.01;
kappa_y = 0.6;
kappa_N = 0.6;
kappa_G = 0.6;

m = 1e5;
k = 1e5;
a = 1e5;

rho_y = 0.1;
rho_N = 0.1;
rho_G = 0.1;

u = 50;

% Equilibrium for reference lines
E_star = find_equilibrium(r, K, mu_y, mu_N, mu_G, eta, lambda, ...
    s_y, s_N, s_G, gamma1, gamma2, gamma3, ...
    kappa_x, kappa_y, kappa_N, kappa_G, ...
    m, k, a, rho_y, rho_N, rho_G, u);

% Time span for simulation
tspan = [0 100];

% Initial conditions
initial_conditions = [
    1e3, 5e5, 5e5, 1e5, 0;
    5e3, 5e5, 5e5, 1e5, 0;
    1e4, 5e5, 5e5, 1e5, 0;
];

% Create figure
figure('Position', [100, 100, 1200, 800]);

colors = lines(size(initial_conditions,1));

for i = 1:size(initial_conditions,1)

    [t, Y] = ode45(@(t,y) model_eqns(t, y, ...
        r, K, mu_y, mu_N, mu_G, eta, lambda, ...
        s_y, s_N, s_G, gamma1, gamma2, gamma3, ...
        kappa_x, kappa_y, kappa_N, kappa_G, ...
        m, k, a, rho_y, rho_N, rho_G, u), ...
        tspan, initial_conditions(i,:)');

    % Plot cancer cells
    subplot(2,2,1);
    semilogy(t, Y(:,1), ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

    % Plot CTLs
    subplot(2,2,2);
    semilogy(t, Y(:,2), ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

    % Phase plane: cancer vs CTLs
    subplot(2,2,3);
    loglog(Y(:,1), Y(:,2), ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

    % Phase plane: NK vs gamma-delta T
    subplot(2,2,4);
    loglog(Y(:,3), Y(:,4), ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

end

% Format plots
subplot(2,2,1);
yline(E_star(1), '--r', 'x*', 'LineWidth', 1.5);
xlabel('Time (days)');
ylabel('Cancer cells');
title('Cancer Cell Dynamics');
grid on;

subplot(2,2,2);
yline(E_star(2), '--r', 'y*', 'LineWidth', 1.5);
xlabel('Time (days)');
ylabel('CTLs');
title('CTL Dynamics');
grid on;

subplot(2,2,3);
plot(E_star(1), E_star(2), ...
    'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Cancer cells');
ylabel('CTLs');
title('Phase Plane: Cancer vs CTLs');
grid on;

subplot(2,2,4);
plot(E_star(3), E_star(4), ...
    'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('NK cells');
ylabel('gamma-delta T cells');
title('Phase Plane: NK vs gamma-delta T');
grid on;

legend_labels = arrayfun(@(i) sprintf('Case %d', i), ...
    1:size(initial_conditions,1), ...
    'UniformOutput', false);

legend(legend_labels, 'Location', 'best');

end


function E_star = find_equilibrium(r, K, mu_y, mu_N, mu_G, ...
    eta, lambda, s_y, s_N, s_G, ...
    gamma1, gamma2, gamma3, ...
    kappa_x, kappa_y, kappa_N, kappa_G, ...
    m, k, a, rho_y, rho_N, rho_G, u)

tspan = [0 2000];

y0 = [
    1e5;
    s_y/mu_y;
    s_N/mu_N;
    s_G/mu_G;
    u/lambda
];

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[~, Y] = ode45(@(t,y) model_eqns(t, y, ...
    r, K, mu_y, mu_N, mu_G, eta, lambda, ...
    s_y, s_N, s_G, gamma1, gamma2, gamma3, ...
    kappa_x, kappa_y, kappa_N, kappa_G, ...
    m, k, a, rho_y, rho_N, rho_G, u), ...
    tspan, y0, options);

E_star = Y(end,:)';

E_star(5) = u/lambda;

end


function dydt = model_eqns(~, y, r, K, mu_y, mu_N, mu_G, ...
    eta, lambda, s_y, s_N, s_G, ...
    gamma1, gamma2, gamma3, ...
    kappa_x, kappa_y, kappa_N, kappa_G, ...
    m, k, a, rho_y, rho_N, rho_G, u)

x = max(y(1), 0);
y_ctl = max(y(2), 0);
N = max(y(3), 0);
G = max(y(4), 0);
D = max(y(5), 0);

dxdt = r*x*(1 - x/K) ...
    - gamma1*x*y_ctl ...
    - gamma2*x*N ...
    - gamma3*x*G ...
    - kappa_x*D*x;

dycdt = s_y ...
    + (rho_y*x*y_ctl)/(m + x) ...
    - mu_y*y_ctl ...
    - eta*x*y_ctl ...
    - kappa_y*D*y_ctl;

dNdt = s_N ...
    + (rho_N*x*N)/(k + x) ...
    - mu_N*N ...
    - kappa_N*D*N;

dGdt = s_G ...
    + (rho_G*x*D)/(a + x) ...
    - mu_G*G ...
    - kappa_G*D*G;

dDdt = u - lambda*D;

dydt = [dxdt; dycdt; dNdt; dGdt; dDdt];

end
