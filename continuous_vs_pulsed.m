clc;
clear;
close all;

%% ================= PARAMETERS =================

params.r = 0.18;
params.K = 1e9;

params.mu_y = 0.0412;
params.mu_N = 0.0412;
params.mu_G = 0.0412;

params.eta = 1e-7;
params.lambda = 4.16;

params.sy = 1.3e4;
params.sN = 1.3e4;
params.sG = 0.1;

params.gamma1 = 1.1e-8;
params.gamma2 = 1e-6;
params.gamma3 = 1e-6;

params.kappa_x = 1.0;
params.kappa_y = 0.2;
params.kappa_N = 0.2;
params.kappa_G = 0.2;

params.m = 1e5;
params.k = 1e5;
params.a = 1e5;

params.rho_y = 1.0;
params.rho_N = 1.0;
params.rho_G = 1.0;

params.u = 1.0;

%% ================= INITIAL CONDITIONS =================

X0 = [10, 1e4, 1e4, 100, 0];

tspan = [0 100];

%% ================= CONTINUOUS CHEMOTHERAPY =================

[t_cont, X_cont] = ode45( ...
    @(t,X) model_continuous(t,X,params), ...
    tspan, X0);

%% ================= PULSED CHEMOTHERAPY =================

[t_pulse, X_pulse] = ode45( ...
    @(t,X) model_pulsed(t,X,params), ...
    tspan, X0);

%% ================= COMPARISON PLOTS =================

figure;

% Tumor dynamics
subplot(2,2,1)

plot(t_cont, X_cont(:,1), 'LineWidth', 2);
hold on;

plot(t_pulse, X_pulse(:,1), 'LineWidth', 2);

xlabel('Time (days)');
ylabel('Cancer Cells');
title('Tumor Dynamics');

legend('Continuous', 'Pulsed');
grid on;

% Immune cell dynamics
subplot(2,2,2)

plot(t_cont, X_cont(:,2), 'LineWidth', 2);
hold on;

plot(t_cont, X_cont(:,3), 'LineWidth', 2);
plot(t_cont, X_cont(:,4), 'LineWidth', 2);

xlabel('Time (days)');
ylabel('Immune Cells');

legend('CTLs', 'NK', '\gamma\delta T');

title('Immune Cell Dynamics');
grid on;

% Drug concentration
subplot(2,2,3)

plot(t_cont, X_cont(:,5), 'LineWidth', 2);
hold on;

plot(t_pulse, X_pulse(:,5), 'LineWidth', 2);

xlabel('Time (days)');
ylabel('Drug');

title('Drug Concentration');

legend('Continuous', 'Pulsed');
grid on;

% Phase plane
subplot(2,2,4)

plot(X_cont(:,1), X_cont(:,2), 'LineWidth', 2);
hold on;

plot(X_pulse(:,1), X_pulse(:,2), 'LineWidth', 2);

xlabel('Tumor');
ylabel('CTLs');

title('Phase Plane (Tumor vs CTL)');

legend('Continuous', 'Pulsed');
grid on;


%% ================= CONTINUOUS MODEL =================

function dXdt = model_continuous(t, X, params)

x = X(1);
y = X(2);
N = X(3);
G = X(4);
D = X(5);

u = params.u;

dxdt = params.r*x*(1 - x/params.K) ...
    - params.gamma1*x*y ...
    - params.gamma2*x*N ...
    - params.gamma3*x*G ...
    - params.kappa_x*D*x;

dydt = params.sy ...
    + (params.rho_y*x*y)/(params.m + x) ...
    - params.mu_y*y ...
    - params.eta*x*y ...
    - params.kappa_y*D*y;

dNdt = params.sN ...
    + (params.rho_N*x*N)/(params.k + x) ...
    - params.mu_N*N ...
    - params.kappa_N*D*N;

dGdt = params.sG ...
    + (params.rho_G*x*G)/(params.a + x) ...
    - params.mu_G*G ...
    - params.kappa_G*D*G;

dDdt = u - params.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end


%% ================= PULSED MODEL =================

function dXdt = model_pulsed(t, X, params)

x = X(1);
y = X(2);
N = X(3);
G = X(4);
D = X(5);

% Drug administered 1 day every 7 days

if mod(t,7) < 1
    u = params.u;
else
    u = 0;
end

dxdt = params.r*x*(1 - x/params.K) ...
    - params.gamma1*x*y ...
    - params.gamma2*x*N ...
    - params.gamma3*x*G ...
    - params.kappa_x*D*x;

dydt = params.sy ...
    + (params.rho_y*x*y)/(params.m + x) ...
    - params.mu_y*y ...
    - params.eta*x*y ...
    - params.kappa_y*D*y;

dNdt = params.sN ...
    + (params.rho_N*x*N)/(params.k + x) ...
    - params.mu_N*N ...
    - params.kappa_N*D*N;

dGdt = params.sG ...
    + (params.rho_G*x*G)/(params.a + x) ...
    - params.mu_G*G ...
    - params.kappa_G*D*G;

dDdt = u - params.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end