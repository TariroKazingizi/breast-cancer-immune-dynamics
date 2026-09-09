%% ============================================================
% TUMOR-IMMUNE DYNAMICS UNDER ADAPTIVE CHEMOTHERAPY
% ============================================================

clc;
clear;
close all;

%% ================= PARAMETERS =================

params.r = 0.32;              % Tumor growth rate
params.K = 1e9;

params.mu_y = 0.0412;
params.mu_N = 0.0412;
params.mu_G = 0.0412;

params.eta = 5e-8;
params.lambda = 4.16;

params.sy = 1.3e4;
params.sN = 1.3e4;
params.sG = 0.1;

params.gamma1 = 4e-9;
params.gamma2 = 4e-7;
params.gamma3 = 4e-7;

params.kappa_x = 0.8;
params.kappa_y = 0.1;
params.kappa_N = 0.1;
params.kappa_G = 0.1;

params.m = 1e5;
params.k = 1e5;
params.a = 1e5;

params.rho_y = 0.8;
params.rho_N = 0.8;
params.rho_G = 0.8;

%% ================= THERAPY PARAMETERS =================

params.u_constant = 1.2;      % Constant drug level
params.umax = 1.6;             % Adaptive maximum
params.h = 5e4;                % Half-saturation for adaptive

%% ================= SIMULATION SETTINGS =================

tspan = [0 300];

X0 = [10, 2e4, 2e4, 200, 0];
% Early detection: small tumor

%% ================= CASE 1: CONSTANT CHEMO =================

[t_const, X_const] = ode45( ...
    @(t,X) model_constant(t,X,params), ...
    tspan, X0);

%% ================= CASE 2: ADAPTIVE CHEMO =================

[t_adapt, X_adapt] = ode45( ...
    @(t,X) model_adaptive(t,X,params), ...
    tspan, X0);

%% ================= CASE 3: PARAMETER VARIATION =================

params_high = params;
params_high.r = 0.40;          % Faster tumor growth

[t_high, X_high] = ode45( ...
    @(t,X) model_adaptive(t,X,params_high), ...
    tspan, X0);

%% ================= PLOTS =================

figure

subplot(2,2,1)

plot(t_const, X_const(:,1),'r','LineWidth',2)
hold on

plot(t_adapt, X_adapt(:,1),'b','LineWidth',2)

xlabel('Time (days)')
ylabel('Tumor Cells')
legend('Constant','Adaptive')
title('Tumor Dynamics')
grid on

subplot(2,2,2)

plot(t_adapt, X_adapt(:,2),'b','LineWidth',2)
hold on

plot(t_adapt, X_adapt(:,3),'g','LineWidth',2)
plot(t_adapt, X_adapt(:,4),'m','LineWidth',2)

xlabel('Time (days)')
ylabel('Immune Cells')
legend('CTL','NK','\gamma\delta T')
title('Immune Dynamics (Adaptive)')
grid on

subplot(2,2,3)

u_adaptive = params.umax .* X_adapt(:,1) ./ ...
    (X_adapt(:,1) + params.h);

plot(t_adapt, u_adaptive,'k','LineWidth',2)

xlabel('Time (days)')
ylabel('u(t)')
title('Adaptive Drug Profile')
grid on

subplot(2,2,4)

plot(t_adapt, X_adapt(:,1),'b','LineWidth',2)
hold on

plot(t_high, X_high(:,1),'r--','LineWidth',2)

xlabel('Time (days)')
ylabel('Tumor Cells')
legend('r = 0.32','r = 0.40')
title('Parameter Variation (Tumor Growth Rate)')
grid on

sgtitle('Tumor–Immune Dynamics under Adaptive Chemotherapy')


%% ================= MODEL: CONSTANT CHEMO =================

function dXdt = model_constant(~, X, p)

x = X(1);
y = X(2);
N = X(3);
G = X(4);
D = X(5);

u = p.u_constant;

dxdt = p.r*x*(1 - x/p.K) ...
    - p.gamma1*x*y ...
    - p.gamma2*x*N ...
    - p.gamma3*x*G ...
    - p.kappa_x*D*x;

dydt = p.sy ...
    + (p.rho_y*x*y)/(p.m+x) ...
    - p.mu_y*y ...
    - p.eta*x*y ...
    - p.kappa_y*D*y;

dNdt = p.sN ...
    + (p.rho_N*x*N)/(p.k+x) ...
    - p.mu_N*N ...
    - p.kappa_N*D*N;

dGdt = p.sG ...
    + (p.rho_G*x*G)/(p.a+x) ...
    - p.mu_G*G ...
    - p.kappa_G*D*G;

dDdt = u - p.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end


%% ================= MODEL: ADAPTIVE CHEMO =================

function dXdt = model_adaptive(~, X, p)

x = X(1);
y = X(2);
N = X(3);
G = X(4);
D = X(5);

u = p.umax * x/(x + p.h);

dxdt = p.r*x*(1 - x/p.K) ...
    - p.gamma1*x*y ...
    - p.gamma2*x*N ...
    - p.gamma3*x*G ...
    - p.kappa_x*D*x;

dydt = p.sy ...
    + (p.rho_y*x*y)/(p.m+x) ...
    - p.mu_y*y ...
    - p.eta*x*y ...
    - p.kappa_y*D*y;

dNdt = p.sN ...
    + (p.rho_N*x*N)/(p.k+x) ...
    - p.mu_N*N ...
    - p.kappa_N*D*N;

dGdt = p.sG ...
    + (p.rho_G*x*G)/(p.a+x) ...
    - p.mu_G*G ...
    - p.kappa_G*D*G;

dDdt = u - p.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end