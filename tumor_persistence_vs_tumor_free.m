clc;
clear;
close all;

%% ================= PARAMETERS =================
params.r = 0.32;
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

tspan = [0 300];

%% ================= EARLY DETECTION INITIAL CONDITION =================
% Small tumor burden
X0 = [10, 2e4, 2e4, 200, 0];

%% ================= CASE 1: TUMOR PERSISTENCE =================
params.u = 0.9; % Weak drug

[t1, X1] = ode45(@(t,X) model(t,X,params), tspan, X0);

tumor_persistent = X1(:,1);

%% ================= CASE 2: TUMOR-FREE =================
params.u = 1.6; % Strong drug

[t2, X2] = ode45(@(t,X) model(t,X,params), tspan, X0);

tumor_free = X2(:,1);

%% ================= PLOTS =================

figure

subplot(1,2,1)

plot(t1, tumor_persistent, 'r', 'LineWidth', 2)

xlabel('Time (days)')
ylabel('Tumor Cells')

title('Early Detection: Tumor Persistence (u = 0.9)')

grid on

subplot(1,2,2)

plot(t2, tumor_free, 'b', 'LineWidth', 2)

xlabel('Time (days)')
ylabel('Tumor Cells')

title('Early Detection: Tumor-Free (u = 1.6)')

grid on

sgtitle('Early Detection: Tumor-Free vs Tumor Persistence')

%% ================= SAME GRAPH COMPARISON =================

figure

plot(t1, tumor_persistent, 'r', 'LineWidth', 2)

hold on

plot(t2, tumor_free, 'b', 'LineWidth', 2)

xlabel('Time (days)')
ylabel('Tumor Cells')

legend('Persistent (u=0.9)', 'Tumor-Free (u=1.6)')

title('Early Detection Comparison')

grid on

%% ================= MODEL FUNCTION =================

function dXdt = model(~, X, params)

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