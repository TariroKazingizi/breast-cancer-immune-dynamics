clc;
clear;
close all;

%% ============================================================
% STABILITY AND BIFURCATION ANALYSIS
% Breast Cancer-Immune Dynamics with Chemotherapy
% ============================================================

%% ================= PARAMETERS =================

params = get_parameters();

%% ================= BIFURCATION ANALYSIS =================

u_values = linspace(0.2, 2.0, 40);
tumor_max = zeros(size(u_values));

X0 = [10, 2e4, 2e4, 200, 0];
tspan = [0 500];

for i = 1:length(u_values)

    params.u = u_values(i);

    [~, X] = ode45( ...
        @(t,X) model_autonomous(t,X,params), ...
        tspan, X0);

    x = X(:,1);

    % Remove transient: keep last 30%
    x_ss = x(round(end*0.7):end);

    tumor_max(i) = max(x_ss);

end

%% ================= BIFURCATION DIAGRAM =================

figure;

plot(u_values, tumor_max, 'LineWidth', 2);

xlabel('Drug amplitude u');
ylabel('Maximum Tumor (after transient)');
title('Bifurcation Diagram (Tumor max vs u)');
grid on;


%% ================= EIGENVALUE ANALYSIS =================

% Choose drug amplitude near transition
params.u = 1.05;

X0 = [10, 2e4, 2e4, 200, 0];

[t, X] = ode45( ...
    @(t,X) model_autonomous(t,X,params), ...
    [0 300], X0);

E = X(end,:)';

J = numerical_jacobian( ...
    @(X) autonomous_vector(X,params), E);

eigvals = eig(J);

disp(' ');
disp('Approximate Equilibrium:');
disp(E);

disp('Eigenvalues of Jacobian:');
disp(eigvals);


%% ================= PARAMETERS FUNCTION =================

function params = get_parameters()

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

params.u = 1.05;

end


%% ================= ODE MODEL =================

function dXdt = model_autonomous(~, X, params)

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


%% ================= AUTONOMOUS VECTOR =================

function F = autonomous_vector(X, params)

F = model_autonomous(0, X, params);

end


%% ================= NUMERICAL JACOBIAN =================

function J = numerical_jacobian(f, X)

n = length(X);
J = zeros(n);

h = 1e-6;

for i = 1:n

    X1 = X;
    X2 = X;

    X1(i) = X1(i) + h;
    X2(i) = X2(i) - h;

    J(:,i) = (f(X1) - f(X2)) / (2*h);

end

end
