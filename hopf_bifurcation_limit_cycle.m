% ============================================================
% BREAST CANCER MODEL
% HOPF BIFURCATION AND LIMIT CYCLE
% ============================================================

clc;
clear;
close all;

%% ================= PARAMETERS (Tuned Near Hopf) =================

params.r = 0.28;          % Slightly higher tumor growth
params.K = 1e9;

params.mu_y = 0.0412;
params.mu_N = 0.0412;
params.mu_G = 0.0412;

params.eta = 5e-8;        % Reduced immune suppression
params.lambda = 4.16;

params.sy = 1.3e4;
params.sN = 1.3e4;
params.sG = 0.1;

% Reduced killing to allow oscillations
params.gamma1 = 5e-9;
params.gamma2 = 5e-7;
params.gamma3 = 5e-7;

% Reduced immune drug toxicity
params.kappa_x = 0.8;
params.kappa_y = 0.2;
params.kappa_N = 0.2;
params.kappa_G = 0.2;

params.m = 1e5;
params.k = 1e5;
params.a = 1e5;

params.rho_y = 0.8;
params.rho_N = 0.8;
params.rho_G = 0.8;

params.u = 1.2;           % Drug amplitude

%% ================= INITIAL CONDITIONS =================

X0 = [10, 2e4, 2e4, 200, 0];

tspan = [0 250];

%% ================= SOLVE SYSTEM =================

[t, X] = ode45(@(t,X) model_pulsed(t,X,params), ...
               tspan, X0);

x = X(:,1);
y = X(:,2);
N = X(:,3);
G = X(:,4);
D = X(:,5);

%% ================= EQUILIBRIUM (APPROXIMATE) =================

E = X(end,:);

fprintf('\nApproximate Equilibrium:\n');
disp(E);

%% ================= PLOTS =================

figure;

% ------------------------------------------------------------
% 1. Tumor oscillations
% ------------------------------------------------------------
subplot(2,2,1)

plot(t, x, 'r', 'LineWidth', 2);
hold on;

yline(E(1), '--k', 'x*', 'LineWidth', 1.5);

xlabel('Time');
ylabel('Tumor');
title('Tumor Oscillations (Near Hopf)');
grid on;


% ------------------------------------------------------------
% 2. Immune cell oscillations
% ------------------------------------------------------------
subplot(2,2,2)

semilogy(t, y, 'b', 'LineWidth', 2);
hold on;

semilogy(t, N, 'g', 'LineWidth', 2);
semilogy(t, G, 'm', 'LineWidth', 2);

yline(E(2), '--b');
yline(E(3), '--g');
yline(E(4), '--m');

xlabel('Time');
ylabel('Immune Cells (log scale)');

legend('CTL', 'NK', '\gamma\delta T');

title('Immune Oscillations (Log Scale)');
grid on;


% ------------------------------------------------------------
% 3. Phase plane: Tumor vs CTL
% ------------------------------------------------------------
subplot(2,2,3)

plot(x, y, 'LineWidth', 2);
hold on;

plot(E(1), E(2), 'ro', 'LineWidth', 2);

xlabel('Tumor');
ylabel('CTL');

title('Spiral Coils → Hopf');
grid on;


% ------------------------------------------------------------
% 4. Drug concentration
% ------------------------------------------------------------
subplot(2,2,4)

plot(t, D, 'k', 'LineWidth', 2);

xlabel('Time');
ylabel('Drug');

title('Pulsed Drug');
grid on;


%% ================= MODEL FUNCTION =================

function dXdt = model_pulsed(t, X, params)

x = X(1);
y = X(2);
N = X(3);
G = X(4);
D = X(5);

% Pulsed 1 day every 7 days
if mod(t,7) < 1
    u = params.u;
else
    u = 0;
end

% Tumor
dxdt = params.r*x*(1 - x/params.K) ...
     - params.gamma1*x*y ...
     - params.gamma2*x*N ...
     - params.gamma3*x*G ...
     - params.kappa_x*D*x;

% CTLs
dydt = params.sy ...
     + (params.rho_y*x*y)/(params.m + x) ...
     - params.mu_y*y ...
     - params.eta*x*y ...
     - params.kappa_y*D*y;

% NK cells
dNdt = params.sN ...
     + (params.rho_N*x*N)/(params.k + x) ...
     - params.mu_N*N ...
     - params.kappa_N*D*N;

% Gamma-delta T cells
dGdt = params.sG ...
     + (params.rho_G*x*G)/(params.a + x) ...
     - params.mu_G*G ...
     - params.kappa_G*D*G;

% Drug
dDdt = u - params.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end