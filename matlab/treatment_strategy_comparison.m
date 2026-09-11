%DIFFERENT TREATMENT SCHEMES

clc;
clear;
close all;

%% ================= PARAMETERS =================
p.r = 0.32;
p.K = 1e9;
p.mu_y = 0.0412;
p.mu_N = 0.0412;
p.mu_G = 0.0412;
p.eta = 5e-8;
p.lambda = 4.16;
p.sy = 1.3e4;
p.sN = 1.3e4;
p.sG = 0.1;
p.gamma1 = 4e-9;
p.gamma2 = 4e-7;
p.gamma3 = 4e-7;
p.kappa_x = 0.8;
p.kappa_y = 0.1;
p.kappa_N = 0.1;
p.kappa_G = 0.1;
p.m = 1e5;
p.k = 1e5;
p.a = 1e5;
p.rho_y = 0.8;
p.rho_N = 0.8;
p.rho_G = 0.8;

%% Treatment parameters
p.u_constant = 1.2;
p.u_high = 2.5;
p.umax = 1.6;
p.h = 5e4;

%% ================= SIMULATION SETTINGS =================
tspan = [0 300];
X0 = [10, 2e4, 2e4, 200, 0]; % Early detection

%% ================= RUN SIMULATIONS =================
[t_none, X_none] = ode45(@(t,X) model_none(t,X,p), tspan, X0);
[t_const, X_const] = ode45(@(t,X) model_constant(t,X,p), tspan, X0);
[t_high, X_high] = ode45(@(t,X) model_high(t,X,p), tspan, X0);
[t_adapt, X_adapt] = ode45(@(t,X) model_adaptive(t,X,p), tspan, X0);

%% ================= TUMOR COMPARISON =================
figure
plot(t_none, X_none(:,1),'k','LineWidth',2)
hold on
plot(t_const, X_const(:,1),'r','LineWidth',2)
plot(t_high, X_high(:,1),'m','LineWidth',2)
plot(t_adapt, X_adapt(:,1),'b','LineWidth',2)

xlabel('Time (days)')
ylabel('Tumor Cells')
legend('No Treatment','Constant','High Constant','Adaptive')
title('Tumor Comparison')
grid on

%% ================= IMMUNE COMPARISON (CTL only) =================
figure
plot(t_none, X_none(:,2),'k','LineWidth',2)
hold on
plot(t_const, X_const(:,2),'r','LineWidth',2)
plot(t_high, X_high(:,2),'m','LineWidth',2)
plot(t_adapt, X_adapt(:,2),'b','LineWidth',2)

xlabel('Time (days)')
ylabel('CTL Cells')
legend('No Treatment','Constant','High Constant','Adaptive')
title('Immune (CTL) Comparison')
grid on

%% ================= DRUG PROFILES =================
figure

u_none = zeros(size(t_none));
u_const = p.u_constant*ones(size(t_const));
u_high = p.u_high*ones(size(t_high));
u_adapt = p.umax .* X_adapt(:,1) ./ (X_adapt(:,1)+p.h);

plot(t_none, u_none,'k','LineWidth',2)
hold on
plot(t_const, u_const,'r','LineWidth',2)
plot(t_high, u_high,'m','LineWidth',2)
plot(t_adapt, u_adapt,'b','LineWidth',2)

xlabel('Time (days)')
ylabel('u(t)')
legend('No Treatment','Constant','High Constant','Adaptive')
title('Drug Administration Profiles')
grid on

%% ================= CUMULATIVE DRUG EXPOSURE =================
U_const = trapz(t_const, u_const);
U_high = trapz(t_high, u_high);
U_adapt = trapz(t_adapt, u_adapt);

fprintf('\nCUMULATIVE DRUG EXPOSURE\n')
fprintf('Constant: %.2f\n', U_const)
fprintf('High Constant: %.2f\n', U_high)
fprintf('Adaptive: %.2f\n', U_adapt)

%% ================= MODEL FUNCTIONS =================

function dXdt = model_none(~,X,p)

x=X(1); y=X(2); N=X(3); G=X(4); D=X(5);

u = 0;

dxdt = p.r*x*(1-x/p.K) ...
 - p.gamma1*x*y - p.gamma2*x*N - p.gamma3*x*G;

dydt = p.sy + (p.rho_y*x*y)/(p.m+x) ...
 - p.mu_y*y - p.eta*x*y;

dNdt = p.sN + (p.rho_N*x*N)/(p.k+x) ...
 - p.mu_N*N;

dGdt = p.sG + (p.rho_G*x*G)/(p.a+x) ...
 - p.mu_G*G;

dDdt = u - p.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end

function dXdt = model_constant(~,X,p)

x=X(1); y=X(2); N=X(3); G=X(4); D=X(5);

u = p.u_constant;

dxdt = p.r*x*(1-x/p.K) ...
 - p.gamma1*x*y - p.gamma2*x*N - p.gamma3*x*G ...
 - p.kappa_x*D*x;

dydt = p.sy + (p.rho_y*x*y)/(p.m+x) ...
 - p.mu_y*y - p.eta*x*y - p.kappa_y*D*y;

dNdt = p.sN + (p.rho_N*x*N)/(p.k+x) ...
 - p.mu_N*N - p.kappa_N*D*N;

dGdt = p.sG + (p.rho_G*x*G)/(p.a+x) ...
 - p.mu_G*G - p.kappa_G*D*G;

dDdt = u - p.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end

function dXdt = model_high(~,X,p)

p.u_constant = p.u_high;
dXdt = model_constant([],X,p);

end

function dXdt = model_adaptive(~,X,p)

x=X(1); y=X(2); N=X(3); G=X(4); D=X(5);

u = p.umax * x/(x + p.h);

dxdt = p.r*x*(1-x/p.K) ...
 - p.gamma1*x*y - p.gamma2*x*N - p.gamma3*x*G ...
 - p.kappa_x*D*x;

dydt = p.sy + (p.rho_y*x*y)/(p.m+x) ...
 - p.mu_y*y - p.eta*x*y - p.kappa_y*D*y;

dNdt = p.sN + (p.rho_N*x*N)/(p.k+x) ...
 - p.mu_N*N - p.kappa_N*D*N;

dGdt = p.sG + (p.rho_G*x*G)/(p.a+x) ...
 - p.mu_G*G - p.kappa_G*D*G;

dDdt = u - p.lambda*D;

dXdt = [dxdt; dydt; dNdt; dGdt; dDdt];

end