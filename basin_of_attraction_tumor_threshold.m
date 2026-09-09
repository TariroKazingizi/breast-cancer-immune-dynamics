clear;
clc;
close all;

%% Parameters

p.r = 0.18;
p.K = 1e9;

p.gamma1 = 5e-8;
p.gamma2 = 1e-6;
p.gamma3 = 1e-6;

p.kappa_x = 0; % no chemo for basin study

p.s_y = 1.3e4;
p.s_N = 1.3e4;
p.s_G = 0.1;

p.mu_y = 0.0412;
p.mu_N = 0.0412;
p.mu_G = 0.0412;

p.eta = 1e-7;

p.m = 1e5;
p.k = 1e5;
p.a = 1e5;

p.lambda = 4.16;

p.rho_y = 0.01;
p.rho_N = 0.01;
p.rho_G = 1e-4;

tspan = [0 300];

%% Initial tumor range

tumor_init = logspace(1,7,40);
final_state = zeros(size(tumor_init));

for i = 1:length(tumor_init)

    X0 = [tumor_init(i); 1e5; 1e5; 1e3; 0];

    [~, X] = ode45(@(t,X) model(t,X,p), tspan, X0);

    final_state(i) = X(end,1);

end

%% Plot

figure;

semilogx(tumor_init, final_state, 'o-', 'LineWidth', 2);

xlabel('Initial Tumor Cells');
ylabel('Final Tumor Cells');

title('Basin of Attraction: Tumor Threshold Dynamics');

grid on;

%% MODEL

function dXdt = model(~,X,p)

x = X(1);
y = X(2);
N = X(3);
G = X(4);
D = X(5);

dx = p.r*x*(1 - x/p.K) ...
    - p.gamma1*x*y ...
    - p.gamma2*x*N ...
    - p.gamma3*x*G;

dy = p.s_y ...
    + (p.rho_y*x)/(p.m+x) ...
    - p.mu_y*y ...
    - p.eta*x*y;

dN = p.s_N ...
    + (p.rho_N*x)/(p.k+x) ...
    - p.mu_N*N;

dG = p.s_G ...
    + (p.rho_G*x*D)/(p.a+x) ...
    - p.mu_G*G;

dD = -p.lambda*D;

dXdt = [dx; dy; dN; dG; dD];

end