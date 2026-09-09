%% ============================================================
% BREAST CANCER MODEL
% EARLY vs LATE DETECTION
% WITH PERIODIC CHEMOTHERAPY u(t)
% ============================================================

clear; close all; clc;

fprintf('Running Early vs Late Detection Simulation...\n');

%% ===============================
% PARAMETERS
%% ===============================

p.r = 0.18;              % Tumor growth rate
p.K = 5e8;               % Carrying capacity

% Immune killing
p.gamma1 = 5e-9;
p.gamma2 = 5e-7;
p.gamma3 = 5e-7;

% Baseline immune production
p.sy = 5e3;
p.sN = 5e3;
p.sG = 200;

% Natural death
p.mu_y = 0.04;
p.mu_N = 0.04;
p.mu_G = 0.04;

% Immune stimulation
p.rho_y = 1e-7;
p.rho_N = 1e-7;
p.rho_G = 1e-7;

p.m = 1e5;
p.k = 1e5;
p.a = 1e5;

p.eta = 1e-8;

% Drug effects
p.kappa_x = 1.2;
p.kappa_y = 0.2;
p.kappa_N = 0.2;
p.kappa_G = 0.2;

p.lambda = 1.0;          % Drug clearance

%% ===============================
% CHEMOTHERAPY SCHEDULE
%% ===============================

p.u0 = 1.0;              % Dose strength when ON
p.cycle_length = 21;     % 21-day cycle
p.treatment_days = 5;    % First 5 days ON

%% ===============================
% INITIAL CONDITIONS
%% ===============================

% Early detection (single transformed cell)
IC_early = [1; 5e4; 5e4; 5e3; 0];

% Late detection (clinically detectable tumor)
IC_late = [1e7; 2e4; 2e4; 2e3; 0];

tspan = [0 100];

%% ===============================
% SOLVE ODE SYSTEM
%% ===============================

[tE, solE] = ode45(@(t,y) tumor_model(t,y,p), ...
                   tspan, IC_early);

[tL, solL] = ode45(@(t,y) tumor_model(t,y,p), ...
                   tspan, IC_late);

fprintf('Simulation complete.\n');

%% ===============================
% PLOTS
%% ===============================

figure('Position',[100 100 1300 800])

% -----------------------------
% Tumor Comparison
% -----------------------------

subplot(3,1,1)

semilogy(tE, solE(:,1),'b','LineWidth',2);
hold on

semilogy(tL, solL(:,1),'r','LineWidth',2);

title('Tumor Dynamics')
ylabel('Tumor Cells (log scale)')
legend('Early Detection','Late Detection')
grid on

% -----------------------------
% Immune Comparison (CTL)
% -----------------------------

subplot(3,1,2)

plot(tE, solE(:,2),'b','LineWidth',2);
hold on

plot(tL, solL(:,2),'r','LineWidth',2);

title('CTL Dynamics')
ylabel('CTL Population')
legend('Early Detection','Late Detection')
grid on

% -----------------------------
% Drug Concentration
% -----------------------------

subplot(3,1,3)

plot(tE, solE(:,5),'k','LineWidth',2);

title('Drug Concentration (Periodic Chemotherapy)')
xlabel('Time (days)')
ylabel('Drug Level')
grid on

sgtitle('Early vs Late Detection Under Periodic Chemotherapy')


%% ============================================================
% MODEL FUNCTION WITH TIME-DEPENDENT u(t)
% ============================================================

function dydt = tumor_model(t,y,p)

x = max(y(1),0);
CTL = max(y(2),0);
NK = max(y(3),0);
Gd = max(y(4),0);
D = max(y(5),0);

% -----------------------------
% PERIODIC CHEMOTHERAPY u(t)
% -----------------------------

if mod(t,p.cycle_length) <= p.treatment_days
    u = p.u0;
else
    u = 0;
end

dydt = zeros(5,1);

% -----------------------------
% Tumor
% -----------------------------

dydt(1) = p.r*x*(1 - x/p.K) ...
    - p.gamma1*x*CTL ...
    - p.gamma2*x*NK ...
    - p.gamma3*x*Gd ...
    - p.kappa_x*x*D;

% -----------------------------
% CTL
% -----------------------------

dydt(2) = p.sy ...
    + (p.rho_y*x*CTL)/(p.m + x) ...
    - p.mu_y*CTL ...
    - p.eta*x*CTL ...
    - p.kappa_y*CTL*D;

% -----------------------------
% NK
% -----------------------------

dydt(3) = p.sN ...
    + (p.rho_N*x*NK)/(p.k + x) ...
    - p.mu_N*NK ...
    - p.kappa_N*NK*D;

% -----------------------------
% Gamma delta T
% -----------------------------

dydt(4) = p.sG ...
    + (p.rho_G*x*Gd)/(p.a + x) ...
    - p.mu_G*Gd ...
    - p.kappa_G*Gd*D;

% -----------------------------
% Drug
% -----------------------------

dydt(5) = u - p.lambda*D;

end