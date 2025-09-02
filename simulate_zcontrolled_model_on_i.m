%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_zcontrolled_model_on_i.m
%
% Description:
%   This script generates Fig. 5 of the manuscript:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%
%   It includes:
%   - Z-control applied to i(t) 
%   - Control functions G(t) for forward and backward scenarios
%   - Single figure with three square panels:
%       - Top: system dynamics
%       - Bottom-left: forward G(t)
%       - Bottom-right: backward G(t)
%
% Author: Angela Martiradonna
% Date: Semptember 2025
%
% License:
%   Creative Commons Attribution 4.0 International License (CC BY 4.0)
%
% Citation:
%   If you use this script or any part of the code, please cite the original article:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% === PARAMETERS ===
mu = 0.05; rho = 0.25; lambda = 0.02; p = 0.7;
q = 0.8; a = 0.5; k = 2; gamma = 0.2;

T = 200;
tspan = [0 T];
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% Colors
col_u = [0, 0.4471, 0.7412];
col_b = [0.8500, 0.3250, 0.0980];
col_i = [0.4660, 0.6740, 0.1880];
col_m = [0.4940, 0.1840, 0.5560];
col_n = [0, 0, 0];
deep_red = [0.6350, 0.0780, 0.1840];

%% === Z-CONTROL ON i(t) ===
zcontrol_scenario = struct('label', 'IIb', 'alpha', 2.0, 'sigma', 0.95, 'zeta', 0.03);
id = 0.01;

% Function for G(t)
compute_Gt_i = @(t, y, alpha, sigma, zeta) ...
    (1 - p) * rho .* y(:,1) + sigma - alpha * p .* y(:,3) + ...
    (zeta * (y(:,3) - id) - (lambda + mu) .* y(:,3) + gamma * (1 - q) .* y(:,4) .* y(:,1)) ./ y(:,2);

% Forward and backward scenarios
forward_scenari = {
    struct('label', 'Ia', 'alpha', 0.4, 'sigma', 0.6, 'zeta', 0.03, 'style', '-'),
    struct('label', 'Ib', 'alpha', 0.4, 'sigma', 0.8, 'zeta', 0.03, 'style', '-.')
};
backward_scenari = {
    struct('label', 'II',  'alpha', 2.0, 'sigma', 0.6,  'zeta', 0.03, 'style', '-'),
    struct('label', 'IIa', 'alpha', 2.0, 'sigma', 0.8,  'zeta', 0.03, 'style', '-.'),
    struct('label', 'IIb', 'alpha', 2.0, 'sigma', 0.95, 'zeta', 0.03, 'style', ':')
};

%% === FIGURE ===
figure(1); clf;
set(gcf, 'Position', [100, -100, 900, 1000]);

% --- Top panel: dynamics of the system ---
ax1 = axes('Position', [0.25, 0.62, 0.50, 0.30]); hold on;
alpha = zcontrol_scenario.alpha;
sigma = zcontrol_scenario.sigma;
zeta  = zcontrol_scenario.zeta;

u0 = 0.4; b0 = 0.5; m0 = k * b0; i0 = 1 - u0 - b0;
y0 = [u0, b0, i0, m0];

ode_system = @(t, y) [
    mu - rho * y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
    rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + zeta*(y(3) - id);
    -zeta*(y(3) - id);
    a * k * y(2) - a * y(4)
];

[t, y] = ode45(ode_system, tspan, y0, options);
n_t = y(:,1) + y(:,2) + y(:,3);

plot(ax1, t, y(:,1), '-', 'Color', col_u, 'LineWidth', 2.5);
plot(ax1, t, y(:,2), '-', 'Color', col_b, 'LineWidth', 2.5);
plot(ax1, t, y(:,3), '-', 'Color', col_i, 'LineWidth', 2.5);
plot(ax1, t, y(:,4), '-', 'Color', col_m, 'LineWidth', 2.5);
plot(ax1, t, n_t, '--', 'Color', col_n, 'LineWidth', 2.5);

title(ax1, 'Solutions of the Z-controlled model', 'FontName', 'Times New Roman');
xlabel(ax1, 't'); grid(ax1, 'on'); box(ax1, 'on');
legend(ax1, {'u(t)', 'b(t)', 'i(t)', 'm(t)', 'n(t)'}, ...
       'Location', 'best', 'Orientation', 'horizontal', 'FontSize', 10);
set(ax1, 'FontSize', 20, 'FontName', 'Times New Roman');

% --- Bottom-left panel: forward G(t) ---
ax2 = axes('Position', [0.10, 0.12, 0.35, 0.35]); hold on;
for i = 1:length(forward_scenari)
    s = forward_scenari{i};
    y0 = [0.4, 0.5, 0.1, k*0.5];

    ode_system = @(t, y) [
        mu - rho * y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + s.zeta*(y(3) - id);
        -s.zeta*(y(3) - id);
        a * k * y(2) - a * y(4)
    ];

    [t, y] = ode45(ode_system, tspan, y0, options);
    G_t = compute_Gt_i(t, y, s.alpha, s.sigma, s.zeta);
    plot(ax2, t, G_t, s.style, 'Color', deep_red, 'LineWidth', 2.5, 'DisplayName', s.label);
end
title(ax2, 'Forward scenarios', 'FontName', 'Times New Roman');
xlabel(ax2, 't'); ylabel(ax2, 'G(t)');
legend(ax2, 'Location', 'southeast');
grid(ax2, 'on'); box(ax2, 'on');
ylim(ax2, [0.4, 1]);
set(ax2, 'FontSize', 20, 'FontName', 'Times New Roman');

% --- Bottom-right panel: backward G(t) ---
ax3 = axes('Position', [0.55, 0.12, 0.35, 0.35]); hold on;
for i = 1:length(backward_scenari)
    s = backward_scenari{i};
    y0 = [0.4, 0.5, 0.1, k*0.5];

    ode_system = @(t, y) [
        mu - rho * y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + s.zeta*(y(3) - id);
        -s.zeta*(y(3) - id);
        a * k * y(2) - a * y(4)
    ];

    [t, y] = ode45(ode_system, tspan, y0, options);
    G_t = compute_Gt_i(t, y, s.alpha, s.sigma, s.zeta);
    plot(ax3, t, G_t, s.style, 'Color', deep_red, 'LineWidth', 2.5, 'DisplayName', s.label);
end
title(ax3, 'Backward scenarios', 'FontName', 'Times New Roman');
xlabel(ax3, 't'); ylabel(ax3, 'G(t)');
legend(ax3, 'Location', 'southeast');
grid(ax3, 'on'); box(ax3, 'on');
ylim(ax3, [0.4, 1]);
set(ax3, 'FontSize', 20, 'FontName', 'Times New Roman');

%% === EXPORT FIGURE ===
print(gcf, 'Fig5_zcontrol_on_i.tiff', '-dtiff', '-r300');
print(gcf, 'Fig5_zcontrol_on_i.eps',  '-depsc');
