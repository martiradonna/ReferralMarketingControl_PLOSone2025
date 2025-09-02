%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_zcontrolled_model_on_b.m
%
% Description:
%   This script generates Fig. 4 of the manuscript:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%
%   It includes:
%   - Z-control applied to b(t)
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

clear all; close all;

%% === PARAMETERS ===
mu = 0.05; rho = 0.25; lambda = 0.02; p = 0.7;
q = 0.8; a = 0.5; k = 2; gamma = 0.2;
b0 = 0.5; m0 = k * b0; u0 = 0.4;
i0 = 1 - u0 - b0;
y0 = [u0, b0, i0, m0];

T = 200;
tspan = [0 T];
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% Colors
col_u = [0, 0.4471, 0.7412];
col_b = [0.8500, 0.3250, 0.0980];
col_i = [0.4660, 0.6740, 0.1880];
col_m = [0.4940, 0.1840, 0.5560];
col_n = [0, 0, 0];
purple = [0.6350 0.0780 0.1840];

%% === Z-CONTROLLED SCENARIOS ===
scenari = {
    struct('label', 'Ia',  'alpha', 0.4, 'sigma', 0.6,  'b_d', 0.7, 'xi', 0.03),  % Forward
    struct('label', 'IIb', 'alpha', 2.0, 'sigma', 0.95, 'b_d', 0.7, 'xi', 0.03)  % Backward
};

forward_scenari = {
    struct('label', 'Ia', 'alpha', 0.4, 'sigma', 0.6, 'b_d', 0.7, 'xi', 0.03, 'style', '-'),
    struct('label', 'Ib', 'alpha', 0.4, 'sigma', 0.8, 'b_d', 0.7, 'xi', 0.03, 'style', '-.')
};
backward_scenari = {
    struct('label', 'II',  'alpha', 2.0, 'sigma', 0.6,  'b_d', 0.7, 'xi', 0.03, 'style', '-'),
    struct('label', 'IIa', 'alpha', 2.0, 'sigma', 0.8,  'b_d', 0.7, 'xi', 0.03, 'style', '-.'),
    struct('label', 'IIb', 'alpha', 2.0, 'sigma', 0.95, 'b_d', 0.7, 'xi', 0.03, 'style', ':')
};

compute_Gt = @(t, y, alpha, sigma, b_d, xi) ...
    -p * rho .* y(:,1) - alpha * p .* y(:,3) + sigma + mu ...
    - (xi * (y(:,2) - b_d) + lambda * y(:,3) + gamma * q .* y(:,4) .* y(:,1)) ./ y(:,2);

%% === FIGURE ===
figure(1); clf;
set(gcf, 'Position', [100, -100, 900, 1000]);

% --- Top panel: dynamics ---
ax1 = axes('Position', [0.25, 0.62, 0.50, 0.30]); hold on;
for i = 1:2
    s = scenari{i};
    ode_system = @(t, y) [
        mu - rho * y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        -s.xi*(y(2) - s.b_d);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + s.xi*(y(2) - s.b_d);
        a * k * y(2) - a * y(4)
    ];
    [t, y] = ode45(ode_system, tspan, y0, options);
    n_t = y(:,1) + y(:,2) + y(:,3);
    plot(ax1, t, y(:,1), '-', 'Color', col_u, 'LineWidth', 2.5);
    plot(ax1, t, y(:,2), '-', 'Color', col_b, 'LineWidth', 2.5);
    plot(ax1, t, y(:,3), '-', 'Color', col_i, 'LineWidth', 2.5);
    plot(ax1, t, y(:,4), '-', 'Color', col_m, 'LineWidth', 2.5);
    plot(ax1, t, n_t , '--', 'Color', col_n, 'LineWidth', 2.5);
end
title(ax1, 'Solutions of the Z-controlled model', 'FontName', 'Times New Roman');
xlabel(ax1, 't'); grid(ax1, 'on'); box(ax1, 'on');
legend(ax1, {'u(t)', 'b(t)', 'i(t)', 'm(t)', 'n(t)'}, ...
       'Location', 'best', 'Orientation', 'horizontal', 'FontSize', 10);
set(ax1, 'FontSize', 20, 'FontName', 'Times New Roman');

% --- Bottom-left panel: forward G(t) ---
ax2 = axes('Position', [0.10, 0.12, 0.35, 0.35]); hold on;
for i = 1:length(forward_scenari)
    s = forward_scenari{i};
    ode_system = @(t, y) [
        mu - rho * y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        -s.xi*(y(2) - s.b_d);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + s.xi*(y(2) - s.b_d);
        a * k * y(2) - a * y(4)
    ];
    [t, y] = ode45(ode_system, tspan, y0, options);
    G_t = compute_Gt(t, y, s.alpha, s.sigma, s.b_d, s.xi);
    plot(ax2, t, G_t, s.style, 'Color', purple, 'LineWidth', 2.5, 'DisplayName', s.label);
end
title(ax2, 'Forward scenarios', 'FontName', 'Times New Roman');
xlabel(ax2, 't'); ylabel(ax2, 'G(t)');
legend(ax2, 'Location', 'southeast');
grid(ax2, 'on'); box(ax2, 'on');
ylim(ax2, [0.1, 0.8]);
set(ax2, 'FontSize', 20, 'FontName', 'Times New Roman');

% --- Bottom-right panel: backward G(t) ---
ax3 = axes('Position', [0.55, 0.12, 0.35, 0.35]); hold on;
for i = 1:length(backward_scenari)
    s = backward_scenari{i};
    ode_system = @(t, y) [
        mu - rho * y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        -s.xi*(y(2) - s.b_d);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + s.xi*(y(2) - s.b_d);
        a * k * y(2) - a * y(4)
    ];
    [t, y] = ode45(ode_system, tspan, y0, options);
    G_t = compute_Gt(t, y, s.alpha, s.sigma, s.b_d, s.xi);
    plot(ax3, t, G_t, s.style, 'Color', purple, 'LineWidth', 2.5, 'DisplayName', s.label);
end
title(ax3, 'Backward scenarios', 'FontName', 'Times New Roman');
xlabel(ax3, 't'); ylabel(ax3, 'G(t)');
legend(ax3, 'Location', 'southeast');
grid(ax3, 'on'); box(ax3, 'on');
ylim(ax3, [0.1, 0.8]);
set(ax3, 'FontSize', 20, 'FontName', 'Times New Roman');

%% === EXPORT ===
print(gcf, 'Fig4_zcontrol_on_b.tiff', '-dtiff', '-r300');
print(gcf, 'Fig4_zcontrol_on_b.eps',  '-depsc');
