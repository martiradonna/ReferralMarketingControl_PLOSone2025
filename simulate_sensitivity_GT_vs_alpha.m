%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_sensitivity_GT_vs_alpha.m
%
% Description:
%   This script generates Fig. 6 of the manuscript:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%
%   It computes G(T) at final time for a range of α and:
%   - relative broadcasters increase (%) under Z-control on b(t)
%   - relative inerts reduction (%) under Z-control on i(t)
%
% Author: Angela Martiradonna
% Date: August 2025
%
% License:
%   Creative Commons Attribution 4.0 International License (CC BY 4.0)
%
% Citation:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% === COMMON PARAMETERS ===
mu = 0.05; rho = 0.25; lambda_ = 0.02; p = 0.7;
q = 0.8; a = 0.5; k = 2; gamma = 0.2;
sigma = 0.95; xi = 0.03; zeta = 0.03;
b0 = 0.5; m0 = k * b0; u0 = 0.4; i0 = 1 - u0 - b0;
y0 = [u0, b0, i0, m0];

T = 1e6;
alpha = linspace(0.01, 3, 100)';
alpha_star = 1.1684;

%% === GRID FOR Z-CONTROL ON b(t) ===
b_d_vec = linspace(b0, 1 - mu/(rho + gamma), 50);
incremento = (b_d_vec - b0) / b0 * 100;
GT_b = zeros(length(alpha), length(b_d_vec));

for j = 1:length(b_d_vec)
    b_d = b_d_vec(j);
    [~, y] = solve_ode_b(mu, rho, gamma, xi, a, k, b_d, y0, [0, T]);
    uT = y(end,1); bT = y(end,2); iT = y(end,3); mT = y(end,4);
    term = (xi*(bT - b_d) + lambda_ * iT + gamma*q*mT*uT)/bT;
    GT_b(:,j) = -p*rho*uT - alpha*p*iT + sigma + mu - term;
end

%% === GRID FOR Z-CONTROL ON i(t) ===
i_d_vec = linspace(i0, 0, 50);
decremento = (i0 - i_d_vec) / i0 * 100;
GT_i = zeros(length(alpha), length(i_d_vec));

for j = 1:length(i_d_vec)
    id = i_d_vec(j);
    [~, y] = solve_ode_i(mu, rho, gamma, zeta, a, k, id, y0, [0, T]);
    uT = y(end,1); bT = y(end,2); iT = y(end,3); mT = y(end,4);
    term = (zeta*(iT - id) - (lambda_ + mu)*iT + gamma*(1 - q)*mT*uT)/bT;
    GT_i(:,j) = (1 - p)*rho*uT + sigma - alpha*p*iT + term;
end

%% === FIGURE ===
figure(1); clf;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Left panel: Z-control on b(t)
ax1 = axes('Position', [0.08, 0.15, 0.40, 0.75]);
[X1, Y1] = meshgrid(incremento, alpha);
contourf(ax1, X1, Y1, GT_b, 100, 'LineStyle', 'none');
colormap(ax1, jet); clim(ax1, [0.1 1]);
%colorbar(ax1, 'Position', [0.49, 0.15, 0.015, 0.75]);
xlabel(ax1, 'Relative broadcasters increase (%)', 'FontName', 'Times New Roman');
ylabel(ax1, '\alpha', 'FontName', 'Times New Roman');
title(ax1, 'G(T), Z-control on b(t), T = 10^6', 'FontName', 'Times New Roman');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 20);
yline(ax1, alpha_star, '--k', 'LineWidth', 1.5);

% Adjust y-ticks
yt = yticks(ax1);
if all(abs(yt - alpha_star) > 1e-6)
    yt = sort(unique([yt, alpha_star]));
end
yt_labels = arrayfun(@(x) sprintf('%.2f', x), yt, 'UniformOutput', false);
idx_star = find(abs(yt - alpha_star) < 1e-6);
yt_labels{idx_star} = '\alpha*';
set(ax1, 'YTick', yt, 'YTickLabel', yt_labels);

% Right panel: Z-control on i(t)
ax2 = axes('Position', [0.55, 0.15, 0.40, 0.75]);
[X2, Y2] = meshgrid(decremento, alpha);
contourf(ax2, X2, Y2, GT_i, 100, 'LineStyle', 'none');
colormap(ax2, jet); clim(ax2, [0.1 1]);
colorbar(ax2, 'Position', [0.96, 0.15, 0.015, 0.75]);
xlabel(ax2, 'Relative inerts reduction (%)', 'FontName', 'Times New Roman');
ylabel(ax2, '\alpha', 'FontName', 'Times New Roman');
title(ax2, 'G(T), Z-control on i(t), T = 10^6', 'FontName', 'Times New Roman');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 20);
yline(ax2, alpha_star, '--k', 'LineWidth', 1.5);

yt = yticks(ax2);
if all(abs(yt - alpha_star) > 1e-6)
    yt = sort(unique([yt, alpha_star]));
end
yt_labels = arrayfun(@(x) sprintf('%.2f', x), yt, 'UniformOutput', false);
idx_star = find(abs(yt - alpha_star) < 1e-6);
yt_labels{idx_star} = '\alpha*';
set(ax2, 'YTick', yt, 'YTickLabel', yt_labels);

%% === EXPORT ===
print(gcf, 'Fig6_GT.tiff', '-dtiff', '-r300');
print(gcf, 'Fig6_GT.eps',  '-depsc');

%% === FUNCTIONS ===
function [T, Y] = solve_ode_b(mu, rho, gamma, xi, a, k, b_d, y0, tspan)
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    ode = @(t, y)[
        mu - rho*y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        -xi*(y(2) - b_d);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + xi*(y(2) - b_d);
        a*k*y(2) - a*y(4)
    ];
    [T, Y] = ode45(ode, tspan, y0, options);
end

function [T, Y] = solve_ode_i(mu, rho, gamma, zeta, a, k, id, y0, tspan)
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    ode = @(t, y)[
        mu - rho*y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
        rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + zeta*(y(3) - id);
        -zeta*(y(3) - id);
        a*k*y(2) - a*y(4)
    ];
    [T, Y] = ode45(ode, tspan, y0, options);
end
