%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_uncontrolled_model.m
%
% Description:
%   This script generates Figs 2 and 3 of the manuscript:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%
%   It simulates the uncontrolled model under:
%   - Forward scenarios (Fig. 2)
%   - Backward scenarios (Fig. 3)
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

%% === MODEL PARAMETERS AND INITIAL CONDITIONS ===
mu = 0.05; rho = 0.25; lambda = 0.02; p = 0.7;
q = 0.8; a = 0.5; k = 2; gamma = 0.2;

y0 = [0.4, 0.5, 0.1, 1];  % [u0, b0, i0, m0]
T = 200;
tspan = [0 T];
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% Colors for each variable
var_colors = [
    0, 0.4471, 0.7412;       % u(t) - blue
    0.8500, 0.3250, 0.0980;  % b(t) - orange
    0.2, 0.6, 0.2;           % i(t) - green
    0.4940, 0.1840, 0.5560;  % m(t) - purple
    0, 0, 0                  % n(t) - black
];

%% === MODEL DEFINITION ===
referral_model_rhs = @(alpha, sigma) @(t, y) [
    mu - rho * y(2) * y(1) - mu * y(1) - gamma * y(4) * y(1);
    p * rho * y(2) * y(1) - sigma * y(2) + alpha * p * y(2) * y(3) ...
        - mu * y(2) + lambda * y(3) + gamma * q * y(4) * y(1);
    (1 - p) * rho * y(2) * y(1) + sigma * y(2) - alpha * p * y(2) * y(3) ...
        - lambda * y(3) - mu * y(3) + gamma * (1 - q) * y(4) * y(1);
    a * k * y(2) - a * y(4)
];

%% === FORWARD SCENARIOS (Fig. 2) ===
fig1 = figure;
clf(fig1);
set(fig1, 'Position', [100, 100, 500, 450]);
tl1 = tiledlayout(fig1, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

forward_cases = {
    struct('alpha', 0.4, 'sigma', 0.6, 'label', '(Ia)'), 
    struct('alpha', 0.4, 'sigma', 0.8, 'label', '(Ib)')
};

for i = 1:length(forward_cases)
    alpha = forward_cases{i}.alpha;
    sigma = forward_cases{i}.sigma;

    [t, y] = ode45(referral_model_rhs(alpha, sigma), tspan, y0, options);
    n_t = y(:,1) + y(:,2) + y(:,3);

    ax = nexttile(i); hold on;

    h = gobjects(1,5);
    h(1) = plot(t, y(:,1), 'Color', var_colors(1,:), 'LineWidth', 2);
    h(2) = plot(t, y(:,2), 'Color', var_colors(2,:), 'LineWidth', 2);
    h(3) = plot(t, y(:,3), 'Color', var_colors(3,:), 'LineWidth', 2);
    h(4) = plot(t, y(:,4), 'Color', var_colors(4,:), 'LineWidth', 2);
    h(5) = plot(t, n_t , '--', 'Color', var_colors(5,:), 'LineWidth', 2);

    if i == 1
        h_forward = h;
    end

    title(['Uncontrolled model - Forward scenario ', forward_cases{i}.label], 'FontName', 'Times New Roman');
    xlabel('t', 'FontName', 'Times New Roman');
    grid on; box on;
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
end

lgf = legend(h_forward, {'u(t)', 'b(t)', 'i(t)', 'm(t)', 'n(t)'}, ...
    'Orientation', 'horizontal', 'Location', 'southoutside');
lgf.Layout.Tile = 'south';

%% === BACKWARD SCENARIOS (Fig. 3) ===
fig2 = figure;
clf(fig2);
set(fig2, 'Position', [950, 100, 500, 700]);
tl2 = tiledlayout(fig2, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

backward_cases = {
    struct('alpha', 2.0, 'sigma', 0.6,  'label', '(II)'), 
    struct('alpha', 2.0, 'sigma', 0.8,  'label', '(IIa)'), 
    struct('alpha', 2.0, 'sigma', 0.95, 'label', '(IIb)')
};

for i = 1:length(backward_cases)
    alpha = backward_cases{i}.alpha;
    sigma = backward_cases{i}.sigma;

    [t, y] = ode45(referral_model_rhs(alpha, sigma), tspan, y0, options);
    n_t = y(:,1) + y(:,2) + y(:,3);

    ax = nexttile(i); hold on;

    h = gobjects(1,5);
    h(1) = plot(t, y(:,1), 'Color', var_colors(1,:), 'LineWidth', 2);
    h(2) = plot(t, y(:,2), 'Color', var_colors(2,:), 'LineWidth', 2);
    h(3) = plot(t, y(:,3), 'Color', var_colors(3,:), 'LineWidth', 2);
    h(4) = plot(t, y(:,4), 'Color', var_colors(4,:), 'LineWidth', 2);
    h(5) = plot(t, n_t , '--', 'Color', var_colors(5,:), 'LineWidth', 2);

    if i == 1
        h_backward = h;
    end

    title(['Uncontrolled model - Backward scenario ', backward_cases{i}.label], 'FontName', 'Times New Roman');
    xlabel('t', 'FontName', 'Times New Roman');
    grid on; box on;
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
end

lgb = legend(h_backward, {'u(t)', 'b(t)', 'i(t)', 'm(t)', 'n(t)'}, ...
    'Orientation', 'horizontal', 'Location', 'southoutside');
lgb.Layout.Tile = 'south';

%% === SAVE FIGURES ===

% Forward scenario – Fig. 2
print(fig1, 'Fig2_uncontrolled_forward.tiff', '-dtiff', '-r300')
print(fig1, 'Fig2_uncontrolled_forward.eps', '-depsc')

% Backward scenario – Fig. 3
print(fig2, 'Fig3_uncontrolled_backward.tiff', '-dtiff', '-r300')
print(fig2, 'Fig3_uncontrolled_backward.eps', '-depsc')