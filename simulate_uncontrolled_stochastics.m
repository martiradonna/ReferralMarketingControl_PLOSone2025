%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_uncontrolled_stochastics.m
%
% Description:
%   This script generates the dynamic plots in Fig. 7 of the manuscript:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%
%   It shows the time evolution of the state variables u(t), b(t), i(t), m(t),
%   and the total population n(t), under stochastic fluctuations of α(t),
%   in two scenarios (forward and backward) for the uncontrolled model.
%
% License:
%   This code is released under the Creative Commons Attribution 4.0
%   International License (CC BY 4.0).
%   You are free to share and adapt the material
%   under the condition that appropriate credit is given.
%
% Citation:
%   If you use this script or any part of the code, please cite the original article:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE.
%
% Author: Angela Martiradonna
% Date: August 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% === PARAMETERS ===
mu = 0.05; rho = 0.25; lambda = 0.02; p = 0.7;
q = 0.8; a = 0.5; k = 2; gamma = 0.2;

%% === TIME DISCRETIZATION ===
T = 200; dt = 0.05;
t = 0:dt:T; Nt = length(t);

%% === INITIAL CONDITIONS ===
y0 = [0.4; 0.5; 0.1; 1];  % [u; b; i; m]

%% === SINE-WIENER PARAMETERS ===
B = 0.1;
tau = 10;
Nreal = 10;  % Number of stochastic realizations

%% === COLORS ===
col_u = [0.00, 0.45, 0.74];
col_b = [0.85, 0.33, 0.10];
col_i = [0.20, 0.60, 0.20];
col_m = [0.49, 0.18, 0.56];
col_n = 'k';
col_det = [0, 0, 0];

%% === DYNAMICS FUNCTION ===
rhs = @(y, alpha, sigma) [
    mu - rho*y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
    p*rho*y(2)*y(1) - sigma*y(2) + alpha*p*y(2)*y(3) - mu*y(2) + lambda*y(3) + gamma*q*y(4)*y(1);
    (1-p)*rho*y(2)*y(1) + sigma*y(2) - alpha*p*y(2)*y(3) - lambda*y(3) - mu*y(3) + gamma*(1-q)*y(4)*y(1);
    a*k*y(2) - a*y(4)
];

%% === SCENARIOS (Forward / Backward) ===
scenarios = {
    struct('alpha_bar', 0.4, 'sigma', 0.6, 'label', 'Uncontrolled model – Forward (Ia) – stochastic'), ...
    struct('alpha_bar', 2.0, 'sigma', 0.6, 'label', 'Uncontrolled model – Backward (II) – stochastic')
};

%% === FIGURE ===
figure();

for i = 1:2
    alpha_bar = scenarios{i}.alpha_bar;
    sigma = scenarios{i}.sigma;
    label = scenarios{i}.label;

    idx_alpha = 2*i - 1;
    idx_dynamics = 2*i;

    subplot(2,2,idx_alpha); hold on; grid on;
    title('\alpha(t) realizations');
    xlabel('t'); ylabel('\alpha(t)');
    set(gca, 'FontSize', 24, 'FontName', 'Times New Roman');
    box on;

    subplot(2,2,idx_dynamics); hold on; grid on;
    title(label);
    xlabel('t');
    set(gca, 'FontSize', 24, 'FontName', 'Times New Roman');
    box on;

    for r = 1:Nreal
        % Generate alpha(t) with sine-Wiener process
        dW = sqrt(dt) * randn(1,Nt);
        W = cumsum(dW);
        alpha_vec = alpha_bar * (1 + B * sin(sqrt(2/tau) * W));

        % Simulate system with RK2
        Y = zeros(4, Nt); Y(:,1) = y0;
        for n = 1:Nt-1
            a1 = alpha_vec(n); a2 = alpha_vec(n+1);
            k1 = rhs(Y(:,n), a1, sigma);
            k2 = rhs(Y(:,n) + dt*k1, a2, sigma);
            Y(:,n+1) = Y(:,n) + dt/2 * (k1 + k2);
        end

        % Plot alpha(t)
        subplot(2,2,idx_alpha)
        plot(t, alpha_vec, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

        % Plot state variables
        subplot(2,2,idx_dynamics)
        h1 = plot(t, Y(1,:), 'Color', col_u, 'LineWidth', 1); % u
        h2 = plot(t, Y(2,:), 'Color', col_b, 'LineWidth', 1); % b
        h3 = plot(t, Y(3,:), 'Color', col_i, 'LineWidth', 1); % i
        h4 = plot(t, Y(4,:), 'Color', col_m, 'LineWidth', 1); % m
        h5 = plot(t, sum(Y(1:3,:)), '--', 'Color', col_n, 'LineWidth', 1.5); % n = u + b + i
    end

    % Add reference line on alpha(t)
    subplot(2,2,idx_alpha)
    yline(alpha_bar, 'k--', 'LineWidth', 2);
    ylim([0,3]);

    % Deterministic dynamics (optional)
    Ydet = zeros(4, Nt); Ydet(:,1) = y0;
    for n = 1:Nt-1
        k1 = rhs(Ydet(:,n), alpha_bar, sigma);
        k2 = rhs(Ydet(:,n) + dt*k1, alpha_bar, sigma);
        Ydet(:,n+1) = Ydet(:,n) + dt/2 * (k1 + k2);
    end

    % Optional overlay of deterministic curves (commented out)
    % subplot(2,2,idx_dynamics)
    % plot(t, Ydet(1,:), '-',  'Color', col_det, 'LineWidth', 2.5); % u
    % plot(t, Ydet(2,:), '-',  'Color', col_det, 'LineWidth', 2.5); % b
    % plot(t, Ydet(3,:), '--', 'Color', col_det, 'LineWidth', 2.5); % i
    % plot(t, Ydet(4,:), ':',  'Color', col_det, 'LineWidth', 2.5); % m

    % Legend
    legend([h1 h2 h3 h4 h5], {'u(t)', 'b(t)', 'i(t)', 'm(t)', 'n(t)'}, ...
        'Orientation', 'horizontal', 'Location', 'best');
end

%% === SAVE FIGURE ===
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % fullscreen
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'Fig7_uncontrolled_stochastics.tiff', '-dtiff', '-r300');
print(gcf, 'Fig7_uncontrolled_stochastics.eps',  '-depsc');
