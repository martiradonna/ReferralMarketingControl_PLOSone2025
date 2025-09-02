%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_G_zcontrol_stochastic.m
%
% Description:
%   This script generates the dynamic plots in Fig. 8 of the manuscript:
%   Lacitignola D., Martiradonna A. (2025, forthcoming). Can we enhance trust in the circular economy through referral marketing control? %   PLOS ONE
%
%   It shows the time evolution of the gain function G(t) under stochastic
%   fluctuations of α(t) in backward scenarios (II, IIa, IIb), comparing:
%   - Z-control on broadcasters b(t)
%   - Z-control on inerts i(t)
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

clear; close all; clc;

%% === PARAMETERS ===
mu = 0.05; rho = 0.25; lambda = 0.02; p = 0.7;
q = 0.8; a = 0.5; k = 2; gamma = 0.2;
xi = 0.03; zeta = 0.03; b_d = 0.7; id = 0.1;

alpha_bar = 2.0; B = 0.1; tau = 10; Nreal = 10;

T = 200; dt = 0.05; t = 0:dt:T; Nt = length(t);

b0 = 0.5; m0 = k * b0; u0 = 0.4; i0 = 1 - u0 - b0;
y0 = [u0; b0; i0; m0];

%% === SYSTEM DYNAMICS ===
rhs_b = @(y) [
    mu - rho*y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
    -xi*(y(2) - b_d);
    rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + xi*(y(2) - b_d);
    a*k*y(2) - a*y(4)
];

rhs_i = @(y) [
    mu - rho*y(2)*y(1) - mu*y(1) - gamma*y(4)*y(1);
    rho*y(2)*y(1) + gamma*y(4)*y(1) - mu*y(2) - mu*y(3) + zeta*(y(3) - id);
    -zeta*(y(3) - id);
    a*k*y(2) - a*y(4)
];

%% === INTEGRATE DYNAMICS (deterministic path for each control) ===
Y_b = zeros(4,Nt); Y_b(:,1) = y0;
for n = 1:Nt-1
    k1 = rhs_b(Y_b(:,n));
    k2 = rhs_b(Y_b(:,n) + dt*k1);
    Y_b(:,n+1) = Y_b(:,n) + dt/2*(k1 + k2);
end

Y_i = zeros(4,Nt); Y_i(:,1) = y0;
for n = 1:Nt-1
    k1 = rhs_i(Y_i(:,n));
    k2 = rhs_i(Y_i(:,n) + dt*k1);
    Y_i(:,n+1) = Y_i(:,n) + dt/2*(k1 + k2);
end

%% === SCENARIOS: backward with different sigma ===
scenarios = {
    struct('label', 'II',  'sigma', 0.6,  'style', '-'), ...
    struct('label', 'IIa', 'sigma', 0.8,  'style', '-.'), ...
    struct('label', 'IIb', 'sigma', 0.95, 'style', ':')
};

purple = [0.5, 0.1, 0.4];

%% === FIGURE: G(t) under stochastic α(t) ===
figure;
set(gcf, 'Position', [100, 100, 1000, 400], 'Color', 'w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% LEFT: Z-control on b(t)
nexttile; hold on;
title('Z-control on b(t) – stochastic', 'FontName','Times New Roman');
xlabel('t'); ylabel('G(t)');
set(gca,'FontSize',20,'FontName','Times New Roman'); grid on; box on;
ylim([0,1])

for s = scenarios
    scen = s{1}; sigma = scen.sigma;

    for r = 1:Nreal
        alpha_vec = alpha_bar * (1 + B * sin(sqrt(2/tau) * cumsum(sqrt(dt)*randn(1,Nt))));
        G_b = -p*rho.*Y_b(1,:) - alpha_vec.*p.*Y_b(3,:) + sigma + mu ...
            - (xi*(Y_b(2,:) - b_d) + lambda*Y_b(3,:) + gamma*q*Y_b(4,:).*Y_b(1,:)) ./ Y_b(2,:);
        plot(t, G_b, 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
    end

    G_b_det = -p*rho.*Y_b(1,:) - alpha_bar*p.*Y_b(3,:) + sigma + mu ...
        - (xi*(Y_b(2,:) - b_d) + lambda*Y_b(3,:) + gamma*q*Y_b(4,:).*Y_b(1,:)) ./ Y_b(2,:);
    plot(t, G_b_det, scen.style, 'Color', purple, 'LineWidth', 2.2, 'DisplayName', scen.label);
end
legend('Location','southeast');

% RIGHT: Z-control on i(t)
nexttile; hold on;
title('Z-control on i(t) – stochastic', 'FontName','Times New Roman');
xlabel('t'); ylabel('G(t)');
set(gca,'FontSize',20,'FontName','Times New Roman'); grid on; box on;
ylim([0,1])

for s = scenarios
    scen = s{1}; sigma = scen.sigma;

    for r = 1:Nreal
        alpha_vec = alpha_bar * (1 + B * sin(sqrt(2/tau) * cumsum(sqrt(dt)*randn(1,Nt))));
        G_i = (1 - p)*rho.*Y_i(1,:) + sigma - alpha_vec.*p.*Y_i(3,:) ...
            + (zeta*(Y_i(3,:) - id) - (lambda + mu)*Y_i(3,:) + gamma*(1 - q)*Y_i(4,:).*Y_i(1,:)) ./ Y_i(2,:);
        plot(t, G_i, 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
    end

    G_i_det = (1 - p)*rho.*Y_i(1,:) + sigma - alpha_bar*p.*Y_i(3,:) ...
        + (zeta*(Y_i(3,:) - id) - (lambda + mu)*Y_i(3,:) + gamma*(1 - q)*Y_i(4,:).*Y_i(1,:)) ./ Y_i(2,:);
    plot(t, G_i_det, scen.style, 'Color', purple, 'LineWidth', 2.2, 'DisplayName', scen.label);
end
legend('Location','southeast');

%% === SAVE FIGURE ===
print(gcf, 'Fig8_zcontrol_G_backward.tiff', '-dtiff', '-r300');
print(gcf, 'Fig8_zcontrol_G_backward.eps',  '-depsc');
