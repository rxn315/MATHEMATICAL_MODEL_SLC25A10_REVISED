% ===========================================
% Local Sensitivity Analysis for J_mal (UPDATED REA model)
% Complex-step derivative (no finite differences)
% TP is FIXED (NOT a parameter)
% Params: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
% ===========================================
clear; clc; close all;

% Load Palmier competition data
load('Palmier_DIC_Competition_Exp_Data.mat');

% -----------------------------
% Fixed phosphate capacity (constant, not estimated)
% -----------------------------
TP_FIXED = 1.0;   % choose consistent with your REA flux units

% -----------------------------
% Parameters (NO TP)
% -----------------------------
param_names = { ...
    'T^{s}_{max}', 'T^{m}_{max}', ...
    'K^{s}_{m}',   'K^{m}_{m}',   'K^{p}_{m}', ...
    '\lambda_{21}','\lambda_{31}'};
n_params = numel(param_names);

% Baseline parameters:
% [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
params0 = [70, 70, 1.17, 0.23, 0.93, 1.0, 1.0];

% -----------------------------
% Complex-step settings
% -----------------------------
h = 0.01;   % complex-step size (can be tiny; no subtractive cancellation)

% How many points to evaluate
n_pts = numel(M_with_S);

% Store sensitivity coefficients at each data point
LSC_MATRIX = zeros(n_pts, n_params);

% -----------------------------
% Sensitivity loop over data points
% -----------------------------
for j = 1:n_pts

    % Experimental conditions (your current choice)
    Pm   = 2;   Sm   = 0.2;   Mm   = 0.2;
    Pims = 0;
    Sims = 0.5;
    Mims = 1 / M_with_S(j);   % keep your original mapping

    y0 = [Mm, Mims, Sm, Sims, Pm, Pims];

    % Baseline flux (net, using your compute_flux_REA_UPDATED)
    [~, J_mal_base, ~] = compute_flux_REA_UPDATED(y0, params0, TP_FIXED);

    % Avoid division by ~0 in normalization
    if abs(J_mal_base) < 1e-12
        LSC_MATRIX(j,:) = NaN;
        continue;
    end

    % Local Sensitivity Coefficients at this point
    LSC_vec = zeros(1, n_params);

    for i = 1:n_params
        % Complex-step perturbation (additive, not multiplicative)
        p_cs = params0;
        p_cs(i) = params0(i) + 1i*h;

        % Complex-evaluated flux
        [~, J_mal_cs, ~] = compute_flux_REA_UPDATED(y0, p_cs, TP_FIXED);

        % Complex-step derivative: dJ/dp_i = imag(J(p+i h))/h
        dJdpi = imag(J_mal_cs) / h;

        % Normalized LSC = (p/J) * dJ/dp
        LSC_vec(i) = (params0(i) / J_mal_base) * dJdpi;
    end

    LSC_MATRIX(j,:) = LSC_vec;
end

% Mean LSC across points (ignore NaNs)
mean_LSCs_mal = mean(LSC_MATRIX, 1, 'omitnan');

save('LSCs_mal_complexstep.mat','mean_LSCs_mal','LSC_MATRIX','params0','TP_FIXED','h');

% -----------------------------
% Plot
% -----------------------------
figure('Color','w');
barh(mean_LSCs_mal);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names, 'FontSize', 14);
xlabel('Normalized Sensitivity (LSC) for J_{mal}', 'FontSize', 16);
ylabel('Parameters', 'FontSize', 16);
title('Mean Local Sensitivity Coefficients for Malate Flux (UPDATED REA, complex-step, TP fixed)', 'FontSize', 18);
grid on;


% ==========================================================
% UPDATED REA flux computation (TP is FIXED input)
% Complex-step SAFE version:
%   - No abs()/sign() branches in denominators
%   - Small analytic regularization added as "+ eps" (constant)
%
% Uses UPDATED denominator + single phi:
%   delta1 = 1 + Sm/Ks + Mm/Km + Pm/Kp
%   delta2 = 1 + Sims/Ks + Mims/Km + Pims/Kp
%   phi    = (Sm + lambda21*Mm + lambda31*Pm) / (Sims + lambda21*Mims + lambda31*Pims)
%   Den    = delta1 + phi*delta2
% ==========================================================
function [J_succ, J_mal, J_pho] = compute_flux_REA_UPDATED(y, params, TP)

    % States: [Mm, Mims, Sm, Sims, Pm, Pims]
    Mm   = y(1);  Mims = y(2);
    Sm   = y(3);  Sims = y(4);
    Pm   = y(5);  Pims = y(6);

    % Parameters
    Ts_max   = params(1);
    Tm_max   = params(2);
    Ks_m     = params(3);
    Km_m     = params(4);
    Kp_m     = params(5);
    lambda21 = params(6);
    lambda31 = params(7);

    % Saturation polynomials (ALL substrates on each side)
    delta1 = 1 + Sm./Ks_m   + Mm./Km_m   + Pm./Kp_m;
    delta2 = 1 + Sims./Ks_m + Mims./Km_m + Pims./Kp_m;

    % Small analytic regularization (constant, keeps holomorphic form)
    epsReg = 1e-30;

    % Phi term (single phi)
    denom_phi = (Sims + lambda21*Mims + lambda31*Pims) + epsReg;
    phi = (Sm + lambda21*Mm + lambda31*Pm) / denom_phi;

    % Common denominator
    Den = (delta1 + phi*delta2) + epsReg;

    % Forward-cycle rates
    v1f_s = Ts_max * (Sm  - phi*Sims) / (Ks_m * Den);
    v1f_m = Tm_max * (Mm  - phi*Mims) / (Km_m * Den);
    v2f_p = TP     * (phi*Pims - Pm)  / (Kp_m * Den);

    % Reverse-cycle rates
    v1r_s = Ts_max * (phi*Sims - Sm)  / (Ks_m * Den);
    v1r_m = Tm_max * (phi*Mims - Mm)  / (Km_m * Den);
    v2r_p = TP     * (Pm - phi*Pims)  / (Kp_m * Den);

    % Net fluxes (as in your current code)
    J_succ = v1r_s - v1f_s;
    J_mal  = v1r_m - v1f_m;
    J_pho  = J_succ + J_mal;

    % If you ever want phosphate net explicitly:
    % J_pho = v2r_p - v2f_p;
end
