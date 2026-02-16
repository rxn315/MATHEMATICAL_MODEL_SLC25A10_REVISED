% ===========================================
% Local Sensitivity Analysis for J_pho (UPDATED REA model)
% TP is FIXED (not a parameter)
% Params: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
% ===========================================

clear; clc; close all;

% Load Palmieri competition data
load('Palmier_DIC_Competition_Exp_Data.mat');  % expects M_with_S

% -----------------------------
% Interpret M_with_S
% -----------------------------
use_inverse_M = true;   % true: Mims = 1/M_with_S(j);  false: Mims = M_with_S(j)

% -----------------------------
% Fixed phosphate capacity (NOT estimated)
% -----------------------------
TP_FIXED = 70;   % choose your baseline TP (same units as Ts_max,Tm_max in this REA algebra)

% ---------- Parameter names (7 parameters; TP removed) ----------
param_names = { ...
    'T^{s}_{max}', 'T^{m}_{max}', ...
    'K^{s}_{m}',   'K^{m}_{m}',   'K^{p}_{m}', ...
    '\lambda_{21}', '\lambda_{31}'};
n_params = numel(param_names);

% ---------- Baseline parameters ----------
% Order: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
params0 = [70, 70, 1.17, 0.23, 0.93, 1.0, 1.0];

% ---------- Finite-difference perturbation ----------
perturbation = 0.01;  % ±1%

% ---------- Preallocate ----------
n_pts = numel(M_with_S);
LSC_MATRIX = zeros(n_pts, n_params);

% ---------- Loop over experimental conditions ----------
for j = 1:n_pts

    % Operating point (choose values consistent with your analysis)
    Pm   = 2;
    Sm   = 0.2;
    Mm   = 0.2;
    Pims = 0;
    Sims = 0.5;

    if use_inverse_M
        Mims = 1 / M_with_S(j);
    else
        Mims = M_with_S(j);
    end

    % State vector: [Mm, Mims, Sm, Sims, Pm, Pims]
    y0 = [Mm, Mims, Sm, Sims, Pm, Pims];

    % Compute baseline phosphate flux (UPDATED REA convention)
    [~, ~, J_pho_base] = compute_flux_REA_UPDATED(y0, params0, TP_FIXED);

    % Avoid division by ~0
    if abs(J_pho_base) < 1e-12
        LSC_MATRIX(j,:) = NaN;
        continue;
    end

    % Local Sensitivity Coefficients for J_pho
    LSC_vec = zeros(1, n_params);

    for i = 1:n_params
        % Perturb parameter i by ±1%
        p_plus  = params0; p_plus(i)  = params0(i) * (1 + perturbation);
        p_minus = params0; p_minus(i) = params0(i) * (1 - perturbation);

        % Perturbed phosphate fluxes
        [~, ~, J_pho_plus]  = compute_flux_REA_UPDATED(y0, p_plus,  TP_FIXED);
        [~, ~, J_pho_minus] = compute_flux_REA_UPDATED(y0, p_minus, TP_FIXED);

        % Central finite difference derivative dJ/dp_i
        dJdpi = (J_pho_plus - J_pho_minus) / (2 * perturbation * params0(i));

        % Normalized local sensitivity coefficient
        LSC_vec(i) = (params0(i) / J_pho_base) * dJdpi;
    end

    % Store this condition’s LSCs
    LSC_MATRIX(j, :) = LSC_vec;
end

% ---------- Mean LSCs across conditions ----------
mean_LSCs_pho = mean(LSC_MATRIX, 1, 'omitnan');
save('LSCs_pho.mat','mean_LSCs_pho');

% ---------- Plot ----------
figure('Color','w');
barh(mean_LSCs_pho);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names, 'FontSize', 14);
xlabel('Normalized Sensitivity (LSC) for J_{pho}', 'FontSize', 16);
ylabel('Parameters', 'FontSize', 16);
title('Mean Local Sensitivity Coefficients for Phosphate Flux (UPDATED REA, TP fixed)', 'FontSize', 18);
grid on;


% ==========================================================
% UPDATED REA flux computation (new model form; TP is fixed)
% Returns also forward/reverse phosphate rates if you want them
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

    % --- Saturation polynomials (ALL substrates on each side) ---
    delta1 = 1 + Sm/Ks_m   + Mm/Km_m   + Pm/Kp_m;
    delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

    % --- Phi term (single phi in updated model) ---
    denom_phi = (Sims + lambda21*Mims + lambda31*Pims);
    if abs(denom_phi) < 1e-12
        denom_phi = 1e-12;
    end
    phi1 = (Sm + lambda21*Mm + lambda31*Pm) / denom_phi;

    % --- Common denominator (same for forward + reverse) ---
    Den = (delta1 + phi1*delta2);
    if abs(Den) < 1e-12
        Den = 1e-12;
    end

    % --- Forward-cycle fluxes ---
    v1f_s = Ts_max * (Sm  - phi1*Sims) / (Ks_m * Den);
    v1f_m = Tm_max * (Mm  - phi1*Mims) / (Km_m * Den);
    v2f_p = TP     * (phi1*Pims - Pm)  / (Kp_m * Den);

    % --- Reverse-cycle fluxes ---
    v1r_s = Ts_max * (phi1*Sims - Sm)  / (Ks_m * Den);
    v1r_m = Tm_max * (phi1*Mims - Mm)  / (Km_m * Den);
    v2r_p = TP     * (Pm - phi1*Pims)  / (Kp_m * Den);

    % --- Output conventions (match your updated Palmieri code) ---
    J_succ = v1r_s-v1f_s;
    J_mal  = v1r_m-v1f_m;
    %J_pho  = v2r_p-v2f_p
    J_pho  = J_succ+J_mal;
end
