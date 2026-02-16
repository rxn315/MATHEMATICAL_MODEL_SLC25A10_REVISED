% ===========================================
% Local Sensitivity Analysis for J_pho (UPDATED REA model)
% Complex-step derivative (no finite differences)
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
TP_FIXED = 70;   % same flux-capacity units as Ts_max,Tm_max in this REA algebra

% -----------------------------
% Parameters (NO TP)
% -----------------------------
param_names = { ...
    'T^{s}_{max}', 'T^{m}_{max}', ...
    'K^{s}_{m}',   'K^{m}_{m}',   'K^{p}_{m}', ...
    '\lambda_{21}', '\lambda_{31}'};
n_params = numel(param_names);

% Baseline parameters:
% [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
params0 = [70, 70, 1.17, 0.23, 0.93, 1.0, 1.0];

% -----------------------------
% Complex-step settings
% -----------------------------
h = 0.01;   % complex-step size (can be extremely small)

% -----------------------------
% Preallocate
% -----------------------------
n_pts = numel(M_with_S);
LSC_MATRIX = zeros(n_pts, n_params);

% -----------------------------
% Loop over operating points
% -----------------------------
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

    % Baseline phosphate flux
    [~, ~, J_pho_base] = compute_flux_REA_UPDATED_cs(y0, params0, TP_FIXED);

    % Avoid division by ~0 in normalization
    if abs(J_pho_base) < 1e-12
        LSC_MATRIX(j,:) = NaN;
        continue;
    end

    % Complex-step local sensitivity coefficients
    LSC_vec = zeros(1, n_params);

    for i = 1:n_params
        p_cs = params0;
        p_cs(i) = params0(i) + 1i*h;

        [~, ~, J_pho_cs] = compute_flux_REA_UPDATED_cs(y0, p_cs, TP_FIXED);

        % Complex-step derivative
        dJdpi = imag(J_pho_cs) / h;

        % Normalized LSC
        LSC_vec(i) = (params0(i) / J_pho_base) * dJdpi;
    end

    LSC_MATRIX(j,:) = LSC_vec;
end

% -----------------------------
% Summary + save
% -----------------------------
mean_LSCs_pho = mean(LSC_MATRIX, 1, 'omitnan');
save('LSCs_pho_complexstep.mat','mean_LSCs_pho','LSC_MATRIX','params0','TP_FIXED','h','use_inverse_M');

% -----------------------------
% Plot
% -----------------------------
figure('Color','w');
barh(mean_LSCs_pho);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names, 'FontSize', 14);
xlabel('Normalized Sensitivity (LSC) for J_{pho}', 'FontSize', 16);
ylabel('Parameters', 'FontSize', 16);
title('Mean Local Sensitivity Coefficients for Phosphate Flux (UPDATED REA, complex-step, TP fixed)', 'FontSize', 18);
grid on;


% ==========================================================
% UPDATED REA flux computation (complex-step SAFE)
% TP is fixed input.
%
% Notes for complex-step:
%   - Avoid abs(), sign(), and "if abs(x)<..." logic on complex values.
%   - Use an analytic regularization "+ epsReg" (constant) in denominators.
% ==========================================================
function [J_succ, J_mal, J_pho] = compute_flux_REA_UPDATED_cs(y, params, TP)

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

    % Small analytic regularization (keeps holomorphic form)
    epsReg = 1e-30;

    % Phi term (single phi in updated model)
    denom_phi = (Sims + lambda21*Mims + lambda31*Pims) + epsReg;
    phi = (Sm + lambda21*Mm + lambda31*Pm) / denom_phi;

    % Common denominator (same for forward + reverse)
    Den = (delta1 + phi*delta2) + epsReg;

    % Forward-cycle rates
    v1f_s = Ts_max * (Sm  - phi*Sims) / (Ks_m * Den);
    v1f_m = Tm_max * (Mm  - phi*Mims) / (Km_m * Den);
    v2f_p = TP     * (phi*Pims - Pm)  / (Kp_m * Den);

    % Reverse-cycle rates
    v1r_s = Ts_max * (phi*Sims - Sm)  / (Ks_m * Den);
    v1r_m = Tm_max * (phi*Mims - Mm)  / (Km_m * Den);
    v2r_p = TP     * (Pm - phi*Pims)  / (Kp_m * Den);

    % Net fluxes (match your current convention)
    J_succ = v1r_s - v1f_s;
    J_mal  = v1r_m - v1f_m;

    % Keep exactly your convention:
    J_pho  = J_succ + J_mal;

    % If you ever want phosphate net explicitly instead:
    % J_pho = v2r_p - v2f_p;
end
