% ===========================================
% Local Sensitivity Analysis (UPDATED REA ping–pong model, 1:1 DIC–Pi)
% Params: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
% TP is FIXED (NOT a parameter)
%
% Fluxes are computed with the UPDATED REA algebra:
%   delta1 = 1 + Sm/Ks + Mm/Km + Pm/Kp
%   delta2 = 1 + Sims/Ks + Mims/Km + Pims/Kp
%   phi    = (Sm + lambda21*Mm + lambda31*Pm) / (Sims + lambda21*Mims + lambda31*Pims)
%   Den    = delta1 + phi*delta2
% Net: J_succ = v1r_s - v1f_s,  J_mal = v1r_m - v1f_m,  J_pho = J_succ + J_mal
% ===========================================

clear; clc; close all;

% Load experimental competition data
load('Palmier_DIC_Competition_Exp_Data.mat');  % provides M_with_S, etc.

% ---------- Fixed phosphate capacity (constant, not estimated) ----------
TP_FIXED = 1.0;   % choose consistent with your flux unit convention

% ---------- Parameter names (7 parameters) ----------
param_names = { ...
    'T^{s}_{max}', ...     % succinate capacity
    'T^{m}_{max}', ...     % malate capacity
    'K^{s}_{m}',   ...     % succinate KD
    'K^{m}_{m}',   ...     % malate KD
    'K^{p}_{m}',   ...     % phosphate KD
    '\lambda_{21}',...     % lambda2/lambda1
    '\lambda_{31}'};       % lambda3/lambda1

n_params = numel(param_names);

% ---------- Baseline parameter vector ----------
% Order: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
params0 = [70, 70, 1.17, 0.23, 0.93, 1.0, 1.0];

% ---------- Finite-difference perturbation ----------
perturbation = 0.01;  % ±1%

% ---------- Loop settings ----------
n_pts = numel(M_with_S);
LSC_MATRIX = zeros(n_pts, n_params);

% ---------- Sensitivity loop ----------
for j = 1:n_pts   % change to "for j = 1" if you only want the first condition

    % Initial conditions for experiment j (match your dataset/interpretation)
    P_m0 = 2;
    S_m0 = 0.2;
    M_m0 = 0.2;
    P_c0 = 0;
    M_c0 = 1 / M_with_S(j);
    S_c0 = 0.5;

    % State vector: y = [Mm, Mims, Sm, Sims, Pm, Pims]
    y0 = [M_m0, M_c0, S_m0, S_c0, P_m0, P_c0];

    % Compute baseline fluxes
    [J_succ0, J_mal0, J_pho0] = compute_flux_REA_UPDATED(y0, params0, TP_FIXED);

    % Choose which flux to analyze:
    J_base = J_succ0;        % succinate
    % J_base = J_mal0;       % malate
    % J_base = J_pho0;       % phosphate (enforced 1:1)

    if abs(J_base) < 1e-12
        LSC_MATRIX(j,:) = NaN;
        continue;
    end

    LSC_vec = zeros(1, n_params);

    for i = 1:n_params
        p_plus  = params0; p_plus(i)  = params0(i) * (1 + perturbation);
        p_minus = params0; p_minus(i) = params0(i) * (1 - perturbation);

        [Jsp, Jmp, Jpp] = compute_flux_REA_UPDATED(y0, p_plus,  TP_FIXED);
        [Jsm, Jmm, Jpm] = compute_flux_REA_UPDATED(y0, p_minus, TP_FIXED);

        % Match the choice used for J_base:
        J_plus  = Jsp;  J_minus = Jsm;     % succinate
        % J_plus = Jmp; J_minus = Jmm;     % malate
        % J_plus = Jpp; J_minus = Jpm;     % phosphate

        dJdpi = (J_plus - J_minus) / (2 * perturbation * params0(i));
        LSC_vec(i) = (params0(i) / J_base) * dJdpi;
    end

    LSC_MATRIX(j,:) = LSC_vec;
end

mean_LSCs = mean(LSC_MATRIX, 1, 'omitnan');
save('LSCs_succ.mat','mean_LSCs');

% ---------- Plot ----------
figure('Color','w');
barh(mean_LSCs);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names, 'FontSize', 14);
xlabel('Normalized Sensitivity (LSC)', 'FontSize', 16);
ylabel('Parameters', 'FontSize', 16);
title('Mean Local Sensitivity Coefficients (UPDATED REA 1:1 DIC--Pi, TP fixed)', 'FontSize', 18);
grid on;


% ==========================================================
% UPDATED REA flux computation (TP fixed; 1:1 enforced)
% ==========================================================
function [J_succ, J_mal, J_pho] = compute_flux_REA_UPDATED(y, params, TP)

    % States: y = [Mm, Mims, Sm, Sims, Pm, Pims]
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

    % Phi term (single phi)
    denom_phi = (Sims + lambda21*Mims + lambda31*Pims);
    epsPhi = 1e-12;
    if abs(denom_phi) < epsPhi
        denom_phi = epsPhi;
    end
    phi = (Sm + lambda21*Mm + lambda31*Pm) / denom_phi;

    % Common denominator
    Den = (delta1 + phi*delta2);
    epsDen = 1e-12;
    if abs(Den) < epsDen
        Den = epsDen;
    end

    % Forward-cycle rates
    v1f_s = Ts_max * (Sm  - phi*Sims) / (Ks_m * Den);
    v1f_m = Tm_max * (Mm  - phi*Mims) / (Km_m * Den);
    v2f_p = TP     * (phi*Pims - Pm)  / (Kp_m * Den);

    % Reverse-cycle rates
    v1r_s = Ts_max * (phi*Sims - Sm)  / (Ks_m * Den);
    v1r_m = Tm_max * (phi*Mims - Mm)  / (Km_m * Den);
    v2r_p = TP     * (Pm - phi*Pims)  / (Kp_m * Den);

    % Net fluxes
    J_succ = v1r_s - v1f_s;
    J_mal  = v1r_m - v1f_m;

    % Strict 1:1 DIC--Pi exchange (enforced)
    J_pho  = J_succ + J_mal;

    

    
    % v2f_p and v2r_p are computed but not returned; keep if you want diagnostics.
end
