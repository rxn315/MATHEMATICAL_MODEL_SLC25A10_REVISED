% ===========================================
% Local Sensitivity Analysis (UPDATED REA ping–pong model, 1:1 DIC–Pi)
% COMPLEX-STEP derivative (no finite differences)
%
% Params: [Ts_max, Tm_max, Ks_m, Km_m, Kp_m, lambda21, lambda31]
% TP is FIXED (NOT a parameter)
%
% UPDATED REA algebra:
%   delta1 = 1 + Sm/Ks + Mm/Km + Pm/Kp
%   delta2 = 1 + Sims/Ks + Mims/Km + Pims/Kp
%   phi    = (Sm + lambda21*Mm + lambda31*Pm) / (Sims + lambda21*Mims + lambda31*Pims)
%   Den    = delta1 + phi*delta2
% Net: J_succ = v1r_s - v1f_s,  J_mal = v1r_m - v1f_m,  J_pho = J_succ + J_mal
% ===========================================

clear; clc; close all;

% Load experimental competition data
load('Palmier_DIC_Competition_Exp_Data.mat');  % provides M_with_S

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

% ---------- Complex-step settings ----------
h = 0.01;  % complex-step size (very small is fine)

% ---------- Loop settings ----------
n_pts = numel(M_with_S);
LSC_MATRIX = zeros(n_pts, n_params);

% ---------- Sensitivity loop ----------
for j = 1:n_pts   % change to "for j = 1" to use only the first condition

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
    [J_succ0, J_mal0, J_pho0] = compute_flux_REA_UPDATED_cs(y0, params0, TP_FIXED);

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
        p_cs = params0;
        p_cs(i) = params0(i) + 1i*h;   % complex-step perturbation

        [J_succ_cs, J_mal_cs, J_pho_cs] = compute_flux_REA_UPDATED_cs(y0, p_cs, TP_FIXED);

        % Pick the perturbed flux consistent with J_base selection
        J_cs = J_succ_cs;       % succinate
        % J_cs = J_mal_cs;      % malate
        % J_cs = J_pho_cs;      % phosphate

        % Complex-step derivative dJ/dp_i
        dJdpi = imag(J_cs) / h;

        % Normalized LSC = (p/J) * dJ/dp
        LSC_vec(i) = (params0(i) / J_base) * dJdpi;
    end

    LSC_MATRIX(j,:) = LSC_vec;
end

mean_LSCs = mean(LSC_MATRIX, 1, 'omitnan');
save('LSCs_succ_complexstep.mat','mean_LSCs','LSC_MATRIX','params0','TP_FIXED','h');

% ---------- Plot ----------
figure('Color','w');
barh(mean_LSCs);
set(gca, 'YTick', 1:n_params, 'YTickLabel', param_names, 'FontSize', 14);
xlabel('Normalized Sensitivity (LSC)', 'FontSize', 16);
ylabel('Parameters', 'FontSize', 16);
title('Mean Local Sensitivity Coefficients (UPDATED REA, complex-step, TP fixed)', 'FontSize', 18);
grid on;


% ==========================================================
% UPDATED REA flux computation (COMPLEX-STEP SAFE)
% TP fixed; 1:1 enforced via J_pho = J_succ + J_mal
%
% Complex-step notes:
%   - Avoid abs(), sign(), and "if abs(x)<..." branching on complex values.
%   - Use analytic regularization "+ epsReg" in denominators.
% ==========================================================
function [J_succ, J_mal, J_pho] = compute_flux_REA_UPDATED_cs(y, params, TP)

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

    % Analytic regularization (constant) to keep holomorphic operations
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

    % Net fluxes
    J_succ = v1r_s - v1f_s;
    J_mal  = v1r_m - v1f_m;

    % Strict 1:1 DIC--Pi exchange (enforced)
    J_pho  = J_succ + J_mal;

     % keep v2f_p,v2r_p for debugging if you want
end
