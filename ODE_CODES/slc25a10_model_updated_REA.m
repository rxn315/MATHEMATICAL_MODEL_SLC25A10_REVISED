%% =====================================================================
% Reduced-parameter intact-mitochondria REA model for SLC25A10.
%
% Parameters are the reduced Bayesian-inference parameters:
%   [Ts_max, Tm_max, Tp_max, K_D_s, K_D_m, K_D_p]
%
% The coupling ratios are not independent ODE inputs. They are derived as:
%   lambda_21 = (Tm_max*K_D_s)/(Ts_max*K_D_m)
%   lambda_31 = (Tp_max*K_D_s)/(Ts_max*K_D_p)
%
% Phosphate flux is constrained by electroneutral DIC-Pi exchange:
%   J_pho = J_suc + J_mal
%% =====================================================================
function dYdt = slc25a10_model_updated_REA( ...
    t, y, Vims_g, Vm_g, Ks_m, Km_m, Kp_m, Ts_max, Tm_max, Tp_max) %#ok<INUSD>
% y = [Mm; Mims; Sm; Sims; Pm; Pims]

    eps_safe = 1.0e-12;

    Mm   = y(1);  Mims = y(2);
    Sm   = y(3);  Sims = y(4);
    Pm   = y(5);  Pims = y(6);

    lambda21 = (Tm_max * Ks_m) / (Ts_max * Km_m);
    lambda31 = (Tp_max * Ks_m) / (Ts_max * Kp_m);

    delta_m   = 1 + Sm/Ks_m   + Mm/Km_m   + Pm/Kp_m;
    delta_ims = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

    phi_num = Sm   + lambda21*Mm   + lambda31*Pm;
    phi_den = Sims + lambda21*Mims + lambda31*Pims;
    if abs(phi_den) < eps_safe
        if phi_den < 0
            phi_den = -eps_safe;
        else
            phi_den = eps_safe;
        end
    end

    phi = phi_num / phi_den;
    den = delta_m + phi*delta_ims;
    if abs(den) < eps_safe
        if den < 0
            den = -eps_safe;
        else
            den = eps_safe;
        end
    end

    v1f_s = Ts_max * (Sm  - phi*Sims) / (Ks_m * den);
    v1f_m = Tm_max * (Mm  - phi*Mims) / (Km_m * den);
    v1r_s = Ts_max * (phi*Sims - Sm)  / (Ks_m * den);
    v1r_m = Tm_max * (phi*Mims - Mm)  / (Km_m * den);

    % Ts/Tm/Tp are in umol/min/g. Divide by 60 and multiply by 1e-3 below
    % to obtain mM/s after volume normalization.
    J_suc = (v1r_s - v1f_s) / 60;
    J_mal = (v1r_m - v1f_m) / 60;
    J_pho = J_suc + J_mal;

    dMm_dt   =  (1/Vm_g)   * J_mal * 1.0e-3;
    dMims_dt = -(1/Vims_g) * J_mal * 1.0e-3;

    dSm_dt   =  (1/Vm_g)   * J_suc * 1.0e-3;
    dSims_dt = -(1/Vims_g) * J_suc * 1.0e-3;

    dPm_dt   = -(1/Vm_g)   * J_pho * 1.0e-3;
    dPims_dt =  (1/Vims_g) * J_pho * 1.0e-3;

    dYdt = [dMm_dt; dMims_dt; dSm_dt; dSims_dt; dPm_dt; dPims_dt];
end
