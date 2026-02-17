%% =====================================================================
% Local function: UPDATED REA MODEL (NEW)
% =====================================================================
function dYdt = slc25a10_model_updated_REA( ...
    t, y, Vims_g, Vm_g, Ks_m, Km_m, Kp_m, Ts_max, Tm_max, Tp_max, lambda21, lambda31)
% y = [Mm; Mims; Sm; Sims; Pm; Pims]
    
    eps_safe = 1e-12;

    Mm   = y(1);  Mims = y(2);
    Sm   = y(3);  Sims = y(4);
    Pm   = y(5);  Pims = y(6);

    % Saturation polynomials (ALL substrates on each side)
    delta_m   = 1 + Sm/Ks_m   + Mm/Km_m   + Pm/Kp_m;
    delta_ims = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

    % phi (single; NO phi2; NO internal volume scaling inside phi)
    phi_num = Sm   + lambda21*Mm   + lambda31*Pm;
    phi_den = Sims + lambda21*Mims + lambda31*Pims;

    phi = phi_num / (phi_den + eps_safe);
    den = delta_m + phi*delta_ims + eps_safe;
    % if abs(phi_den) < eps_safe
    %     phi_den = sign(phi_den + eps_safe) * eps_safe;
    % end
    % phi = phi_num / phi_den;
    % 
    % % Single denominator
    % den = delta_m + phi * delta_ims;
    % if abs(den) < eps_safe
    %     den = sign(den + eps_safe) * eps_safe;
    % end

    % Forward / reverse cycle rates
    v1f_s = Ts_max * (Sm  - phi*Sims) / (Ks_m * den);
    v1f_m = Tm_max * (Mm  - phi*Mims) / (Km_m * den);
    %v2f_p = Tp_max * (phi*Pims - Pm)  / (Kp_m * den);

    v1r_s = Ts_max * (phi*Sims - Sm)  / (Ks_m * den);
    v1r_m = Tm_max * (phi*Mims - Mm)  / (Km_m * den);
    %v2r_p = Tp_max * (Pm - phi*Pims)  / (Kp_m * den);

    % Net fluxes
    J_succ = (v1r_s - v1f_s)/60;
    J_mal  = (v1r_m - v1f_m)/60;
    J_pho  = J_mal+J_succ;

    % ODEs (volume-aware)
    dMm_dt   =  (1/Vm_g)   * J_mal * 10^-3;
    dMims_dt = -(1/Vims_g) * J_mal * 10^-3;

    dSm_dt   =  (1/Vm_g)   * J_succ * 10^-3;
    dSims_dt = -(1/Vims_g) * J_succ * 10^-3;

    dPm_dt   = -(1/Vm_g)   * J_pho * 10^-3;
    dPims_dt =  (1/Vims_g) * J_pho * 10^-3;


    dYdt = [dMm_dt; dMims_dt; dSm_dt; dSims_dt; dPm_dt; dPims_dt];
end
