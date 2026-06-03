%% =====================================================================
% Reduced-parameter proteoliposome REA model for SLC25A10.
%
% Parameters are the reduced Bayesian-inference parameters:
%   [Tmax_m, Tmax_p, K_D_m, K_D_p]
%
% Succinate is absent in the proteoliposome datasets used here. The
% identifiable coupling ratio is derived as:
%   rho_pm = lambda_31/lambda_21
%          = (Tmax_p*K_D_m)/(Tmax_m*K_D_p)
%
% This model uses the direct J_pho form used in the reduced Figure 5
% proteoliposome calibration.
%% =====================================================================
function dYdt = slc25a10_model_liposome_reduced_Jpho( ...
    t, y, Vext_g, Vint_g, K_D_m, K_D_p, Tmax_m, Tmax_p) %#ok<INUSD>
% y = [Mm; Mext; Pm; Pext]

    eps_safe = 1.0e-12;

    Mm   = y(1);  Mext = y(2);
    Pm   = y(3);  Pext = y(4);

    rho_pm = (Tmax_p * K_D_m) / (Tmax_m * K_D_p);

    delta_int = 1 + Mm/K_D_m   + Pm/K_D_p;
    delta_ext = 1 + Mext/K_D_m + Pext/K_D_p;

    phi_den = Mext + rho_pm*Pext;
    if abs(phi_den) < eps_safe
        if phi_den < 0
            phi_den = -eps_safe;
        else
            phi_den = eps_safe;
        end
    end

    phi = (Mm + rho_pm*Pm) / phi_den;
    den = delta_int + phi*delta_ext;
    if abs(den) < eps_safe
        if den < 0
            den = -eps_safe;
        else
            den = eps_safe;
        end
    end

    J_mal = Tmax_m * (phi*Mext - Mm) / (K_D_m * den);
    J_pho = Tmax_p * (phi*Pext - Pm) / (K_D_p * den);

    % Fluxes are in mmol/min/g and volumes in L/g, so dC/dt is mM/min.
    dMm_dt   =  (1/Vint_g) * J_mal;
    dMext_dt = -(1/Vext_g) * J_mal;
    dPm_dt   =  (1/Vint_g) * J_pho;
    dPext_dt = -(1/Vext_g) * J_pho;

    dYdt = [dMm_dt; dMext_dt; dPm_dt; dPext_dt];
end
