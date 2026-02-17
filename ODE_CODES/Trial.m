clear; clc;

% ===================== Parameters =====================
load('MCMC_Result_Palmier_modes2.mat')

% IMPORTANT:
% This updated model includes Tp_max and uses the NEW REA denominator:
%   den = delta_m + phi*delta_ims
%
% If your posterior_mean ordering is different, adjust indices here.

Ts_max   = modes(1);   % succinate Tmax
Tm_max   = modes(2);   % malate Tmax

Ks_m     = modes(3);   % succinate Km (mM)
Km_m     = modes(4);   % malate Km (mM)
Kp_m     = modes(5);   % phosphate Km (mM)

Tp_max   = 0;
lambda21 = modes(6);
lambda31 = modes(7);

% ---------------- Volumes per Grame of Protein ----------------
Vm_g   = 7*10^-4; % L/g
Vims_g = 0.50; % L/g

% ===================== Data / ICs =====================
load('Palmier_DIC_competition_Exp_Data.mat'); % expects M_with_S
assert(exist('M_with_S','var')==1, 'M_with_S not found in MAT file.');

P_m0 = 10;
S_m0 = 0.2;
M_m0 = 0.2;

Pims0 = 10^-12;
Mims0 = (1./M_with_S(1));
Sims0 = 0.9;

y0 = [M_m0, Mims0, S_m0, Sims0, P_m0, Pims0]; % [Mm, Mims, Sm, Sims, Pm, Pims]

% ===================== Time span =====================
tspan = [0 30];  % minutes (your axis labels)

% ===================== ODE solve (UPDATED MODEL) =====================
odefun = @(t,y) slc25a10_model_updated_REA( ...
    t, y, Vims_g, Vm_g, Ks_m, Km_m, Kp_m, Ts_max, Tm_max, Tp_max, lambda21, lambda31);

opts = odeset('RelTol',1e-12,'AbsTol',1e-16);
[Time, y] = ode15s(odefun, tspan, y0, opts);