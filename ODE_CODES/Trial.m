%% =====================================================================
% Intact-mitochondria ODE trial using reduced Bayesian parameters.
%
% This script loads the reduced uniform-prior posterior result when it is
% available. If no result file is found, it falls back to the posterior means
% reported in the revised intact-mitochondria table.
%% =====================================================================

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
repo_dir = fileparts(script_dir);

theta_point = load_mito_reduced_point(repo_dir, script_dir);

Ts_max = theta_point(1);
Tm_max = theta_point(2);
Tp_max = theta_point(3);
Ks_m   = theta_point(4);
Km_m   = theta_point(5);
Kp_m   = theta_point(6);

lambda21 = (Tm_max * Ks_m) / (Ts_max * Km_m);
lambda31 = (Tp_max * Ks_m) / (Ts_max * Kp_m);

fprintf('Using reduced intact-mitochondria parameters:\n');
fprintf('Ts=%.6g, Tm=%.6g, Tp=%.6g, Ks=%.6g, Km=%.6g, Kp=%.6g\n', ...
    Ts_max, Tm_max, Tp_max, Ks_m, Km_m, Kp_m);
fprintf('Derived lambda21=%.6g, lambda31=%.6g\n', lambda21, lambda31);

% Volumes per gram of protein.
Vm_g   = 7.0e-4;  % L/g
Vims_g = 0.49;    % L/g

data_file = fullfile(script_dir, 'Palmier_DIC_Competition_Exp_Data.mat');
if exist(data_file, 'file') ~= 2
    error('Missing data file: %s', data_file);
end
data = load(data_file);
if ~isfield(data, 'M_with_S')
    error('Palmier_DIC_Competition_Exp_Data.mat must contain M_with_S.');
end

Mm0   = 0.20;
Mims0 = 1 ./ data.M_with_S(1);
Sm0   = 0.20;
Sims0 = 0.50;
Pm0   = 2.0;
Pims0 = 0.0;

y0 = [Mm0; Mims0; Sm0; Sims0; Pm0; Pims0];
tspan_sec = [0, 50];

odefun = @(t, y) slc25a10_model_updated_REA( ...
    t, y, Vims_g, Vm_g, Ks_m, Km_m, Kp_m, Ts_max, Tm_max, Tp_max);

opts = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-13);
[Time_sec, Y] = ode15s(odefun, tspan_sec, y0, opts);

figure('Color', 'w', 'Name', 'Reduced intact-mitochondria ODE trial');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

labels = {'M_m', 'M_c', 'S_m', 'S_c', 'P_m', 'P_c'};
for k = 1:6
    nexttile;
    plot(Time_sec, Y(:, k), 'LineWidth', 1.8);
    xlabel('Time (s)');
    ylabel([labels{k}, ' (mM)']);
    grid on;
end

%% =============================== Local functions ==========================
function theta_point = load_mito_reduced_point(repo_dir, script_dir)
    candidates = { ...
        fullfile(script_dir, 'MCMC_Result_Palmier2_reduced_uniform.mat'), ...
        fullfile(repo_dir, 'Bayesian_MCMC_MH_CODES', 'Mito_calibration', 'MCMC_Result_Palmier2_reduced_uniform.mat'), ...
        fullfile(repo_dir, 'FIGURE_4', 'reduced_parameter_inference', 'MCMC_Result_Palmier2_reduced_uniform.mat'), ...
        fullfile(repo_dir, 'FIGURE_4', 'reduced_parameter_inference', 'New Folder', 'MCMC_Result_Palmier2_reduced_uniform.mat')};

    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2
            result = load(candidates{k});
            if isfield(result, 'posterior_mean') && numel(result.posterior_mean) >= 6
                theta_point = result.posterior_mean(1:6);
                return;
            elseif isfield(result, 'theta_samples') && size(result.theta_samples, 2) >= 6
                theta_point = mean(result.theta_samples(:, 1:6), 1);
                return;
            end
        end
    end

    warning('Reduced mitochondrial posterior result was not found. Using table posterior means.');
    theta_point = [145.07, 212.85, 285.45, 2.4242, 0.37291, 2.3154];
end
