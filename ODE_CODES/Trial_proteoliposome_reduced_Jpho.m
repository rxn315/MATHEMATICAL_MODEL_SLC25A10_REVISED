%% =====================================================================
% Proteoliposome ODE trial using reduced Bayesian parameters.
%
% This script loads the reduced Gaussian-prior/inverse-gamma posterior
% result when available. If no result file is found, it falls back to the
% posterior means reported in the revised proteoliposome table.
%% =====================================================================

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
repo_dir = fileparts(script_dir);

theta_point = load_liposome_reduced_point(repo_dir, script_dir);

Tmax_m = theta_point(1);
Tmax_p = theta_point(2);
K_D_m  = theta_point(3);
K_D_p  = theta_point(4);

rho_pm = (Tmax_p * K_D_m) / (Tmax_m * K_D_p);

fprintf('Using reduced proteoliposome parameters:\n');
fprintf('Tmax_m=%.6g, Tmax_p=%.6g, K_D_m=%.6g, K_D_p=%.6g\n', ...
    Tmax_m, Tmax_p, K_D_m, K_D_p);
fprintf('Derived rho_pm=%.6g\n', rho_pm);

% Example compartment volumes per gram of protein.
Vint_g = 2.436;  % L/g, proteoliposome internal aqueous volume
Vext_g = 0.49;   % L/g, external assay medium

% Example malate uptake against internal phosphate.
Mm0   = 0.0;
Mext0 = 1.0;
Pm0   = 15.0;
Pext0 = 0.0;

y0 = [Mm0; Mext0; Pm0; Pext0];
tspan_min = [0, 5];

odefun = @(t, y) slc25a10_model_liposome_reduced_Jpho( ...
    t, y, Vext_g, Vint_g, K_D_m, K_D_p, Tmax_m, Tmax_p);

opts = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-13);
[Time_min, Y] = ode15s(odefun, tspan_min, y0, opts);

figure('Color', 'w', 'Name', 'Reduced proteoliposome ODE trial');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Time_min, Y(:, 1), 'LineWidth', 1.8); hold on;
plot(Time_min, Y(:, 2), 'LineWidth', 1.8);
xlabel('Time (min)');
ylabel('Malate (mM)');
legend({'M_{in}', 'M_{out}'}, 'Location', 'best');
grid on;

nexttile;
plot(Time_min, Y(:, 3), 'LineWidth', 1.8); hold on;
plot(Time_min, Y(:, 4), 'LineWidth', 1.8);
xlabel('Time (min)');
ylabel('Phosphate (mM)');
legend({'P_{in}', 'P_{out}'}, 'Location', 'best');
grid on;

%% =============================== Local functions ==========================
function theta_point = load_liposome_reduced_point(repo_dir, script_dir)
    candidates = { ...
        fullfile(script_dir, 'MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat'), ...
        fullfile(repo_dir, 'Bayesian_MCMC_MH_CODES', 'proteoliposome_calibration', 'MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat'), ...
        fullfile(repo_dir, 'FIGURE_5', 'reduced_parameter_Jpho', 'MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat'), ...
        fullfile(repo_dir, 'FIGURE_5', 'reduced_parameter_Jpho', 'New Folder', 'MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat')};

    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2
            result = load(candidates{k});
            if isfield(result, 'posterior_mean') && numel(result.posterior_mean) >= 4
                theta_point = result.posterior_mean(1:4);
                return;
            elseif isfield(result, 'theta_samples') && size(result.theta_samples, 2) >= 4
                theta_point = mean(result.theta_samples(:, 1:4), 1);
                return;
            end
        end
    end

    warning('Reduced proteoliposome posterior result was not found. Using table posterior means.');
    theta_point = [9.0682, 8.3548, 0.47927, 1.3709];
end
