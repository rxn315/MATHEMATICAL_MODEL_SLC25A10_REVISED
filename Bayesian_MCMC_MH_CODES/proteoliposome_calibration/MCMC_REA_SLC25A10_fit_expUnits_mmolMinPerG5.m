% ============================================================
% Reduced-parameter Metropolis-Hastings MCMC for SLC25A10
% proteoliposome Figure 5 data, using direct phosphate flux J_pho.
% This version uses Gaussian priors for the kinetic parameters and an
% inverse-gamma prior for sigma^2.
%
% This no-exMal1mM version excludes:
%   Proteoliposome_1993_exMal_1mM_inPho.mat
%
% This liposome dataset has Sm = Sims = 0, so lambda_21 and lambda_31
% are not separately identifiable. Their identifiable ratio is
%
%   rho_pm = lambda_31/lambda_21
%          = (Tmax_p*K_D_m)/(Tmax_m*K_D_p).
%
% Parameters sampled:
%   [1] Tmax_m   mmol/min/g protein
%   [2] Tmax_p   mmol/min/g protein
%   [3] K_D_m    mM
%   [4] K_D_p    mM
%   [5] sigma2   (mmol/min/g protein)^2
%
% Outputs written in this folder:
%   MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
%   posterior_summary_table_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.tex
%   Liposom_FIGURE4_gaussian_IGsigma_no_exMal1mM.png
%   Liposom_FIGURE4_gaussian_IGsigma_no_exMal1mM.fig
%   FIGURE_5_reduced_Jpho_gaussian_IGsigma_no_exMal1mM_fit_main.png
%   FIGURE_5_reduced_Jpho_gaussian_IGsigma_no_exMal1mM_fit_exMal_inPho.png
%
% No Statistics and Machine Learning Toolbox functions are used.
% ============================================================

clear; clc; close all;
rng(42);

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

result_file = fullfile(script_dir, ...
    'MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat');

% Set this to false if you only want to regenerate figures from result_file.
run_mcmc = true;

%% ----------------- Step 1: Load Experimental Data -----------------
data = load_liposome_data(script_dir);
n_total = numel(data.J_mal_exp) + numel(data.J_pho_exp);

%% ----------------- Step 2: Prior and MCMC Settings ----------------
% theta = [Tmax_m, Tmax_p, K_D_m, K_D_p, sigma2]
kinetic_prior_mu = [ ...
    6;  ... % Tmax_m, mmol/min/g
    6;  ... % Tmax_p, mmol/min/g
    0.49 ;   ... % K_D_m, mM
    1.41];     % K_D_p, mM

kinetic_prior_sigma = [ ...
    2.4;  ... % Tmax_m
    2.4;  ... % Tmax_p
    0.05; ... % K_D_m
    0.35];    % K_D_p

% Match the inverse-gamma prior to the sigma2 mean/std used previously.
sigma2_prior_mean = 3;
sigma2_prior_std  = 1;
alpha_0 = 2 + (sigma2_prior_mean / sigma2_prior_std)^2;
beta_0  = sigma2_prior_mean * (alpha_0 - 1);

% Positivity is enforced because all five parameters are physical rates,
% affinities, or variances. No upper truncation is imposed.
theta_lb = [1.0e-8; 1.0e-8; 1.0e-8; 1.0e-8; 1.0e-8];
theta_ub = [Inf; Inf; Inf; Inf; Inf];

theta0 = [kinetic_prior_mu; sigma2_prior_mean];

n_samples    = 1000000;
burn_in      = 10000;
thin         = 1;
z_step       = 0.1;
report_every = 2000;

if run_mcmc || exist(result_file, 'file') ~= 2
    fprintf('Running reduced Figure 5 MCMC with Gaussian kinetic priors, inverse-gamma sigma2, and no exMal=1mM/inPho dataset...\n');

    theta = theta0(:);
    p_all = numel(theta);
    theta_samples_raw = zeros(n_samples, p_all);
    acc = false(n_samples, 1);

    loglik_curr = log_likelihood(theta(1:4), theta(5), data, n_total);
    logprior_curr = log_prior_gaussian_IGsigma(theta, theta_lb, theta_ub, ...
        kinetic_prior_mu, kinetic_prior_sigma, alpha_0, beta_0);

    for i = 1:n_samples
        log_theta_prop = log(theta) + z_step * randn(p_all, 1);
        theta_prop = exp(log_theta_prop);

        logprior_prop = log_prior_gaussian_IGsigma(theta_prop, theta_lb, theta_ub, ...
            kinetic_prior_mu, kinetic_prior_sigma, alpha_0, beta_0);
        if isfinite(logprior_prop)
            loglik_prop = log_likelihood(theta_prop(1:4), theta_prop(5), data, n_total);

            % Correct Hastings factor for log-normal random-walk proposals.
            log_hastings = sum(log(theta_prop) - log(theta));

            log_alpha = (loglik_prop - loglik_curr) + ...
                        (logprior_prop - logprior_curr) + log_hastings;

            if log(rand) < log_alpha
                theta = theta_prop;
                loglik_curr = loglik_prop;
                logprior_curr = logprior_prop;
                acc(i) = true;
            end
        end

        theta_samples_raw(i, :) = theta.';

        if mod(i, report_every) == 0
            fprintf('Iter %7d | acc(last %d)=%.3f | acc(all)=%.3f | sigma2=%.4g\n', ...
                i, report_every, mean(acc(i-report_every+1:i)), mean(acc(1:i)), theta(5));
        end
    end

    keep_idx = (burn_in+1):thin:n_samples;
    theta_samples = theta_samples_raw(keep_idx, :);

    rho_pm_samples_raw = derived_rho_pm(theta_samples_raw);
    rho_pm_samples = derived_rho_pm(theta_samples);

    posterior_samples_plot = [theta_samples(:, 1:4), rho_pm_samples, theta_samples(:, 5)];
    posterior_mean = mean(posterior_samples_plot, 1);
    posterior_std = std(posterior_samples_plot, 0, 1);
    posterior_ci = local_percentile(posterior_samples_plot, [2.5 97.5], 1);

    save(result_file, ...
        'theta_samples_raw', 'theta_samples', ...
        'rho_pm_samples_raw', 'rho_pm_samples', ...
        'posterior_samples_plot', 'posterior_mean', 'posterior_std', ...
        'posterior_ci', 'acc', 'z_step', 'burn_in', 'thin', ...
        'theta_lb', 'theta_ub', 'theta0', ...
        'kinetic_prior_mu', 'kinetic_prior_sigma', ...
        'sigma2_prior_mean', 'sigma2_prior_std', 'alpha_0', 'beta_0', ...
        'data', 'n_total');

    fprintf('\nSaved MCMC result to:\n%s\n', result_file);
else
    fprintf('Loading existing MCMC result:\n%s\n', result_file);
    result = load(result_file);
    theta_samples_raw = result.theta_samples_raw;
    theta_samples = result.theta_samples;
    acc = result.acc;
    z_step = result.z_step;
    burn_in = result.burn_in;
    thin = result.thin;
    theta_lb = result.theta_lb;
    theta_ub = result.theta_ub;
    kinetic_prior_mu = result.kinetic_prior_mu;
    kinetic_prior_sigma = result.kinetic_prior_sigma;
    sigma2_prior_mean = result.sigma2_prior_mean;
    sigma2_prior_std = result.sigma2_prior_std;
    alpha_0 = result.alpha_0;
    beta_0 = result.beta_0;
    data = result.data;
    n_total = result.n_total;

    rho_pm_samples_raw = derived_rho_pm(theta_samples_raw);
    rho_pm_samples = derived_rho_pm(theta_samples);
    posterior_samples_plot = [theta_samples(:, 1:4), rho_pm_samples, theta_samples(:, 5)];
    posterior_mean = mean(posterior_samples_plot, 1);
    posterior_std = std(posterior_samples_plot, 0, 1);
    posterior_ci = local_percentile(posterior_samples_plot, [2.5 97.5], 1);
end

%% ----------------- Step 3: Posterior Predictions -----------------
n_post = size(theta_samples, 1);
n_mal = numel(data.J_mal_exp);
n_pho = numel(data.J_pho_exp);

J_M_pred = zeros(n_post, n_mal);
J_P_pred = zeros(n_post, n_pho);

for i = 1:n_post
    [Jm, Jp] = model_fluxes_reduced_Jpho(theta_samples(i, 1:4), data);
    J_M_pred(i, :) = Jm;
    J_P_pred(i, :) = Jp;
end

mean_J_M = mean(J_M_pred, 1);
CI_J_M = local_percentile(J_M_pred, [2.5 97.5], 1);

mean_J_P = mean(J_P_pred, 1);
CI_J_P = local_percentile(J_P_pred, [2.5 97.5], 1);

%% ----------------- Step 4: Save Summary and Figures ---------------
write_summary_table(fullfile(script_dir, ...
    'posterior_summary_table_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.tex'), ...
    posterior_samples_plot, acc, z_step);

fprintf('\nPosterior summary, reduced Figure 5 model with Gaussian kinetic priors, inverse-gamma sigma2, and no exMal=1mM/inPho dataset:\n');
print_summary_to_console(posterior_mean, posterior_std, posterior_ci, acc);

plot_liposom_posterior(theta_samples_raw, theta_samples);
save_current_figure(script_dir, 'Liposom_FIGURE4_gaussian_IGsigma_no_exMal1mM');

plot_fit_main(data, mean_J_M, CI_J_M, mean_J_P, CI_J_P);
save_current_figure(script_dir, 'FIGURE_5_reduced_Jpho_gaussian_IGsigma_no_exMal1mM_fit_main');

plot_fit_exmal_inpho(data, mean_J_M, CI_J_M);
save_current_figure(script_dir, 'FIGURE_5_reduced_Jpho_gaussian_IGsigma_no_exMal1mM_fit_exMal_inPho');

fprintf('\nSaved reduced Figure 5 outputs in:\n%s\n', script_dir);

%% ========================= Helper Functions =========================
function data = load_liposome_data(script_dir)
% Keep experimental rates in mmol/min/g protein.

S = load(fullfile(script_dir, 'Proteoliposome_intPho_exMAL.mat'));
n_m0 = numel(S.intPho_exMal);
Pm_v_mal0 = 15 * ones(1, n_m0);
Mc_v_mal0 = S.intPho_exMal(:).';
J_mal_exp0 = S.intPho_exMal_Rate(:).';
Mm_v_mal0 = zeros(1, n_m0);
Pc_v_mal0 = zeros(1, n_m0);

S = load(fullfile(script_dir, 'Proteoliposome_intMAl_exPHO.mat'));
n_p1 = numel(S.intMal_exPho);
Pm_v_pho1 = zeros(1, n_p1);
Mm_v_pho1 = 15 * ones(1, n_p1);
Pc_v_pho1 = S.intMal_exPho(:).';
Mc_v_pho1 = zeros(1, n_p1);
J_pho_exp1 = S.intMal_exPho_Rate(:).';

S = load(fullfile(script_dir, 'Proteoliposome_1993_intPho_10mM_exMAL.mat'));
n_m2 = numel(S.intPho_10mM_exMa);
Pm_v_mal2 = 1.0 * ones(1, n_m2);
Mc_v_mal2 = S.intPho_10mM_exMa(:).';
J_mal_exp2 = S.intPho_10mM_exMa_Rate(:).';
Mm_v_mal2 = zeros(1, n_m2);
Pc_v_mal2 = zeros(1, n_m2);

S = load(fullfile(script_dir, 'Proteoliposome_1993_intPho_25mM_exMAL.mat'));
n_m3 = numel(S.intPho_25mM_exMa);
Pm_v_mal3 = 2.5 * ones(1, n_m3);
Mc_v_mal3 = S.intPho_25mM_exMa(:).';
J_mal_exp3 = S.intPho_25mM_exMa_Rate(:).';
Mm_v_mal3 = zeros(1, n_m3);
Pc_v_mal3 = zeros(1, n_m3);

S = load(fullfile(script_dir, 'Proteoliposome_1993_exMal_020mM_inPho.mat'));
n_m6 = numel(S.exMal_020mM_inPho);
Pm_v_mal6 = S.exMal_020mM_inPho(:).';
Mc_v_mal6 = 0.2 * ones(1, n_m6);
J_mal_exp6 = S.exMal_020mM_inPho_Rate(:).';
Mm_v_mal6 = zeros(1, n_m6);
Pc_v_mal6 = zeros(1, n_m6);

Mc_v_mal = [Mc_v_mal0, Mc_v_mal2, Mc_v_mal3, Mc_v_mal6];
Mm_v_mal = [Mm_v_mal0, Mm_v_mal2, Mm_v_mal3, Mm_v_mal6];
Pm_v_mal = [Pm_v_mal0, Pm_v_mal2, Pm_v_mal3, Pm_v_mal6];
Pc_v_mal = [Pc_v_mal0, Pc_v_mal2, Pc_v_mal3, Pc_v_mal6];
J_mal_exp = [J_mal_exp0, J_mal_exp2, J_mal_exp3, J_mal_exp6];

idx_m0 = 1:n_m0;
idx_m2 = idx_m0(end) + (1:n_m2);
idx_m3 = idx_m2(end) + (1:n_m3);
idx_m6 = idx_m3(end) + (1:n_m6);

Mc_v_pho = Mc_v_pho1;
Mm_v_pho = Mm_v_pho1;
Pm_v_pho = Pm_v_pho1;
Pc_v_pho = Pc_v_pho1;
J_pho_exp = J_pho_exp1;
idx_p1 = 1:numel(J_pho_exp);

data = struct( ...
    'Mc_v_mal', Mc_v_mal, 'Mm_v_mal', Mm_v_mal, ...
    'Pm_v_mal', Pm_v_mal, 'Pc_v_mal', Pc_v_mal, ...
    'J_mal_exp', J_mal_exp, ...
    'Mc_v_pho', Mc_v_pho, 'Mm_v_pho', Mm_v_pho, ...
    'Pm_v_pho', Pm_v_pho, 'Pc_v_pho', Pc_v_pho, ...
    'J_pho_exp', J_pho_exp, ...
    'idx_m0', idx_m0, 'idx_m2', idx_m2, 'idx_m3', idx_m3, ...
    'idx_m6', idx_m6, 'idx_p1', idx_p1);
end

function logp = log_prior_gaussian_IGsigma(theta, theta_lb, theta_ub, ...
    kinetic_prior_mu, kinetic_prior_sigma, alpha_0, beta_0)
if any(theta <= 0) || any(theta < theta_lb) || any(theta > theta_ub)
    logp = -inf;
    return;
end

if any(kinetic_prior_sigma <= 0)
    logp = -inf;
    return;
end

z = (theta(1:4) - kinetic_prior_mu(:)) ./ kinetic_prior_sigma(:);
logp_kinetic = -0.5 * sum(z.^2) - ...
    sum(log(kinetic_prior_sigma(:))) - numel(z) * 0.5 * log(2*pi);

logp_sigma2 = log_inv_gamma(theta(5), alpha_0, beta_0);
logp = logp_kinetic + logp_sigma2;
end

function v = log_inv_gamma(sigma2, alpha, beta)
if sigma2 <= 0
    v = -inf;
    return;
end
v = (-alpha - 1) * log(sigma2) - beta / sigma2 - gammaln(alpha) + alpha * log(beta);
end

function loglik = log_likelihood(params4, sigma2, data, n_total)
if sigma2 <= 0
    loglik = -inf;
    return;
end
sse = obj_sse(params4, data);
loglik = -n_total/2 * log(2*pi*sigma2) - 0.5/sigma2 * sse;
end

function sse = obj_sse(params4, data)
[Jm, Jp] = model_fluxes_reduced_Jpho(params4, data);
sse = sum((Jm - data.J_mal_exp).^2) + sum((Jp - data.J_pho_exp).^2);
end

function [Jm, Jp] = model_fluxes_reduced_Jpho(params4, data)
Tmax_m = params4(1);
Tmax_p = params4(2);
K_D_m = params4(3);
K_D_p = params4(4);

n_mal = numel(data.J_mal_exp);
Jm = zeros(1, n_mal);
for i = 1:n_mal
    Mm = data.Mm_v_mal(i);
    Mims = data.Mc_v_mal(i);
    Pm = data.Pm_v_mal(i);
    Pims = data.Pc_v_mal(i);
    [~, Jm_i, ~] = compute_flux_condition_reduced( ...
        Tmax_m, Tmax_p, K_D_m, K_D_p, Mm, Mims, Pm, Pims);
    Jm(i) = Jm_i;
end

n_pho = numel(data.J_pho_exp);
Jp = zeros(1, n_pho);
for i = 1:n_pho
    Mm = data.Mm_v_pho(i);
    Mims = data.Mc_v_pho(i);
    Pm = data.Pm_v_pho(i);
    Pims = data.Pc_v_pho(i);
    [~, ~, Jp_i] = compute_flux_condition_reduced( ...
        Tmax_m, Tmax_p, K_D_m, K_D_p, Mm, Mims, Pm, Pims);
    Jp(i) = Jp_i;
end
end

function [J_succ, J_mal, J_pho] = compute_flux_condition_reduced( ...
    Tmax_m, Tmax_p, K_D_m, K_D_p, Mm, Mims, Pm, Pims)

rho_pm = (Tmax_p * K_D_m) / (Tmax_m * K_D_p);

delta1 = 1 + Mm/K_D_m + Pm/K_D_p;
delta2 = 1 + Mims/K_D_m + Pims/K_D_p;

phi_den = Mims + rho_pm * Pims;
if abs(phi_den) < 1.0e-14
    phi_den = 1.0e-14;
end

phi = (Mm + rho_pm * Pm) / phi_den;
den = delta1 + phi * delta2;

J_succ = 0;
J_mal = Tmax_m * (phi*Mims - Mm) / (K_D_m * den);
J_pho = Tmax_p * (phi*Pims - Pm) / (K_D_p * den);
end

function rho_pm = derived_rho_pm(theta_samples)
% theta_samples columns: [Tmax_m, Tmax_p, K_D_m, K_D_p, sigma2]
rho_pm = (theta_samples(:, 2) .* theta_samples(:, 3)) ./ ...
         (theta_samples(:, 1) .* theta_samples(:, 4));
end

function write_summary_table(out_file, samples6, acc, z_step)
% samples6 columns: [Tmax_m, Tmax_p, K_D_m, K_D_p, rho_pm, sigma2]
param_names = { ...
    'T^{m}_{\\max}', 'T^{p}_{\\max}', 'K^{m}_{D}', ...
    'K^{p}_{D}', '\\rho_{pm}', '\\sigma^2'};

means = mean(samples6, 1);
medians = median(samples6, 1);
stds = std(samples6, 0, 1);
ci = local_percentile(samples6, [2.5 97.5], 1);
modes = zeros(1, size(samples6, 2));
for j = 1:size(samples6, 2)
    [xi, f] = local_gaussian_kde(samples6(:, j), 600);
    [~, ix] = max(f);
    modes(j) = xi(ix);
end

accepted_count = sum(acc(:));
total_steps = numel(acc);
a_str = sprintf('%.3g', z_step);
fmt = @(x) sprintf('%.5g', x);

fid = fopen(out_file, 'w');
if fid < 0
    error('Could not open table file for writing: %s', out_file);
end

fprintf(fid, '\\begin{table}[ht]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Reduced-parameter SLC25A10 proteoliposome posterior summary using direct $J_{pho}$, Gaussian kinetic priors, inverse-gamma $\\sigma^2$, and excluding the $M_{out}=1.0$ mM internal-phosphate series.}\n');
fprintf(fid, '\\label{tab:liposome_reduced_jpho_gaussian_IGsigma_no_exMal1mM}\n');
fprintf(fid, '\\small\n');
fprintf(fid, '\\begin{tabular}{l c c c c c c c c c}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, '\\textbf{Parameter} & Mean & Median & Mode (KDE) & Std & 2.5\\%% & 97.5\\%% & $a$ & Accepted & Total \\\\\n');
fprintf(fid, '\\hline\n');
for j = 1:numel(param_names)
    fprintf(fid, '%s & %s & %s & %s & %s & %s & %s & %s & %d & %d \\\\\n', ...
        param_names{j}, fmt(means(j)), fmt(medians(j)), fmt(modes(j)), ...
        fmt(stds(j)), fmt(ci(1, j)), fmt(ci(2, j)), ...
        a_str, accepted_count, total_steps);
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
end

function print_summary_to_console(means, stds, ci, acc)
names = {'Tmax_m', 'Tmax_p', 'K_D_m', 'K_D_p', 'rho_pm', 'sigma2'};
units = {'mmol/min/g', 'mmol/min/g', 'mM', 'mM', 'dimensionless', '(mmol/min/g)^2'};
for j = 1:numel(names)
    fprintf('%-8s = %.6g  (std %.3g)  CI[%.6g, %.6g]  %s\n', ...
        names{j}, means(j), stds(j), ci(1, j), ci(2, j), units{j});
end
fprintf('Acceptance rate: %.3f\n', mean(acc(:)));
end

function plot_liposom_posterior(theta_raw, theta_post)
plot_raw = [theta_raw(:, 1:4), derived_rho_pm(theta_raw), theta_raw(:, 5)];
plot_post = [theta_post(:, 1:4), derived_rho_pm(theta_post), theta_post(:, 5)];
labels = {'T^{m}_{max}', 'T^{p}_{max}', 'K^{m}_{D}', ...
          'K^{p}_{D}', '\rho_{pm}', '\sigma^2'};

figure('Color', 'w', 'Position', [80 80 1920 940], ...
    'Name', 'Proteoliposome reduced posterior traces and PDFs');

n_plot = numel(labels);
n_raw = size(plot_raw, 1);
plot_step = max(1, floor(n_raw / 50000));
trace_idx = 1:plot_step:n_raw;

for k = 1:n_plot
    ax = subplot(2, 6, k);
    plot(ax, trace_idx, plot_raw(trace_idx, k), 'k-', 'LineWidth', 0.8);
    xlabel(ax, {'MCMC', 'iteration'}, 'FontWeight', 'bold');
    ylabel(ax, labels{k}, 'FontWeight', 'bold', 'Interpreter', 'tex');
    grid(ax, 'on');
    set(ax, 'FontSize', 14, 'FontWeight', 'bold');
end

for k = 1:n_plot
    ax = subplot(2, 6, 6 + k);
    hold(ax, 'on');
    v = plot_post(:, k);
    histogram(ax, v, 50, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
        'FaceColor', [0.43921568627451 0.6 0.929411764705882], ...
        'DisplayName', 'Posterior Samples');
    [xi, f] = local_gaussian_kde(v, 600);
    plot(ax, xi, f, 'r-', 'LineWidth', 1.8, 'DisplayName', 'PDF');
    xlabel(ax, labels{k}, 'FontWeight', 'bold', 'Interpreter', 'tex');
    ylabel(ax, 'pdf', 'FontWeight', 'bold');
    grid(ax, 'on');
    set(ax, 'FontSize', 14, 'FontWeight', 'bold');
    if k == 1
        lg = legend(ax, 'show');
        set(lg, 'FontSize', 13, 'EdgeColor', [1 1 1]);
    end
    hold(ax, 'off');
end

for k = 1:12
    ax = subplot(2, 6, k);
    letter = sprintf('%c)', char('a' + k - 1));
    text(ax, -0.18, 1.08, letter, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 22, 'Color', [0 0 0]);
end
end

function plot_fit_main(data, mean_J_M, CI_J_M, mean_J_P, CI_J_P)
figure('Color', 'w', 'Position', [80 80 1200 850], ...
    'Name', 'Reduced Figure 5 posterior prediction bands');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = nexttile;
idx = data.idx_m0;
plot_prediction_band(ax, data.Mc_v_mal(idx), mean_J_M(idx), CI_J_M(:, idx), ...
    data.J_mal_exp(idx), 'External malate (mM)', ...
    'J_M (mmol min^{-1} g^{-1})', 'Malate, P_{in}=15 mM');

ax = nexttile;
idx = data.idx_m2;
plot_prediction_band(ax, data.Mc_v_mal(idx), mean_J_M(idx), CI_J_M(:, idx), ...
    data.J_mal_exp(idx), 'External malate (mM)', ...
    'J_M (mmol min^{-1} g^{-1})', 'Malate, P_{in}=1.0 mM');

ax = nexttile;
idx = data.idx_m3;
plot_prediction_band(ax, data.Mc_v_mal(idx), mean_J_M(idx), CI_J_M(:, idx), ...
    data.J_mal_exp(idx), 'External malate (mM)', ...
    'J_M (mmol min^{-1} g^{-1})', 'Malate, P_{in}=2.5 mM');

ax = nexttile;
idx = data.idx_p1;
plot_prediction_band(ax, data.Pc_v_pho(idx), mean_J_P(idx), CI_J_P(:, idx), ...
    data.J_pho_exp(idx), 'External phosphate (mM)', ...
    'J_P (mmol min^{-1} g^{-1})', 'Phosphate, M_{in}=15 mM');
end

function plot_fit_exmal_inpho(data, mean_J_M, CI_J_M)
figure('Color', 'w', 'Position', [120 120 1150 460], ...
    'Name', 'Reduced Figure 5 exMal inPho prediction band, no exMal=1mM dataset');
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = nexttile;
idx = data.idx_m6;
plot_prediction_band(ax, data.Pm_v_mal(idx), mean_J_M(idx), CI_J_M(:, idx), ...
    data.J_mal_exp(idx), 'Internal phosphate (mM)', ...
    'J_M (mmol min^{-1} g^{-1})', 'Malate, M_{out}=0.2 mM');
end

function plot_prediction_band(ax, x, y_mean, y_ci, y_exp, xlab, ylab, ttl)
x = x(:).';
y_mean = y_mean(:).';
y_exp = y_exp(:).';
y_lo = y_ci(1, :);
y_hi = y_ci(2, :);
[x_sort, ord] = sort(x);

hold(ax, 'on');
fill(ax, [x_sort, fliplr(x_sort)], ...
    [y_lo(ord), fliplr(y_hi(ord))], ...
    'r', 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'DisplayName', '95% CI');
plot(ax, x_sort, y_mean(ord), 'r--', 'LineWidth', 2, 'DisplayName', 'Posterior mean');
plot(ax, x_sort, y_exp(ord), 'b-', 'LineWidth', 2, 'DisplayName', 'Experimental');
xlabel(ax, xlab, 'FontWeight', 'bold');
ylabel(ax, ylab, 'FontWeight', 'bold');
title(ax, ttl, 'FontWeight', 'bold', 'Interpreter', 'tex');
grid(ax, 'on');
set(ax, 'FontSize', 13, 'FontWeight', 'bold');
legend(ax, 'Location', 'best');
hold(ax, 'off');
end

function save_current_figure(script_dir, base_name)
png_file = fullfile(script_dir, [base_name, '.png']);
fig_file = fullfile(script_dir, [base_name, '.fig']);

if exist('exportgraphics', 'file') == 2
    exportgraphics(gcf, png_file, 'Resolution', 300);
else
    print(gcf, png_file, '-dpng', '-r300');
end

if exist('savefig', 'file') == 2
    savefig(gcf, fig_file);
end
fprintf('Saved %s\n', png_file);
end

function q = local_percentile(x, p, dim)
if nargin < 3
    dim = 1;
end

if isvector(x)
    q = local_percentile_vector(x, p);
    return;
end

if dim == 1
    q = zeros(numel(p), size(x, 2));
    for j = 1:size(x, 2)
        q(:, j) = local_percentile_vector(x(:, j), p);
    end
elseif dim == 2
    q = zeros(size(x, 1), numel(p));
    for i = 1:size(x, 1)
        q(i, :) = local_percentile_vector(x(i, :), p).';
    end
else
    error('local_percentile supports only dim 1 or dim 2.');
end
end

function q = local_percentile_vector(v, p)
v = v(:);
v = v(isfinite(v));
v = sort(v);
n = numel(v);
q = zeros(numel(p), 1);

if n == 0
    q(:) = NaN;
    return;
end
if n == 1
    q(:) = v;
    return;
end

for i = 1:numel(p)
    pos = 1 + (n - 1) * p(i) / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q(i) = v(lo);
    else
        w = pos - lo;
        q(i) = (1 - w) * v(lo) + w * v(hi);
    end
end
end

function [xi, f] = local_gaussian_kde(v, n_grid)
% Toolbox-free Gaussian KDE using Silverman's bandwidth rule.
if nargin < 2
    n_grid = 600;
end

v = v(:);
v = v(isfinite(v));
if isempty(v)
    xi = 0;
    f = 0;
    return;
end

max_points = 50000;
if numel(v) > max_points
    idx = round(linspace(1, numel(v), max_points));
    v = v(idx);
end

lo = min(v);
hi = max(v);
span = hi - lo;
if span <= 0
    span = max(abs(lo), 1) * 0.02;
    lo = lo - span;
    hi = hi + span;
end

pad = 0.08 * span;
if lo >= 0
    x_lo = max(0, lo - pad);
else
    x_lo = lo - pad;
end
xi = linspace(x_lo, hi + pad, n_grid);

n = numel(v);
bw = 1.06 * std(v, 0) * n^(-1/5);
if ~isfinite(bw) || bw <= 0
    bw = max(abs(mean(v)), 1) * 1.0e-3;
end

f = zeros(size(xi));
norm_const = 1 / (sqrt(2*pi) * bw * n);
chunk_size = 50;
for i = 1:chunk_size:numel(xi)
    idx = i:min(i + chunk_size - 1, numel(xi));
    z = (xi(idx).' - v.') / bw;
    f(idx) = norm_const * sum(exp(-0.5 * z.^2), 2).';
end
end
