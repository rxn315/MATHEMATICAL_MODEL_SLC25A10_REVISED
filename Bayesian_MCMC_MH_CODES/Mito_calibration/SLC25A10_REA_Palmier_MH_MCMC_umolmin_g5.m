% ============================================================
% Reduced-parameter MCMC for Palmieri DIC data using uniform priors.
%
% This version estimates six kinetic parameters:
%   [1] Ts_max_exp   (umol/min/g)
%   [2] Tm_max_exp   (umol/min/g)
%   [3] Tp_max_exp   (umol/min/g)
%   [4] Ks_D         (mM)
%   [5] Km_D         (mM)
%   [6] Kp_D         (mM)
%
% The coupling ratios are not inferred and have no prior/support checks.
% They are only derived after sampling for reporting:
%   lambda21 = (Tm_max_exp * Ks_D) / (Ts_max_exp * Km_D)
%   lambda31 = (Tp_max_exp * Ks_D) / (Ts_max_exp * Kp_D)
%
% theta also includes sigma2 for the Gaussian likelihood:
%   theta = [Ts_max_exp, Tm_max_exp, Tp_max_exp, Ks_D, Km_D, Kp_D, sigma2]
%
% Uniform priors are used for every inferred quantity in theta. Edit
% theta_lb and theta_ub below if you want wider or narrower support.
%
% Output:
%   MCMC_Result_Palmier2_reduced_uniform.mat
%   posterior_summary_table_P6_lambdas_derived_uniform.tex
%   Mito_FIGURE4_uniform.png
% ============================================================

clear; clc; close all;
rng(42);

%% ----------------- Fixed matrix concentrations ----------------------
Pm_mat_fix = 2;     % mM
Sm_mat_fix = 0.2;   % mM
Mm_mat_fix = 0.2;   % mM

%% ----------------- Load experimental data ---------------------------
load('Palmier_DIC_Exp_Data.mat');             % S, T_S
Sc_v_succ  = S(:)';
Mc_v_succ  = zeros(size(Sc_v_succ));
Pc_v_succ  = zeros(size(Sc_v_succ));
J_succ_exp = T_S(:)';

load('Palmier_DIC_Competition_Exp_Data.mat'); % M_with_S, T_M_with_S, M_control, T_M_control
Mc_v_mal   = (1 ./ M_with_S(:))';
Sc_v_mal   = 0.5 * ones(size(Mc_v_mal));
Pc_v_mal   = zeros(size(Mc_v_mal));
J_mal_exp  = (1 ./ T_M_with_S(:))';

% Keep the original Figure 4 convention: the phosphate-associated panel
% uses the malate-control data from Palmier_DIC_Competition_Exp_Data.mat.


load('Palmier_DIC_Exp_Data.mat'); 
Pc_v_pho   = (1 ./ P(:))' ;
Mc_v_pho   = zeros(size(Pc_v_pho));
Sc_v_pho   = zeros(size(Pc_v_pho));

J_pho_exp  = (1 ./ T_P(:))';

 
data = struct( ...
    'Mc_v_succ',Mc_v_succ,'Sc_v_succ',Sc_v_succ,'Pc_v_succ',Pc_v_succ,'J_succ_exp',J_succ_exp, ...
    'Mc_v_mal', Mc_v_mal, 'Sc_v_mal', Sc_v_mal, 'Pc_v_mal', Pc_v_mal, 'J_mal_exp', J_mal_exp, ...
    'Mc_v_pho', Mc_v_pho, 'Sc_v_pho', Sc_v_pho, 'Pc_v_pho', Pc_v_pho, 'J_pho_exp', J_pho_exp);

n_total = numel(J_succ_exp) + numel(J_mal_exp) + numel(J_pho_exp);

%% ----------------- Uniform priors -----------------------------------
% theta = [Ts_max, Tm_max, Tp_max, Ks_D, Km_D, Kp_D, sigma2]
% Units: Tmax values are umol/min/g, K_D values are mM, sigma2 is
% (umol/min/g)^2. Bounds are intentionally broad and can be edited.
theta_lb = [1.0e-6; 1.0e-6; 1.0e-6; 1.0e-6; 1.0e-6; 1.0e-6; 1.0e-8];
theta_ub = [1000;    1000;    1000; 10;     10;      10;     1.0e3];

log_prior_all = @(theta7) log_prior_uniform(theta7, theta_lb, theta_ub);

log_likelihood = @(params6, sigma2) ...
    -n_total/2 * log(2*pi*sigma2) ...
    - 0.5/sigma2 * obj_sse(params6, data, Pm_mat_fix, Sm_mat_fix, Mm_mat_fix);

%% ----------------- Metropolis-Hastings -------------------------------
% theta = [Ts_max, Tm_max, Tp_max, Ks_D, Km_D, Kp_D, sigma2]
theta = [69.5885; 62.9399; 62.9399; 1.1434; 0.2393; 1.8511; 1.9377];

n_samples    = 1000000;
burn_in      = 10000;
thin         = 1;
z_step       = 0.1;
report_every = 2000;

p_all = numel(theta);
theta_samples_raw = zeros(n_samples, p_all);
acc = false(n_samples, 1);

loglik_curr   = log_likelihood(theta(1:6), theta(7));
logprior_curr = log_prior_all(theta);

for i = 1:n_samples
    log_theta_prop = log(theta) + z_step * randn(p_all, 1);
    theta_prop = exp(log_theta_prop);

    logprior_prop = log_prior_all(theta_prop);
    if ~isfinite(logprior_prop)
        theta_samples_raw(i,:) = theta.';
        continue;
    end

    loglik_prop  = log_likelihood(theta_prop(1:6), theta_prop(7));
    log_hastings = sum(log(theta) - log(theta_prop));

    log_alpha = (loglik_prop - loglik_curr) + ...
                (logprior_prop - logprior_curr) + ...
                log_hastings;

    if log(rand) < log_alpha
        theta         = theta_prop;
        loglik_curr   = loglik_prop;
        logprior_curr = logprior_prop;
        acc(i)        = true;
    end

    theta_samples_raw(i,:) = theta.';
    if mod(i, report_every) == 0
        fprintf('Iter %7d | acc(last %d)=%.3f | acc(all)=%.3f | sigma2=%.4g\n', ...
            i, report_every, mean(acc(i-report_every+1:i)), mean(acc(1:i)), theta(7));
    end
end
fprintf('\nOverall acceptance rate: %.3f\n', mean(acc));

keep_idx = (burn_in+1):thin:n_samples;
theta_samples = theta_samples_raw(keep_idx, :);

lambda21_samples_raw = derived_lambda21(theta_samples_raw(:,1:6));
lambda31_samples_raw = derived_lambda31(theta_samples_raw(:,1:6));
lambda21_samples = lambda21_samples_raw(keep_idx);
lambda31_samples = lambda31_samples_raw(keep_idx);

posterior_samples_plot = [theta_samples(:,1:6), lambda21_samples, lambda31_samples];
posterior_mean = mean(theta_samples, 1);
posterior_std  = std(theta_samples, 0, 1);
lambda21_mean = mean(lambda21_samples);
lambda31_mean = mean(lambda31_samples);

fprintf('\nPosterior means:\n');
fprintf('Ts_max=%.6g, Tm_max=%.6g, Tp_max=%.6g, Ks_D=%.6g, Km_D=%.6g, Kp_D=%.6g, sigma2=%.6g\n', ...
    posterior_mean(1), posterior_mean(2), posterior_mean(3), posterior_mean(4), ...
    posterior_mean(5), posterior_mean(6), posterior_mean(7));
fprintf('Derived lambda21=%.6g, lambda31=%.6g\n', lambda21_mean, lambda31_mean);

%% ----------------- Save results -------------------------------------
save('MCMC_Result_Palmier2_reduced_uniform.mat', ...
    'theta_samples', 'theta_samples_raw', 'posterior_samples_plot', ...
    'lambda21_samples', 'lambda31_samples', ...
    'lambda21_samples_raw', 'lambda31_samples_raw', ...
    'lambda21_mean', 'lambda31_mean', 'posterior_mean', 'posterior_std', ...
    'acc', 'z_step', 'burn_in', 'thin', 'n_samples', ...
    'Pm_mat_fix', 'Sm_mat_fix', 'Mm_mat_fix', ...
    'theta_lb', 'theta_ub');

write_summary_tex(theta_samples, lambda21_samples, lambda31_samples, acc, z_step);
plot_reduced_posterior(theta_samples_raw, theta_samples, ...
    lambda21_samples_raw, lambda31_samples_raw, lambda21_samples, lambda31_samples);

if exist('exportgraphics', 'file') == 2
    exportgraphics(gcf, 'Mito_FIGURE4_uniform.png', 'Resolution', 300);
else
    print(gcf, 'Mito_FIGURE4_uniform.png', '-dpng', '-r300');
end
savefig(gcf, 'Mito_FIGURE4_uniform.fig');

%% ============================================================
% Helper functions
% ============================================================

function sse = obj_sse(params6, data, Pm_fix, Sm_fix, Mm_fix)
fluxes = obj_func_palmier(params6, data, Pm_fix, Sm_fix, Mm_fix);
J_S = fluxes{1};
J_M = fluxes{2};
J_P = fluxes{3};
sse = sum((J_S - data.J_succ_exp).^2) + ...
      sum((J_M - data.J_mal_exp).^2) + ...
      sum((J_P - data.J_pho_exp).^2);
end

function lp = log_prior_uniform(theta7, theta_lb, theta_ub)
theta7 = theta7(:);
if numel(theta7) ~= numel(theta_lb)
    error('theta must have %d elements.', numel(theta_lb));
end

if any(theta7 < theta_lb) || any(theta7 > theta_ub)
    lp = -inf; return;
end

% Uniform prior density is constant inside support. Including the constant
% makes the definition explicit; it cancels in the MH ratio.
lp = -sum(log(theta_ub - theta_lb));
end

function out = obj_func_palmier(params6, data, Pm_fix, Sm_fix, Mm_fix)
Ts_max = params6(1);
Tm_max = params6(2);
Tp_max = params6(3);
K_D_s  = params6(4);
K_D_m  = params6(5);
K_D_p  = params6(6);

Pm = Pm_fix; Sm = Sm_fix; Mm = Mm_fix;

n_succ = numel(data.Sc_v_succ);
J_S = zeros(1, n_succ);
for i = 1:n_succ
    [Js,~,~] = compute_flux_condition_REA( ...
        Ts_max, Tm_max, Tp_max, K_D_s, K_D_m, K_D_p, ...
        Sm, data.Sc_v_succ(i), Mm, data.Mc_v_succ(i), Pm, data.Pc_v_succ(i));
    J_S(i) = Js;
end

n_mal = numel(data.Mc_v_mal);
J_M = zeros(1, n_mal);
for i = 1:n_mal
    [~,Jm,~] = compute_flux_condition_REA( ...
        Ts_max, Tm_max, Tp_max, K_D_s, K_D_m, K_D_p, ...
        Sm, data.Sc_v_mal(i), Mm, data.Mc_v_mal(i), Pm, data.Pc_v_mal(i));
    J_M(i) = Jm;
end

n_pho = numel(data.Mc_v_pho);
J_P = zeros(1, n_pho);
for i = 1:n_pho
    [~,~,Jp] = compute_flux_condition_REA( ...
        Ts_max, Tm_max, Tp_max, K_D_s, K_D_m, K_D_p, ...
        Sm, data.Sc_v_pho(i), Mm, data.Mc_v_pho(i), Pm, data.Pc_v_pho(i));
    J_P(i) = Jp;
end

out = {J_S, J_M, J_P};
end

function [J_succ, J_mal, J_pho] = compute_flux_condition_REA( ...
    Ts_max, Tm_max, Tp_max, K_D_s, K_D_m, K_D_p, ...
    Sm, Sims, Mm, Mims, Pm, Pims)

delta1 = 1 + Sm/K_D_s   + Mm/K_D_m   + Pm/K_D_p;
delta2 = 1 + Sims/K_D_s + Mims/K_D_m + Pims/K_D_p;

ratio_m_to_s = (Tm_max * K_D_s) / (Ts_max * K_D_m);
ratio_p_to_s = (Tp_max * K_D_s) / (Ts_max * K_D_p);

phi_den = Sims + ratio_m_to_s*Mims + ratio_p_to_s*Pims;
if abs(phi_den) < 1.0e-14
    phi_den = sign(phi_den + 1.0e-14) * 1.0e-14;
end

phi = (Sm + ratio_m_to_s*Mm + ratio_p_to_s*Pm) / phi_den;
den = delta1 + phi*delta2;

J_succ = Ts_max * (phi*Sims - Sm) / (K_D_s * den);
J_mal  = Tm_max * (phi*Mims - Mm) / (K_D_m * den);

% Original Figure 4 convention.
%J_pho  = J_mal;
J_pho  = Tp_max * (phi*Pims - Pm) / (K_D_p * den);
end

function lambda21 = derived_lambda21(params6)
if isvector(params6)
    p = params6(:).';
    lambda21 = (p(2) * p(4)) / (p(1) * p(5));
else
    lambda21 = (params6(:,2) .* params6(:,4)) ./ (params6(:,1) .* params6(:,5));
end
end

function lambda31 = derived_lambda31(params6)
if isvector(params6)
    p = params6(:).';
    lambda31 = (p(3) * p(4)) / (p(1) * p(6));
else
    lambda31 = (params6(:,3) .* params6(:,4)) ./ (params6(:,1) .* params6(:,6));
end
end

function write_summary_tex(theta_samples, lambda21_samples, lambda31_samples, acc, z_step)
samples = [theta_samples(:,1:6), lambda21_samples, lambda31_samples, theta_samples(:,7)];
names = {'$T^s_{\max}$', '$T^m_{\max}$', '$T^p_{\max}$', '$K^s_D$', '$K^m_D$', '$K^p_D$', ...
    '$\lambda_{21}$ (derived)', '$\lambda_{31}$ (derived)', '$\sigma^2$'};

means = mean(samples, 1);
medians = median(samples, 1);
stds = std(samples, 0, 1);
cis = local_percentile(samples, [2.5 97.5]);
fmt = @(x) sprintf('%.5g', x);

fid = fopen('posterior_summary_table_P6_lambdas_derived_uniform.tex', 'w');
fprintf(fid, '\\begin{table}[ht]\n\\centering\n');
fprintf(fid, '\\caption{Reduced-parameter MCMC posterior summary with bounded uniform priors. $\\lambda_{21}$ and $\\lambda_{31}$ are derived from each posterior sample.}\n');
fprintf(fid, '\\label{tab:post_reduced_lambda_derived_uniform}\n\\small\n');
fprintf(fid, '\\begin{tabular}{l c c c c c c c}\\hline\n');
fprintf(fid, '\\textbf{Quantity} & Mean & Median & Std & 2.5\\%% & 97.5\\%% & Accepted & Total \\\\\n');
fprintf(fid, '\\hline\n');
for j = 1:numel(names)
    fprintf(fid, '%s & %s & %s & %s & %s & %s & %d & %d \\\\\n', ...
        names{j}, fmt(means(j)), fmt(medians(j)), fmt(stds(j)), ...
        fmt(cis(1,j)), fmt(cis(2,j)), sum(acc(:)), numel(acc));
end
fprintf(fid, '\\hline\n\\end{tabular}\n\\end{table}\n');
fclose(fid);

fprintf('Wrote posterior_summary_table_P6_lambdas_derived_uniform.tex (z_step=%g)\n', z_step);
end

function plot_reduced_posterior(theta_raw, theta_post, lambda21_raw, lambda31_raw, lambda21_post, lambda31_post)
plot_raw = [theta_raw(:,1:6), lambda21_raw, lambda31_raw];
plot_post = [theta_post(:,1:6), lambda21_post, lambda31_post];
labels = {'T^{s}_{max}', 'T^{m}_{max}', 'T^{p}_{max}', ...
          'K^{s}_{D}', 'K^{m}_{D}', 'K^{p}_{D}', ...
          '\lambda_{21}', '\lambda_{31}'};

figure('Color', 'w', 'Position', [80 80 1365 768], ...
    'Name', 'Reduced Palmieri posterior traces and PDFs, uniform priors');

for k = 1:8
    ax = subplot(4,4,k);
    plot(ax, plot_raw(:,k), 'k-', 'LineWidth', 0.8);
    xlabel(ax, {'MCMC', 'iteration'}, 'FontWeight', 'bold');
    ylabel(ax, labels{k}, 'FontWeight', 'bold', 'Interpreter', 'tex');
    grid(ax, 'on'); %box(ax, 'on');
    set(ax, 'FontSize', 16, 'FontWeight', 'bold');
end

for k = 1:8
    ax = subplot(4,4,8+k);
    hold(ax, 'on');
    v = plot_post(:,k);
    histogram(ax, v, 50, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
        'FaceColor', [0.43921568627451 0.6 0.929411764705882], ...
        'DisplayName', 'Posterior Samples');
    [xi, f] = local_gaussian_kde(v, 600);
    plot(ax, xi, f, 'r-', 'LineWidth', 1.8, 'DisplayName', 'PDF');
    xlabel(ax, labels{k}, 'FontWeight', 'bold', 'Interpreter', 'tex');
    ylabel(ax, 'pdf', 'FontWeight', 'bold');
    grid(ax, 'on'); %box(ax, 'on');
    set(ax, 'FontSize', 16, 'FontWeight', 'bold');
    if k == 1
        lg = legend(ax, 'show');
        set(lg, 'FontSize', 20, 'Orientation', 'horizontal', 'EdgeColor', [1 1 1]);
    end
    hold(ax, 'off');
end

letters = arrayfun(@(i) sprintf('%c)', char('a'+i-1)), 1:16, 'UniformOutput', false);
for i = 1:16
    ax = subplot(4,4,i);
    text(ax, -0.18, 1.08, letters{i}, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 24, 'Color', [0 0 0]);
end
end

function p = local_percentile(x, probs)
% Toolbox-free percentile calculation along columns.
x = sort(x, 1);
n = size(x, 1);
n_cols = size(x, 2);
p = zeros(numel(probs), n_cols);

for i = 1:numel(probs)
    pos = 1 + (n - 1) * probs(i) / 100;
    lo = floor(pos);
    hi = ceil(pos);
    w = pos - lo;

    lo = max(1, min(n, lo));
    hi = max(1, min(n, hi));

    p(i,:) = (1 - w) * x(lo,:) + w * x(hi,:);
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

% Keep plotting fast for long MCMC chains without changing saved samples.
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
