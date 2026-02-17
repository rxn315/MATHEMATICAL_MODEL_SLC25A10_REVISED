% ============================================================
% MCMC for Palmieri DIC data (succinate + malate competition + phosphate)
% UPDATED REA-based ping–pong mechanism (surface-form algebra)
%
% THIS VERSION (CHANGE REQUEST):
%   - Fix (do NOT estimate):
%       Pm_mat, Sm_mat, Mm_mat (set below)
%   - Estimate:
%       [Ts_max_exp, Tm_max_exp, Ks_m, Km_m]  : Gaussian priors (positive)
%       [Kp_m]                                : Gaussian prior (positive)  <-- CHANGED
%       [lambda21, lambda31]                  : Uniform priors
%       [sigma^2]                             : Inv-Gamma prior
%   - TP_max_exp remains REMOVED
%
% Units:
%   Concentrations: mM
%   Model flux internal: mol/s/g
%   Experimental flux (and SSE/likelihood): µmol/min/g
%
% Estimated parameters (theta):
%   [1]  Ts_max_exp   (µmol/min/g)   (Gaussian, >0)
%   [2]  Tm_max_exp   (µmol/min/g)   (Gaussian, >0)
%   [3]  Ks_m         (mM)           (Gaussian, >0)
%   [4]  Km_m         (mM)           (Gaussian, >0)
%   [5]  Kp_m         (mM)           (Gaussian, >0)  mean=0.93, sd=0.40
%   [6]  lambda21     (dimensionless) (UNIFORM)
%   [7]  lambda31     (dimensionless) (UNIFORM)
%   [8]  sigma^2      ((µmol/min/g)^2) (Inv-Gamma)
% ============================================================

clear; clc; close all;
rng(42);

%% ----------------- Step 0: Geometry conversion ---------------------
fact_model_to_exp = 60 * 1e6;   % mol/s/g -> µmol/min/g

%% ----------------- Fixed matrix concentrations ----------------------
Pm_mat_fix = 2;     % mM
Sm_mat_fix = 0.2;   % mM
Mm_mat_fix = 0.2;   % mM

%% ----------------- Step 1: Load Experimental Data (NO conversion) ---
load('Palmier_DIC_Exp_Data.mat');           % provides S, T_S
Sc_v_succ  = S(:)';                          % [S_ims] varied (mM)
Mc_v_succ  = zeros(size(Sc_v_succ));         % [M_ims] = 0
Pc_v_succ  = zeros(size(Sc_v_succ));         % [P_ims] = 0
J_succ_exp = T_S(:)';                        % µmol/min/g

load('Palmier_DIC_Competition_Exp_Data.mat'); % M_with_S, T_M_with_S
Mc_v_mal   = (1./M_with_S(:))';               % [M_ims] varied (mM)
Sc_v_mal   = 0.5 * ones(size(Mc_v_mal));      % [S_ims] fixed (mM)
Pc_v_mal   = zeros(size(Mc_v_mal));           % [P_ims] = 0
J_mal_exp  = (1./T_M_with_S(:))';             % µmol/min/g

load('Palmier_DIC_Pho_Exp_Data.mat');         % M_control, T_M_control
Mc_v_pho   = (1./M_control(:))';              % [M_ims] varied (mM)
Sc_v_pho   = zeros(size(Mc_v_pho));           % [S_ims] = 0
Pc_v_pho   = zeros(size(Mc_v_pho));           % [P_ims] = 0 (replace if varied)
J_pho_exp  = (1./T_M_control(:))';            % µmol/min/g

data = struct( ...
   'Mc_v_succ',Mc_v_succ,'Sc_v_succ',Sc_v_succ,'Pc_v_succ',Pc_v_succ,'J_succ_exp',J_succ_exp, ...
   'Mc_v_mal', Mc_v_mal, 'Sc_v_mal', Sc_v_mal, 'Pc_v_mal', Pc_v_mal, 'J_mal_exp', J_mal_exp, ...
   'Mc_v_pho', Mc_v_pho, 'Sc_v_pho', Sc_v_pho, 'Pc_v_pho', Pc_v_pho, 'J_pho_exp', J_pho_exp);

n_succ  = numel(J_succ_exp);
n_mal   = numel(J_mal_exp);
n_pho   = numel(J_pho_exp);
n_total = n_succ + n_mal + n_pho;

%% ----------------- Step 2: Priors -----------------------------------
% Gaussian priors for 1..4 (positive)
mu_vec    = [64.2; 69.04; 1.17; 0.23];
sigma_vec = [ 5.1;  5.60; 0.10; 0.02];

% Gaussian prior for Kp_m (positive)  <-- CHANGED
mu_kp    = 0.93;
sigma_kp = 0.40;

% Uniform priors for lambda21, lambda31
a_c = 0.01;  b_c = 100;

% Inv-Gamma for sigma^2 in EXP units
alpha_0 = 1;
beta_0  = 1;

log_prior_all = @(params7, sigma2) ...
    log_prior_gauss_gauss_uniform_fixedM(params7, sigma2, mu_vec, sigma_vec, mu_kp, sigma_kp, a_c, b_c, alpha_0, beta_0);

log_likelihood = @(params7, sigma2) ...
    -n_total/2*log(2*pi*sigma2) - 0.5/sigma2 * obj_sse(params7, data, fact_model_to_exp, Pm_mat_fix, Sm_mat_fix, Mm_mat_fix);

%% ----------------- Step 3: Metropolis–Hastings ----------------------
% theta = [params(1:7); sigma2]
% params(1:7) = [Ts_max_exp, Tm_max_exp, Ks_m, Km_m, Kp_m, lambda21, lambda31]
theta = [60; 60; 1.0; 0.3; 1.0; 0.1; 0.1; 4e-2];

n_samples    = 1000000;
burn_in      = 10000;
thin         = 1;
z_step       = 0.01;
report_every = 2000;

p_all = numel(theta);  % 8
theta_samples_raw = zeros(n_samples, p_all);
acc = false(n_samples,1);

loglik_curr   = log_likelihood(theta(1:7), theta(8));
logprior_curr = log_prior_all(theta(1:7), theta(8));

for i = 1:n_samples
    log_theta_prop = log(theta) + z_step * randn(p_all,1);
    theta_prop = exp(log_theta_prop);  % positivity

    % Support checks:
    % 1..5 >0 ; 6..7 in [a_c,b_c] ; 8>0
    if any(theta_prop(1:5) <= 0) || ...
       any(theta_prop(6:7) < a_c | theta_prop(6:7) > b_c) || ...
       theta_prop(8) <= 0
        theta_samples_raw(i,:) = theta.';
        continue;
    end

    loglik_prop   = log_likelihood(theta_prop(1:7), theta_prop(8));
    logprior_prop = log_prior_all(theta_prop(1:7), theta_prop(8));
    log_hastings  = sum(log(theta) - log(theta_prop));

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
        fprintf('Iter %7d | acc(last %d)=%.3f | acc(all)=%.3f | sigma^2=%.4g\n', ...
            i, report_every, mean(acc(i-report_every+1:i)), mean(acc(1:i)), theta(8));
    end
end
fprintf('\nOverall acceptance rate: %.3f\n', mean(acc));

keep_idx      = (burn_in+1):thin:n_samples;
theta_samples = theta_samples_raw(keep_idx,:);

%% ----------------- Step 4: Posterior Flux Predictions (EXP units) ---
n_post = size(theta_samples, 1);

J_S_pred  = zeros(n_post, n_succ);
J_M_pred  = zeros(n_post, n_mal);
J_P_pred  = zeros(n_post, n_pho);

for i = 1:n_post
    params7 = theta_samples(i, 1:7);
    fluxes  = forward_fluxes(params7, data, fact_model_to_exp, Pm_mat_fix, Sm_mat_fix, Mm_mat_fix);
    J_S_pred(i,:) = fluxes{1};
    J_M_pred(i,:) = fluxes{2};
    J_P_pred(i,:) = fluxes{3};
end

mean_J_S = mean(J_S_pred,1);
mean_J_M = mean(J_M_pred,1);
mean_J_P = mean(J_P_pred,1);

CI_J_S   = prctile(J_S_pred,[2.5,97.5],1);
CI_J_M   = prctile(J_M_pred,[2.5,97.5],1);
CI_J_P   = prctile(J_P_pred,[2.5,97.5],1);

%% ----------------- Step 5: Flux Plots (EXP units) -------------------
figure('Color','w','Name','Palmieri DIC: Posterior credible bands (EXP units)');

subplot(1,3,1); hold on;
x = data.Sc_v_succ;
fill([x, fliplr(x)], [CI_J_S(1,:), fliplr(CI_J_S(2,:))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
plot(x, mean_J_S, 'r--','LineWidth',2,'DisplayName','Posterior mean');
plot(x, data.J_succ_exp, 'b-','LineWidth',2,'DisplayName','Experimental');
xlabel('[S_{ims}] (mM)'); ylabel('J_S (\mumol min^{-1} g^{-1})');
title('Succinate flux'); legend('Location','best'); grid on; hold off;

subplot(1,3,2); hold on;
x = data.Mc_v_mal;
fill([x, fliplr(x)], [CI_J_M(1,:), fliplr(CI_J_M(2,:))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
plot(x, mean_J_M, 'r--','LineWidth',2,'DisplayName','Posterior mean');
plot(x, data.J_mal_exp, 'b-','LineWidth',2,'DisplayName','Experimental');
xlabel('[M_{ims}] (mM)'); ylabel('J_M (\mumol min^{-1} g^{-1})');
title('Malate flux (S_{ims}=0.5 mM)'); legend('Location','best'); grid on; hold off;

subplot(1,3,3); hold on;
x = data.Mc_v_pho;
fill([x, fliplr(x)], [CI_J_P(1,:), fliplr(CI_J_P(2,:))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
plot(x, mean_J_P, 'r--','LineWidth',2,'DisplayName','Posterior mean');
plot(x, data.J_pho_exp, 'b-','LineWidth',2,'DisplayName','Experimental');
xlabel('[M_{ims}] (mM)'); ylabel('J_P (\mumol min^{-1} g^{-1})');
title('Phosphate-associated flux'); legend('Location','best'); grid on; hold off;

%% ----------------- Step 6: Prior draws (for plotted params) ---------
n_prior = 10000;

% Truncated-at-zero Gaussians for 1..4
Ts_pd = truncate(makedist('Normal','mu',mu_vec(1),'sigma',sigma_vec(1)), 0, Inf);
Tm_pd = truncate(makedist('Normal','mu',mu_vec(2),'sigma',sigma_vec(2)), 0, Inf);
Ks_pd = truncate(makedist('Normal','mu',mu_vec(3),'sigma',sigma_vec(3)), 0, Inf);
Km_pd = truncate(makedist('Normal','mu',mu_vec(4),'sigma',sigma_vec(4)), 0, Inf);

Ts_max_prior = random(Ts_pd, n_prior, 1);
Tm_max_prior = random(Tm_pd, n_prior, 1);
Ks_m_prior   = random(Ks_pd, n_prior, 1);
Km_m_prior   = random(Km_pd, n_prior, 1);

% Truncated-at-zero Gaussian for Kp_m  <-- CHANGED
Kp_pd      = truncate(makedist('Normal','mu',mu_kp,'sigma',sigma_kp), 0, Inf);
Kp_m_prior = random(Kp_pd, n_prior, 1);

lambda21_prior = unifrnd(a_c, b_c, n_prior, 1);
lambda31_prior = unifrnd(a_c, b_c, n_prior, 1);
sigma2_prior   = 1 ./ gamrnd(alpha_0, 1/beta_0, n_prior, 1);

%% ----------------- Step 7: Prior vs Posterior -----------------------
post   = theta_samples;
labels = {'T^{s}_{max,exp}','T^{m}_{max,exp}','K^{s}_m','K^{m}_m','K^{p}_m', ...
          '\lambda_{21}','\lambda_{31}','\sigma^2'};

priors = {Ts_max_prior, Tm_max_prior, Ks_m_prior, Km_m_prior, Kp_m_prior, ...
          lambda21_prior, lambda31_prior, sigma2_prior};

figure('Color','w','Name','Prior vs Posterior (Pm,Sm,Mm fixed; Kp Gaussian)');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
for k = 1:8
    nexttile; hold on;
    histogram(priors{k}, 60,'Normalization','pdf', ...
        'FaceColor',[0.2 0.5 1],'EdgeColor','none','DisplayName','Prior');
    histogram(post(:,k), 60,'Normalization','pdf', ...
        'FaceColor',[1 0.3 0.3],'EdgeColor','none', ...
        'FaceAlpha',0.7,'DisplayName','Posterior');
    title(labels{k},'Interpreter','tex','FontWeight','bold');
    legend('Location','best'); grid on; hold off;
end

%% ----------------- Step 8: Posterior summaries ----------------------
posterior_mean = mean(theta_samples,1);
posterior_std  = std(theta_samples,0,1);

fprintf(['\nFixed: Pm_mat=%.3f mM, Sm_mat=%.3f mM, Mm_mat=%.3f mM\n'], ...
    Pm_mat_fix, Sm_mat_fix, Mm_mat_fix);

fprintf(['\nPosterior means (EXP sampling; capacities in µmol/min/g):\n', ...
         'Ts_max=%.3f, Tm_max=%.3f, Ks_m=%.3f, Km_m=%.3f, Kp_m=%.3f,\n', ...
         'lambda21=%.3f, lambda31=%.3f, sigma2=%.3g\n'], ...
         posterior_mean(1),posterior_mean(2),posterior_mean(3),posterior_mean(4), ...
         posterior_mean(5),posterior_mean(6),posterior_mean(7),posterior_mean(8));

fprintf(['Posterior stds:\n', ...
         'Ts_max=%.3f, Tm_max=%.3f, Ks_m=%.3f, Km_m=%.3f, Kp_m=%.3f,\n', ...
         'lambda21=%.3f, lambda31=%.3f, sigma2=%.3g\n'], ...
         posterior_std(1),posterior_std(2),posterior_std(3),posterior_std(4), ...
         posterior_std(5),posterior_std(6),posterior_std(7),posterior_std(8));

%% ----------------- Step 9: Traces ----------------------------------
figure('Color','w','Name','MCMC traces (Pm,Sm,Mm fixed; Kp Gaussian)');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
for k = 1:8
    nexttile;
    plot(theta_samples_raw(:,k),'LineWidth',0.8);
    xlabel('iteration');
    ylabel(labels{k},'Interpreter','tex');
    grid on;
end

% ============================================================
% Helper functions
% ============================================================

function sse = obj_sse(params7, data, fact_model_to_exp, Pm_fix, Sm_fix, Mm_fix)
fluxes_exp = obj_func_palmier(params7, data, fact_model_to_exp, Pm_fix, Sm_fix, Mm_fix, true);
J_S = fluxes_exp{1};
J_M = fluxes_exp{2};
J_P = fluxes_exp{3};
sse = sum((J_S - data.J_succ_exp).^2) + sum((J_M - data.J_mal_exp).^2) + sum((J_P - data.J_pho_exp).^2);
end

function fluxes_exp = forward_fluxes(params7, data, fact_model_to_exp, Pm_fix, Sm_fix, Mm_fix)
fluxes_exp = obj_func_palmier(params7, data, fact_model_to_exp, Pm_fix, Sm_fix, Mm_fix, true);
end

function lp = log_prior_gauss_gauss_uniform_fixedM(params7, sigma2, mu_vec, sigma_vec, mu_kp, sigma_kp, a_c, b_c, alpha_0, beta_0)
% params7 = [Ts_max_exp, Tm_max_exp, Ks_m, Km_m, Kp_m, lambda21, lambda31]
if sigma2 <= 0
    lp = -inf; return;
end

% Support checks
if any(params7(1:5) <= 0)
    lp = -inf; return;
end
if any(params7(6:7) < a_c | params7(6:7) > b_c)
    lp = -inf; return;
end

% Gaussian log-prior for 1..4 (truncation-at-zero normalizing constants omitted)
x = params7(1:4);
lp_gauss_14 = -0.5 * sum(((x - mu_vec) ./ sigma_vec).^2) ...
             - sum(log(sigma_vec)) - numel(x)/2 * log(2*pi);

% Gaussian log-prior for Kp_m (param 5)  <-- CHANGED
kp = params7(5);
lp_kp = -0.5 * ((kp - mu_kp)/sigma_kp)^2 - log(sigma_kp) - 0.5*log(2*pi);

% Uniform lambdas: constant on support
lp_sig2 = log_inv_gamma(sigma2, alpha_0, beta_0);

lp = lp_gauss_14 + lp_kp + lp_sig2;
end

function v = log_inv_gamma(s2, alpha, beta)
if s2 <= 0
    v = -inf; return;
end
v = (-alpha-1)*log(s2) - beta/s2 - gammaln(alpha) + alpha*log(beta);
end

function out = obj_func_palmier(params7, data, fact_model_to_exp, Pm_fix, Sm_fix, Mm_fix, return_fluxes)
if nargin < 7
    return_fluxes = true;
end

Ts_max_exp = params7(1);    % µmol/min/g
Tm_max_exp = params7(2);    % µmol/min/g
Ks_m       = params7(3);    % mM
Km_m       = params7(4);    % mM
Kp_m       = params7(5);    % mM
lambda21   = params7(6);    % λ_ms
lambda31   = params7(7);    % λ_ps

Pm_mat = Pm_fix; Sm_mat = Sm_fix; Mm_mat = Mm_fix;

Ts_max_surf = Ts_max_exp / fact_model_to_exp; % mol/s/g
Tm_max_surf = Tm_max_exp / fact_model_to_exp; % mol/s/g

% Succinate dataset
n_succ = numel(data.Sc_v_succ);
J_S_surf = zeros(1,n_succ);
for i = 1:n_succ
    Sm   = Sm_mat;             Sims = data.Sc_v_succ(i);
    Mm   = Mm_mat;             Mims = data.Mc_v_succ(i);
    Pm   = Pm_mat;             Pims = data.Pc_v_succ(i);
    [Js,~,~] = compute_flux_condition_REA_SURF( ...
        Ts_max_surf,Tm_max_surf,Ks_m,Km_m,Kp_m,lambda21,lambda31, ...
        Sm,Sims,Mm,Mims,Pm,Pims);
    J_S_surf(i) = Js;
end

% Malate dataset
n_mal = numel(data.Mc_v_mal);
J_M_surf = zeros(1,n_mal);
for i = 1:n_mal
    Sm   = Sm_mat;             Sims = data.Sc_v_mal(i);
    Mm   = Mm_mat;             Mims = data.Mc_v_mal(i);
    Pm   = Pm_mat;             Pims = data.Pc_v_mal(i);
    [~,Jm,~] = compute_flux_condition_REA_SURF( ...
        Ts_max_surf,Tm_max_surf,Ks_m,Km_m,Kp_m,lambda21,lambda31, ...
        Sm,Sims,Mm,Mims,Pm,Pims);
    J_M_surf(i) = Jm;
end

% Phosphate dataset
n_pho = numel(data.Mc_v_pho);
J_P_surf = zeros(1,n_pho);
for i = 1:n_pho
    Sm   = Sm_mat;             Sims = data.Sc_v_pho(i);
    Mm   = Mm_mat;             Mims = data.Mc_v_pho(i);
    Pm   = Pm_mat;             Pims = data.Pc_v_pho(i);
    [~,~,Jp] = compute_flux_condition_REA_SURF( ...
        Ts_max_surf,Tm_max_surf,Ks_m,Km_m,Kp_m,lambda21,lambda31, ...
        Sm,Sims,Mm,Mims,Pm,Pims);
    J_P_surf(i) = Jp;
end

J_S_exp = J_S_surf * fact_model_to_exp;
J_M_exp = J_M_surf * fact_model_to_exp;
J_P_exp = J_P_surf * fact_model_to_exp;

if return_fluxes
    out = {J_S_exp, J_M_exp, J_P_exp};
else
    err_s = J_S_exp - data.J_succ_exp;
    err_m = J_M_exp - data.J_mal_exp;
    err_p = J_P_exp - data.J_pho_exp;
    out   = sum(err_s.^2 + err_m.^2 + err_p.^2);
end
end

function [J_succ,J_mal,J_pho] = compute_flux_condition_REA_SURF( ...
    Ts_max,Tm_max,Ks_m,Km_m,Kp_m,lambda21,lambda31, ...
    Sm,Sims,Mm,Mims,Pm,Pims)

delta1 = 1 + Sm/Ks_m   + Mm/Km_m   + Pm/Kp_m;
delta2 = 1 + Sims/Ks_m + Mims/Km_m + Pims/Kp_m;

phi1 = (Sm + lambda21*Mm + lambda31*Pm) / ...
       (Sims + lambda21*Mims + lambda31*Pims);

den = (delta1 + phi1*delta2);

v1r_s = Ts_max * (phi1*Sims - Sm)  / (Ks_m * den);
v1r_m = Tm_max * (phi1*Mims - Mm)  / (Km_m * den);

J_succ = v1r_s;
J_mal  = v1r_m;

% Keep your original convention:
J_pho  = v1r_m;
end
