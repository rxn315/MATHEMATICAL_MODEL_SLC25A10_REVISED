% ============================================================
% Metropolis–Hastings MCMC for UPDATED SLC25A10 REA ping–pong model
% (TP REMOVED) + (Ks_m REMOVED)
%
% Sampling directly in EXPERIMENTAL flux units:
%   Flux units: mmol/min/g protein  (EXPERIMENTAL + likelihood space)
%   Concentrations: mM
%
% Updated model changes:
%   - Unified total transporter conservation across all states
%   - Saturation polynomials include ALL substrates on each side
%   - Denominator uses (delta1 + phi1*delta2)
%
% THIS VERSION:
%   - Removes Ks_m from parameter vector and from model algebra.
%   - Assumes Sm = Sims = 0 in these datasets, so succinate terms drop out.
%
% Parameters sampled (positive):
%   [1] Tm_max_exp   (mmol/min/g)   malate Tmax in exp units
%   [2] Km_m         (mM)
%   [3] Kp_m         (mM)
%   [4] lambda21     (dimensionless) lambda_ms
%   [5] lambda31     (dimensionless) lambda_ps
%   [6] sigma2       (mmol/min/g)^2
% ============================================================
clear; clc; close all;
rng(42);

%% ----------------- Constants ---------------------

fact_model_to_exp = 60 * 1e3;  % mol/s/g -> mmol/min/g

%% ----------------- Step 1: Load Experimental Data -----------------
% Keep experimental rates in mmol/min/g.

% ---------- Older MALATE flux dataset: intPho_exMAL (M0) ----------
load('Proteoliposome_intPho_exMAL.mat');   % intPho_exMal, intPho_exMal_Rate
n_m0        = numel(intPho_exMal);
Pm_v_mal0   = 15 * ones(1, n_m0);          % internal phosphate 15 mM
Mc_v_mal0   = intPho_exMal(:)';            % external malate (mM)
J_mal_exp0  = intPho_exMal_Rate(:)';       % [mmol/min/g]
Sm_v_mal0   = 0 * ones(1, n_m0);
Mm_v_mal0   = 0 * ones(1, n_m0);
Pc_v_mal0   = 0 * ones(1, n_m0);
Sc_v_mal0   = 0 * ones(1, n_m0);

% ---------- Older PHOSPHATE-associated dataset: intMAl_exPHO (P1) ----------
load('Proteoliposome_intMAl_exPHO.mat');   % intMal_exPho, intMal_exPho_Rate
n_p1        = numel(intMal_exPho);
Pm_v_pho1   = 0  * ones(1, n_p1);          % matrix phosphate (mM)
Sm_v_pho1   = 0  * ones(1, n_p1);
Mm_v_pho1   = 15 * ones(1, n_p1);          % internal malate 15 mM
Pc_v_pho1   = intMal_exPho(:)';            % IMS phosphate varied (mM)
Mc_v_pho1   = 0  * ones(1, n_p1);
Sc_v_pho1   = 0  * ones(1, n_p1);
J_pho_exp1  = intMal_exPho_Rate(:)';       % [mmol/min/g]

% ---------- 1993 MALATE flux: P_in = 1 mM (M2) ----------
load('Proteoliposome_1993_intPho_10mM_exMAL.mat');   % intPho_10mM_exMa, intPho_10mM_exMa_Rate
n_m2        = numel(intPho_10mM_exMa);
Pm_v_mal2   = 1.0 * ones(1, n_m2);
Mc_v_mal2   = intPho_10mM_exMa(:)';
J_mal_exp2  = intPho_10mM_exMa_Rate(:)';   % [mmol/min/g]
Sm_v_mal2   = 0 * ones(1, n_m2);
Mm_v_mal2   = 0 * ones(1, n_m2);
Pc_v_mal2   = 0 * ones(1, n_m2);
Sc_v_mal2   = 0 * ones(1, n_m2);

% ---------- 1993 MALATE flux: P_in = 2.5 mM (M3) ----------
load('Proteoliposome_1993_intPho_25mM_exMAL.mat');   % intPho_25mM_exMa, intPho_25mM_exMa_Rate
n_m3        = numel(intPho_25mM_exMa);
Pm_v_mal3   = 2.5 * ones(1, n_m3);
Mc_v_mal3   = intPho_25mM_exMa(:)';
J_mal_exp3  = intPho_25mM_exMa_Rate(:)';   % [mmol/min/g]
Sm_v_mal3   = 0 * ones(1, n_m3);
Mm_v_mal3   = 0 * ones(1, n_m3);
Pc_v_mal3   = 0 * ones(1, n_m3);
Sc_v_mal3   = 0 * ones(1, n_m3);

% ---------- NEW 1993 MALATE flux: M_out=0.2 mM, P_in varies (M6) ----------
load('Proteoliposome_1993_exMal_020mM_inPho.mat');   % exMal_020mM_inPho, exMal_020mM_inPho_Rate
n_m6        = numel(exMal_020mM_inPho);
Pm_v_mal6   = exMal_020mM_inPho(:)';       % internal phosphate varies
Mc_v_mal6   = 0.2 * ones(1, n_m6);
J_mal_exp6  = exMal_020mM_inPho_Rate(:)';  % [mmol/min/g]
Sm_v_mal6   = 0 * ones(1, n_m6);
Mm_v_mal6   = 0 * ones(1, n_m6);
Pc_v_mal6   = 0 * ones(1, n_m6);
Sc_v_mal6   = 0 * ones(1, n_m6);

% ---------- NEW 1993 MALATE flux: M_out=1.0 mM, P_in varies (M8) ----------
load('Proteoliposome_1993_exMal_1mM_inPho.mat');     % exMal_1mM_inPho, exMal_1mM_inPho_Rate
n_m8        = numel(exMal_1mM_inPho);
Pm_v_mal8   = exMal_1mM_inPho(:)';         % internal phosphate varies
Mc_v_mal8   = 1.0 * ones(1, n_m8);
J_mal_exp8  = exMal_1mM_inPho_Rate(:)';    % [mmol/min/g]
Sm_v_mal8   = 0 * ones(1, n_m8);
Mm_v_mal8   = 0 * ones(1, n_m8);
Pc_v_mal8   = 0 * ones(1, n_m8);
Sc_v_mal8   = 0 * ones(1, n_m8);

% ---------- Concatenate MALATE datasets ----------
Mc_v_mal  = [Mc_v_mal0,  Mc_v_mal2,  Mc_v_mal3,  Mc_v_mal6,  Mc_v_mal8];
Sm_v_mal  = [Sm_v_mal0,  Sm_v_mal2,  Sm_v_mal3,  Sm_v_mal6,  Sm_v_mal8];
Sc_v_mal  = [Sc_v_mal0,  Sc_v_mal2,  Sc_v_mal3,  Sc_v_mal6,  Sc_v_mal8];
Mm_v_mal  = [Mm_v_mal0,  Mm_v_mal2,  Mm_v_mal3,  Mm_v_mal6,  Mm_v_mal8];
Pm_v_mal  = [Pm_v_mal0,  Pm_v_mal2,  Pm_v_mal3,  Pm_v_mal6,  Pm_v_mal8];
Pc_v_mal  = [Pc_v_mal0,  Pc_v_mal2,  Pc_v_mal3,  Pc_v_mal6,  Pc_v_mal8];
J_mal_exp = [J_mal_exp0, J_mal_exp2, J_mal_exp3, J_mal_exp6, J_mal_exp8];

idx_m0 = 1:n_m0;
idx_m2 = idx_m0(end) + (1:n_m2);
idx_m3 = idx_m2(end) + (1:n_m3);
idx_m6 = idx_m3(end) + (1:n_m6);
idx_m8 = idx_m6(end) + (1:n_m8);

% ---------- Phosphate dataset ----------
Mc_v_pho  = Mc_v_pho1;
Sm_v_pho  = Sm_v_pho1;
Sc_v_pho  = Sc_v_pho1;
Mm_v_pho  = Mm_v_pho1;
Pm_v_pho  = Pm_v_pho1;
Pc_v_pho  = Pc_v_pho1;
J_pho_exp = J_pho_exp1;
idx_p1    = 1:numel(J_pho_exp);

data = struct( ...
   'Mc_v_mal',Mc_v_mal,'Sm_v_mal',Sm_v_mal,'Sc_v_mal',Sc_v_mal,'Mm_v_mal',Mm_v_mal,'Pm_v_mal',Pm_v_mal,'Pc_v_mal',Pc_v_mal,'J_mal_exp',J_mal_exp, ...
   'Mc_v_pho',Mc_v_pho,'Sm_v_pho',Sm_v_pho,'Sc_v_pho',Sc_v_pho,'Mm_v_pho',Mm_v_pho,'Pm_v_pho',Pm_v_pho,'Pc_v_pho',Pc_v_pho,'J_pho_exp',J_pho_exp, ...
   'idx_m0',idx_m0,'idx_m2',idx_m2,'idx_m3',idx_m3,'idx_m6',idx_m6,'idx_m8',idx_m8,'idx_p1',idx_p1);

n_total = numel(J_mal_exp) + numel(J_pho_exp);

%% ----------------- Step 2: Priors -----------------
% params = [Tm_max, Km_m, Kp_m, lambda21, lambda31]
lb_Tmax = 1e-4;
ub_Tmax = 2e1;

lb_vec = [lb_Tmax; 0; 0];   % [Tm_max, Km, Kp]
ub_vec = [ub_Tmax; 1.5; 30];

a_c = 0.01;  b_c = 0.15;
alpha_0 = 1; beta_0 = 1;

log_prior_all = @(params, sigma2) log_prior_mix(params, sigma2, lb_vec, ub_vec, a_c, b_c, alpha_0, beta_0);

log_likelihood = @(params, sigma2) ...
    -n_total/2*log(2*pi*sigma2) - 0.5/sigma2 * obj_sse_exp(params, data, fact_model_to_exp);

%% ----------------- Step 3: Metropolis-Hastings --------------------
theta = [ ...
    1;     ... % Tm_max_exp (mmol/min/g)
    1;     ... % Km_m (mM)
    1;     ... % Kp_m (mM)
    0.15;  ... % lambda21
    0.15;   ... % lambda31
    4e-2]; ... % sigma2

n_samples    = 1000000;
burn_in      = 10000;
thin         = 1;
z_step       = 0.01;
report_every = 2000;

p_all = numel(theta); % 6
theta_samples_raw = zeros(n_samples, p_all);
acc = false(n_samples,1);

loglik_curr   = log_likelihood(theta(1:5), theta(6));
logprior_curr = log_prior_all(theta(1:5), theta(6));

for i = 1:n_samples
    log_theta_prop = log(theta) + z_step * randn(p_all,1);
    theta_prop = exp(log_theta_prop);

    if any(theta_prop(1:3) < lb_vec | theta_prop(1:3) > ub_vec) || ...
       any(theta_prop(4:5) < a_c   | theta_prop(4:5) > b_c)   || ...
       theta_prop(6) <= 0
        theta_samples_raw(i,:) = theta.';
        continue;
    end

    loglik_prop   = log_likelihood(theta_prop(1:5), theta_prop(6));
    logprior_prop = log_prior_all(theta_prop(1:5), theta_prop(6));
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
        fprintf('Iter %7d | acc(last %d)=%.3f | acc(all)=%.3f | sigma2=%.4g\n', ...
            i, report_every, mean(acc(i-report_every+1:i)), mean(acc(1:i)), theta(6));
    end
end
fprintf('\nOverall acceptance rate: %.3f\n', mean(acc));

keep_idx      = (burn_in+1):thin:n_samples;
theta_samples = theta_samples_raw(keep_idx,:);

%% ----------------- Step 4: Posterior Predictions -----------------
n_post = size(theta_samples,1);
n_mal  = numel(data.J_mal_exp);
n_pho  = numel(data.J_pho_exp);

J_M_pred = zeros(n_post, n_mal);
J_P_pred = zeros(n_post, n_pho);

for i = 1:n_post
    params = theta_samples(i,1:5);
    [Jm_exp, Jp_exp] = model_fluxes_exp_UPDATED(params, data, fact_model_to_exp);
    J_M_pred(i,:) = Jm_exp;
    J_P_pred(i,:) = Jp_exp;
end

mean_J_M = mean(J_M_pred,1);
CI_J_M   = prctile(J_M_pred,[2.5,97.5],1);

mean_J_P = mean(J_P_pred,1);
CI_J_P   = prctile(J_P_pred,[2.5,97.5],1);

%% ----------------- Step 5: Plots (same layout as before) -----------------
figure('Color','w','Name','Posterior bands (Ks removed; TP removed)');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% M0
nexttile; hold on;
idx = data.idx_m0;
fill([data.Mc_v_mal(idx), fliplr(data.Mc_v_mal(idx))], ...
     [CI_J_M(1,idx),      fliplr(CI_J_M(2,idx))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none');
plot(data.Mc_v_mal(idx), mean_J_M(idx), 'r--','LineWidth',2);
plot(data.Mc_v_mal(idx), data.J_mal_exp(idx), 'b-','LineWidth',2);
xlabel('External malate (mM)'); ylabel('J_M (mmol min^{-1} g^{-1})');
title('Malate: intPho\_exMAL (P_{in}=15 mM)');
grid on; hold off;

% M2
nexttile; hold on;
idx = data.idx_m2;
fill([data.Mc_v_mal(idx), fliplr(data.Mc_v_mal(idx))], ...
     [CI_J_M(1,idx),      fliplr(CI_J_M(2,idx))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none');
plot(data.Mc_v_mal(idx), mean_J_M(idx), 'r--','LineWidth',2);
plot(data.Mc_v_mal(idx), data.J_mal_exp(idx), 'b-','LineWidth',2);
xlabel('External malate (mM)'); ylabel('J_M (mmol min^{-1} g^{-1})');
title('Malate: 1993 (P_{in}=1.0 mM)');
grid on; hold off;

% M3
nexttile; hold on;
idx = data.idx_m3;
fill([data.Mc_v_mal(idx), fliplr(data.Mc_v_mal(idx))], ...
     [CI_J_M(1,idx),      fliplr(CI_J_M(2,idx))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none');
plot(data.Mc_v_mal(idx), mean_J_M(idx), 'r--','LineWidth',2);
plot(data.Mc_v_mal(idx), data.J_mal_exp(idx), 'b-','LineWidth',2);
xlabel('External malate (mM)'); ylabel('J_M (mmol min^{-1} g^{-1})');
title('Malate: 1993 (P_{in}=2.5 mM)');
grid on; hold off;

% Phosphate P1
nexttile; hold on;
idx = data.idx_p1;
fill([data.Pc_v_pho(idx), fliplr(data.Pc_v_pho(idx))], ...
     [CI_J_P(1,idx),      fliplr(CI_J_P(2,idx))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none');
plot(data.Pc_v_pho(idx), mean_J_P(idx), 'r--','LineWidth',2);
plot(data.Pc_v_pho(idx), data.J_pho_exp(idx), 'b-','LineWidth',2);
xlabel('IMS phosphate (mM)'); ylabel('J_P (mmol min^{-1} g^{-1})');
title('Phosphate-associated: intMAl\_exPHO');
grid on; hold off;

figure('Color','w','Name','Posterior bands - exMal_inPho (Ks removed; TP removed)');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% M6: M_out = 0.2 mM, P_in varies
nexttile; hold on;
idx = data.idx_m6;
pin = data.Pm_v_mal(idx);
fill([pin, fliplr(pin)], ...
     [CI_J_M(1,idx), fliplr(CI_J_M(2,idx))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
plot(pin, mean_J_M(idx), 'r--','LineWidth',2,'DisplayName','Posterior mean');
plot(pin, data.J_mal_exp(idx), 'b-','LineWidth',2,'DisplayName','Experimental');
xlabel('Internal phosphate (mM)');
ylabel('J_M (mmol min^{-1} g^{-1})');
title('Malate: exMal=0.2 mM, inPho varying');
legend('Location','best'); grid on; hold off;

% M8: M_out = 1.0 mM, P_in varies
nexttile; hold on;
idx = data.idx_m8;
pin = data.Pm_v_mal(idx);
fill([pin, fliplr(pin)], ...
     [CI_J_M(1,idx), fliplr(CI_J_M(2,idx))], ...
     'r','FaceAlpha',0.2,'EdgeColor','none','DisplayName','95% CI');
plot(pin, mean_J_M(idx), 'r--','LineWidth',2,'DisplayName','Posterior mean');
plot(pin, data.J_mal_exp(idx), 'b-','LineWidth',2,'DisplayName','Experimental');
xlabel('Internal phosphate (mM)');
ylabel('J_M (mmol min^{-1} g^{-1})');
title('Malate: exMal=1.0 mM, inPho varying');
legend('Location','best'); grid on; hold off;

%% ----------------- Step 6: Posterior summary (with 95% CI) -----------------
posterior_mean = mean(theta_samples,1);
posterior_std  = std(theta_samples,0,1);
posterior_ci   = prctile(theta_samples,[2.5 97.5],1);

fprintf('\nPosterior summary (Ks removed; TP removed):\n');
fprintf('Tm_max_exp  = %.6g  (std %.3g)  CI[%.6g, %.6g]   mmol/min/g\n', posterior_mean(1),posterior_std(1),posterior_ci(1,1),posterior_ci(2,1));
fprintf('Km_m        = %.6g  (std %.3g)  CI[%.6g, %.6g]   mM\n',        posterior_mean(2),posterior_std(2),posterior_ci(1,2),posterior_ci(2,2));
fprintf('Kp_m        = %.6g  (std %.3g)  CI[%.6g, %.6g]   mM\n',        posterior_mean(3),posterior_std(3),posterior_ci(1,3),posterior_ci(2,3));
fprintf('lambda21    = %.6g  (std %.3g)  CI[%.6g, %.6g]\n',             posterior_mean(4),posterior_std(4),posterior_ci(1,4),posterior_ci(2,4));
fprintf('lambda31    = %.6g  (std %.3g)  CI[%.6g, %.6g]\n',             posterior_mean(5),posterior_std(5),posterior_ci(1,5),posterior_ci(2,5));
fprintf('sigma2      = %.6g  (std %.3g)  CI[%.6g, %.6g]   (mmol/min/g)^2\n', posterior_mean(6),posterior_std(6),posterior_ci(1,6),posterior_ci(2,6));

%% ----------------- Step 7: Prior vs Posterior ----------------------
n_prior = 10000;

Tm_max_prior   = unifrnd(lb_vec(1), ub_vec(1), n_prior, 1);
Km_m_prior     = unifrnd(lb_vec(2), ub_vec(2), n_prior, 1);
Kp_m_prior     = unifrnd(lb_vec(3), ub_vec(3), n_prior, 1);
lambda21_prior = unifrnd(a_c, b_c, n_prior, 1);
lambda31_prior = unifrnd(a_c, b_c, n_prior, 1);
sigma2_prior   = 1 ./ gamrnd(alpha_0, 1/beta_0, n_prior, 1);

post = theta_samples;

labels = { ...
    'T^{m}_{max} (mmol min^{-1} g^{-1})', ...
    'K^{m}_m (mM)', ...
    'K^{p}_m (mM)', ...
    '\lambda_{21}', ...
    '\lambda_{31}', ...
    '\sigma^2 ((mmol min^{-1} g^{-1})^2)'};

priors = {Tm_max_prior, Km_m_prior, Kp_m_prior, lambda21_prior, lambda31_prior, sigma2_prior};

figure('Color','w','Name','Prior vs Posterior (Ks removed; TP removed)');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for k = 1:6
    nexttile; hold on;
    histogram(priors{k}, 60, 'Normalization','pdf', ...
        'FaceColor',[0.2 0.5 1], 'EdgeColor','none', 'DisplayName','Prior');
    histogram(post(:,k), 60, 'Normalization','pdf', ...
        'FaceColor',[1 0.3 0.3], 'EdgeColor','none', 'FaceAlpha',0.7, ...
        'DisplayName','Posterior');
    title(labels{k},'Interpreter','tex','FontWeight','bold');
    legend('Location','best'); grid on; hold off;
end

% ========================= Helper Functions =========================
function sse = obj_sse_exp(params, data, fact_model_to_exp)
[Jm_exp, Jp_exp] = model_fluxes_exp_UPDATED(params, data, fact_model_to_exp);
sse = sum((Jm_exp - data.J_mal_exp).^2) + sum((Jp_exp - data.J_pho_exp).^2);
end

function lp = log_prior_mix(params, sigma2, lb_vec, ub_vec, a_c, b_c, alpha_0, beta_0)
if sigma2 <= 0
    lp = -inf; return;
end
% params = [Tm_max, Km, Kp, lambda21, lambda31]
if any(params(1:3) < lb_vec | params(1:3) > ub_vec)
    lp = -inf; return;
end
if any(params(4:5) < a_c | params(4:5) > b_c)
    lp = -inf; return;
end
lp = log_inv_gamma(sigma2, alpha_0, beta_0);
end

function v = log_inv_gamma(s2, alpha, beta)
if s2 <= 0
    v = -inf; return;
end
v = (-alpha-1)*log(s2) - beta/s2 - gammaln(alpha) + alpha*log(beta);
end

function [Jm_exp, Jp_exp] = model_fluxes_exp_UPDATED(params, data, fact_model_to_exp)
Tm_max_exp = params(1);
Km_m       = params(2);
Kp_m       = params(3);
lambda21   = params(4);
lambda31   = params(5);

Tm_max_surface = Tm_max_exp / fact_model_to_exp;  % mol/s/g

% ---- Malate dataset ----
n_mal = numel(data.J_mal_exp);
Jm_surface = zeros(1,n_mal);
for i = 1:n_mal
    Mm   = data.Mm_v_mal(i);  Mims = data.Mc_v_mal(i);
    Pm   = data.Pm_v_mal(i);  Pims = data.Pc_v_mal(i);
    [~,Jm,~] = compute_flux_condition_REA_UPDATED( ...
        Tm_max_surface,Km_m,Kp_m,lambda21,lambda31, ...
        Mm,Mims,Pm,Pims);
    Jm_surface(i) = Jm;
end

% ---- Phosphate dataset ----
n_pho = numel(data.J_pho_exp);
Jp_surface = zeros(1,n_pho);
for i = 1:n_pho
    Mm   = data.Mm_v_pho(i);  Mims = data.Mc_v_pho(i);
    Pm   = data.Pm_v_pho(i);  Pims = data.Pc_v_pho(i);
    [~,~,Jp] = compute_flux_condition_REA_UPDATED( ...
        Tm_max_surface,Km_m,Kp_m,lambda21,lambda31, ...
        Mm,Mims,Pm,Pims);
    Jp_surface(i) = Jp;
end

Jm_exp = Jm_surface * fact_model_to_exp;
Jp_exp = Jp_surface * fact_model_to_exp;
end

function [J_succ,J_mal,J_pho] = compute_flux_condition_REA_UPDATED( ...
    Tm_max,Km_m,Kp_m,lambda21,lambda31, ...
    Mm,Mims,Pm,Pims)

delta1 = 1 + Mm/Km_m + Pm/Kp_m;
delta2 = 1 + Mims/Km_m + Pims/Kp_m;

% NOTE: with Sm=Sims=0, phi1 uses only (lambda21*Mm + lambda31*Pm) terms
phi1 = (lambda21*Mm + lambda31*Pm) / (lambda21*Mims + lambda31*Pims);

Den = (delta1 + phi1*delta2);



v1r_m = Tm_max * (phi1*Mims - Mm) / (Km_m * Den);
v1f_m = Tm_max * (Mm - phi1*Mims) / (Km_m * Den);

J_succ = 0;
J_mal  = v1r_m;
J_pho  = v1f_m;
end
