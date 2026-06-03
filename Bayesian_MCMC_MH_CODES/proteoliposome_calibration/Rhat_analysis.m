%% ========================================================================
% Rhat analysis for reduced proteoliposome MCMC chains.
%
% Expected reduced chain files:
%   chain1_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
%   chain2_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
%   chain3_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
%   chain4_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
%
% theta_samples columns:
%   [Tmax_m, Tmax_p, K_D_m, K_D_p, sigma2]
%
% rho_pm is derived for diagnostics.
%% ========================================================================

clear; clc; close all;

chain_pattern = 'chain%d_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat';
n_chains = 4;

quantity_names = {'Tmax_m', 'Tmax_p', 'K_D_m', 'K_D_p', 'rho_pm_derived', 'sigma2'};

chains_cell = cell(n_chains, 1);
for c = 1:n_chains
    chain_file = sprintf(chain_pattern, c);
    if exist(chain_file, 'file') ~= 2
        error('Missing reduced chain file: %s', chain_file);
    end

    S = load(chain_file);
    if isfield(S, 'diagnostic_samples')
        diagnostic_samples = S.diagnostic_samples;
    elseif isfield(S, 'theta_samples')
        theta_samples = S.theta_samples;
        rho_pm = derived_rho_pm(theta_samples);
        diagnostic_samples = [theta_samples(:, 1:4), rho_pm, theta_samples(:, 5)];
    else
        error('%s must contain diagnostic_samples or theta_samples.', chain_file);
    end
    chains_cell{c} = diagnostic_samples;
end

N = min(cellfun(@(x) size(x, 1), chains_cell));
P = numel(quantity_names);
chains = zeros(n_chains, N, P);
for c = 1:n_chains
    chains(c, :, :) = chains_cell{c}(1:N, 1:P);
end

Rhat = compute_rhat(chains);

fprintf('Reduced proteoliposome Rhat values:\n');
for p = 1:P
    fprintf('%-20s Rhat = %.4f\n', quantity_names{p}, Rhat(p));
end

figure('Color', 'w');
bar(Rhat);
hold on;
plot([0, P+1], [1.1, 1.1], 'r--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:P, 'XTickLabel', quantity_names, 'XTickLabelRotation', 45);
ylabel('Rhat');
title('Reduced proteoliposome Rhat diagnostics');
grid on;

save('Rhat_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat', ...
    'Rhat', 'quantity_names', 'N', 'n_chains');

%% =============================== Local functions ==========================
function Rhat = compute_rhat(chains)
    [M, N, P] = size(chains);
    Rhat = zeros(1, P);
    for p = 1:P
        theta_p = squeeze(chains(:, :, p));
        chain_means = mean(theta_p, 2);
        B = N * var(chain_means, 0, 1);
        W = mean(var(theta_p, 0, 2));
        V = ((N - 1) / N) * W + (1 / N) * B;
        Rhat(p) = sqrt(V / W);
    end
end

function rho_pm = derived_rho_pm(theta_samples)
    rho_pm = (theta_samples(:, 2) .* theta_samples(:, 3)) ./ ...
             (theta_samples(:, 1) .* theta_samples(:, 4));
end
