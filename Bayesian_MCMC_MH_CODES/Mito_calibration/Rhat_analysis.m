%% ========================================================================
% Rhat analysis for reduced intact-mitochondria MCMC chains.
%
% Expected reduced chain files:
%   chain1_result_palmier_reduced_uniform.mat
%   chain2_result_palmier_reduced_uniform.mat
%   chain3_result_palmier_reduced_uniform.mat
%   chain4_result_palmier_reduced_uniform.mat
%
% theta_samples columns:
%   [Ts_max, Tm_max, Tp_max, K_D_s, K_D_m, K_D_p, sigma2]
%
% lambda21 and lambda31 are derived for diagnostics.
%% ========================================================================

clear; clc; close all;

chain_pattern = 'chain%d_result_palmier_reduced_uniform.mat';
n_chains = 4;

quantity_names = {'Ts_max', 'Tm_max', 'Tp_max', ...
    'K_D_s', 'K_D_m', 'K_D_p', ...
    'lambda21_derived', 'lambda31_derived', 'sigma2'};

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
        lambda21 = derived_lambda21(theta_samples(:, 1:6));
        lambda31 = derived_lambda31(theta_samples(:, 1:6));
        diagnostic_samples = [theta_samples(:, 1:6), lambda21, lambda31, theta_samples(:, 7)];
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

fprintf('Reduced intact-mitochondria Rhat values:\n');
for p = 1:P
    fprintf('%-20s Rhat = %.4f\n', quantity_names{p}, Rhat(p));
end

figure('Color', 'w');
bar(Rhat);
hold on;
plot([0, P+1], [1.1, 1.1], 'r--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:P, 'XTickLabel', quantity_names, 'XTickLabelRotation', 45);
ylabel('Rhat');
title('Reduced intact-mitochondria Rhat diagnostics');
grid on;

save('Rhat_reduced_uniform.mat', 'Rhat', 'quantity_names', 'N', 'n_chains');

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

function lambda21 = derived_lambda21(params6)
    lambda21 = (params6(:, 2) .* params6(:, 4)) ./ ...
               (params6(:, 1) .* params6(:, 5));
end

function lambda31 = derived_lambda31(params6)
    lambda31 = (params6(:, 3) .* params6(:, 4)) ./ ...
               (params6(:, 1) .* params6(:, 6));
end
