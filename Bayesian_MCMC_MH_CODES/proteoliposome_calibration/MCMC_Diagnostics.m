%% ========================================================================
% Single-chain diagnostics for reduced proteoliposome MCMC output.
%
% This script uses theta_samples from:
%   MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
% or from:
%   chain1_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat
%
% Diagnostics include trace plots, toolbox-free ACF estimates, ESS, and MCSE.
%% ========================================================================

clear; clc; close all;

S = load_first_existing({ ...
    'MCMC_Result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat', ...
    'chain1_result_liposome_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat'});

if ~isfield(S, 'theta_samples')
    error('The selected result file must contain theta_samples.');
end

theta_samples = S.theta_samples;
rho_pm = derived_rho_pm(theta_samples);
samples = [theta_samples(:, 1:4), rho_pm, theta_samples(:, 5)];

labels = {'Tmax_m', 'Tmax_p', 'K_D_m', 'K_D_p', 'rho_pm_derived', 'sigma2'};

[N, P] = size(samples);
max_lag = min(100, floor(N / 5));

figure('Color', 'w', 'Name', 'Reduced proteoliposome trace plots');
tiledlayout(P, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for j = 1:P
    nexttile;
    plot(1:N, samples(:, j), 'k-');
    ylabel(labels{j}, 'Interpreter', 'none');
    grid on;
end
xlabel('Post-burn-in iteration');

acf = zeros(max_lag + 1, P);
lags = (0:max_lag).';
figure('Color', 'w', 'Name', 'Reduced proteoliposome ACF');
tiledlayout(P, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for j = 1:P
    acf(:, j) = local_acf(samples(:, j), max_lag);
    nexttile;
    stem(lags, acf(:, j), 'filled');
    xlim([0 max_lag]);
    ylabel(labels{j}, 'Interpreter', 'none');
    grid on;
end
xlabel('Lag');

ESS = zeros(P, 1);
for j = 1:P
    rho = acf(2:end, j);
    gsum = 0;
    for k = 1:2:(numel(rho)-1)
        pair = rho(k) + rho(k+1);
        if pair <= 0
            break;
        end
        gsum = gsum + pair;
    end
    ESS(j) = N / (1 + 2*gsum);
end

B = 50;
m = floor(N / B);
idx = 1:(m * B);
Sbm = reshape(samples(idx, :), m, B, P);
bm = squeeze(mean(Sbm, 1));
mcse = sqrt(var(bm, 0, 1) / B).';
postSD = std(samples, 0, 1).';
ratio = mcse ./ postSD;

diagnostics = table(labels(:), ESS, mcse, postSD, ratio, ...
    'VariableNames', {'Quantity', 'ESS', 'MCSE', 'PosteriorSD', 'MCSE_over_SD'});
disp(diagnostics);

save('MCMC_Diagnostics_reduced_Jpho_gaussian_IGsigma_no_exMal1mM.mat', ...
    'diagnostics', 'labels', 'ESS', 'mcse', 'postSD', 'ratio', 'acf', 'lags');

%% =============================== Local functions ==========================
function S = load_first_existing(files)
    for k = 1:numel(files)
        if exist(files{k}, 'file') == 2
            fprintf('Loading %s\n', files{k});
            S = load(files{k});
            return;
        end
    end
    error('None of the expected result files were found.');
end

function acf = local_acf(x, max_lag)
    x = x(:);
    x = x - mean(x, 'omitnan');
    denom = sum(x.^2);
    acf = zeros(max_lag + 1, 1);
    if denom <= 0
        acf(1) = 1;
        return;
    end
    for lag = 0:max_lag
        acf(lag + 1) = sum(x(1:end-lag) .* x(1+lag:end)) / denom;
    end
end

function rho_pm = derived_rho_pm(theta_samples)
    rho_pm = (theta_samples(:, 2) .* theta_samples(:, 3)) ./ ...
             (theta_samples(:, 1) .* theta_samples(:, 4));
end
