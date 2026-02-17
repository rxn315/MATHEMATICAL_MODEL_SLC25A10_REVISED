% MATLAB code to compute Rhat statistic for 4 chains, each with 9000 samples and 8 parameters
% Input: chains are stored as a 3D array 'chains' of size M x N x P
%   M = number of chains (4)
%   N = samples per chain (9000)
%   P = number of parameters (8)
% Output: Rhat for each parameter

clear; clc;

% Parameters
M = 4;        % Number of chains
N = 990000;     % Samples per chain
P = 8;        % Number of parameters

% Load your chain data here
% Expected format: chains is a 3D array of size M x N x P
% Example: Generate dummy data (replace with your actual data)
chains = randn(M, N, P);  % Dummy data: M chains, N samples, P parameters
% Replace with your actual data, e.g.:
load('chain1_result_Palmier.mat');
chain1=theta_samples;
load('chain2_result_palmier.mat');
chain2=theta_samples;
load('chain3_result_palmier.mat');
chain3=theta_samples;
load('chain4_result_palmier.mat');
chain4=theta_samples;


chains(1, :, :) = chain1; % where your_chain1 is 9000 x 8
chains(2, :, :) = chain2;
chains(3, :, :) = chain3;
chains(4, :, :) = chain4;

% Initialize Rhat array for each parameter
Rhat = zeros(1, P);

for p = 1:P
    % Extract data for parameter p across all chains
    theta_p = squeeze(chains(:, :, p)); % M x N matrix for parameter p

    % Compute chain means (mean over samples for each chain)
    chain_means = mean(theta_p, 2); % M x 1 vector

    % Grand mean across all chains for parameter p
    grand_mean = mean(chain_means);

    % Between-chain variance B (sample variance of chain means, scaled by N)
    B = N * var(chain_means); % var uses 1/(M-1) denominator

    % Within-chain variance W
    within_vars = var(theta_p, 0, 2); % M x 1, variance for each chain
    W = mean(within_vars); % Average over chains

    % Pooled variance estimate V
    V = ((N - 1) / N) * W + (1 / N) * B;

    % Rhat statistic for parameter p
    Rhat(p) = sqrt(V / W);
end

% Display results
fprintf('Rhat values for each parameter:\n');
for p = 1:P
    fprintf('Parameter %d: Rhat = %.4f\n', p, Rhat(p));
    if Rhat(p) < 1.1
        fprintf('  - Convergence looks good (Rhat < 1.1).\n');
    else
        fprintf('  - Potential convergence issues (Rhat >= 1.1).\n');
    end
end

% Optional: Plot Rhat values
figure;
bar(1:P, Rhat);
hold on;
plot([0, P+1], [1.1, 1.1], 'r--', 'LineWidth', 1.5); % Threshold line
xlabel('Parameter');
ylabel('Rhat');
title('Rhat Statistic for Convergence Diagnosis');
grid on;