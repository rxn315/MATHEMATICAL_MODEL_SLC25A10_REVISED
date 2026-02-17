% MATLAB function to calculate optimal batch size for MCMC chains (multivariate case)
% Input: chain_matrix - a matrix of size N x P (9000 samples x 8 parameters)
% Output: b_opt - optimal batch size (scalar, same for all parameters)
%         mcse_matrix - P x 1 vector of MCSE for each parameter using batch means
%         ess_matrix - P x 1 vector of approximate ESS for each parameter
% Heuristic: b_opt ≈ floor(N^(1/3)) for reasonable mixing chains
% Reference: Liu, Y., Vats, D., & Flegal, J. M. (2021). Batch size selection for variance estimators in MCMC.

function [b_opt, mcse_matrix, ess_matrix] = optimalBatchSizeMultivariate(chain_matrix)
    % Validate input
    [N, P] = size(chain_matrix);
    if N ~= 9000
        warning('Number of samples is %d, not 9000. Adjusting accordingly.', N);
    end
    if P ~= 8
        warning('Number of parameters is %d, not 8. Proceeding.', P);
    end
    if ~ismatrix(chain_matrix)
        error('Input "chain_matrix" must be a matrix of size N x P (samples x parameters).');
    end
    
    % Heuristic optimal batch size: floor(N^(1/3))
    % Ensures at least 30 batches for reliable variance estimation
    b_opt = floor(N^(1/3));
    min_batches = 30;
    if N / b_opt < min_batches
        b_opt = floor(N / min_batches);
    end
    b_opt = max(1, b_opt);  % Minimum batch size 1
    
    % Number of batches
    m = floor(N / b_opt);
    
    % Initialize outputs
    mcse_matrix = zeros(P, 1);
    ess_matrix = zeros(P, 1);
    
    % For each parameter
    for p = 1:P
        theta = chain_matrix(:, p);  % Extract column for parameter p (N x 1)
        
        % Truncate to m * b_opt samples to avoid partial batches
        theta_trunc = theta(1 : m * b_opt);
        
        % Reshape into m batches of size b_opt
        batches = reshape(theta_trunc, b_opt, m);  % b_opt rows x m columns
        batch_means = mean(batches, 1);  % Mean of each batch (1 x m)
        
        % Sample variance of batch means (1/(m-1) denominator)
        var_batch_means = var(batch_means, 0, 2);  % 1 x 1 (along dimension 2)
        
        % MCSE for the overall mean using batch means
        mcse_matrix(p) = sqrt(var_batch_means * (b_opt / N));
        
        % Approximate Effective Sample Size (ESS) from batch means variance
        % ESS ≈ N * (var_sample / var_batch_means) / b_opt, but more precisely:
        % Asymptotic variance σ² ≈ b_opt * var_batch_means
        % Then ESS = N * (sample_var / asymptotic_var)
        sample_var = var(theta, 0);  % Population variance for large N
        asymptotic_var = b_opt * var_batch_means;
        if asymptotic_var > 0
            ess_matrix(p) = N * (sample_var / asymptotic_var);
        else
            ess_matrix(p) = N;  % Fallback if no variance
        end
    end
    
    % Display results
    fprintf('Chain dimensions: %d samples x %d parameters\n', N, P);
    fprintf('Optimal batch size b_opt = %d (heuristic: floor(N^(1/3)) ≈ %.1f)\n', b_opt, N^(1/3));
    fprintf('Number of batches m = %d\n', m);
    fprintf('\nMCSE and ESS for each parameter:\n');
    fprintf('Param\tMCSE\t\tESS\n');
    for p = 1:P
        fprintf('%d\t%.6f\t%.1f\n', p, mcse_matrix(p), ess_matrix(p));
    end
end

% Example usage (with dummy data for N=9000, P=8)
% chain_matrix = randn(9000, 8);  % Replace with your actual MCMC chain matrix
% [b_opt, mcse_matrix, ess_matrix] = optimalBatchSizeMultivariate(chain_matrix);