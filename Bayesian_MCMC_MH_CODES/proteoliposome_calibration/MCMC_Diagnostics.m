load('chain1_result.mat')
samples=theta_samples;
chain=theta_samples;
[b_opt, mcse_matrix, ess_matrix] = optimalBatchSizeMultivariate(chain);

% INPUT: samples (N x P), N=9000
[N,P] = size(samples);

% ===== TRACE PLOTS =====
figure; 
for j = 1:P
    subplot(P,1,j); plot(1:N, samples(:,j), '-'); 
    xlabel('iteration'); ylabel(sprintf('\\theta_%d',j)); 
    title(sprintf('Trace: \\theta_%d', j));
end

% ===== AUTOCORRELATION (ACF) =====
% max lag ~ min(1000, N/5) is typical
maxLag = min(100, floor(N/5));
acf = cell(P,1);
figure;
for j = 1:P
    [acf{j}, lags] = autocorr(samples(:,j), 'NumLags', maxLag); % Stats TB
    subplot(P,1,j); stem(lags, acf{j}, 'filled'); 
    xlim([0 maxLag]); xlabel('lag'); ylabel('ACF'); 
    title(sprintf('ACF: \\theta_%d', j));
end

% ===== EFFECTIVE SAMPLE SIZE (ESS) =====
% ESS = N / (1 + 2 * sum_{k>=1} rho_k), with sum truncated at first even-odd
% positive pair or until acf becomes negative for many lags.
ESS = zeros(P,1);
for j = 1:P
    rho = acf{j}(2:end);  % drop lag 0
    % Geyer initial positive sequence trick (simple version)
    gsum = 0; 
    for k = 1:2:length(rho)-1
        pair = rho(k) + (k+1<=length(rho))*rho(k+1);
        if pair <= 0, break; end
        gsum = gsum + pair;
    end
    ESS(j) = N / (1 + 2*gsum);
end

% ===== MCSE (Monte Carlo Standard Error) via batch means =====
% Choose ~50–100 batches; batch size >= 30 is common.
B = 50;                          % number of batches
m = floor(N / B);                % batch size
idx = (1:(m*B))';
Sbm = reshape(samples(idx,:), m, B, P);  % (m x B x P)
bm  = squeeze(mean(Sbm,1));              % (B x P) batch means
mcse = sqrt(var(bm, 0, 1) / B);          % (1 x P) MCSE of posterior mean
mcse = mcse(:);                          % (P x 1)

% Posterior SD (for MCSE/SD ratio):
postSD = std(samples, 0, 1)';            % (P x 1)
ratio  = mcse ./ postSD;                 % want < 0.10 ideally

% ===== RANK PLOT (needs multiple chains to be meaningful) =====
% With a single chain, skip rank plots. See multi-chain code below.
