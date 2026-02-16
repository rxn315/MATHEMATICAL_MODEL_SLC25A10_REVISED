% ============================================================
% Compare finite-difference vs complex-step LSA results
% Uses parameter NAMES (not indices) in tables + plots
% ============================================================
clear; clc; close all;

% -----------------------
% FILES (edit if needed)
% -----------------------
pairs = { ...
    'succinate', 'LSCs_succ.mat', 'LSCs_succ_complexstep.mat'; ...
    'malate',    'LSCs_mal.mat',  'LSCs_mal_complexstep.mat'; ...
    'phosphate', 'LSCs_pho.mat',  'LSCs_pho_complexstep.mat'};

% If your scripts used different parameter order/names, set them here:
% (leave empty to auto-generate p1, p2, ...)
param_names = { ...
    'T^{s}_{max}','T^{m}_{max}','K^{s}_{m}','K^{m}_{m}','K^{p}_{m}','\lambda_{21}','\lambda_{31}'};

epsRel = 1e-12;  % for safe relative error

% -----------------------
% Helper to extract mean LSC vector robustly
% -----------------------
extract_mean = @(S, fluxname) local_extract_mean(S, fluxname);

% -----------------------
% Helper to build parameter labels (names preferred)
% -----------------------
make_param_labels = @(n) local_param_labels(n, param_names);

% -----------------------
% Loop over fluxes
% -----------------------
all_rows = [];
all_metrics = struct();

for r = 1:size(pairs,1)
    fluxname = pairs{r,1};
    fd_file  = pairs{r,2};
    cs_file  = pairs{r,3};

    if ~isfile(fd_file), error('Missing FD file: %s', fd_file); end
    if ~isfile(cs_file), error('Missing CS file: %s', cs_file); end

    FD = load(fd_file);
    CS = load(cs_file);

    fd = extract_mean(FD, fluxname);
    cs = extract_mean(CS, fluxname);

    fd = fd(:); cs = cs(:);

    if numel(fd) ~= numel(cs)
        error('Length mismatch for %s: FD=%d, CS=%d', fluxname, numel(fd), numel(cs));
    end

    % Parameter labels (names if provided, else p1..pn)
    pLabels = make_param_labels(numel(cs));           % string array (n x 1)
    xCats   = categorical(pLabels, pLabels, 'Ordinal', true);  % preserve order

    % ---------- Errors per parameter ----------
    err_abs = abs(cs - fd);
    err_rel = err_abs ./ max(abs(cs), epsRel);   % relative to CS

    % ---------- Summary metrics ----------
    metrics = struct();
    metrics.n         = numel(cs);
    metrics.MAE_abs   = mean(err_abs,'omitnan');
    metrics.RMSE_abs  = sqrt(mean((cs - fd).^2,'omitnan'));
    metrics.Max_abs   = max(err_abs,[],'omitnan');

    metrics.MAE_rel   = mean(err_rel,'omitnan');
    metrics.Max_rel   = max(err_rel,[],'omitnan');

    metrics.L2_rel    = norm(cs - fd, 2) / max(norm(cs,2), epsRel);
    metrics.Linf_rel  = norm(cs - fd, inf) / max(norm(cs,inf), epsRel);

    metrics.corr      = corr(fd, cs, 'Rows','complete');

    % sign mismatch count (excluding tiny values)
    tiny = 1e-14;
    s_fd = sign(fd); s_cs = sign(cs);
    mask = (abs(fd) > tiny) & (abs(cs) > tiny);
    metrics.sign_mismatch = sum(s_fd(mask) ~= s_cs(mask));
    metrics.sign_total    = sum(mask);

    all_metrics.(fluxname) = metrics;

    % ---------- Build a per-parameter table-like matrix ----------
    row = table();
    row.Flux      = repmat(string(fluxname), numel(cs), 1);
    row.ParamName = pLabels;                 % <-- name replaces index
    row.FD        = fd;
    row.CS        = cs;
    row.AbsError  = err_abs;
    row.RelError  = err_rel;

    all_rows = [all_rows; row]; %#ok<AGROW>

    % ---------- Plot per-flux errors ----------
    figure('Color','w','Name',sprintf('FD vs CS errors: %s', fluxname));
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    bar(xCats, err_abs);
    grid on;
    xlabel('Parameter');
    ylabel('|CS - FD|');
    title(sprintf('%s: Absolute error', fluxname));
    xtickangle(45);

    nexttile;
    bar(xCats, err_rel);
    grid on;
    xlabel('Parameter');
    ylabel('|CS - FD| / |CS|');
    title(sprintf('%s: Relative error (vs CS)', fluxname));
    xtickangle(45);
end

% -----------------------
% Print summary metrics
% -----------------------
fprintf('\n================ FD vs CS Summary ================\n');
fluxes = fieldnames(all_metrics);
for i = 1:numel(fluxes)
    f = fluxes{i};
    m = all_metrics.(f);
    fprintf('\n[%s]\n', upper(f));
    fprintf('  MAE_abs   = %.3e\n', m.MAE_abs);
    fprintf('  RMSE_abs  = %.3e\n', m.RMSE_abs);
    fprintf('  Max_abs   = %.3e\n', m.Max_abs);
    fprintf('  MAE_rel   = %.3e\n', m.MAE_rel);
    fprintf('  Max_rel   = %.3e\n', m.Max_rel);
    fprintf('  L2_rel    = %.3e\n', m.L2_rel);
    fprintf('  Linf_rel  = %.3e\n', m.Linf_rel);
    fprintf('  corr(FD,CS)= %.6f\n', m.corr);
    fprintf('  sign mismatch = %d / %d\n', m.sign_mismatch, m.sign_total);
end
fprintf('\n==================================================\n');

% -----------------------
% Save outputs
% -----------------------
save('FD_vs_CS_error_metrics.mat','all_metrics','all_rows');
writetable(all_rows,'FD_vs_CS_error_per_parameter.csv');

disp('Saved: FD_vs_CS_error_metrics.mat');
disp('Saved: FD_vs_CS_error_per_parameter.csv');

% ============================================================
% Local helper: extract mean LSC vector from a loaded struct
% ============================================================
function v = local_extract_mean(S, fluxname)
    if isfield(S,'mean_LSCs')
        v = S.mean_LSCs; return;
    end

    switch lower(fluxname)
        case 'succinate'
            candidates = {'mean_LSCs_succ','mean_LSCs_succinate'};
        case 'malate'
            candidates = {'mean_LSCs_mal','mean_LSCs_malate'};
        case 'phosphate'
            candidates = {'mean_LSCs_pho','mean_LSCs_phosphate'};
        otherwise
            candidates = {};
    end

    for k = 1:numel(candidates)
        if isfield(S, candidates{k})
            v = S.(candidates{k});
            return;
        end
    end

    fn = fieldnames(S);
    for k = 1:numel(fn)
        x = S.(fn{k});
        if isnumeric(x) && isvector(x) && numel(x) > 1
            v = x;
            warning('Using fallback field "%s" as mean LSC vector.', fn{k});
            return;
        end
    end

    error('Could not find mean LSC vector in the provided .mat file.');
end

% ============================================================
% Local helper: build parameter labels (names if provided)
% ============================================================
function pLabels = local_param_labels(n, param_names)
    if ~isempty(param_names) && numel(param_names) == n
        pLabels = string(param_names(:));
    else
        % fallback to p1..pn
        pLabels = compose("p%d", (1:n)');
    end
end
