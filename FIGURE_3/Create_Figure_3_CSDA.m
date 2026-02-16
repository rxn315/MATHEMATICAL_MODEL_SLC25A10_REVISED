%% Run LSA scripts (COMPLEX-STEP versions)
% These scripts should save mean_LSCs_* (or mean_LSCs) computed via complex-step.
run('LSA_Jsucc_REA_complexstep.m');   % e.g., saves LSCs_succ_complexstep.mat
run('LSA_Jmal_REA_complexstep.m');    % e.g., saves LSCs_mal_complexstep.mat
run('LSA_Jpho_REA_complexstep.m');    % e.g., saves LSCs_pho_complexstep.mat

%% Load LSC results (complex-step)
data_succ = load('LSCs_succ_complexstep.mat');
data_mal  = load('LSCs_mal_complexstep.mat');
data_pho  = load('LSCs_pho_complexstep.mat');

% ---------- helper to extract mean LSCs robustly ----------
get_mean = @(S,preferred,alt) ...
    (isfield(S,preferred) * S.(preferred)) + (~isfield(S,preferred) * S.(alt));

% Succinate
if isfield(data_succ,'mean_LSCs_succ')
    mean_LSCs_succ = data_succ.mean_LSCs_succ;
elseif isfield(data_succ,'mean_LSCs')
    mean_LSCs_succ = data_succ.mean_LSCs;
else
    error('LSCs_succ_complexstep.mat does not contain mean_LSCs_succ or mean_LSCs');
end

% Malate
if isfield(data_mal,'mean_LSCs_mal')
    mean_LSCs_mal = data_mal.mean_LSCs_mal;
elseif isfield(data_mal,'mean_LSCs')
    mean_LSCs_mal = data_mal.mean_LSCs;
else
    error('LSCs_mal_complexstep.mat does not contain mean_LSCs_mal or mean_LSCs');
end

% Phosphate
if isfield(data_pho,'mean_LSCs_pho')
    mean_LSCs_pho = data_pho.mean_LSCs_pho;
elseif isfield(data_pho,'mean_LSCs')
    mean_LSCs_pho = data_pho.mean_LSCs;
else
    error('LSCs_pho_complexstep.mat does not contain mean_LSCs_pho or mean_LSCs');
end

%% Create combined LSC figure (J_pho, J_mal, J_succ)
% (reuse your existing createfigure; it should just plot the vectors)
createfigure(mean_LSCs_pho, mean_LSCs_mal, mean_LSCs_succ);

%% (Optional) Save figure with a name that reflects complex-step
set(gcf,'Color','w');
exportgraphics(gcf,'LSC_combined_complexstep.png','Resolution',300);
