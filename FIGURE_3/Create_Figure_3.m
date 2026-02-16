%% Run LSA scripts for each flux (REA model)
% Make sure these scripts use the REA-based compute_flux_REA
run('LSA_Jsucc_REA.m');   % saves LSCs_succ.mat
run('LSA_Jmal_REA.m');    % saves LSCs_mal.mat
run('LSA_Jpho_REA.m');    % saves LSCs_pho.mat

%% Load LSC results
data_succ = load('LSCs_succ.mat');  % succinate
data_mal  = load('LSCs_mal.mat');   % malate
data_pho  = load('LSCs_pho.mat');   % phosphate

% Handle possible variable name differences inside .mat files
if isfield(data_succ,'mean_LSCs_succ')
    mean_LSCs_succ = data_succ.mean_LSCs_succ;
elseif isfield(data_succ,'mean_LSCs')
    mean_LSCs_succ = data_succ.mean_LSCs;
else
    error('LSCs_succ.mat does not contain mean_LSCs_succ or mean_LSCs');
end

if isfield(data_mal,'mean_LSCs_mal')
    mean_LSCs_mal = data_mal.mean_LSCs_mal;
elseif isfield(data_mal,'mean_LSCs')
    mean_LSCs_mal = data_mal.mean_LSCs;
else
    error('LSCs_mal.mat does not contain mean_LSCs_mal or mean_LSCs');
end

if isfield(data_pho,'mean_LSCs_pho')
    mean_LSCs_pho = data_pho.mean_LSCs_pho;
elseif isfield(data_pho,'mean_LSCs')
    mean_LSCs_pho = data_pho.mean_LSCs;
else
    error('LSCs_pho.mat does not contain mean_LSCs_pho or mean_LSCs');
end

%% Create combined LSC figure (J_pho, J_mal, J_succ)
createfigure(mean_LSCs_pho, mean_LSCs_mal, mean_LSCs_succ);
