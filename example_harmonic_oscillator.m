%EXAMPLE_HARMONIC_OSCILLATOR 
% Thermal relaxation of a harmonic oscillator.
% The ensemble starts from U = V = 0, not from thermal equilibrium.


%% Input parameters in consistent reduced units
cfg = struct();
cfg.M = 1.0;
cfg.K = 1.0;
cfg.Z = 1.0;
cfg.kBT = 1.0;
cfg.dt = 0.1;
cfg.t_end = 20.0;
cfg.n_trajectories = 20000;
cfg.n_paths_to_store = 4;
cfg.seed = 2026;
cfg.png_resolution = 300;


%% Create an independent output directory and preserve the source files
scriptPath = [mfilename('fullpath'), '.m'];
[sourceDirectory, scriptName] = fileparts(scriptPath);
addpath(sourceDirectory);
outputRoot = fullfile(sourceDirectory, scriptName);
if exist(outputRoot,'dir') ~= 7
    [ok, message] = mkdir(outputRoot);
    if ~ok, error('PGLGV:OutputDirectory','%s',message); end
end
while true
    runDirectory = fullfile(outputRoot, datestr(now,'yyyymmdd_HHMMSS'));
    if exist(runDirectory,'dir') ~= 7, break; end
    pause(0.1);
end
[ok, message] = mkdir(runDirectory);
if ~ok, error('PGLGV:OutputDirectory','%s',message); end
sourceNames = {[scriptName,'.m'], 'integrator_PGLGV_setup.m', 'integrator_PGLGV.m'};
for sourceIndex = 1:numel(sourceNames)
    [ok, message] = copyfile(fullfile(sourceDirectory,sourceNames{sourceIndex}), ...
        fullfile(runDirectory,sourceNames{sourceIndex}));
    if ~ok, error('PGLGV:SourceSnapshot','%s',message); end
end
diary(fullfile(runDirectory,'run.log'));
rng_before = rng;
rng(cfg.seed,'twister');
rng_initial = rng;
completed_steps = 0;
run_status = 'initializing';
matlab_version = version;
save(fullfile(runDirectory,'inputs.mat'), 'cfg', 'rng_before', 'rng_initial', 'matlab_version');

try
    validateattributes(cfg.n_trajectories, {'double'}, ...
        {'scalar','integer','finite','>=',2});
    validateattributes(cfg.n_paths_to_store, {'double'}, ...
        {'scalar','integer','positive','<=',cfg.n_trajectories});
    validateattributes(cfg.t_end, {'double'}, {'scalar','finite','positive'});
    nSteps = round(cfg.t_end/cfg.dt);
    if nSteps < 1 || abs(nSteps*cfg.dt-cfg.t_end) > 100*eps*max(1,cfg.t_end)
        error('PGLGV:TimeGrid','t_end must be an integer multiple of dt.');
    end
    time = (0:nSteps)'*cfg.dt;
    op = integrator_PGLGV_setup(cfg.M,cfg.K,cfg.Z,cfg.kBT,cfg.dt);
    fprintf('Harmonic oscillator with %d independent trajectories\n',cfg.n_trajectories);
    fprintf('  M = %.6g, K = %.6g, Z = %.6g, kBT = %.6g\n', ...
        cfg.M,cfg.K,cfg.Z,cfg.kBT);
    fprintf('  dt = %.6g, t_end = %.6g, steps = %d, RNG seed = %d\n', ...
        cfg.dt,cfg.t_end,nSteps,cfg.seed);
    fprintf('  Scalar admissible step limit = %.9g\n',op.dt_limit);
    fprintf('  Scaled minimum covariance eigenvalue = %.6e\n', ...
        op.diagnostics.scaled_minimum_eigenvalue);
    fprintf('  Normalized covariance-balance residual = %.6e\n', ...
        op.diagnostics.covariance_balance_residual);
    fprintf('  Output directory = %s\n',runDirectory);

    Ac = [0,1; -cfg.K/cfg.M,-cfg.Z/cfg.M];
    E = expm(cfg.dt*Ac);
    Qc = op.Sigma_eq-E*op.Sigma_eq*E';
    Qc = (Qc+Qc')/2;
    Sigma_discrete = zeros(2,2);
    Sigma_continuous = zeros(2,2);
    U = zeros(1,cfg.n_trajectories);
    V = zeros(1,cfg.n_trajectories);
    paths_U = nan(nSteps+1,cfg.n_paths_to_store);
    paths_V = nan(nSteps+1,cfg.n_paths_to_store);
    paths_U(1,:) = 0;
    paths_V(1,:) = 0;
    mean_U = zeros(nSteps+1,1);
    mean_V = zeros(nSteps+1,1);
    moments_ensemble = zeros(nSteps+1,3);
    moments_discrete = zeros(nSteps+1,3);
    moments_continuous = zeros(nSteps+1,3);
    covariance_sampling_error = zeros(nSteps+1,1);
    covariance_time_error = zeros(nSteps+1,1);
    progressStride = max(1,ceil(nSteps/5));
    simulationClock = tic;
    run_status = 'running';

    for stepIndex = 1:nSteps
        [U,V] = integrator_PGLGV(op,U,V);
        Sigma_discrete = op.A*Sigma_discrete*op.A'+op.Q;
        Sigma_continuous = E*Sigma_continuous*E'+Qc;
        X = [U;V];
        meanX = mean(X,2);
        X_centered = bsxfun(@minus,X,meanX);
        Sigma_sample = (X_centered*X_centered')/(cfg.n_trajectories-1);
        row = stepIndex+1;
        mean_U(row) = meanX(1);
        mean_V(row) = meanX(2);
        moments_ensemble(row,:) = [Sigma_sample(1,1),Sigma_sample(2,2),Sigma_sample(1,2)];
        moments_discrete(row,:) = [Sigma_discrete(1,1),Sigma_discrete(2,2),Sigma_discrete(1,2)];
        moments_continuous(row,:) = [Sigma_continuous(1,1),Sigma_continuous(2,2),Sigma_continuous(1,2)];
        covariance_sampling_error(row) = norm(op.S_eq \ ...
            ((Sigma_sample-Sigma_discrete)/op.S_eq'),'fro')/sqrt(2);
        covariance_time_error(row) = norm(op.S_eq \ ...
            ((Sigma_discrete-Sigma_continuous)/op.S_eq'),'fro')/sqrt(2);
        paths_U(row,:) = U(1:cfg.n_paths_to_store);
        paths_V(row,:) = V(1:cfg.n_paths_to_store);
        completed_steps = stepIndex;
        if mod(stepIndex,progressStride) == 0 || stepIndex == nSteps
            fprintf('  Completed %d / %d steps in %.2f s\n', ...
                stepIndex,nSteps,toc(simulationClock));
        end
    end
    elapsed_seconds = toc(simulationClock);
    rng_final = rng;
    run_status = 'simulation_complete';
    fprintf('Final ensemble variance ratios (canonical value 1)\n');
    fprintf('  Displacement = %.6f\n',moments_ensemble(end,1)/op.Sigma_eq(1,1));
    fprintf('  Velocity     = %.6f\n',moments_ensemble(end,2)/op.Sigma_eq(2,2));
    fprintf('  Normalized displacement-velocity covariance = %.6e\n', ...
        moments_ensemble(end,3)/sqrt(op.Sigma_eq(1,1)*op.Sigma_eq(2,2)));
    fprintf('Maximum normalized covariance time-discretization error = %.6e\n',max(covariance_time_error));
    fprintf('Finite ensemble variances retain sampling error.\n');

    save(fullfile(runDirectory,'results.mat'), 'cfg','op','time','paths_U','paths_V', ...
        'U','V','mean_U','mean_V','moments_ensemble','moments_discrete', ...
        'moments_continuous','covariance_sampling_error','covariance_time_error', ...
        'rng_initial','rng_final','rng_before','completed_steps','elapsed_seconds', ...
        'run_status','matlab_version','-v7');
    momentTable = table(time,mean_U,mean_V, ...
        moments_ensemble(:,1),moments_ensemble(:,2),moments_ensemble(:,3), ...
        moments_discrete(:,1),moments_discrete(:,2),moments_discrete(:,3), ...
        moments_continuous(:,1),moments_continuous(:,2),moments_continuous(:,3), ...
        covariance_sampling_error,covariance_time_error, ...
        'VariableNames',{'time','mean_U','mean_V', ...
        'var_U_ensemble','var_V_ensemble','cov_UV_ensemble', ...
        'var_U_discrete','var_V_discrete','cov_UV_discrete', ...
        'var_U_continuous','var_V_continuous','cov_UV_continuous', ...
        'covariance_sampling_error','covariance_time_error'});
    writetable(momentTable,fullfile(runDirectory,'harmonic_moments.csv'));

    fig = figure('Name','Harmonic displacement trajectories');
    plot(time,paths_U,'LineWidth',1.0);
    xlabel('Time'); ylabel('Displacement u'); grid on; box on;
    local_save_figure(fig,runDirectory,'fig_harmonic_displacement',cfg.png_resolution);

    fig = figure('Name','Harmonic velocity trajectories');
    plot(time,paths_V,'LineWidth',1.0);
    xlabel('Time'); ylabel('Stored endpoint velocity v'); grid on; box on;
    local_save_figure(fig,runDirectory,'fig_harmonic_velocity',cfg.png_resolution);

    names = {'displacement','velocity'};
    labels = {'Displacement variance / (k_B*T/K)', ...
        'Velocity variance / (k_B*T/M)'};
    for component = 1:2
        fig = figure('Name',['Harmonic ',names{component},' variance']);
        targetVariance = op.Sigma_eq(component,component);
        plot(time,moments_ensemble(:,component)/targetVariance,'-','LineWidth',2.0);
        hold on;
        plot(time,moments_discrete(:,component)/targetVariance,'--','LineWidth',1.5);
        plot(time,moments_continuous(:,component)/targetVariance,':','LineWidth',1.5);
        plot([time(1),time(end)],[1,1],'-.','LineWidth',1.5);
        xlabel('Time'); ylabel(labels{component}); grid on; box on;
        legend('Ensemble','Discrete covariance recursion','Continuous reference', ...
            'Canonical value','Location','southeast');
        local_save_figure(fig,runDirectory, ...
            ['fig_harmonic_',names{component},'_variance'],cfg.png_resolution);
    end
    run_status = 'completed';
    save(fullfile(runDirectory,'results.mat'),'run_status','-append');
    fprintf('Completed. MAT, CSV, source snapshots and FIG/EPS/PNG figures are saved.\n');
    diary off;
catch ME
    failure = struct('identifier',ME.identifier,'message',ME.message, ...
        'stack',ME.stack,'completed_steps',completed_steps,'run_status',run_status);
    rng_at_failure = rng;
    save(fullfile(runDirectory,'failure.mat'),'failure','cfg','rng_initial','rng_at_failure');
    fprintf(2,'Calculation stopped after %d steps. %s\n',completed_steps,ME.message);
    diary off;
    rethrow(ME);
end

%% helper
function local_save_figure(fig,folder,baseName,resolution)
    set(fig,'PaperPositionMode','auto');
    savefig(fig,fullfile(folder,[baseName,'.fig']));
    print(fig,fullfile(folder,[baseName,'.eps']),'-depsc','-painters');
    print(fig,fullfile(folder,[baseName,'.png']),'-dpng',sprintf('-r%d',resolution));
end
