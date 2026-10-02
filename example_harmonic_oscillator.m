% EXAMPLE_HARMONIC_OSCILLATOR
% Thermal relaxation of a harmonic oscillator from U = V = 0.

clc;
close all;

M = 1.0;
K = 1.0;
Z = 1.0;
kBT = 1.0;
dt = 0.1;
tEnd = 10.0;
nTraj = 20000;
nPath = 10;
seed = 2026;
pngResolution = 300;

scriptPath = [mfilename('fullpath'),'.m'];
[sourceDirectory,scriptName] = fileparts(scriptPath);
addpath(sourceDirectory);
outputRoot = fullfile(sourceDirectory,scriptName);
mkdir(outputRoot);
runDirectory = fullfile(outputRoot,datestr(now,'yyyymmdd_HHMMSS'));
while exist(runDirectory,'dir')
    pause(1);
    runDirectory = fullfile(outputRoot,datestr(now,'yyyymmdd_HHMMSS'));
end
mkdir(runDirectory);

rng(seed,'twister');
nSteps = round(tEnd/dt);
time = (0:nSteps)'*dt;
op = integrator_COVLD_setup(M,K,Z,kBT,dt);

Ac = [0,1;
      -K/M,-Z/M];
E = expm(dt*Ac);
Qc = op.Sigma_eq-E*op.Sigma_eq*E';
Qc = (Qc+Qc')/2;
SigmaD = zeros(2);
SigmaC = zeros(2);
U = zeros(1,nTraj);
V = zeros(1,nTraj);
pathsU = zeros(nSteps+1,nPath);
pathsV = zeros(nSteps+1,nPath);
varE = zeros(nSteps+1,2);
varD = zeros(nSteps+1,2);
varC = zeros(nSteps+1,2);

for n = 1:nSteps
    [U,V] = integrator_COVLD(op,U,V);
    SigmaD = op.A*SigmaD*op.A'+op.Q;
    SigmaC = E*SigmaC*E'+Qc;
    varE(n+1,:) = var([U;V],0,2)';
    varD(n+1,:) = diag(SigmaD)';
    varC(n+1,:) = diag(SigmaC)';
    pathsU(n+1,:) = U(1:nPath);
    pathsV(n+1,:) = V(1:nPath);
end

save(fullfile(runDirectory,'results.mat'),...
    'M','K','Z','kBT','dt','tEnd',...
    'nTraj','nPath','seed','op','time',...
    'U','V','pathsU','pathsV','varE','varD','varC');

names = {'displacement','velocity'};
pathData = {pathsU,pathsV};
pathLabels = {'Displacement','Velocity'};
varianceLabels = {'Displacement variance / (k_BT/K)','Velocity variance / (k_BT/M)'};

for j = 1:2
    fig = figure('Color','w');
    plot(time,pathData{j},'LineWidth',1.0);
    xlabel('Time');
    ylabel(pathLabels{j});
    grid on;
    box on;
    local_save_figure(fig,runDirectory,['fig_harmonic_',names{j}],pngResolution);
end

for j = 1:2
    target = op.Sigma_eq(j,j);
    fig = figure('Color','w');
    hold on;
    plot([time(1),time(end)],[1,1],'-','LineWidth',1.5,'Color','k');
    plot(time,varE(:,j)/target,'-','LineWidth',1.5,'Color','b');
    plot(time,varD(:,j)/target,'--','LineWidth',2.5,'Color','r');
    plot(time,varC(:,j)/target,':','LineWidth',2.0,'Color','g');
    xlabel('Time');
    ylabel(varianceLabels{j});
    ylim([0,1.2]);
    grid on;
    box on;
    legend('Canonical value','Ensemble','Discrete covariance recursion','Continuous reference','Location','southeast');
    local_save_figure(fig,runDirectory,['fig_harmonic_',names{j},'_variance'],pngResolution);
end

fprintf('Completed. Results saved in %s\n',runDirectory);

function local_save_figure(fig,folder,name,resolution)
savefig(fig,fullfile(folder,[name,'.fig']));
print(fig,fullfile(folder,[name,'.eps']),'-depsc','-painters');
print(fig,fullfile(folder,[name,'.png']),'-dpng',sprintf('-r%d',resolution));
end
