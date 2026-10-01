function [U_next, V_next] = integrator_COVLD(op, U, V, R_n, R_next, eta)
% INTEGRATOR_COVLD
% Advance the covariance-matched Langevin update.
% U and V are N-by-L arrays with one trajectory per column.
% R_n and R_next are physical forces of size N-by-1 or N-by-L.
% An omitted R_n is zero, and an omitted R_next equals R_n.

narginchk(3,6);
if ~isstruct(op) || ~isscalar(op) || ~isfield(op,'type') || ...
        ~strcmp(op.type,'covariance-matched-langevin') || ...
        ~isfield(op,'admissibility_status') || ...
        ~strcmp(op.admissibility_status,'strict_interior_resolved')
    error('PGLGV:InvalidOperator', 'Use op returned by integrator_PGLGV_setup.');
end
validateattributes(U, {'double'}, {'real','finite','nonempty','2d'}, mfilename, 'U');
validateattributes(V, {'double'}, {'real','finite','nonempty','2d'}, mfilename, 'V');
if size(U,1) ~= op.N || ~isequal(size(U),size(V))
    error('PGLGV:StateSize', 'U and V must have the same N-by-Ntraj size.');
end
nTraj = size(U,2);
if nargin < 4 || isempty(R_n)
    R_n = zeros(op.N,1);
end
if nargin < 5 || isempty(R_next)
    R_next = R_n;
end
local_force_size(R_n, op.N, nTraj, 'R_n');
local_force_size(R_next, op.N, nTraj, 'R_next');
if nargin < 6 || isempty(eta)
    eta = randn(2*op.N,nTraj);
else
    validateattributes(eta, {'double'}, {'real','finite','nonempty','2d'}, mfilename, 'eta');
    if ~isequal(size(eta), [2*op.N,nTraj])
        error('PGLGV:NoiseSize', 'eta must have size 2*N-by-Ntraj.');
    end
end

h = op.dt;
dU = op.Gu*eta;
dV = op.Gv*eta;
F_n = op.K*U;
predictorRhs = bsxfun(@plus, op.Bp*V-3*h*F_n, 3*h*R_n);
V_predictor = op.Lp' \ (op.Lp \ predictorRhs) + (2/h)*dU;
U_next = U+(h/2)*(V+V_predictor);

F_next = op.K*U_next;
correctorRhs = op.Bq*V-h*(F_n+F_next)+h*(op.K*dU);
correctorRhs = bsxfun(@plus, correctorRhs, h*R_n);
correctorRhs = bsxfun(@plus, correctorRhs, h*R_next);
V_next = op.Lq' \ (op.Lq \ correctorRhs) + dV;
if any(~isfinite(U_next(:))) || any(~isfinite(V_next(:)))
    error('PGLGV:NonfiniteStep', 'The step returned a nonfinite displacement or velocity.');
end
end

function local_force_size(R,N,nTraj,name)
validateattributes(R, {'double'}, {'real','finite','nonempty','2d'}, mfilename, name);
if size(R,1) ~= N || ~(size(R,2) == 1 || size(R,2) == nTraj)
    error('PGLGV:ForceSize', '%s must be N-by-1 or N-by-Ntraj.', name);
end
end
