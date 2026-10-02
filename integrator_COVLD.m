function [U_np1, V_np1] = integrator_COVLD(op, U_n, V_n, R_n, R_np1, eta)
% INTEGRATOR_COVLD
% Advance the covariance-matched Langevin update.
% U and V are N-by-L arrays with one trajectory per column.
% R_n and R_next are physical forces of size N-by-1 or N-by-L.
% An omitted R_n is zero, and an omitted R_next equals R_n.

narginchk(3,6);
nTraj = size(U_n,2);
if nargin < 4 || isempty(R_n)
    R_n = zeros(op.N,1);
end
if nargin < 5 || isempty(R_np1)
    R_np1 = R_n;
end
if nargin < 6 || isempty(eta)
    eta = randn(2*op.N,nTraj);
end
dt = op.dt;
dU = op.Gu*eta;
dV = op.Gv*eta;
F_n = op.K*U_n;
R_pred = bsxfun(@plus, op.Bp*V_n-3*dt*F_n, 3*dt*R_n);
V_pred = op.Lp' \ (op.Lp \ R_pred) + (2/dt)*dU;
U_np1  = U_n+(dt/2)*(V_n+V_pred);
F_next = op.K*U_np1;
R_corr = op.Bq*V_n-dt*(F_n+F_next)+dt*(op.K*dU);
R_corr = bsxfun(@plus, R_corr, dt*R_n);
R_corr = bsxfun(@plus, R_corr, dt*R_np1);
V_np1  = op.Lq' \ (op.Lq \ R_corr) + dV;
end
