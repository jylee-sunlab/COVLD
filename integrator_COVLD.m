function [U_np1,V_np1] = integrator_COVLD(op,U_n,V_n,R_n,R_np1,eta)
% INTEGRATOR_COVLD
% Advance one covariance-matched Langevin step.
% U and V are N-by-L arrays with one trajectory per column.

if nargin < 4 || isempty(R_n)
    R_n = zeros(op.N,1);
end
if nargin < 5 || isempty(R_np1)
    R_np1 = R_n;
end
if nargin < 6 || isempty(eta)
    eta = randn(2*op.N,size(U_n,2));
end

dt = op.dt;
dU = op.Gu*eta;
F_n = op.K*U_n;
R_pred = bsxfun(@plus,op.Bp*V_n-3*dt*F_n,3*dt*R_n);
V_pred = op.Lp' \ (op.Lp \ R_pred)+(2/dt)*dU;
U_np1 = U_n+(dt/2)*(V_n+V_pred);
F_np1 = op.K*U_np1;
R_corr = op.Bq*V_n-dt*(F_n+F_np1)+dt*(op.K*dU);
R_corr = bsxfun(@plus,R_corr,dt*R_n);
R_corr = bsxfun(@plus,R_corr,dt*R_np1);
V_np1 = op.Lq' \ (op.Lq \ R_corr)+op.Gv*eta;
end
