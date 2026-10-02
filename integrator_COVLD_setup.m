function op = integrator_COVLD_setup(M, K, Z, kBT, dt)
% INTEGRATOR_COVLD_SETUP
% Prepare the covariance-matched Langevin integrator.
% M, K, and Z are symmetric positive-definite physical matrices.

if isvector(K)
    N = numel(K);
else
    N = size(K,1);
end

M = local_matrix(M,N);
K = local_matrix(K,N);
Z = local_matrix(Z,N);
LM = chol(M,'lower');
chol(K);
chol(Z);

I = eye(N);
Kh = LM \ (K / LM');
Ch = LM \ (Z / LM');
Kh = (Kh+Kh')/2;
Ch = (Ch+Ch')/2;
[Ah,F] = local_scaled_blocks(Kh,Ch,dt);

dh = [dt^(3/2)*ones(N,1);
      sqrt(dt)*ones(N,1)];
d = sqrt(diag(Ch));
d = [d;
     d];
H = F./(d*d');
H = (H+H')/2;
LH = chol(H,'lower');
Gh = bsxfun(@times,sqrt(kBT)*(dh.*d),LH);

A = [LM' \ (Ah(1:N,1:N)*LM'), LM' \ (Ah(1:N,N+1:end)*LM');
     LM' \ (Ah(N+1:end,1:N)*LM'), LM' \ (Ah(N+1:end,N+1:end)*LM')];
G = [LM' \ Gh(1:N,:);
     LM' \ Gh(N+1:end,:)];
Q = G*G';
Q = (Q+Q')/2;
Sigma_eq = kBT*blkdiag(K\I,M\I);

Lp = chol(3*M+dt*Z,'lower');
Lq = chol(2*M+dt*Z,'lower');

op.N = N;
op.M = M;
op.K = K;
op.Z = Z;
op.kBT = kBT;
op.dt = dt;
op.Lp = Lp;
op.Lq = Lq;
op.Bp = 3*M-2*dt*Z;
op.Bq = 2*M-dt*Z;
op.A = A;
op.Q = Q;
op.G = G;
op.Gu = G(1:N,:);
op.Gv = G(N+1:end,:);
op.Sigma_eq = Sigma_eq;
end

function A = local_matrix(A,N)
if isscalar(A)
    A = A*eye(N);
elseif isvector(A)
    A = diag(A(:));
end
A = full((A+A')/2);
end

function [A,F] = local_scaled_blocks(K,C,h)
N = size(K,1);
I = eye(N);
R = (3*I+h*C) \ I;
T = (2*I+h*C) \ I;
Auu = I-(3*h^2/2)*R*K;
Auv = (h/2)*R*(6*I-h*C);
Avu = -h*T*K*(I+Auu);
Avv = T*(2*I-h*C-h*K*Auv);
A = [Auu,Auv;
     Avu,Avv];
B = R*C*(12*I-h*C);
F11 = 6*R*C*R-(h/4)*R*(C*C+9*K)*R;
F12 = (B/2)*T+(h/4)*R*R*(h^2*C*C-30*h*C-18*I)*K*T+(9*h^3/4)*R*K*R*K*T;
F22 = 8*T*C*T-(h^2/2)*T*(B*K+K*B)*T+(h^3/4)*T*K*R*R*(36*I+36*h*C-h^2*C*C)*K*T-(9*h^5/4)*T*K*R*K*R*K*T;
F = [F11,F12;
     F12',F22];
F = (F+F')/2;
end
