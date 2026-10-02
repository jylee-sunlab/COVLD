function op = integrator_COVLD_setup(M, K, Z, kBT, dt)
% INTEGRATOR_COVLD_SETUP
% Prepare a covariance-matched Langevin integrator.
% Inputs are physical M, K and Z, thermal energy kBT, and time step dt.
% M, K and Z must be symmetric positive definite, with positive kBT and dt.

narginchk(5, 5);
if isvector(K)
    N = numel(K);
else
    N = size(K, 1);
end

[M, ~] = local_matrix(M, N);
[K, ~] = local_matrix(K, N);
[Z, ~] = local_matrix(Z, N);
[LM, ~] = chol(M, 'lower');

I = eye(N);
Kh = LM \ (K / LM');
Ch = LM \ (Z / LM');
Kh = (Kh + Kh') / 2;
Ch = (Ch + Ch') / 2;
[Ah, F, ~] = local_scaled_blocks(Kh, Ch, dt);

dh = [dt^(3/2)*ones(N,1);
      sqrt(dt)*ones(N,1)];

frictionDiagonal = diag(Ch);
scaleFloor = sqrt(eps) * max(frictionDiagonal);
d = sqrt(max(frictionDiagonal, scaleFloor));
d = [d;
     d];
H = F ./ (d*d');
H = (H + H') / 2;

[LH, ~] = chol(H, 'lower');
noiseScale = sqrt(kBT) * (dh .* d);
Gh = bsxfun(@times, noiseScale, LH);
Qh = kBT * ((dh*dh') .* F);
Qh = (Qh + Qh') / 2;

A = [LM' \ (Ah(1:N,1:N)*LM'), LM' \ (Ah(1:N,N+1:end)*LM');
     LM' \ (Ah(N+1:end,1:N)*LM'), LM' \ (Ah(N+1:end,N+1:end)*LM')];
G = [LM' \ Gh(1:N,:);
     LM' \ Gh(N+1:end,:)];
Q = [LM' \ (Qh(1:N,1:N)/LM), LM' \ (Qh(1:N,N+1:end)/LM);
     LM' \ (Qh(N+1:end,1:N)/LM), LM' \ (Qh(N+1:end,N+1:end)/LM)];
Q = (Q + Q') / 2;
RK = chol(Kh);
S_eq = sqrt(kBT) * blkdiag(LM' \ (RK \ I), LM' \ I);
Sigma_eq = S_eq * S_eq';

[Lp, ~] = chol(3*M + dt*Z, 'lower');
[Lq, ~] = chol(2*M + dt*Z, 'lower');
dtLimit = NaN;
if N == 1
    omega = sqrt(K/M);
    chi = (Z/M)/omega;
    if chi <= 1
        y = chi*chi;
        omegaLimit = 4*chi*(2*y + sqrt(7*y*y + 18*y + 27))/(y + 3)^2;
    else
        y = (1/chi)^2;
        omegaLimit = (4/chi)*(2 + sqrt(7 + 18*y + 27*y*y))/(1 + 3*y)^2;
    end
    dtLimit = omegaLimit/omega;
end

op = struct();
op.N = N;
op.M = M;
op.K = K;
op.Z = Z;
op.kBT = kBT;
op.dt = dt;
op.Lp = Lp;
op.Lq = Lq;
op.Bp = 3*M - 2*dt*Z;
op.Bq = 2*M - dt*Z;
op.A = A;
op.Q = Q;
op.G = G;
op.Gu = G(1:N,:);
op.Gv = G(N+1:end,:);
op.Sigma_eq = Sigma_eq;
op.S_eq = S_eq;
op.dt_limit = dtLimit;
end

function [A, residual] = local_matrix(A, N)
if isscalar(A)
    A = A*eye(N);
elseif isvector(A)
    A = diag(A(:));
end
A = full(A);
residual = norm(A-A','fro') / max(norm(A,'fro'), realmin);
A = (A+A')/2;
end

function [A,F,termBound] = local_scaled_blocks(K,C,h)
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
f11a = 6*R*C*R;
f11b = -(h/4)*R*(C*C+9*K)*R;
f12a = (B/2)*T;
f12b = (h/4)*R*R*(h^2*C*C-30*h*C-18*I)*K*T;
f12c = (9*h^3/4)*R*K*R*K*T;
f22a = 8*T*C*T;
f22b = -(h^2/2)*T*(B*K+K*B)*T;
f22c = (h^3/4)*T*K*R*R*(36*I+36*h*C-h^2*C*C)*K*T;
f22d = -(9*h^5/4)*T*K*R*K*R*K*T;
F11 = f11a+f11b;
F12 = f12a+f12b+f12c;
F22 = f22a+f22b+f22c+f22d;
F = [F11,F12;
     F12',F22];
F = (F+F')/2;
termBound = [abs(f11a)+abs(f11b), abs(f12a)+abs(f12b)+abs(f12c);
    (abs(f12a)+abs(f12b)+abs(f12c))', abs(f22a)+abs(f22b)+abs(f22c)+abs(f22d)];
end
