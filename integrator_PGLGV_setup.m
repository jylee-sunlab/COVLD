function op = integrator_PGLGV_setup(M, K, Z, kBT, dt)
%INTEGRATOR_PGLGV_SETUP Prepare a covariance-matched Langevin integrator.

    narginchk(5, 5);
    validateattributes(K, {'double'}, {'real','finite','nonempty','2d'}, ...
        mfilename, 'K');
    if isvector(K)
        N = numel(K);
    else
        N = size(K, 1);
    end
    validateattributes(kBT, {'double'}, {'real','finite','scalar','positive'}, ...
        mfilename, 'kBT');
    validateattributes(dt, {'double'}, {'real','finite','scalar','positive'}, ...
        mfilename, 'dt');

    symmetryTolerance = 64 * max(1,N) * eps;
    [M, symmetryM] = local_matrix(M, N, 'M', symmetryTolerance);
    [K, symmetryK] = local_matrix(K, N, 'K', symmetryTolerance);
    [Z, symmetryZ] = local_matrix(Z, N, 'Z', symmetryTolerance);
    [LM, flagM] = chol(M, 'lower');
    [~, flagK] = chol(K, 'lower');
    [~, flagZ] = chol(Z, 'lower');
    if flagM ~= 0
        error('PGLGV:MassNotPD', 'M must be positive definite.');
    end
    if flagK ~= 0
        error('PGLGV:StiffnessNotPD', ...
            'K must be positive definite. Unconfined/free-particle modes are not supported by this setup.');
    end
    if flagZ ~= 0
        error('PGLGV:FrictionNotPD', ...
            ['Z must be positive definite for this positive-temperature, confined implementation. ' ...
             'Singular friction admits no positive step for the retained update.']);
    end
    matrixRcond = [rcond(M), rcond(K), rcond(Z)];
    if min(matrixRcond) < 1e-12
        warning('PGLGV:PoorScaling', ...
            'M, K or Z is poorly conditioned (minimum rcond %.3e). Consider rescaling the physical coordinates.', ...
            min(matrixRcond));
    end

    I = eye(N);
    Kh = LM \ (K / LM');
    Ch = LM \ (Z / LM');
    Kh = (Kh + Kh') / 2;
    Ch = (Ch + Ch') / 2;
    [Ah, F, termBound] = local_scaled_blocks(Kh, Ch, dt);
    if any(~isfinite(Ah(:))) || any(~isfinite(F(:))) || any(~isfinite(termBound(:)))
        error('PGLGV:AssemblyRange', 'Matrix assembly overflowed. Rescale the problem or reduce dt.');
    end

    dh = [dt^(3/2)*ones(N,1); sqrt(dt)*ones(N,1)];
    if any(dh == 0) || any(~isfinite(dh))
        error('PGLGV:TimeScalingRange', 'The time-step scaling underflows or overflows.');
    end

    frictionDiagonal = diag(Ch);
    if any(frictionDiagonal <= 0)
        error('PGLGV:MassTransform', 'Positive friction was lost in the mass-coordinate transformation.');
    end
    scaleFloor = sqrt(eps) * max(frictionDiagonal);
    d = sqrt(max(frictionDiagonal, scaleFloor));
    d = [d; d];
    H = F ./ (d*d');
    H = (H + H') / 2;
    scaledTermNorm = norm(termBound ./ (d*d'), 'fro');
    if any(~isfinite(H(:))) || ~isfinite(scaledTermNorm)
        error('PGLGV:CovarianceScalingRange', 'Scaled covariance assembly is not finite.');
    end

    roundoffMultiplier = 64;
    positivityTolerance = roundoffMultiplier*(N+4)*eps * ...
        max([1, norm(H,'fro'), scaledTermNorm]);
    eigenvalues = eig(H);
    minimumEigenvalue = min(eigenvalues);
    if minimumEigenvalue < -positivityTolerance
        error('PGLGV:Inadmissible', ...
            ['The required increment covariance is indefinite (scaled minimum eigenvalue %.6e). ' ...
             'Choose a smaller dt. No covariance correction is applied.'], minimumEigenvalue);
    elseif minimumEigenvalue <= positivityTolerance
        error('PGLGV:AdmissibilityUnresolved', ...
            ['Strict admissibility is unresolved at this step (minimum %.6e, threshold %.6e). ' ...
             'Boundary factors are not implemented. Reduce dt or rescale the problem.'], ...
            minimumEigenvalue, positivityTolerance);
    end
    [LH, flagH] = chol(H, 'lower');
    if flagH ~= 0
        error('PGLGV:NoiseFactorization', 'Cholesky factorization of the scaled covariance failed.');
    end
    noiseScale = sqrt(kBT) * (dh .* d);
    if any(noiseScale == 0) || any(~isfinite(noiseScale))
        error('PGLGV:NoiseScalingRange', 'The physical noise scaling underflows or overflows.');
    end
    Gh = bsxfun(@times, noiseScale, LH);
    Qh = kBT * ((dh*dh') .* F);
    Qh = (Qh + Qh') / 2;

    A = [LM' \ (Ah(1:N,1:N)*LM'), LM' \ (Ah(1:N,N+1:end)*LM'); ...
         LM' \ (Ah(N+1:end,1:N)*LM'), LM' \ (Ah(N+1:end,N+1:end)*LM')];
    G = [LM' \ Gh(1:N,:); LM' \ Gh(N+1:end,:)];
    Q = [LM' \ (Qh(1:N,1:N)/LM), LM' \ (Qh(1:N,N+1:end)/LM); ...
         LM' \ (Qh(N+1:end,1:N)/LM), LM' \ (Qh(N+1:end,N+1:end)/LM)];
    Q = (Q + Q') / 2;
    RK = chol(Kh);
    S0 = blkdiag(RK \ I, I);
    S_eq = sqrt(kBT) * blkdiag(LM' \ (RK \ I), LM' \ I);
    Sigma_eq = S_eq * S_eq';

    residualTolerance = 5e-11;
    factorResidual = norm(LH*LH' - H, 'fro') / max(norm(H,'fro'), realmin);
    W = S0 \ (Ah*S0);
    J = S0 \ (Gh/sqrt(kBT));
    balanceResidual = norm(eye(2*N) - W*W' - J*J', 'fro') / sqrt(2*N);
    physicalFactorResidual = norm(G*G' - Q, 'fro') / max(norm(Q,'fro'), realmin);
    if any(~isfinite([A(:); G(:); Q(:); Sigma_eq(:); factorResidual; balanceResidual; physicalFactorResidual]))
        error('PGLGV:NonfiniteBuild', 'The physical coefficients or residuals are not finite.');
    end
    if max([factorResidual, balanceResidual, physicalFactorResidual]) > residualTolerance
        error('PGLGV:CovarianceResidual', ...
            'The noise-factor or covariance-balance residual exceeds %.3e. Rescale the problem.', ...
            residualTolerance);
    end

    [Lp, flagP] = chol(3*M + dt*Z, 'lower');
    [Lq, flagQ] = chol(2*M + dt*Z, 'lower');
    if flagP ~= 0 || flagQ ~= 0
        error('PGLGV:MechanicalFactorization', 'A mass-plus-friction factorization failed.');
    end
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
    op.type = 'covariance-matched-langevin';
    op.version = '001';
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
    op.admissibility_status = 'strict_interior_resolved';
    op.diagnostics = struct( ...
        'symmetry_tolerance', symmetryTolerance, ...
        'input_symmetry_residuals', [symmetryM, symmetryK, symmetryZ], ...
        'matrix_rcond', matrixRcond, ...
        'friction_scaling_floor', scaleFloor, ...
        'roundoff_multiplier', roundoffMultiplier, ...
        'scaled_minimum_eigenvalue', minimumEigenvalue, ...
        'positivity_tolerance', positivityTolerance, ...
        'scaled_factor_residual', factorResidual, ...
        'physical_factor_residual', physicalFactorResidual, ...
        'covariance_balance_residual', balanceResidual, ...
        'residual_tolerance', residualTolerance, ...
        'spectral_radius', max(abs(eig(Ah))), ...
        'energy_norm', norm(W,2));
end

function [A, residual] = local_matrix(A, N, name, tolerance)
    validateattributes(A, {'double'}, {'real','finite','nonempty','2d'}, mfilename, name);
    if isscalar(A)
        A = A*eye(N);
    elseif isvector(A)
        if numel(A) ~= N
            error('PGLGV:MatrixSize', '%s must have %d diagonal entries.', name, N);
        end
        A = diag(A(:));
    elseif ~isequal(size(A), [N,N])
        error('PGLGV:MatrixSize', '%s must be %d-by-%d.', name, N, N);
    end
    A = full(A);
    residual = norm(A-A','fro') / max(norm(A,'fro'), realmin);
    if residual > tolerance
        error('PGLGV:NonSymmetric', '%s has a nonsymmetric part above the roundoff threshold.', name);
    end
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
    A = [Auu,Auv; Avu,Avv];
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
    F = [F11,F12; F12',F22];
    F = (F+F')/2;
    termBound = [abs(f11a)+abs(f11b), abs(f12a)+abs(f12b)+abs(f12c); ...
        (abs(f12a)+abs(f12b)+abs(f12c))', abs(f22a)+abs(f22b)+abs(f22c)+abs(f22d)];
end
