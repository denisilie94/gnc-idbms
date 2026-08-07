function [F_best, mu_best, N_tot] = taap(F0, m, n, field, opts)
%TAAP  Grassmannian frame computation via Targeted Accelerated Alternating
%  Projections (Massion & Massart, UCLouvain, SampTA 2025). MATLAB port of
%  the authors' reference implementation (taap-main/taap.py, v1.0, 30.01.2025).
%
%  [F, mu, N_tot] = taap(F0, m, n, field, opts)
%    F0    - initial frame (m x n), unit-norm columns
%    field - 'real' | 'complex'
%    opts  - struct with fields (defaults = the authors' run-script values):
%      .beta    = 2       target-update aggressiveness
%      .N_budg  = 100000  total budget of inner AAP iterations
%      .tau     = 1e-6    stop when delta_t = mu_best - t falls below this
%      .N_p     = 100     inner patience: iterations without accepted progress
%      .eps_p   = 1e-3    progress acceptance threshold (fraction of delta_t)
%      .eps_s   = 1e-1    success threshold (fraction of delta_t)
%      .accel   = true    Nesterov-style acceleration
%      .verbose = false
%      .bound   = []      coherence lower bound; default taap_bound(m,n,field)
%  Output:
%    F_best  - designed frame (m x n)
%    mu_best - its mutual coherence
%    N_tot   - total inner iterations used
%
%  Algorithm: alternate (i) clipping the Gram off-diagonal to the target t
%  (sig_proj_convex) and (ii) projecting onto rank-m PSD matrices via the
%  top-m eigenpairs (spec_proj_positive_truncated), with acceleration; an
%  outer loop adapts t toward/away from the lower bound depending on whether
%  the inner AAP succeeded (mu close to t) or stalled.
%
%  Port note: the Python taap() signature defaults are eps_p=1e-1, eps_s=1e-3,
%  but BOTH author run scripts (run_taap.py, run_taap_torch.py) pass
%  eps_p=1e-3, eps_s=1e-1; we default to the run-script values.

    if nargin < 5, opts = struct(); end
    beta_   = def(opts,'beta',2);
    N_budg  = def(opts,'N_budg',100000);
    tau     = def(opts,'tau',1e-6);
    N_p     = def(opts,'N_p',100);
    eps_p   = def(opts,'eps_p',1e-3);
    eps_s   = def(opts,'eps_s',1e-1);
    accel   = def(opts,'accel',true);
    verbose = def(opts,'verbose',false);
    lb      = def(opts,'bound',[]);
    if isempty(lb), lb = taap_bound(m, n, field); end

    if verbose
        fprintf('mu_0_AAP  \tmu_AAP    \ttarget    \tdelta_t   \tN_AAP \tN_tot\n');
    end

    N_tot   = 0;
    G_best  = F0' * F0;                       % ' = ctranspose: both fields
    mu_best = coh_gram(G_best, n);
    t       = lb;
    delta_t = mu_best - t;

    while ~(delta_t < tau || N_tot > N_budg)
        G_AAP  = G_best;
        mu_AAP = mu_best;
        k_AAP  = 0;

        c_k1 = 1.0;          % c_{k-1}
        G_k2 = G_best;       % G_{k-2}
        G_k1 = G_best;       % G_{k-1}
        k = 1;

        while ~(mu_AAP - t < eps_s*delta_t || k - k_AAP > N_p)
            if accel
                c_k = sqrt(4*c_k1^2 + 1)/2 + 1/2;
                Y_k = G_k1 + (c_k1 - 1)/c_k * (G_k1 - G_k2);
                G_k = spec_proj(sig_proj(Y_k, n, t), m, n);
            else
                c_k = c_k1;
                G_k = spec_proj(sig_proj(G_k1, n, t), m, n);
            end

            mu_k = coh_gram(normalize_gram(G_k), n);

            if mu_AAP - mu_k > eps_p*delta_t   % accepted progress
                G_AAP  = G_k;
                mu_AAP = mu_k;
                k_AAP  = k;
            end

            k = k + 1;
            G_k2 = G_k1;
            G_k1 = G_k;
            c_k1 = c_k;
            N_tot = N_tot + 1;
        end

        if verbose
            fprintf('%.6f \t%.6f \t%.6f \t%.6f \t%-6d \t%-6d\n', ...
                    mu_best, mu_AAP, t, delta_t, k-1, N_tot);
        end

        if mu_AAP - t < eps_s*delta_t          % success: push target down
            t = max(mu_AAP - beta_*delta_t, lb);
        else                                   % stall: relax target
            t = max(mu_AAP - delta_t/beta_, lb);
        end

        G_best  = normalize_gram(G_AAP);
        mu_best = mu_AAP;
        delta_t = mu_best - t;
    end

    F_best = reconstruct_frame(G_best, m, n);
end

function v = def(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function G = sig_proj(G, n, t)
%SIG_PROJ  Clip off-diagonal Gram magnitudes to t (keep sign/phase), diag = 1.
    aG  = abs(G);
    idx = aG >= t;
    G(idx) = t * (G(idx) ./ aG(idx));
    G(1:n+1:end) = 1;
end

function Gt = spec_proj(G, m, n)
%SPEC_PROJ  Rank-m projection from the top-m eigenpairs (spectrum NOT clipped
%  to >=0 here, exactly as in the reference implementation).
    G = (G + G')/2;
    if m <= n/4 && n >= 800     % partial decomposition pays off
        [Q, D] = eigs(G, m, 'largestreal');
        lam = diag(D);
    else
        [Q, D] = eig(G, 'vector');
        [~, ord] = sort(D, 'ascend');
        Q = Q(:, ord(n-m+1:end));
        lam = D(ord(n-m+1:end));
    end
    Gt = Q * (lam .* Q');
end

function G = normalize_gram(G)
    xinv = 1 ./ sqrt(real(diag(G)));
    G = G .* (xinv * xinv.');
end

function mu = coh_gram(G, n)
%COH_GRAM  Mutual coherence of a Gram matrix with unit diagonal.
    mu = max(max(abs(G - eye(n))));
end

function F = reconstruct_frame(G, m, n)
    [Q, D] = eig((G + G')/2, 'vector');
    [~, ord] = sort(D, 'ascend');
    lam = max(D(ord(n-m+1:end)), 0);
    F = diag(sqrt(lam)) * Q(:, ord(n-m+1:end))';
end
