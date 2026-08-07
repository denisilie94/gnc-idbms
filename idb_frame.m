function [F, cohv] = idb_frame(m, n, desired_coh, K, gamma0, rho, search_it, tol_stop, btype, p)
%IDB_FRAME  Incoherent frame design by block-coordinate descent on a barrier.
%  One sweep updates every atom in random order; each atom takes a gradient
%  step accepted by a backtracking (step-halving) rule, and the mutual
%  coherence is recomputed once per sweep.
%
%  Input:
%    m, n         - frame size (m x n, m < n)
%    desired_coh  - target mutual coherence mu (drives the barrier and the stop)
%    K            - maximum number of sweeps
%    gamma0       - initial gradient step size
%    rho          - per-sweep step-size decrease factor (< 1)
%    search_it    - maximum number of halvings in the line search
%    tol_stop     - stopping tolerance
%    btype        - 'IDB'   : weighted quadratic + hinge barrier
%                   'IDBMS' : weighted quadratic + multi-similarity barrier
%    p            - parameter struct:
%                     .lambda        barrier weight
%                     .beta          barrier temperature (IDBMS)
%                     .bsched        [beta_lo beta_hi] geometric annealing
%                                    over the K sweeps (IDBMS); omit for a
%                                    fixed temperature
%                     .momentum      heavy-ball coefficient eta (0 = off)
%                     .mom_restart   true to zero the velocities when the
%                                    coherence rises above 1.02*best
%  Output:
%    F            - frame with the best coherence found over the sweeps
%    cohv         - mutual coherence after each sweep

    p.mu = desired_coh;              % keep the barrier target in sync

    F       = normcols(randn(m,n));
    cohv    = zeros(1,K);
    coh_min = 1;
    F_best  = F;
    diag_idx = 1:n+1:n*n;            % linear indices of the Gram diagonal

    % geometric annealing of the barrier temperature (graduated non-convexity)
    anneal = isfield(p,'bsched') && ~isempty(p.bsched);
    if anneal
        b_lo = p.bsched(1);
        b_hi = p.bsched(2);
    end

    % per-atom heavy-ball momentum
    mom = 0;
    if isfield(p,'momentum'), mom = p.momentum; end
    mom_restart = isfield(p,'mom_restart') && p.mom_restart;
    if mom > 0, Pv = zeros(m,n); end     % accepted increment per atom

    for k = 1 : K
        if anneal
            t = (k-1)/max(K-1,1);
            p.beta = b_lo * (b_hi/b_lo)^t;
        end

        perm = randperm(n-1) + 1;    % a sweep of the atoms in random order
        for j = perm
            d = F(:,j);
            v = F'*d;
            v(j) = 0;

            % ---- barrier gradient (the only method-specific part) ----
            switch upper(btype)
                case 'IDB',   g = frame_grad_idb(v, d, F, j, p);
                case 'IDBMS', g = frame_grad_idbms(v, d, F, j, p);
                otherwise
                    error('idb_frame:unknownBarrier', ...
                          'unknown barrier "%s" (IDB | IDBMS)', btype);
            end

            % ---- backtracking line search --------------------------------
            % with dn = (d - gamma*g)/||d - gamma*g||, we have
            % F'*dn = (v - gamma*F'*g)/||d - gamma*g||, so each trial is O(n)
            Fg = F'*g;
            Fg(j) = 0;
            c_max = max(abs(v));
            g_now = gamma0;

            if mom > 0
                pj = Pv(:,j);
                Fpj = F'*pj;
                Fpj(j) = 0;
                for i = 1 : search_it
                    dvec = d - g_now*g + mom*pj;
                    nrm  = norm(dvec);
                    vn   = (v - g_now*Fg + mom*Fpj)/nrm;
                    if max(abs(vn)) < c_max, break; end
                    g_now = g_now/2;
                end
                newd = dvec/nrm;
                Pv(:,j) = newd - d;          % store the accepted increment
                F(:,j)  = newd;
            else
                for i = 1 : search_it
                    dvec = d - g_now*g;
                    nrm  = norm(dvec);
                    vn   = (v - g_now*Fg)/nrm;
                    if max(abs(vn)) < c_max, break; end
                    g_now = g_now/2;
                end
                F(:,j) = dvec/nrm;
            end
        end
        gamma0 = gamma0*rho;

        % mutual coherence of the sweep (unit-norm columns => diag(G) = 1)
        G = F'*F;
        G(diag_idx) = 0;
        cohv(k) = max(abs(G(:)));
        if cohv(k) < coh_min
            coh_min = cohv(k);
            F_best  = F;
        end
        if mom > 0 && mom_restart && cohv(k) > 1.02*coh_min
            Pv = zeros(m,n);         % adaptive restart: the velocity overshot
        end
        if coh_min - desired_coh < tol_stop     % close enough to the target
            cohv = cohv(1:k);
            break
        end
    end

    F = F_best;
end
