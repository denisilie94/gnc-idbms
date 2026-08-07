function [Fbest, coh_best, iter_total_count] = bisection_frame(mu_min, mu_max, ...
                                     nr_it_bis, m, n, K, gamma0, rho, search_it, ...
                                     tol_stop, btype, p, verbose)
%BISECTION_FRAME  Best-coherence frame design by bisection on the target.
%  A solve (idb_frame) is run at the middle of the current coherence interval.
%  The design counts as a success when the achieved coherence comes within
%  10*tol_stop of the target, in which case the upper end of the interval is
%  lowered; otherwise the lower end is raised.
%
%  Input:
%    mu_min, mu_max - initial search interval for the target coherence
%    nr_it_bis      - maximum number of bisection steps
%    m, n           - frame size
%    K, gamma0, rho, search_it, tol_stop - idb_frame parameters
%    btype, p       - barrier name and parameters (see idb_frame)
%    verbose        - print per-bisection progress (default false)
%  Output:
%    Fbest, coh_best  - best frame found and its mutual coherence
%    iter_total_count - total number of sweeps over all solves

    if nargin < 13 || isempty(verbose), verbose = false; end

    iter_total_count = 0;
    coh_best = 1;
    Fbest    = [];

    for i_bis = 1 : nr_it_bis
        mu = (mu_min + mu_max)/2;

        [F, coh] = idb_frame(m, n, mu, K, gamma0, rho, search_it, tol_stop, btype, p);
        iter_total_count = iter_total_count + length(coh);
        coh_crt = min(coh);

        if coh_crt < coh_best
            coh_best = coh_crt;
            Fbest    = F;
        end

        if coh_crt < mu + 10*tol_stop        % successful design
            if coh_crt > mu_max              % no further improvement possible
                break
            end
            mu_max = coh_best;
        else
            mu_min = mu;
            mu_max = min(coh_crt, mu_max);   % the upper end may also improve
        end

        if verbose
            fprintf('  bis %2d: mu = %.5f  coh = %.5f  best = %.5f  [%.5f, %.5f]\n', ...
                    i_bis, mu, coh_crt, coh_best, mu_min, mu_max);
        end

        if mu_max - mu_min < tol_stop        % interval exhausted
            break
        end
    end
end
