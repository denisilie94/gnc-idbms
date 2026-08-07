function [F, mu, info] = design_frame(method, m, n, varargin)
%DESIGN_FRAME  Design an incoherent m x n frame with one of the four methods.
%  This is the single place where the tuned recipe of each method is defined.
%
%  [F, mu]       = design_frame(method, m, n)
%  [F, mu, info] = design_frame(method, m, n, 'seed', 0, 'K', 2000, ...)
%
%  method (case-insensitive):
%    'idb'        Distance-barrier baseline: weighted quadratic + hinge barrier,
%                 block-coordinate descent under a bisection on the target
%                 coherence. lambda follows the published n/m rule.
%    'idbms'      GNC-IDBMS: the IDB quadratic with the hinge replaced by an
%                 annealed multi-similarity (log-sum-exp) barrier, c: 2 -> 500,
%                 lambda = 12.                       [best for n/m <~ 3]
%    'idbms_mom'  GNC-IDBMS plus per-atom heavy-ball momentum 0.6 with adaptive
%                 restart.                           [best for n/m >~ 4]
%    'taap'       Targeted Accelerated Alternating Projections, authors'
%                 settings.
%  Aliases: 'gnc' = 'idbms', 'gnc_mom' = 'idbms_mom'.
%
%  Name/value options (defaults):
%    'seed'    (0)      RNG seed for the random initial frame
%    'K'       (2000)   sweeps per bisection step (IDB / IDBMS family)
%    'Kbis'    (15)     number of bisection steps
%    'gamma0'  ([])     initial step size; [] uses the method's default
%    'field'   ('real') TAAP field: 'real' | 'complex'
%    'verbose' (false)  print per-bisection / per-TAAP-round progress
%
%  Output:
%    F     - designed m x n frame with unit-norm columns
%    mu    - its mutual coherence, max_{i~=j} |f_i' f_j|
%    info  - struct with the parameters actually used

    ip = inputParser;
    ip.addParameter('seed',    0);
    ip.addParameter('K',       2000);
    ip.addParameter('Kbis',    15);
    ip.addParameter('gamma0',  []);
    ip.addParameter('field',   'real');
    ip.addParameter('verbose', false);
    ip.parse(varargin{:});
    o = ip.Results;

    switch lower(method)                  % accept the paper's short names
        case 'gnc',     method = 'idbms';
        case 'gnc_mom', method = 'idbms_mom';
    end

    bound = getBound(m, n);
    nm    = n/m;
    info  = struct('method',lower(method),'m',m,'n',n,'seed',o.seed, ...
                   'bound',bound,'nm',nm,'K',o.K,'Kbis',o.Kbis);

    % ---------------------------- TAAP ----------------------------------
    if strcmpi(method, 'taap')
        rng(o.seed);
        F0   = normcols(randn(m,n));
        opts = struct('beta',2,'N_budg',100000,'tau',1e-6,'N_p',100, ...
                      'eps_p',1e-3,'eps_s',1e-1,'accel',true, ...
                      'verbose',o.verbose);
        [F, mu, N] = taap(F0, m, n, o.field, opts);
        info.N_tot = N;
        info.opts  = opts;
        return
    end

    % ------------------- IDB / IDBMS family recipes ----------------------
    switch lower(method)
        case 'idb'
            % published trade-off rule for the barrier weight
            if     nm <  2,  lam = 0.2;
            elseif nm < 10,  lam = 0.5;
            elseif nm <= 20, lam = 1;
            else,            lam = 2;
            end
            p = struct('lambda', lam);
            btype = 'IDB';  g0 = 0.1;
            info.lambda = lam;

        case 'idbms'
            p = struct('lambda', 12, 'bsched', gnc_schedule(bound));
            p.beta = p.bsched(1);
            btype = 'IDBMS'; g0 = 0.1;
            info.lambda = 12; info.bsched = p.bsched;

        case 'idbms_mom'
            p = struct('lambda', 12, 'bsched', gnc_schedule(bound), ...
                       'momentum', 0.6, 'mom_restart', true);
            p.beta = p.bsched(1);
            btype = 'IDBMS'; g0 = 0.1;
            info.lambda = 12; info.bsched = p.bsched; info.momentum = 0.6;

        otherwise
            error('design_frame:unknownMethod', ...
                  'unknown method "%s" (idb | idbms | idbms_mom | taap)', method);
    end

    % the published special case for the extreme 20 x 5000 configuration
    if m == 20 && n == 5000 && strcmpi(method,'idb'), g0 = 0.01; end
    if ~isempty(o.gamma0), g0 = o.gamma0; end
    info.gamma0 = g0;
    info.btype  = btype;

    % search interval: Welch bound below, coherence of a random frame above
    rng(o.seed);
    Fr     = normcols(randn(m,n));
    mu_max = max(max(abs(Fr'*Fr - eye(n))));
    rng(o.seed);                     % re-seed: the solver draws its own frame
    [F, mu] = bisection_frame(bound, mu_max, o.Kbis, m, n, o.K, g0, ...
                              0.999, 5, 1e-4, btype, p, o.verbose);
end

function bs = gnc_schedule(bound)
%GNC_SCHEDULE  Annealing range for the barrier temperature beta.
%  The schedule is specified through the dimensionless product c = beta*mu,
%  run from c_0 = 2 to c_K = 500; here mu is approximated by the coherence
%  lower bound, which fixes the range once per frame size.
    bs = [max(round(2/bound), 2), round(500/bound)];
end
