%RUN_EXPERIMENT  Single entry point for the incoherent-frame experiments.
%
%  Edit the CONFIGURATION block below (frame sizes and methods), then run this
%  script. It designs a frame with each selected method for each frame size and
%  prints a table of the resulting mutual coherences, together with the
%  theoretical lower bound and the run time.
%
%  Methods:
%    'idb'        IDB, the distance-barrier baseline (hinge barrier)
%    'idbms'      GNC-IDBMS: IDB quadratic + annealed multi-similarity barrier
%    'idbms_mom'  GNC-IDBMS with per-atom heavy-ball momentum
%    'taap'       Targeted Accelerated Alternating Projections
%  Aliases: 'gnc' = 'idbms', 'gnc_mom' = 'idbms_mom'.
%
%  Requires only base MATLAB (no toolboxes).

clear; clc;
addpath(fileparts(mfilename('fullpath')));

%% ============================ CONFIGURATION ============================

% Frame sizes, one [m n] pair per row. Add or remove rows freely.
% The default is a small, fast sanity check: it takes about 25 s in total and
% reproduces the published FLIP row 50 x 60 (idb 0.0624, idbms 0.06182,
% idbms_mom 0.06280, taap 0.06303).
frames = [ 50   60 ];   % FLIP, near square (n/m = 1.2)

% Further published sizes; uncomment to reproduce more of the tables. Run time
% grows quickly with n, from seconds at 64 x 128 to hours in the last rows.
% frames = [ 50   60      % FLIP,  near square       -> idbms wins
%            64  128      % CPM,   n/m = 2           -> idbms wins
%            64  256      % CPM,   n/m = 4           -> idbms_mom wins
%            23  500 ];   % TELET, n/m = 21.7        -> idbms_mom wins

% Methods to run (any subset, in any order).
methods = {'idb', 'idbms', 'idbms_mom', 'taap'};

seed    = 0;        % RNG seed for the random initial frame
K       = 2000;     % sweeps per bisection step (IDB / IDBMS family)
Kbis    = 15;       % number of bisection steps on the target coherence
verbose = false;    % true prints per-bisection (or per-TAAP-round) progress
savecsv = '';       % e.g. 'results.csv' to append the results; '' disables

%% =======================================================================

nF = size(frames, 1);
nM = numel(methods);
coh = nan(nF, nM);
sec = nan(nF, nM);

fprintf('\nIncoherent frame design | seed %d | K = %d | Kbis = %d\n', ...
        seed, K, Kbis);
fprintf('%s\n', repmat('-', 1, 13 + 11 + 12*nM));
fprintf('%5s %5s %10s', 'm', 'n', 'bound');
for jm = 1:nM, fprintf(' %11s', methods{jm}); end
fprintf('\n%s\n', repmat('-', 1, 13 + 11 + 12*nM));

for i = 1:nF
    m = frames(i,1);
    n = frames(i,2);
    fprintf('%5d %5d %10.5f', m, n, getBound(m, n));
    for jm = 1:nM
        tic;
        [~, coh(i,jm)] = design_frame(methods{jm}, m, n, ...
                            'seed', seed, 'K', K, 'Kbis', Kbis, ...
                            'verbose', verbose);
        sec(i,jm) = toc;
        fprintf(' %11.5f', coh(i,jm));
    end
    fprintf('\n');
end

fprintf('%s\n', repmat('-', 1, 13 + 11 + 12*nM));

% --- run times -----------------------------------------------------------
fprintf('\nRun time (seconds)\n');
fprintf('%5s %5s %10s', 'm', 'n', '');
for jm = 1:nM, fprintf(' %11s', methods{jm}); end
fprintf('\n');
for i = 1:nF
    fprintf('%5d %5d %10s', frames(i,1), frames(i,2), '');
    for jm = 1:nM, fprintf(' %11.1f', sec(i,jm)); end
    fprintf('\n');
end

% --- best method per size ------------------------------------------------
fprintf('\nLowest coherence per frame size\n');
for i = 1:nF
    [best, ib] = min(coh(i,:));
    fprintf('  %4d x %-5d  %-12s %.5f   (n/m = %.1f)\n', ...
            frames(i,1), frames(i,2), methods{ib}, best, ...
            frames(i,2)/frames(i,1));
end
fprintf('\n');

% --- optional CSV --------------------------------------------------------
if ~isempty(savecsv)
    fid = fopen(savecsv, 'a');
    for i = 1:nF
        for jm = 1:nM
            fprintf(fid, '%d,%d,%s,%d,%.6f,%.2f\n', frames(i,1), frames(i,2), ...
                    methods{jm}, seed, coh(i,jm), sec(i,jm));
        end
    end
    fclose(fid);
    fprintf('Appended %d rows to %s\n\n', nF*nM, savecsv);
end
