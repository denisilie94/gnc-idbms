# gnc-idbms

MATLAB source for *Incoherent frame design via graduated non-convexity of a
multi-similarity barrier*.

Self-contained implementation of the three methods compared in the paper:
**IDB**, **GNC-IDBMS** (with and without momentum) and **TAAP**. Requires base
MATLAB only — no toolboxes.

## Quick start

Open `run_experiment.m`, edit the configuration block, and run it:

```matlab
frames  = [ 90 100        % one [m n] pair per row
            64 256
            23 500 ];
methods = {'idb', 'idbms', 'idbms_mom', 'taap'};
seed    = 0;
```

It designs a frame with each method for each size and prints the mutual
coherence, the lower bound and the run time. That single file is the only one
you need to edit to change the frame sizes or the methods.

To call a single design directly:

```matlab
[F, mu] = design_frame('idbms', 90, 100, 'seed', 0);
```

## Methods

| name | description |
|---|---|
| `idb` | Distance-barrier baseline: weighted quadratic + hinge barrier, block-coordinate descent under a bisection on the target coherence; `lambda` follows the published `n/m` rule. |
| `idbms` | **GNC-IDBMS**: the hinge is replaced by an annealed multi-similarity (log-sum-exp) barrier, `c = beta*mu` run from 2 to 500, `lambda = 12`. Best for `n/m <~ 3`. |
| `idbms_mom` | GNC-IDBMS plus per-atom heavy-ball momentum `0.6` with adaptive restart. Best for `n/m >~ 4`. |
| `taap` | Targeted Accelerated Alternating Projections, with the authors' settings. |

Aliases `gnc` and `gnc_mom` map to `idbms` and `idbms_mom`.

## Files

```
run_experiment.m     single driver: set frame sizes + methods, then run
design_frame.m       method dispatcher; the tuned recipes live here
idb_frame.m          block-coordinate solver (barrier, GNC annealing, momentum)
bisection_frame.m    bisection on the target coherence
frame_grad_idb.m     IDB gradient: weighted quadratic + hinge barrier
frame_grad_idbms.m   GNC-IDBMS gradient: quadratic + multi-similarity barrier
stable_softmax.m     numerically stable softmax used by the barrier
taap.m               TAAP (port of the authors' reference implementation)
taap_bound.m         composite coherence lower bound used by TAAP
getBound.m           Welch-type lower bound (Tahir et al., see file header)
normcols.m           unit-norm columns (replaces the toolbox function normc)
```

## Protocol

The published tables use one run per cell at seed 0, with a common budget for
every method: `K = 2000` sweeps, 5 step halvings, `gamma0 = 0.1`, `rho = 0.999`,
tolerance `1e-4`, and 15 bisection steps. These are the defaults of
`design_frame`, so the numbers in the paper are reproduced by

```matlab
[~, mu] = design_frame(method, m, n, 'seed', 0);
```

Note that results can move by a few units in the fifth decimal between machines
or thread counts: the solver is chaotic, so the different summation order of a
different BLAS reduction is amplified by the bisection. The ordering of the
methods is unaffected.

## Citing

If you use this code, or build on it, please cite our work:

> D. C. Ilie-Ablachim and B. Dumitrescu, "Incoherent frame design via graduated
> non-convexity of a multi-similarity barrier", National University of Science
> and Technology Politehnica Bucharest, Department of Automatic Control and
> Computers, 313 Spl. Independenței, 060042 Bucharest, Romania.

```bibtex
@article{IlieAblachim_gncidbms,
  author  = {Ilie-Ablachim, Denis C. and Dumitrescu, Bogdan},
  title   = {Incoherent frame design via graduated non-convexity of a
             multi-similarity barrier},
  journal = {},
  year    = {},
  note    = {Fill in the venue and year once the paper appears}
}
```

The method built on here was introduced in

> D. C. Ilie-Ablachim and B. Dumitrescu, "Incoherent frames design and dictionary
> learning using a distance barrier", *Signal Processing*, 2023.

and the TAAP baseline is due to Massion and Massart, "Targeted accelerated
alternating projections", SampTA 2025.

## Funding

This work is supported by a grant of the Ministry of Research, Innovation and
Digitization, CNCS – UEFISCDI, project number PN-III-P4-PCE-2021-0154, within
PNCDI III.
