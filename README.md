# When Data Augmentation breaks Uncertainty Quantification: Finite-Sample Inference for Data-Augmented Estimators

## Reproducing the experiments

Run both experiment suites from the repository root:

```sh
Rscript run_all_experiments.R
```

Alternatively, uncomment
```
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```
in the file `run_all_experiments.R`
and execute the code in `run_all_experiments.R` in `R` directly.

The script for the coverage expample in Section 4.3 of the paper  (`coverage_example_section4.R`) uses `tidyverse` and `latex2exp`.
The experiment scripts use base R; `plot_results.R` additionally requires
`ggplot2`. The scripts write CSV files to `results/` and regenerate
the appendix diagnostics in `plots/`. Random seeds are fixed in each
experiment script.

The file `coverage_example_section4.R` can be used to reproduce the results in Section 4.3.

For the experiments in Section 6 of the paper, the generated design summaries are
`results/classical_design_summary.csv` and
`results/high_dim_design_summary.csv`; the latter also records the sparsity,
correlation parameter, compatibility constant, Lasso penalty, and localization
radius used in each high-dimensional benchmark.

The interval coverage and relative width of the intervals for each design/augmentation mechanism/direction-combination are reported in `results/classical_results_by_direction.csv` and `results/high_dim_results_by_direction.csv`, respectively. 
The aggregate CSVs `results/classical_results.csv` and `results/high_dim_results.csv` contain the results averaged over the three directions (as in Section D 2.1 and D 2.2 of the paper).
The plot in the folder `plots/` show
pointwise 95% confidence intervals .

The classical script uses seed `1000`, 1,000 repetitions, and one fixed
50/50 row split per design (4 designs in total). The high-dimensional script uses seed `1000`
and 500 repetitions. Its four controlled designs have
`n_cal = n_inf = 120` and `p = 300`. The same fixed designs,
signals, target directions, and augmentation matrices are used throughout all repetitions.

## Probability budgets

All reported procedures target total failure probability `0.025`.

- Credible and Oracle use `alpha = 0.025` directly.
- The classical OLS split interval uses `delta = alpha_inf = 0.0125`.
- The certified high-dimensional interval uses
  `delta = eta = alpha_inf = 0.025 / 3`.
- The zero-certificate diagnostic uses `delta = 0.025 / 3` and
  `alpha_inf = 0.025 - delta` (no `eta` term required).

## High-dimensional compatibility constant

The high-dimensional fixed designs are block-orthogonal between the true
support and nuisance coordinates. Therefore, the smallest eigenvalue of the standardized signal-block Gram matrix

$$\phi_0^2 = \lambda_{\min}(\mathbf{X}_{\mathrm{inf},\mathcal{S}}^\top \mathbf{X}_{\mathrm{inf},\mathcal{S}} / n_{\mathrm{cal}})$$

is a certified compatibility constant where $\lambda_{\min}$ denotes the smallest eigenvalue. The Lasso localization radius is the
paper's theoretical value $12  s  \lambda / \phi_0^2$. No empirical radius
multiplier is used.

## Selection of bias certificate
For every high-dimensional design/augmentation/direction cell, the certificate
is chosen before observing any inference response. Candidate vectors are zero
and 25 ridge solutions

$\mathbf{a} (\kappa) = (\mathbf{X}_{\mathrm{inf}} \mathbf{X}_{\mathrm{inf}}^\top + \kappa \mathbf{I}_{n_{\mathrm{inf}}})^{-1} \mathbf{X}_{\mathrm{inf}} \mathbf{g}_v$

$\kappa = tr(\mathbf{X}_{\mathrm{inf}} \mathbf{X}_{\mathrm{inf}}^\top) / n_{\mathrm{inf}} \cdot 10^t$

with $t$ equally spaced from -5 to 3.

The selected candidate minimizes

$$z_{1 - \eta/2}  \sigma  \Vert \mathbf{a} \Vert_2 + R  \Vert \mathbf{g}_v - \mathbf{X}_{\mathrm{inf}}^\top \mathbf{a} \Vert_\infty.$$

 This search depends on
the fixed design but not on inference noise.

## Data

The classical benchmarks use the original files downloaded from the UCI
Machine Learning Repository and can be found in the `data/` folder:

- `wine.data` (SHA-256 `6be6b1203f3d51df0b553a70e57b8a723cd405683958204f96d23d7cd6aea659`)
- `wdbc.data` (SHA-256 `d606af411f3e5be8a317a5a8b652b425aaf0ff38ca683d5327ffff94c3695f4a`)
- `auto-mpg.data` (SHA-256 `48b830e11feee5572525f8f1691ddb9d38d3d7b7063edcd8fca672c2a5e17d8d`)
