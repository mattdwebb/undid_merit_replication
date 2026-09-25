# MERIT scholarship replication

Code, data, and verified estimates for Panel C of the MERIT scholarship example
in the UN-DID paper. The revised paper adopts the computed results below,
replacing the four UN-DID entries in the earlier draft. CSDID and DID-INT retain
their verified values. The replication uses Stata for CSDID/CSDIDJACK and Julia
for UN-DID and DID-INT.

## Revised results - 22 September 2026

| Aggregation | UN-DID ATT | UN-DID SE | CSDID ATT | CSDID SE | DID-INT ATT | DID-INT SE |
|---|---:|---:|---:|---:|---:|---:|
| simple | 0.0466 | 0.0113 | 0.0464 | 0.0133 | 0.0464 | 0.0102 |
| group | 0.0458 | 0.0133 | 0.0339 | 0.0211 | 0.0458 | 0.0084 |

Download the [results CSV](results/panel_c_results.csv) or the
[replacement LaTeX discussion and table](paper/merit_results.tex).
[LaTeX insertion instructions](paper/README.md) and the
[validation record](results/VALIDATION.md) give additional detail.

### Specifications and standard errors

| Method | Covariates | Reported standard error |
|---|---|---|
| UN-DID | None (`covariates=false` at aggregation) | Final-stage delete-one-subgroup jackknife |
| CSDID | `male`, `black`, `asian` | State-cluster jackknife via `csdidjack` |
| DID-INT | `male`, `black`, `asian`; state-varying coefficients (`ccc="state"`) | Historical final-stage delete-one-subgroup jackknife |

UN-DID and DID-INT explicitly use `both` weighting. Their reported jackknives
delete one of 33 cohort-time estimates for simple aggregation, or one of seven
cohort estimates for group aggregation, and renormalize the remaining weights.
**These are not state-cluster standard errors.** Comparisons therefore involve
different covariate specifications and inference procedures. The replacement
paper text states these distinctions explicitly.

The earlier draft reported UN-DID ATT/SE pairs of 0.0485/0.0110 and
0.0459/0.0188. Those four entries were not reproduced from the raw data and are
superseded by the results above. Their provenance remains documented in
[REPLICATION_FINDINGS.md](REPLICATION_FINDINGS.md); adopting a new baseline
does not resolve the historical discrepancies.

## Requirements and setup

- Stata 16 or later, with `drdid`, `csdid`, and `csdidjack` installed.
- Julia 1.11.7 (the version used for verification), available as `julia` on PATH.
- The repository's Julia environment in `julia/`. Its manifest records the exact
  package trees: DiDInt.jl 0.7.5 and Undid.jl 0.5.2, plus their dependencies.

From the repository directory, initialize the environment once:

```sh
julia --project=julia -e "using Pkg; Pkg.instantiate()"
```

If the Stata commands are missing, install them:

```stata
ssc install drdid, replace
ssc install csdid, replace
net install csdidjack, from("https://raw.githubusercontent.com/liu-yunhan/csdidjack/main/") replace
```

The verified Stata commands are `csdid` 1.81 and `csdidjack` 0.5.2. The DO file
records the commands found on the local machine. Updating packages may change
results; preserve the Julia manifest when reproducing this run.

## Run

Set Stata's working directory to this repository and execute:

```stata
do replicate_merit.do
```

If Julia is not on PATH, first specify its executable:

```stata
global JULIA_EXE "C:/path/to/julia.exe"
do replicate_merit.do
```

The DO file exports one CSV per state and the Julia input, calls
`replicate_julia.jl` as an external process, then runs both CSDID cluster
jackknives. CSV exchange avoids the local `jl save` data-transfer failure;
the `jl`, `undidjl`, and `didintjl` wrappers are not required for this route.
The 51-cluster CSDID jackknives can take a substantial amount of time.
The default command recomputes both in full.

To explicitly reuse the completed CSDID calculations from 22 September 2026:

```stata
do replicate_merit.do reuse_csdid
```

This is the mode used for the integration check. It still recomputes Julia
estimates from raw data, checks the exact `merit.dta` SHA-256, and imports
CSDID values from `diagnostics/verified_csdid_20260922.csv`. The original
completed CSDID log is preserved alongside that file. These are computed
values recorded to nine decimals, not four-decimal target constants. The log's
subsequent failure was in the old Julia bridge, after both CSDID calculations
and their results had completed. `output/replication_status.txt` records the
mode and the revised baseline being checked.

The main outputs, under the Git-ignored `output/` directory, are:

- `panel_c_results.csv` and `panel_c_results.tex`: computed estimates.
- `replicate_merit.log` and `julia.log`: execution logs.
- `replication_status.txt`: agreement with the revised paper at four decimals.
- `didint_inference_comparison.csv`: historical subgroup jackknife, HC3, and
  current state-deletion jackknife shown separately.
- `undid_specifications.csv`: all four weighting choices, with and without
  covariates, for both aggregations.
- `undid_stages/`: freshly generated silo differences.
- `julia_versions.txt`: runtime and loaded package versions.

The script writes the results and exits with Stata `r(459)` if any entry differs
from the revised paper baseline. The baseline is a regression check of the
adopted results, not independent evidence that the earlier draft was correct.
The log separately displays differences from the earlier draft. Targets are
never used to construct estimates.

Randomization inference is not validated or reported by this replication;
the Julia packages receive two permutations only because their API requires a
positive count. RI columns in detailed subgroup exports are diagnostics, not
usable inferential results.

## Independent checks

After the Stata run, use Python 3 (no third-party packages):

```sh
python diagnostics/audit_results.py
```

The audit checks all 1,394 UN-DID estimation differences against raw state-year
means, independently reconstructs the 16 UN-DID specification results,
verifies the historical DID-INT jackknife formula, checks CSDID against its
preserved calculation log, and audits the archived UN-DID summaries. It also
checks agreement with the committed results CSV and reports the earlier draft
comparison separately. It writes `output/replication_audit.json` and
`output/historical_undid_component_comparison.csv`.

## Appendix: CCC variations

Run the five DID-INT specifications in the appendix independently of Panel C:

```stata
do replicate_appendix.do
```

This uses the same pinned Julia environment and raw `merit.dta`. It exports a
temporary CSV from Stata, then calls `replicate_appendix.jl` through an external
Julia process; the Stata `didintjl` wrapper is not needed. If Julia is not on
PATH, set `global JULIA_EXE "C:/path/to/julia.exe"` first. The script generates
`output/appendix_ccc_results.csv`, `output/appendix_ccc_table.tex`, and logs.
The verified numerical snapshot is in [results/appendix_ccc_results.csv](results/appendix_ccc_results.csv),
and the four-decimal LaTeX table is in [paper/appendix_ccc_table.tex](paper/appendix_ccc_table.tex).

All five models use cohort aggregation, `both` weighting, and the covariates
`asian`, `male`, and `black`. The table's “Region” setting is DiDInt.jl's
`ccc="state"`; the other settings are `hom`, `time`, `add` (Two One-Way), and
`int` (Two Way). The displayed SE is the historical final-stage
delete-one-cohort jackknife, not the current package's state-deletion
jackknife and not a state-cluster SE. The script checks every ATT and SE
against the supplied appendix table at four decimals before writing results.
Randomization-inference p-values are not part of this check.

## Data provenance and scope

`merit.dta` contains 42,161 observations, 51 states including DC, and years
1989-2000. Its SHA-256 is
`1509b32bf680bf34783c5f27d58027e67931c85eead8f58c235c004b8887abdc`.

`diagnostics/historical/` contains copies of two old UN-DID summary CSVs from
`Unpool-RA/examples/merit`, retained as evidence. They are not inputs to the
replication and are not treated as verified estimates from the raw data.
The repository covers Panel C and the five appendix CCC specifications above;
it does not re-estimate Panels A and B.
