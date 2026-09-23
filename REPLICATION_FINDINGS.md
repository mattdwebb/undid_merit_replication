# Replication findings - 22 September 2026

## Paper revision

The revised paper adopts the freshly computed UN-DID ATT/SE pairs
0.0466/0.0113 (simple) and 0.0458/0.0133 (group). CSDID and DID-INT retain
their verified values. See the [current table](README.md#revised-results---22-september-2026)
and [replacement LaTeX](paper/merit_results.tex).

The pipeline now checks the revised paper baseline and preserves comparison
with the earlier draft separately. The findings below document the earlier
draft investigation; its four UN-DID discrepancies remain historical facts,
not failures to reproduce the newly adopted table. The new baseline is not
an independent validation of those older values.

The original Stata-Julia transfer error has been bypassed with an external
Julia process and CSV exchange. The investigation computed the estimates from
`merit.dta` and explicitly checked them against the earlier-draft table.

## What reproduces

CSDID's simple and group ATTs and its state-cluster jackknife standard errors
match all four earlier-draft values. The numerical values from the completed run recorded earlier on
22 September 2026 are:

| Aggregation | ATT | State-cluster jackknife SE |
|---|---:|---:|
| simple | 0.046375882 | 0.013347467 |
| group | 0.033862650 | 0.021056207 |

DID-INT's point estimates and the **historical final-stage subgroup jackknife**
also reproduce all four earlier-draft values:

| Aggregation | ATT | Historical subgroup JK | HC3 | Current state-deletion JK |
|---|---:|---:|---:|---:|
| simple | 0.046406233 | 0.010160542 | 0.010318080 | 0.012290995 |
| group | 0.045817214 | 0.008394653 | 0.009067262 | 0.012456524 |

The original DO file incorrectly tried to reproduce the earlier-draft DID-INT
standard errors using HC3. The historical calculation deletes one of the 33
cohort-time estimates (simple), or one of the seven cohort estimates (group),
then recomputes their weighted average. It centers the deletions at the
full-sample estimate. With normalized weights w and subgroup estimates a:

    theta = sum(w[i] * a[i])
    theta_minus_i = (theta - w[i] * a[i]) / (1 - w[i])
    SE = sqrt((m-1)/m * sum((theta_minus_i - theta)^2))

For the final weighted intercept-only regression this equals
`sqrt((m-1)/m) * HC3`. Both the explicit deletion loop and a separate weighted
least-squares implementation verify that identity. Python independently
checks the results. Executing the locally installed DiDInt.jl 0.5.7 source
also produced ATTs 0.04640623337906369 / 0.04581721370543262 and jackknife
SEs 0.010160542435611488 / 0.008394653230371286, independently confirming
the historical-package match (see `output/historical_didint.log`).
This historical calculation appears in DiDInt.jl's
`compute_jknife_se` and `final_regression_results` before its inference
reorganization; for example, see
[historical helpers.jl](https://github.com/ebjamieson97/DiDInt.jl/blob/e81ae94/src/helpers.jl).

**The historical subgroup jackknife is not state-cluster jackknife inference.**
Reproducing those values does not validate the earlier manuscript's claim that all
reported standard errors cluster by state. The current state-deletion values
are provided separately for a substantive decision about inference.

## What remains unresolved

No tested UN-DID specification reproduces the four earlier-draft UN-DID entries.
The original DO file's explicit no-covariate aggregation, with the Stata
wrapper's default `both` weighting made explicit, produces:

| Aggregation | Computed ATT | Published ATT | Computed subgroup JK SE | Published SE |
|---|---:|---:|---:|---:|
| simple | 0.046564232 | 0.0485 | 0.011313697 | 0.0110 |
| group | 0.045822523 | 0.0459 | 0.013331859 | 0.0188 |

Undid.jl's direct Julia default is `att`, while the Stata wrapper defaults to
`both`. The repaired code sets `both` at initialization and aggregation so the
choice is reproducible. All four supported weighting options (`none`, `att`,
`diff`, `both`), both covariate-adjustment choices, and both aggregations are
saved in `output/undid_specifications.csv`. None produces the earlier-draft UN-DID
pair of ATT and SE.

The earlier manuscript said that Panel C includes covariates. The historical DO file
instead called stage three with its default `covariates=false`; the original
repository explicitly preserved that setting. Turning it on with `both`
weighting gives ATTs 0.042989949 and 0.042313575, not the earlier-draft values.
These alternatives are diagnostics, not a recommendation to select a
specification just because it comes closer to a target.

Running the older installed Undid.jl 0.2.7 source on the historical silo
differences also gives unweighted, no-covariate point estimates 0.048171441
and 0.038700404, matching the corresponding fresh calculations rather than
the archived targets. Thus that version downgrade alone does not recover
the earlier-draft point estimates (`output/historical_undid.log`).

## Historical UN-DID files explain the targets, but do not verify them

The two `diagnostics/historical/UNDID_results_*_new.csv` files were copied from
`Unpool-RA/examples/merit`. They contain aggregate ATTs 0.048487265 and
0.045868055, which round to the earlier-draft 0.0485 and 0.0459.

However, the saved component estimates do not all agree with fresh estimates
from the raw data. For example, the archived 1991/1991 estimate is 0.059869261;
the fresh unweighted cohort-time estimate is 0.112702840. The archived number
instead equals the fresh *entire 1991 cohort* average. Seven cohort-time
entries differ; the other 26 agree to the stored precision. This pattern
identifies where to investigate the historical export/editing process; it
does not establish how those changes arose.

The archived standard errors also do not equal a delete-one-subgroup
jackknife of the components saved alongside them:

| Aggregation | Saved SE | JK from saved components | Saved SE times sqrt(50) |
|---|---:|---:|---:|
| simple | 0.001557634 | 0.008811307 | 0.011014136 |
| group | 0.002659423 | 0.006514230 | 0.018804960 |

The last column reproduces the earlier-draft rounded SEs. The old DO-file notes
mention correcting the SE scale. This explains an arithmetic route to the
earlier-draft numbers, but is **not** a reproduction of them through a validated
state-deletion calculation. The original leave-one-out estimates or code
that generated those historical summary files would be needed to resolve
that provenance.

## Validation and scope

The audit checks the data checksum, 42,161 observations, 51 states, and ten
treated states. It verifies all 1,394 raw UN-DID difference cells and their
sample counts, reconstructs 16 UN-DID specification/aggregation combinations,
and independently checks the DID-INT deletion calculation. The audit compares
computed outputs with the revised paper baseline and separately reports the
four historical discrepancies against the earlier draft.

The final integration test uses `do replicate_merit.do reuse_csdid`: Julia
estimates are freshly computed; the completed CSDID 51-state jackknives are
read from the preserved same-day calculation log. Both Julia and the Python
audit check the exact source-data SHA-256. Python also checks the reused
CSDID values against the original logged results. A redundant full CSDID
recheck was stopped after a timed single fit took 33.496 seconds, implying
roughly an hour for 102 repeated fits. Its partial log is retained as
`output/full_recheck_partial.log`; it is not presented as a completed run.
Full CSDID recomputation remains the default when no argument is supplied.

Julia packages are pinned by the repository manifest; the runs use Julia
1.11.7, DiDInt.jl 0.7.5, and Undid.jl 0.5.2. The source data were not modified,
and no results were substituted from the archived tables. The paper revision
adopts the computed UN-DID values and corrects the descriptions of covariates
and inference.
