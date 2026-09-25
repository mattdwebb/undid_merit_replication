# MERIT replication

Data, code, and table outputs for the full-sample MERIT comparison and the
five DID-INT CCC specifications in the appendix.

## Files

- `merit.dta`: source data (42,161 observations, 51 states, 1989–2000).
- `replicate_merit.do` and `replicate_julia.jl`: full-sample comparison.
- `replicate_appendix.do` and `replicate_appendix.jl`: CCC specifications.
- `julia/Project.toml` and `julia/Manifest.toml`: pinned Julia environment.
- `results/*.csv` and `results/*.tex`: computed table values and raw LaTeX
  tables. Run outputs are written under the Git-ignored `output/` directory.

## Requirements

- Stata 16 or later with `drdid`, `csdid`, and `csdidjack` installed.
- Julia 1.11.7 on PATH, or a Julia executable specified in Stata as
  `global JULIA_EXE "C:/path/to/julia.exe"`.

From the repository directory, install the pinned Julia dependencies once:

```sh
julia --project=julia -e "using Pkg; Pkg.instantiate()"
```

If needed, install the Stata commands:

```stata
ssc install drdid, replace
ssc install csdid, replace
net install csdidjack, from("https://raw.githubusercontent.com/liu-yunhan/csdidjack/main/") replace
```

## Run

Set Stata's working directory to the repository root, then run:

```stata
do replicate_merit.do
do replicate_appendix.do
```

The full-sample script exports one CSV per state, runs UN-DID and DID-INT in
Julia through CSV exchange, and computes CSDID state-cluster jackknife SEs in
Stata. The CSDID step can take substantial time. To use the saved CSDID
calculation in `results/csdid_results.csv` instead, run
`do replicate_merit.do reuse_csdid`.

UN-DID and DID-INT report final-stage delete-one-subgroup jackknife SEs;
CSDID reports a state-cluster jackknife SE. The appendix uses cohort
aggregation, `both` weighting, and the `asian`, `male`, and `black` covariates.
Its “Region” row is DiDInt.jl's `ccc="state"` setting. Randomization-inference
p-values are not reported.
