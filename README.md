# MERIT scholarship replication

This repository reproduces Panel C of the MERIT scholarship example with four
Stata-facing estimators: `undidjl`, `csdid`, `csdidjack`, and `didintjl`.

## Published target

| Aggregation | UN-DID ATT | UN-DID SE | CSDID ATT | CSDID SE | DID-INT ATT | DID-INT SE |
|---|---:|---:|---:|---:|---:|---:|
| simple | 0.0485 | 0.0110 | 0.0464 | 0.0133 | 0.0464 | 0.0102 |
| group | 0.0459 | 0.0188 | 0.0339 | 0.0211 | 0.0458 | 0.0084 |

## Files

- `merit.dta`: CPS microdata, 42,161 observations, 51 states (including DC),
  1989–2000.
- `replicate_merit.do`: complete replication.

The DO file creates `output/silos/state_XX.csv` for every state before running
UN-DID. Generated silo files, intermediate UN-DID files, logs, and final table
files stay under `output/`, which is ignored by Git.

## Requirements

- Stata 16 or later (the current `didintjl` wrapper requires Stata 16).
- Julia and Juliaup available from the command line.
- Stata packages `drdid`, `csdid`, `csdidjack`, `didintjl`, and `undidjl`.
- Julia packages `DiDInt.jl` and `Undid.jl` in the Julia environment used by
  Stata's `jl` bridge.

The replication was prepared against these installed command versions:
`csdid` 1.81, `csdidjack` 0.5.2, `didintjl` 0.7.6, and
`undidjl_stage_three` 0.6.4. Package updates can change numerical results, so
the DO file prints the computed table, the published target, and their
difference.

If the Stata packages are missing, run:

```stata
ssc install drdid, replace
ssc install csdid, replace
net install csdidjack, from("https://raw.githubusercontent.com/liu-yunhan/csdidjack/main/") replace
net install didintjl, from("https://raw.githubusercontent.com/ebjamieson97/didintjl/main/") replace
net install undidjl, from("https://raw.githubusercontent.com/ebjamieson97/undidjl/main/") replace
updateundid
```

## Run

Change Stata's working directory to the repository root and run:

```stata
do replicate_merit.do
```

The full run is computationally intensive because `csdidjack` re-estimates
CSDID after omitting each of 51 state clusters. The final outputs are:

- `output/panel_c_results.csv`
- `output/panel_c_results.tex`
- `output/replicate_merit.log`

## Verification status

The copied `merit.dta` matches the source file byte for byte (SHA-256
`1509b32bf680bf34783c5f27d58027e67931c85eead8f58c235c004b8887abdc`).
The current installed CSDID and `csdidjack` commands reproduce all four
published CSDID entries at four decimals. Direct execution of DiDInt.jl with
the historical treated-state mapping also reproduces both DID-INT ATTs. With
the current DiDInt.jl package and HC3 inference, its standard errors are
`0.0103` (simple) and `0.0091` (cohort), versus the published `0.0102` and
`0.0084`. The source of this difference remains unresolved; those two
reported SEs have not been independently reproduced.

On the machine used to prepare the repository, Stata's `jl` 2.0 bridge starts
Julia successfully. The DO file explicitly loads the bundled Stata interface,
which was missing from Julia's load path during batch testing. Even with that
module loaded, `jl save df` returns `r(999)` for a two-row toy dataset in
batch mode. This prevents an end-to-end run of `undidjl` or `didintjl` from
being claimed as verified here. The DO file reports every computed value and
flags deviations from the published table.
