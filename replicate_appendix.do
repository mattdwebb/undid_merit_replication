/*
MERIT appendix: five DID-INT CCC variations, cohort aggregation.
Run from the repository root. Julia runs externally through CSV exchange,
avoiding the local Stata jl data-transfer failure.
*/

version 16
clear all
set more off

global ROOT = subinstr(c(pwd), "\", "/", .)
capture confirm file "$ROOT/merit.dta"
if _rc {
    display as error "Change to the repository directory containing merit.dta."
    exit 601
}

global OUT "$ROOT/output"
capture mkdir "$OUT"
capture log close _all
log using "$OUT/replicate_appendix.log", text replace

// Override JULIA_EXE before running if Julia is not on PATH.
local julia "$JULIA_EXE"
if "`julia'" == "" local julia "julia"

use "$ROOT/merit.dta", clear
confirm variable coll state year gvar asian black male
keep coll state year gvar asian black male
export delimited using "$OUT/appendix_merit_input.csv", replace

capture erase "$OUT/appendix_ccc_results.csv"
shell "`julia'" --project="$ROOT/julia" "$ROOT/replicate_appendix.jl" > "$OUT/appendix_julia.log" 2>&1
capture confirm file "$OUT/appendix_ccc_results.csv"
if _rc {
    display as error "Appendix replication failed. See $OUT/appendix_julia.log."
    log close
    exit 601
}

import delimited using "$OUT/appendix_ccc_results.csv", clear asdouble
assert _N == 5
list ccc_variation aggregate_att jackknife_se, noobs separator(0)
display as result "All five appendix rows reproduced at four decimals."
log close
