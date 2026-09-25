/*
MERIT scholarship replication: Panel C, full sample

Run this file from the repository root, so that merit.dta is in c(pwd).
The script creates one CSV per state for the UN-DID silo stage, estimates
UN-DID, CSDID/CSDIDJACK, and DID-INT, and writes the combined results.
Julia runs externally; DID-INT uses the final-stage subgroup jackknife for
the reported standard error.
*/

version 16
clear all
set more off
set linesize 255
args csdid_mode
if !inlist("`csdid_mode'", "", "reuse_csdid") exit 198

// -----------------------------------------------------------------------------
// Paths and dependency checks
// -----------------------------------------------------------------------------

global ROOT = subinstr(c(pwd), "\", "/", .)
capture confirm file "$ROOT/merit.dta"
if _rc {
    display as error "merit.dta was not found in $ROOT"
    display as error "Change to the repository directory and rerun this DO file."
    exit 601
}

global OUT   "$ROOT/output"
global SILOS "$OUT/silos"
capture mkdir "$OUT"
capture mkdir "$SILOS"

capture log close _all
log using "$OUT/replicate_merit.log", text replace

foreach command in csdid csdidjack {
    capture noisily which `command'
    if _rc {
        display as error "Required command `command' is not installed. See README.md."
        exit 199
    }
}

// Override JULIA_EXE before running if Julia is not on PATH.
local julia "$JULIA_EXE"
if "`julia'" == "" local julia "julia"

// Reference values are a regression check, never inputs to estimation.
matrix target = (0.0466, 0.0113, 0.0464, 0.0133, 0.0464, 0.0102 \ ///
                 0.0458, 0.0133, 0.0339, 0.0211, 0.0458, 0.0084)
matrix rownames target = simple group
matrix colnames target = UNDID_ATT UNDID_SE CSDID_ATT CSDID_SE DIDINT_ATT DIDINT_SE

matrix results = J(2, 6, .)
matrix rownames results = simple group
matrix colnames results = UNDID_ATT UNDID_SE CSDID_ATT CSDID_SE DIDINT_ATT DIDINT_SE

// -----------------------------------------------------------------------------
// 1. Create one physical CSV per state (the artificial UN-DID silos)
// -----------------------------------------------------------------------------

use "$ROOT/merit.dta", clear
confirm variable coll merit male black asian year state gvar
assert inrange(year, 1989, 2000)
assert gvar == 0 | inlist(gvar, 1991, 1993, 1996, 1997, 1998, 1999, 2000)

levelsof state, local(states)
local starts
local ends
local treatments

foreach s of local states {
    quietly summarize year if state == `s', meanonly
    local ymin = r(min)
    local ymax = r(max)
    local starts "`starts' `ymin'"
    local ends   "`ends' `ymax'"

    quietly summarize gvar if state == `s', meanonly
    if r(max) == 0 {
        local treatments "`treatments' control"
    }
    else {
        local gtext : display %9.0f r(max)
        local gtext = trim("`gtext'")
        local treatments "`treatments' `gtext'"
    }

    preserve
        keep if state == `s'
        sort year
        export delimited using "$SILOS/state_`s'.csv", replace
    restore
}

// Run the Julia packages using CSV exchange. This avoids jl save r(999).
use "$ROOT/merit.dta", clear
keep coll state year gvar asian black male
export delimited using "$OUT/merit_input.csv", replace
capture erase "$OUT/julia_results.csv"
shell "`julia'" --project="$ROOT/julia" "$ROOT/replicate_julia.jl" > "$OUT/julia.log" 2>&1
capture confirm file "$OUT/julia_results.csv"
if _rc {
    display as error "Julia replication failed. See $OUT/julia.log."
    log close
    exit 601
}
import delimited using "$OUT/julia_results.csv", clear asdouble
assert _N == 2
assert aggregation[1] == "simple" & aggregation[2] == "group"
forvalues i = 1/2 {
    matrix results[`i',1] = undid_att[`i']
    matrix results[`i',2] = undid_se[`i']
    matrix results[`i',5] = didint_att[`i']
    matrix results[`i',6] = didint_se[`i']
}

// -----------------------------------------------------------------------------
// 3. CSDID point estimates with CSDIDJACK state-cluster jackknife SEs
// -----------------------------------------------------------------------------

if "`csdid_mode'" == "reuse_csdid" {
    // Explicit reuse of a completed calculation, never of the target constants.
    // Julia has verified the exact merit.dta SHA-256 before this branch.
    import delimited using "$ROOT/results/csdid_results.csv", clear asdouble
    assert _N == 2
    assert aggregation[1] == "simple" & aggregation[2] == "group"
    assert data_sha256 == "1509b32bf680bf34783c5f27d58027e67931c85eead8f58c235c004b8887abdc"
    forvalues i = 1/2 {
        matrix results[`i',3] = csdid_att[`i']
        matrix results[`i',4] = csdid_se[`i']
    }
    display as text "Reused saved CSDID calculation from results/csdid_results.csv."
}
else {
    use "$ROOT/merit.dta", clear

    csdid coll male black asian, gvar(gvar) time(year) ///
        agg(simple) reg cluster(state)
    csdidjack
    matrix results[1,3] = r(ATT)
    matrix results[1,4] = r(cv3se)

    csdid coll male black asian, gvar(gvar) time(year) ///
        agg(group) reg cluster(state)
    csdidjack
    matrix results[2,3] = r(ATT)
    matrix results[2,4] = r(cv3se)
}

// -----------------------------------------------------------------------------
// 5. Display, save, and check Panel C
// -----------------------------------------------------------------------------

display as text _newline "Panel C: Full Sample (computed)"
matlist results, format(%9.4f) rowtitle("Agg.")

display as text _newline "Reference values"
matlist target, format(%9.4f) rowtitle("Agg.")

matrix difference = results - target
display as text _newline "Computed minus reference values"
matlist difference, format(%10.6f) rowtitle("Agg.")

file open csv using "$OUT/panel_c_results.csv", write text replace
file write csv "aggregation,undid_att,undid_se,csdid_att,csdid_se,didint_att,didint_se" _n
file write csv "simple," %12.9f (results[1,1]) "," %12.9f (results[1,2]) "," ///
    %12.9f (results[1,3]) "," %12.9f (results[1,4]) "," ///
    %12.9f (results[1,5]) "," %12.9f (results[1,6]) _n
file write csv "group," %12.9f (results[2,1]) "," %12.9f (results[2,2]) "," ///
    %12.9f (results[2,3]) "," %12.9f (results[2,4]) "," ///
    %12.9f (results[2,5]) "," %12.9f (results[2,6]) _n
file close csv

file open tex using "$OUT/panel_c_results.tex", write text replace
file write tex "\begin{tabular}{lrrrrrr}" _n
file write tex "\hline" _n
file write tex " & UN-DID & UN-DID & CSDID & CSDID & DID-INT & DID-INT \\" _n
file write tex "Agg. & ATT & SE & ATT & SE & ATT & SE \\" _n
file write tex "\hline" _n
file write tex "simple & " %6.4f (results[1,1]) " & " %6.4f (results[1,2]) " & " ///
    %6.4f (results[1,3]) " & " %6.4f (results[1,4]) " & " ///
    %6.4f (results[1,5]) " & " %6.4f (results[1,6]) " \\" _n
file write tex "group & " %6.4f (results[2,1]) " & " %6.4f (results[2,2]) " & " ///
    %6.4f (results[2,3]) " & " %6.4f (results[2,4]) " & " ///
    %6.4f (results[2,5]) " & " %6.4f (results[2,6]) " \\" _n
file write tex "\hline" _n
file write tex "\end{tabular}" _n
file close tex

// A value passes when it rounds to the four-decimal reference entry.
local mismatch = 0
forvalues i = 1/2 {
    forvalues j = 1/6 {
        if missing(results[`i',`j']) | abs(round(results[`i',`j'], .0001) - target[`i',`j']) > 1e-10 {
            local mismatch = `mismatch' + 1
        }
    }
}

if `mismatch' == 0 {
    display as result "All 12 entries agree with the reference values at four decimals."
}
else {
    display as error "`mismatch' of 12 entries differ from the reference values at four decimals."
    display as error "See output/undid_specifications.csv for UN-DID alternatives."
}

file open status using "$OUT/replication_status.txt", write text replace
file write status "csdid_mode=`csdid_mode' (empty means full recomputation)" _n
file write status "reference=results/panel_c_results.csv" _n
file write status "mismatches=`mismatch' of 12" _n
file close status
log close
cd "$ROOT"
if `mismatch' > 0 exit 459
