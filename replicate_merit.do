/*
MERIT scholarship replication: Panel C, full sample

Run this file from the repository root, so that merit.dta is in c(pwd).
The script creates one CSV per state for the UN-DID silo stage, estimates
UN-DID, CSDID/CSDIDJACK, and DID-INT, and writes the combined results.
*/

version 16
clear all
set more off
set linesize 255

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

foreach command in csdid csdidjack didintjl create_init_csv create_diff_df undidjl_stage_two undidjl_stage_three {
    capture noisily which `command'
    if _rc {
        display as error "Required command `command' is not installed. See README.md."
        exit 199
    }
}

// Installation/update commands (run manually only when needed):
// ssc install drdid, replace
// ssc install csdid, replace
// net install csdidjack, from("https://raw.githubusercontent.com/liu-yunhan/csdidjack/main/") replace
// net install didintjl, from("https://raw.githubusercontent.com/ebjamieson97/didintjl/main/") replace
// net install undidjl, from("https://raw.githubusercontent.com/ebjamieson97/undidjl/main/") replace
// updateundid

// Keep the expected published values separate from the estimates. They are
// used only for the end-of-file reproduction check.
matrix target = (0.0485, 0.0110, 0.0464, 0.0133, 0.0464, 0.0102 \ ///
                 0.0459, 0.0188, 0.0339, 0.0211, 0.0458, 0.0084)
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

// -----------------------------------------------------------------------------
// 2. UN-DID: initialize, estimate within each silo, then aggregate
// -----------------------------------------------------------------------------

cd "$OUT"

create_init_csv, silo_names("`states'") start_times("`starts'") ///
    end_times("`ends'") treatment_times("`treatments'")

create_diff_df, filepath("$OUT/init.csv") date_format("yyyy") ///
    freq("yearly") covariates("asian black male")

foreach s of local states {
    display as text "UN-DID stage two: state `s'"
    import delimited using "$SILOS/state_`s'.csv", clear case(preserve)
    tostring year, replace format(%9.0f) force
    undidjl_stage_two, filepath("$OUT/empty_diff_df.csv") ///
        local_silo_name("`s'") time_column("year") ///
        outcome_column("coll") local_date_format("yyyy")
}

// These two aggregations correspond to the published simple and group rows.
// The attached historical replication used the default (unadjusted) stage-
// three column, so covariates(false) is explicit here to make that choice stable.
undidjl_stage_three, folder("$OUT/") agg("gt") ///
    covariates("false") nperm(1) seed(1234)
matrix results[1,1] = r(att)
matrix results[1,2] = r(jkse)

undidjl_stage_three, folder("$OUT/") agg("g") ///
    covariates("false") nperm(1) seed(1234)
matrix results[2,1] = r(att)
matrix results[2,2] = r(jkse)

// -----------------------------------------------------------------------------
// 3. CSDID point estimates with CSDIDJACK state-cluster jackknife SEs
// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------
// 4. DID-INT with state-varying CCC and jackknife inference
// -----------------------------------------------------------------------------

use "$ROOT/merit.dta", clear

didintjl, outcome(coll) state(state) time(year) gvar(gvar) ///
    covariates("asian black male") ccc("state") agg("simple") ///
    weighting("both") nperm(1) seed(1234)
matrix results[1,5] = r(att)
matrix results[1,6] = r(jkse)

didintjl, outcome(coll) state(state) time(year) gvar(gvar) ///
    covariates("asian black male") ccc("state") agg("cohort") ///
    weighting("both") nperm(1) seed(1234)
matrix results[2,5] = r(att)
matrix results[2,6] = r(jkse)

// -----------------------------------------------------------------------------
// 5. Display, save, and check Panel C
// -----------------------------------------------------------------------------

display as text _newline "Panel C: Full Sample (computed)"
matlist results, format(%9.4f) rowtitle("Agg.")

display as text _newline "Published target"
matlist target, format(%9.4f) rowtitle("Agg.")

matrix difference = results - target
display as text _newline "Computed minus published target"
matlist difference, format(%10.6f) rowtitle("Agg.")

file open csv using "$OUT/panel_c_results.csv", write text replace
file write csv "aggregation,undid_att,undid_se,csdid_att,csdid_se,didint_att,didint_se" _n
file write csv "simple," %9.6f (results[1,1]) "," %9.6f (results[1,2]) "," ///
    %9.6f (results[1,3]) "," %9.6f (results[1,4]) "," ///
    %9.6f (results[1,5]) "," %9.6f (results[1,6]) _n
file write csv "group," %9.6f (results[2,1]) "," %9.6f (results[2,2]) "," ///
    %9.6f (results[2,3]) "," %9.6f (results[2,4]) "," ///
    %9.6f (results[2,5]) "," %9.6f (results[2,6]) _n
file close csv

file open tex using "$OUT/panel_c_results.tex", write text replace
file write tex "\\begin{tabular}{lrrrrrr}" _n
file write tex "\\hline" _n
file write tex " & UN-DID & UN-DID & CSDID & CSDID & DID-INT & DID-INT \\\\" _n
file write tex "Agg. & ATT & SE & ATT & SE & ATT & SE \\\\" _n
file write tex "\\hline" _n
file write tex "simple & " %6.4f (results[1,1]) " & " %6.4f (results[1,2]) " & " ///
    %6.4f (results[1,3]) " & " %6.4f (results[1,4]) " & " ///
    %6.4f (results[1,5]) " & " %6.4f (results[1,6]) " \\\\" _n
file write tex "group & " %6.4f (results[2,1]) " & " %6.4f (results[2,2]) " & " ///
    %6.4f (results[2,3]) " & " %6.4f (results[2,4]) " & " ///
    %6.4f (results[2,5]) " & " %6.4f (results[2,6]) " \\\\" _n
file write tex "\\hline" _n
file write tex "\\end{tabular}" _n
file close tex

// A value passes when it rounds to the published four-decimal entry.
local mismatch = 0
forvalues i = 1/2 {
    forvalues j = 1/6 {
        if round(results[`i',`j'], .0001) != target[`i',`j'] {
            local mismatch = `mismatch' + 1
        }
    }
}

if `mismatch' == 0 {
    display as result "All 12 estimates reproduce the published table at four decimals."
}
else {
    display as error "`mismatch' of 12 entries differ from the published table at four decimals."
    display as error "Check the package versions recorded near the top of the log."
}

log close
cd "$ROOT"
