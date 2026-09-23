# Called by replicate_merit.do. No Stata/Julia in-process data bridge is needed.
using DataFrames, DiDInt, Undid, DelimitedFiles, LinearAlgebra, SHA

root = @__DIR__
out = joinpath(root, "output")
mkpath(out)
@assert bytes2hex(sha256(read(joinpath(root, "merit.dta")))) ==
    "1509b32bf680bf34783c5f27d58027e67931c85eead8f58c235c004b8887abdc"

function write_csv(path, df)
    # All output fields here are numeric or simple identifiers (no commas).
    open(path, "w") do io
        println(io, join(names(df), ','))
        for row in eachrow(df)
            println(io, join(row, ','))
        end
    end
end

"""Historical final-stage jackknife: delete one subgroup and renormalize weights.

This reproduces DiDInt.jl's former compute_jknife_se on the weighted
intercept-only final regression. It does NOT delete a state from the data.
The historical formula centers deletions at the full-sample estimate.
"""
function subgroup_jackknife(y, weights)
    w = weights ./ sum(weights)
    theta = sum(w .* y)
    deleted = [(sum(w[j] * y[j] for j in eachindex(y) if j != i) /
                sum(w[j] for j in eachindex(y) if j != i)) for i in eachindex(y)]
    se = sqrt((length(y)-1)/length(y) * sum((deleted .- theta).^2))
    # Independent regression calculation, matching the historical package code.
    x = reshape(sqrt.(w), :, 1)
    yw = sqrt.(w) .* y
    regressions = [(x[setdiff(eachindex(y), [i]), :] \
                    yw[setdiff(eachindex(y), [i])])[1] for i in eachindex(y)]
    @assert isapprox(deleted, regressions; atol=1e-12, rtol=1e-10)
    return theta, se, deleted
end

raw, header = readdlm(joinpath(out, "merit_input.csv"), ',', header=true)
df = DataFrame(raw, vec(String.(header)))
@assert nrow(df) == 42161
df.state = string.(Int.(df.state))
df.year = string.(Int.(df.year))
states = sort(unique(df.state))
@assert length(states) == 51
times = [maximum(Int.(df.gvar[df.state .== s])) for s in states]
treated = findall(!=(0), times)
@assert length(treated) == 10

open(joinpath(out, "julia_versions.txt"), "w") do io
    println(io, "Julia=", VERSION)
    println(io, "DiDInt=", pkgversion(DiDInt), " path=", pathof(DiDInt))
    println(io, "Undid=", pkgversion(Undid), " path=", pathof(Undid))
end

# The wrapper's default is both, whereas the Julia function's default is att.
# Set it explicitly at both stages to preserve the original DO specification.
silodir = joinpath(out, "undid_stages")
mkpath(silodir)
cd(silodir) do
    create_init_csv(states, ["1989"], ["2000"],
                    [g == 0 ? "control" : string(g) for g in times])
    create_diff_df(joinpath(silodir, "init.csv"), "yyyy", "yearly";
                   covariates=["asian", "black", "male"], weights="both")
    for s in states
        println("UN-DID silo ", s); flush(stdout)
        undid_stage_two(joinpath(silodir, "empty_diff_df.csv"), s,
                        df[df.state .== s, :], "year", "coll", "yyyy")
    end
end

results = DataFrame(aggregation=["simple", "group"],
    undid_att=zeros(2), undid_se=zeros(2), didint_att=zeros(2), didint_se=zeros(2))
variants = DataFrame(aggregation=String[], covariates=Bool[], weighting=String[],
                     att=Float64[], subgroup_jackknife_se=Float64[])
for (i, agg) in enumerate(["gt", "g"]), cov in [false, true], weighting in ["none", "att", "diff", "both"]
    r = undid_stage_three(silodir; agg=agg, covariates=cov, weighting=weighting,
                          nperm=2, seed=1234, verbose=false)
    push!(variants, (results.aggregation[i], cov, weighting, r.agg_att[1], r.jackknife_se[1]))
    if !cov && weighting == "both"
        results.undid_att[i] = r.agg_att[1]
        results.undid_se[i] = r.jackknife_se[1]
        write_csv(joinpath(out, "undid_$(results.aggregation[i])_subgroups.csv"), r)
    end
end
write_csv(joinpath(out, "undid_specifications.csv"), variants)

inference = DataFrame(aggregation=String[], att=Float64[], historical_subgroup_jackknife_se=Float64[],
    hc3_se=Float64[], current_state_jackknife_se=Float64[], n_subgroups=Int[])
for (i, agg) in enumerate(["simple", "cohort"])
    r = didint("coll", "state", "year", df;
        treated_states=states[treated], treatment_times=string.(times[treated]),
        date_format="yyyy", covariates=["asian", "black", "male"], ccc="state",
        agg=agg, weighting="both", hc=3, nperm=2, seed=1234, verbose=false)
    y = Float64.(r[!, agg == "simple" ? :att_gt : :att_cohort])
    theta, se, deleted = subgroup_jackknife(y, Float64.(r.weights))
    @assert isapprox(theta, r.agg_att[1]; atol=1e-12)
    # For a weighted intercept-only regression the historical SE also equals
    # sqrt((m-1)/m) times HC3. Check this independently of the deletion loop.
    @assert isapprox(se, sqrt((length(y)-1)/length(y))*r.se_agg_att[1]; atol=1e-12)
    results.didint_att[i] = theta
    results.didint_se[i] = se
    push!(inference, (results.aggregation[i], theta, se, r.se_agg_att[1],
                      r.jknifese_agg_att[1], length(y)))
    write_csv(joinpath(out, "didint_$(results.aggregation[i])_subgroups.csv"), r)
    write_csv(joinpath(out, "didint_$(results.aggregation[i])_deletions.csv"),
              DataFrame(deleted_subgroup=1:length(y), att=deleted))
end
write_csv(joinpath(out, "didint_inference_comparison.csv"), inference)
# Write the completion file last; Stata removes any old copy before invoking us.
write_csv(joinpath(out, "julia_results.csv"), results)
show(stdout, MIME("text/plain"), results); println()
