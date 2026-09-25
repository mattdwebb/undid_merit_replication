# Called by replicate_appendix.do; uses the repository's pinned Julia project.
using DataFrames, DiDInt, DelimitedFiles, LinearAlgebra, Printf, SHA

root = @__DIR__
out = joinpath(root, "output")
mkpath(out)
@assert bytes2hex(sha256(read(joinpath(root, "merit.dta")))) ==
    "1509b32bf680bf34783c5f27d58027e67931c85eead8f58c235c004b8887abdc"

raw, header = readdlm(joinpath(out, "appendix_merit_input.csv"), ',', header=true)
df = DataFrame(raw, vec(String.(header)))
@assert nrow(df) == 42161
df.state = string.(Int.(df.state))
df.year = string.(Int.(df.year))
states = sort(unique(df.state))
@assert length(states) == 51
times = [maximum(Int.(df.gvar[df.state .== s])) for s in states]
treated = findall(!=(0), times)
@assert length(treated) == 10

# Final-stage jackknife: delete one cohort estimate, renormalize
# its weight, and center deletions on the full-sample weighted ATT. This is
# not a state-cluster jackknife and differs from current DiDInt.jl's jkse.
function subgroup_jackknife(y, weights)
    w = weights ./ sum(weights)
    theta = sum(w .* y)
    deleted = [(theta - w[i] * y[i]) / (1 - w[i]) for i in eachindex(y)]
    se = sqrt((length(y) - 1) / length(y) * sum((deleted .- theta).^2))
    return theta, se
end

specs = [
    ("Homogeneous", "hom", 0.0450, 0.0092),
    ("Time", "time", 0.0412, 0.0101),
    ("Region", "state", 0.0458, 0.0084),
    ("Two One-Way", "add", 0.0420, 0.0091),
    ("Two Way", "int", 0.0511, 0.0209),
]
results = DataFrame(ccc_variation=String[], ccc_option=String[],
    aggregate_att=Float64[], jackknife_se=Float64[], n_cohorts=Int[])

for (label, ccc, expected_att, expected_se) in specs
    println("Estimating ", label, " (ccc=", ccc, ")"); flush(stdout)
    r = didint("coll", "state", "year", df;
        treated_states=states[treated], treatment_times=string.(times[treated]),
        date_format="yyyy", covariates=["asian", "male", "black"],
        ccc=ccc, agg="cohort", weighting="both", hc=3,
        nperm=2, seed=1234, verbose=false)
    y = Float64.(r.att_cohort)
    att, se = subgroup_jackknife(y, Float64.(r.weights))
    @assert length(y) == 7
    @assert isapprox(att, r.agg_att[1]; atol=1e-12)
    @assert isapprox(se, sqrt(6 / 7) * r.se_agg_att[1]; atol=1e-12)
    @printf("  ATT %.9f, subgroup JK SE %.9f\n", att, se)
    @assert round(att; digits=4) == expected_att "ATT differs from appendix target for $label"
    @assert round(se; digits=4) == expected_se "SE differs from appendix target for $label"
    push!(results, (label, ccc, att, se, length(y)))
end

# Write the completion CSV last, after every model and the LaTeX output pass.
open(joinpath(out, "appendix_ccc_table.tex"), "w") do io
    rowend = repeat("\\", 2)
    println(io, raw"\begin{tabular}{lcc}")
    println(io, raw"\hline")
    println(io, "CCC variation & Aggregate ATT & Jackknife SE ", rowend)
    println(io, raw"\hline")
    for row in eachrow(results)
        @printf(io, "%s & %.4f & %.4f %s\n", row.ccc_variation,
            row.aggregate_att, row.jackknife_se, rowend)
    end
    println(io, raw"\hline\hline")
    println(io, raw"\end{tabular}")
end
open(joinpath(out, "appendix_ccc_results.csv"), "w") do io
    println(io, "ccc_variation,ccc_option,aggregate_att,jackknife_se,n_cohorts")
    for row in eachrow(results)
        @printf(io, "%s,%s,%.9f,%.9f,%d\n", row.ccc_variation,
            row.ccc_option, row.aggregate_att, row.jackknife_se, row.n_cohorts)
    end
end
println("All five appendix rows match the supplied table at four decimals.")
