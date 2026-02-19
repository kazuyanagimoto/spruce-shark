include("model.jl")

using .My
using ProjectRoot
using BenchmarkTools
using CSV
using DataFrames

dir = @projectroot("output")
m = My.Model()

# --- Grid sizes to benchmark -------------------------------------------------
grid_sizes = [25, 50, 100, 200, 400, 800]
n_std = 4.0
years = 1950:2000

results = DataFrame(
    n_b = Int[],
    method = String[],
    time_median_ms = Float64[],
    time_mean_ms = Float64[],
    time_min_ms = Float64[],
    time_max_ms = Float64[],
    allocs = Int[],
    memory_mb = Float64[],
)

for nb in grid_sizes
    # --- DT benchmark ---------------------------------------------------------
    @info "Benchmarking DT  n_b = $nb"
    bm_dt = @benchmark My.simulate(; m = $m, years = $years, method = :dt,
                                     n_b = $nb, n_std = $n_std) seconds = 10
    push!(results, (
        n_b        = nb,
        method     = "DT",
        time_median_ms = median(bm_dt).time / 1e6,
        time_mean_ms   = mean(bm_dt).time   / 1e6,
        time_min_ms    = minimum(bm_dt).time / 1e6,
        time_max_ms    = maximum(bm_dt).time / 1e6,
        allocs     = median(bm_dt).allocs,
        memory_mb  = median(bm_dt).memory / 1024^2,
    ))

    # --- CT benchmark ---------------------------------------------------------
    @info "Benchmarking CT  n_b = $nb"
    bm_ct = @benchmark My.simulate(; m = $m, years = $years, method = :ct,
                                     n_b = $nb, n_std = $n_std) seconds = 10
    push!(results, (
        n_b        = nb,
        method     = "CT",
        time_median_ms = median(bm_ct).time / 1e6,
        time_mean_ms   = mean(bm_ct).time   / 1e6,
        time_min_ms    = minimum(bm_ct).time / 1e6,
        time_max_ms    = maximum(bm_ct).time / 1e6,
        allocs     = median(bm_ct).allocs,
        memory_mb  = median(bm_ct).memory / 1024^2,
    ))
end

# --- Save to CSV --------------------------------------------------------------
CSV.write("$dir/benchmark.csv", results)
@info "Benchmark saved to $dir/benchmark.csv"
println(results)
