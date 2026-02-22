include("model.jl")

using .My
using ProjectRoot
using YAML
using Optim
using CSV
using DataFrames

dir = @projectroot("output")

# Target moments from Greenwood & Guner (2009), Table 3
const TARGETS = (
    share_married_1950 = 0.816,
    share_married_2000 = 0.625,
    prob_divorce_1950 = 0.011,
    prob_divorce_2000 = 0.023,
    prob_marriage_1950 = 0.211,
    prob_marriage_2000 = 0.082,
)

# =============================================================================
# Save DT baseline results
# =============================================================================
m_dt = My.Model()
df_dt = My.simulate(; m = m_dt, years = [1950, 2000], method = :dt)

# =============================================================================
# CT re-estimation: (μₘ, σₘ, η) via Nelder-Mead
# =============================================================================
function loss(x)
    μₘ, log_σₘ, log_η = x
    σₘ = exp(log_σₘ)
    η = exp(log_η)

    m = My.Model(; μₘ, σₘ, η)

    try
        df = My.simulate(; m, years = [1950, 2000], method = :ct)
        r1950 = df[df.t.==1950, :][1, :]
        r2000 = df[df.t.==2000, :][1, :]

        dev = [
            (1 - r1950.s - TARGETS.share_married_1950) / TARGETS.share_married_1950,
            (1 - r2000.s - TARGETS.share_married_2000) / TARGETS.share_married_2000,
            (r1950.pd - TARGETS.prob_divorce_1950) / TARGETS.prob_divorce_1950,
            (r2000.pd - TARGETS.prob_divorce_2000) / TARGETS.prob_divorce_2000,
            (r1950.pm - TARGETS.prob_marriage_1950) / TARGETS.prob_marriage_1950,
            (r2000.pm - TARGETS.prob_marriage_2000) / TARGETS.prob_marriage_2000,
        ]

        return sum(dev .^ 2)
    catch e
        @warn "Solver failed" x e
        return 1e6
    end
end

# --- Initial guess: DT values (naive mapping) --------------------------------
x0 = [m_dt.μₘ, log(m_dt.σₘ), log(m_dt.η)]

@info "Starting CT parameter estimation..."
@info "Initial guess" μₘ = m_dt.μₘ σₘ = m_dt.σₘ η = m_dt.η
@info "Initial objective" f = loss(x0)

result = optimize(
    loss,
    x0,
    NelderMead(),
    Optim.Options(;
        iterations = 5000,
        g_tol = 1e-10,
        show_trace = true,
        show_every = 50,
    ),
)

# --- Extract results ----------------------------------------------------------
x_opt = Optim.minimizer(result)
μₘ_opt = x_opt[1]
σₘ_opt = exp(x_opt[2])
η_opt = exp(x_opt[3])

@info "Estimation complete"
@info "Optimal parameters" μₘ = μₘ_opt σₘ = σₘ_opt η = η_opt
@info "Minimum objective" f = Optim.minimum(result)
@info "Converged" converged = Optim.converged(result)

# --- Simulate at optimal and compare with targets -----------------------------
m_ct = My.Model(; μₘ = μₘ_opt, σₘ = σₘ_opt, η = η_opt)
df_ct = My.simulate(; m = m_ct, years = [1950, 2000], method = :ct, n_b = 200)

r1950 = df_ct[df_ct.t.==1950, :][1, :]
r2000 = df_ct[df_ct.t.==2000, :][1, :]

comparison = DataFrame(
    moment = [
        "Frac. married 1950", "Frac. married 2000",
        "Prob. divorce 1950", "Prob. divorce 2000",
        "Prob. marriage 1950", "Prob. marriage 2000",
    ],
    data = [
        TARGETS.share_married_1950, TARGETS.share_married_2000,
        TARGETS.prob_divorce_1950, TARGETS.prob_divorce_2000,
        TARGETS.prob_marriage_1950, TARGETS.prob_marriage_2000,
    ],
    model_ct = [
        1 - r1950.s, 1 - r2000.s,
        r1950.pd, r2000.pd,
        r1950.pm, r2000.pm,
    ],
)

# DT baseline for comparison
r1950_dt = df_dt[df_dt.t.==1950, :][1, :]
r2000_dt = df_dt[df_dt.t.==2000, :][1, :]
comparison.model_dt = [
    1 - r1950_dt.s, 1 - r2000_dt.s,
    r1950_dt.pd, r2000_dt.pd,
    r1950_dt.pm, r2000_dt.pm,
]

println("\n=== Moment Comparison ===")
println(comparison)

# --- Save results -------------------------------------------------------------
# CT estimated parameters
d_ct = Dict(
    "mu_m" => μₘ_opt,
    "sigma_m" => σₘ_opt,
    "eta" => η_opt,
    "objective" => Optim.minimum(result),
    "converged" => Optim.converged(result),
)
YAML.write_file("$dir/ct_estimated.yaml", d_ct)

CSV.write("$dir/moment.csv", comparison)
@info "Results saved to $dir/ct_estimated.yaml and $dir/moment.csv"

# =============================================================================
# Full simulation 1950–2000 with DT and estimated CT, saved as CSV
# =============================================================================
years_full = 1950:2000

df_dt_full = My.simulate(; m = m_dt, years = years_full, method = :dt)
df_ct_full = My.simulate(; m = m_ct, years = years_full, method = :ct)

df_sim = vcat(df_dt_full, df_ct_full)
CSV.write("$dir/simulation.csv", df_sim)
@info "Full simulation saved to $dir/simulation.csv"
