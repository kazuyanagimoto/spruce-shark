module My

using Roots
using Distributions
using DataFrames
using QuantEcon
using LinearAlgebra

struct Single{TF<:AbstractFloat}
    c::TF
    h::TF
    l::TF
    d::TF
    n::TF
    v::TF
end

struct Married{TF<:AbstractFloat}
    c::TF
    h::TF
    l::TF
    d::TF
    n::TF
    v::TF
end

@kwdef struct Model{TF<:AbstractFloat,TI<:Integer}
    ## Tastes
    β̃::TF = 0.96
    δ::TF = 1 / 47
    β::TF = β̃ * (1 - δ)
    ϕ::TF = log2(1 + 0.7)
    α::TF = 0.278
    ζ::TF = -1.901
    c̄::TF = 0.131
    θ::TF = 0.206
    κ::TF = 0.189

    ## Prices
    p_1950::TF = 9.959
    γ::TF = 0.059
    w_1950::TF = 1.0
    Δw::TF = 0.022

    ## Match Quality
    ϱ::TF = 0.896
    μₘ::TF = 0.521
    σₘ::TF = sqrt(0.680)
    μₛ::TF = -4.252
    σₛ::TF = sqrt(8.063)

    # Continuous-time
    Δt::TF = 1e-3
    λ::TF = 1 / Δt
    ρ̃::TF = -log(β) / Δt
    ν::TF = -log(1 - δ) / Δt
    ρ::TF = ρ̃ + ν
    η::TF = -log(ϱ)

    # Grid
    n_b::TI = 100
    mc::MarkovChain = tauchen(n_b, ϱ, sqrt(1 - ϱ^2) * σₘ, (1 - ϱ) * μₘ, 4.)
    b_grid::Vector{TF} = collect(mc.state_values)
    G::Matrix{TF} = mc.p
    F::Vector{TF} = construct_F(μₛ, σₛ, b_grid)

end

crra(x, ζ) = iszero(ζ) ? log(x) : x^ζ / ζ
u(c, n; m, z) = m.α * log((c - m.c̄) / z^m.ϕ) + (1 - m.α) * crra(n / z^m.ϕ, m.ζ)

fn_p(t; m) = m.p_1950 * exp(-m.γ * (t - 1950))
fn_w(t; m) = m.w_1950 * exp(m.Δw * (t - 1950))

function construct_F(μₛ, σₛ, b_grid)
    Δb = b_grid[2] - b_grid[1]
    F = zeros(length(b_grid))
    F[begin] = cdf(Normal(μₛ, σₛ), b_grid[begin] + 0.5Δb)
    for i in 2:(length(b_grid)-1)
        F[i] = cdf(Normal(μₛ, σₛ), b_grid[i] + 0.5Δb) - cdf(Normal(μₛ, σₛ), b_grid[i] - 0.5Δb)
    end
    F[end] = 1 - cdf(Normal(μₛ, σₛ), b_grid[end] - 0.5Δb)
    return F
end


function Single(p, w; m::Model)
    (; c̄, α, ζ, θ, κ) = m
    C₁ = (p * (1 - θ) / θ)^(1 / (κ - 1))

    h = find_zero(h -> α * (θ * C₁^κ + 1 - θ)^(1 - ζ / κ) * h^(1 - ζ) -
                       (1 - α) * (1 - θ) * (1 - c̄ / w - h - p * C₁ * h), 0.5)
    d = C₁ * h
    n = (θ * d^κ + (1 - θ) * h^κ)^(1 / κ)
    c = w * (1 - h) - w * p * d
    v = u(c, n; m, z=1)
    return Single(c, h, 1 - h, d, n, v)
end

function Married(p, w; m::Model)
    (; c̄, α, ζ, θ, κ, ϕ) = m
    C₁ = (p * (1 - θ) / θ)^(1 / (κ - 1))

    h = find_zero(h -> α * (θ * C₁^κ + 1 - θ)^(1 - ζ / κ) * h^(1 - ζ) -
                       (1 - α) * (1 - θ) * 2^(-ϕ * ζ) * (2 - c̄ / w - h - p * C₁ * h), 0.5)
    d = C₁ * h
    n = (θ * d^κ + (1 - θ) * h^κ)^(1 / κ)
    c = w * (2 - h) - w * p * d
    v = u(c, n; m, z=2)
    return Married(c, h, 2 - h, d, n, v)
end

function solve_dt(uₘ, uₛ; m::Model, tol=1e-6, max_iter=1000)
    (; n_b, b_grid, F, G, β, δ) = m

     # VFI ----------------------------------------------------------------------
    dist = Inf
    iter = 0
    V = zeros(n_b)
    W = 0.0
    V′, W′ = similar(V), W
    while dist > tol && iter < max_iter
        for i = 1:n_b
            V′[i] = uₘ + b_grid[i]
            for j in 1:n_b
                V′[i] += β * max(V[j], W) * G[i, j]
            end
        end
        W′ = uₛ
        for i in 1:n_b
            W′ += β * max(V′[i], W) * F[i]
        end

        dist = maximum(abs, (V′ .- V)) + abs(W′ - W)
        V .= V′
        W = W′
        iter += 1
    end
    
    # Steady State Distributions -----------------------------------------------
    ι = searchsortedfirst(V, W) # V[ι-1] < W ≤ V[ι]
    ω = (W - V[ι-1]) / (V[ι] - V[ι-1])
    P = zeros(n_b + 1, n_b + 1)

    ## Married (column j) → Married/Single
    for j in 1:n_b
        for i in (ι+1):n_b
            P[i, j] = G[j, i]
        end
        P[ι, j] = G[j, ι] * (1 - ω)
        P[n_b+1, j] = sum(G[j, k] for k in 1:(ι-1)) + G[j, ι] * ω
    end

    ## Single (column n_b+1) → Married/Single
    for i in (ι+1):n_b
        P[i, n_b+1] = F[i]
    end
    P[ι, n_b+1] = F[ι] * (1 - ω)
    P[n_b+1, n_b+1] = sum(F[k] for k in 1:(ι-1)) + F[ι] * ω

    prob_marriage = 1 - P[n_b+1, n_b+1]


    ## Solve M̃ = (1-δ)P M̃ + d where d = (0,...,0,δ)
    M̃ = (I - (1 - δ) * P) \ vcat(zeros(n_b), δ)
    s = M̃[n_b+1]

    prob_divorce = sum(P[n_b+1, i] * M̃[i] for i in 1:n_b) / (1 - s)

    return (s=s, pm=prob_marriage, pd=prob_divorce)
end

function simulate_year(t; m::Model)
    AS = Single(fn_p(t; m), fn_w(t; m); m)
    AM = Married(fn_p(t; m), fn_w(t; m); m)
    s, pm, pd = solve_dt(AM.v, AS.v; m)
    return (t=t, s=s, pm=pm, pd=pd)
end

function simulate(; m=Model(), years=1950:2000)
    rows = simulate_year.(years; m)
    return DataFrame(rows)
end


end # module

using Revise
using .My
using ProjectRoot
using YAML

dir = @projectroot("output")
m = My.Model()
d = Dict("phi" => m.ϕ,
    "alpha" => m.α,
    "zeta" => m.ζ,
    "c_bar" => m.c̄,
    "theta" => m.θ,
    "kappa" => m.κ,
    "mu_s" => m.μₛ,
    "sigma_s" => m.σₛ,
    "mu_m" => m.μₘ,
    "sigma_m" => m.σₘ,
    "beta_tilde" => m.β̃,
    "delta" => m.δ,
    "varrho" => m.ϱ,
    "lambda" => m.λ,
    "rho_tilde" => m.ρ̃,
    "nu" => m.ν,
    "eta" => m.η,
    "Delta_t" => m.Δt)
YAML.write_file("$dir/gg09.yaml", d)

df = My.simulate()

using Plots

p1 = plot(df.t, df.pm, label=false, title="Marriages", ylims=(0.05, 0.2))
p2 = plot(df.t, df.pd, label=false, title="Divorces", ylims=(0.01, 0.025))

plot(p1, p2)
