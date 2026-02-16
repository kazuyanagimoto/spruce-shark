module My

using Roots
using Distributions
using DataFrames
using QuantEcon

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

@kwdef struct Model{TF <: AbstractFloat, TI <: Integer}
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
    ρ̃::TF = -log(β) / Δt
    ν::TF = -log(1-δ) / Δt
    ρ::TF = ρ̃ + ν
    η::TF = -log(ϱ)
    λ::TF = 1 / Δt

    # Grid
    n_b::TI = 100
    mc::MarkovChain = tauchen(n_b, ϱ, sqrt(1 - ϱ^2) * σₘ, (1 - ϱ) * μₘ)
    b_grid::Vector{TF} = collect(mc.state_values)
    G::Matrix{TF} = mc.p
    F::Vector{TF} = 1 .- cdf.(Normal(μₛ, σₛ), b_grid)
end

crra(x, ζ) = iszero(ζ) ? log(x) : x^ζ / ζ
u(c, n; m, z) = m.α * log((c - m.c̄)) / z^m.ϕ + (1 - m.α) * crra(n / z^m.ϕ, m.ζ)

fn_p(t; m) = m.p_1950 * exp(-m.γ * (t - 1950))
fn_w(t; m) = m.w_1950 * exp(m.Δw * (t - 1950))


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

function solve(t; m::Model, tol=1e-6, max_iter=1000)
    (; n_b, b_grid, F, G, β, δ, μₛ, σₛ, μₘ, σₘ, ϱ) = m

    AS = Single(fn_p(t; m), fn_w(t; m); m)
    AM = Married(fn_p(t; m), fn_w(t; m); m)

    # VFI
    dist = Inf
    iter = 0
    V = zeros(n_b)
    V′ = similar(V)
    W, W′ = 0., 0.
    while dist > tol && iter < max_iter
        for i = 1:n_b
            V′[i] = AM.v + b_grid[i]
            for j in 1:n_b
                V′[i] += β * max(V[j], W) * G[i, j]
            end
        end
        W′ = AS.v
        for i in 1:n_b
            W′ += β * max(V[i], W) * F[i]
        end

        dist = maximum(abs, (V′ .- V)) + abs(W′ - W)
        V .= V′
        W = W′
        iter += 1
    end

    # Steady State Distributions
    ĩ = searchsortedfirst(V, W) # V[ĩ-1] < W ≤ V[ĩ]
    ω = (W - V[ĩ-1]) / (V[ĩ] - V[ĩ-1]) # weight
    b̃ = (1-ω) * b_grid[ĩ-1] + ω * b_grid[ĩ] # V(b̃) = W
    s = 0.
    M = zeros(n_b)
    M[ĩ:end] .= 1 / (n_b - ĩ + 1)
    M′ = copy(M)
    dist = Inf
    iter = 0
    while dist > tol && iter < max_iter
        for j = ĩ:n_b
            M′[j] = (δ + (1 - δ) * s) * F[j]
            for i in eachindex(b_grid)
                M′[j] += (1 - δ) * G[i, j] * M[i]
            end
        end
        dist = maximum(abs, (M′ .- M))
        M .= M′
        s = 1 - sum(M)
        iter += 1
    end

    prob_marriage = 1 - cdf(Normal(μₛ, σₛ), b̃)
    prob_divorce = sum(G[i, j] * M[i] for i in ĩ:n_b, j = 1:ĩ-1) / (1 - s)

    return (t=t, s=s, us=AS.v, um=AM.v, pm=prob_marriage, pd=prob_divorce)
end

function simulate(;m=Model(), years=1950:2000)
    rows = solve.(years; m)
    return DataFrame(rows)
end


end # module

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

My.solve(1950; m)