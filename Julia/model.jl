module My

using Roots
using Distributions
using DataFrames
using QuantEcon
using LinearAlgebra
using SparseArrays

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

@kwdef struct Model{TF<:AbstractFloat}
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

    # Continuous-time (model period = 1 year)
    Δt::TF = 1.0                      # model period length (years)
    λ::TF = 1 / Δt                    # meeting rate (per year)
    ρ̃::TF = -log(β) / Δt              # discount rate
    ν::TF = -log(1 - δ) / Δt          # exit rate
    ρ::TF = ρ̃ + ν                     # effective discount rate
    η::TF = -log(ϱ) / Δt              # OU mean-reversion speed

    # Numerical
    Δ::TF = 1000.0                    # implicit scheme pseudo-time step
end

crra(x, ζ) = iszero(ζ) ? log(x) : x^ζ / ζ
u(c, n; m, z) = m.α * log((c - m.c̄) / z^m.ϕ) + (1 - m.α) * crra(n / z^m.ϕ, m.ζ)

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

function construct_grid_dt(m::Model; n_b, n_std)
    (; μₛ, σₛ, ϱ, μₘ, σₘ) = m

    mc = tauchen(n_b, ϱ, sqrt(1 - ϱ^2) * σₘ, (1 - ϱ) * μₘ, n_std)
    b_grid = collect(mc.state_values)
    G = mc.p

    Δb = b_grid[2] - b_grid[1]
    F = zeros(length(b_grid))
    F[begin] = cdf(Normal(μₛ, σₛ), b_grid[begin] + 0.5Δb)
    for i in 2:(length(b_grid)-1)
        F[i] = cdf(Normal(μₛ, σₛ), b_grid[i] + 0.5Δb) - cdf(Normal(μₛ, σₛ), b_grid[i] - 0.5Δb)
    end
    F[end] = 1 - cdf(Normal(μₛ, σₛ), b_grid[end] - 0.5Δb)
    return b_grid, F, G
end

function construct_grid_ct(m::Model; n_b, n_std)
    (; μₘ, σₘ, μₛ, σₛ, η) = m
    b_grid = collect(range(μₘ - n_std * σₘ, μₘ + n_std * σₘ, length=n_b))
    Δb = b_grid[2] - b_grid[1]
    diff = η * σₘ^2  # diffusion coefficient

    # OU generator A: upwind drift + central diffusion
    A = spzeros(n_b, n_b)
    for (i, b) in enumerate(b_grid)
        d = η * (μₘ - b)  # drift at grid point i

        if 1 < i < n_b  # interior
            if d >= 0
                A[i, i] -= d / Δb
                A[i, i+1] += d / Δb
            else
                A[i, i] += d / Δb
                A[i, i-1] -= d / Δb
            end
            A[i, i-1] += diff / Δb^2
            A[i, i] -= 2diff / Δb^2
            A[i, i+1] += diff / Δb^2
        elseif i == 1  # lower reflecting boundary
            if d >= 0
                A[i, i] -= d / Δb
                A[i, i+1] += d / Δb
            end
            A[i, i] -= diff / Δb^2
            A[i, i+1] += diff / Δb^2
        else  # i == n_b, upper reflecting boundary
            if d < 0
                A[i, i] += d / Δb
                A[i, i-1] -= d / Δb
            end
            A[i, i-1] += diff / Δb^2
            A[i, i] -= diff / Δb^2
        end
    end

    f = [pdf(Normal(μₛ, σₛ), b) for b in b_grid]

    # Pre-factored implicit matrix
    B = spdiagm(0 => fill(1 / m.Δ + m.ρ, n_b)) - A
    B_fact = lu(B)
    return n_b, b_grid, Δb, A, f, B_fact
end

function solve_dt(vₛ, vₘ, n_b, b_grid, F, G; m::Model, tol=1e-6, max_iter=1000)
    (; β, δ) = m

    # VFI ----------------------------------------------------------------------
    dist = Inf
    iter = 0
    V = zeros(n_b)
    W = 0.0
    V′, W′ = similar(V), W
    while dist > tol && iter < max_iter
        for i = 1:n_b
            V′[i] = vₘ + b_grid[i]
            for j in 1:n_b
                V′[i] += β * max(V[j], W) * G[i, j]
            end
        end
        W′ = vₛ
        for i in 1:n_b
            W′ += β * max(V′[i], W) * F[i]
        end

        dist = maximum(abs, (V′ .- V)) + abs(W′ - W)
        V .= V′
        W = W′
        iter += 1
    end

    # Steady State Distributions -----------------------------------------------
    ι = searchsortedfirst(V, W) # V[ι-1] < W < V[ι]
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

function solve_ct(vₛ, vₘ, n_b, b_grid, Δb, A, f, B_fact; m::Model, tol=1e-6, max_inner=500, max_outer=200)
    (; ρ, λ, ν, Δ, Δt, μₛ, σₛ) = m

    𝐮 = vₘ .+ b_grid
    W = vₛ / ρ
    Ṽ = 𝐮 ./ ρ
    V = max.(Ṽ, W)
    W′ = W
    V′ = similar(V)
    dist_s = Normal(μₛ, σₛ)

    # HJB: nested loop with variational inequality ----------------------------
    for _ in 1:max_outer
        for _ in 1:max_inner
            Ṽ .= B_fact \ (𝐮 .+ V ./ Δ)
            V′ .= max.(Ṽ, W)
            if maximum(abs, V′ .- V) < tol
                break
            end
            V .= V′
        end

        # Update W using smooth threshold from unclamped V --------------------
        ι = searchsortedfirst(Ṽ, W)  # Ṽ[ι-1] < W < Ṽ[ι]
        if 1 < ι < n_b
            ω = clamp((W - Ṽ[ι-1]) / (Ṽ[ι] - Ṽ[ι-1]), 0.0, 1.0)
            b̃ = b_grid[ι-1] + ω * Δb

            P_acc = 1 - cdf(dist_s, b̃)
            EV = (1 - ω) * V[ι-1] * f[ι-1] * Δb
            for k in ι:n_b
                EV += V[k] * f[k] * Δb
            end

            W′ = (vₛ + λ * EV) / (ρ + λ * P_acc)
            W′ = 0.5 * W + 0.5 * W′
        elseif ι == 1
            EV = sum(V[k] * f[k] for k in 1:n_b) * Δb
            W′ = (vₛ + λ * EV) / (ρ + λ)
            W′ = 0.5 * W + 0.5 * W′
        else # ι == n_b + 1
            W′ = vₛ / ρ
        end

        if abs(W - W′) < tol
            break
        end
        W = W′
    end

    # KFE with weighted boundary (single solve, O(N)) --------------------------
    ι = searchsortedfirst(Ṽ, W)  # Ṽ[ι-1] < W < Ṽ[ι]
    if ι <= 1
        return (s=0.01, pm=1.0, pd=0.0)
    elseif ι >= n_b
        return (s=1.0, pm=0.0, pd=0.0)
    end
    ω = clamp((W - Ṽ[ι-1]) / (Ṽ[ι] - Ṽ[ι-1]), 0.0, 1.0)
    b̃ = b_grid[ι-1] + ω * Δb
    P_acc = 1 - cdf(dist_s, b̃)

    ## Solve KFE at two clean absorbing boundaries and interpolate s
    s_lo = _kfe_single(A, f, ι - 1, n_b, λ, ν, Δb)
    s_hi = _kfe_single(A, f, ι,     n_b, λ, ν, Δb)
    s = clamp((1 - ω) * s_lo + ω * s_hi, 0.01, 0.99)

    # Convert CT rates → DT per-period probabilities ---------------------------
    prob_marriage = 1 - exp(-λ * P_acc * Δt)
    divorce_rate = max(λ * P_acc * s / (1 - s) - ν, 0.0)
    prob_divorce = 1 - exp(-divorce_rate * Δt)

    return (s=s, pm=prob_marriage, pd=prob_divorce)
end

function _kfe_single(A, f, i_start, n_b, λ, ν, Δb)
    i_start = max(i_start, 1)
    idx = i_start:n_b
    n_ι = length(idx)

    dd = [A[idx[i], idx[i]] - ν for i in 1:n_ι]
    dl = [A[idx[i], idx[i+1]] for i in 1:n_ι-1]
    du = [A[idx[i+1], idx[i]] for i in 1:n_ι-1]

    T̃ = Tridiagonal(dl, dd, du)
    f̃ = f[idx]

    𝐳 = T̃ \ f̃
    return clamp(1.0 / (1 - λ * Δb * sum(𝐳)), 0.01, 0.99)
end


function simulate_year_dt(t, n_b, b_grid, F, G; m::Model)
    AS = Single(fn_p(t; m), fn_w(t; m); m)
    AM = Married(fn_p(t; m), fn_w(t; m); m)
    s, pm, pd = solve_dt(AS.v, AM.v, n_b, b_grid, F, G; m)
    return (t=t, vs=AS.v, vm=AM.v, s=s, pm=pm, pd=pd, method="DT")
end

function simulate_year_ct(t, n_b, b_grid, Δb, A, f_ct, B_fact; m::Model)
    AS = Single(fn_p(t; m), fn_w(t; m); m)
    AM = Married(fn_p(t; m), fn_w(t; m); m)
    s, pm, pd = solve_ct(AS.v, AM.v, n_b, b_grid, Δb, A, f_ct, B_fact; m)
    return (t=t, vs=AS.v, vm=AM.v, s=s, pm=pm, pd=pd, method="CT")
end

function simulate(; m=Model(), years=1950:2000, method=:dt, n_b=100, n_std=4.0)
    if method == :dt
        b_grid, F, G = construct_grid_dt(m, n_b=n_b, n_std=n_std)
        n_b = length(b_grid)
        rows = [simulate_year_dt(t, n_b, b_grid, F, G; m) for t in years]
        return DataFrame(rows)
    elseif method == :ct
        n_b, b_grid, Δb, A, f_ct, B_fact = construct_grid_ct(m, n_b=n_b, n_std=n_std)
        rows = [simulate_year_ct(t, n_b, b_grid, Δb, A, f_ct, B_fact; m) for t in years]
        return DataFrame(rows)
    else
        error("Invalid method: choose :dt or :ct")
    end
end


end # module
