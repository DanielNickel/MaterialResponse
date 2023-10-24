using CairoMakie
using ForwardDiff
using Newton

struct SlidingSpring
    k::Float64
    Fy::Float64
end

function sliding_spring(m::SlidingSpring, u, up_old)
    F_trial = m.k * u
    print(u)

    Φ_trial = F_trial - m.Fy

    if Φ_trial < 0
        return F_trial, up_old
    else
        rf!(r, x) = residual!(r, x, m, u, up_old)  # Fixed variable name typo
        x0 = 0
        x, ∂r∂x, converged = newtonsolve(x0, rf!)
        Δλ = x
        up = Δλ * ∂r∂x
        F = m.k * u
        return F, up
    end
end

function residual!(r::Float64, x::Float64, m::SlidingSpring, u, up_old)
    # Translate x[1:6] in Mandel notation to SymmetricTensor
    F_trial = m.k * u
    ϕ(x) = abs(x) - m.Fy  # Fixed variable name typo
    Δλ = x
    ∂Φ∂F = gradient(ϕ, F_trial)
    r = (u - up_old) - Δλ * ∂Φ∂F
end

m = SlidingSpring(10.0, 1.0)
u = append!(collect(0:0.01:0.14), collect(0.15:-0.01:-0.15))  # Fixed syntax error in array concatenation
F = zeros(length(u))
up = 0.0
for i in eachindex(u)
    global up
    F[i], up = sliding_spring(m, u[i], up)
end
lines(u, F)
