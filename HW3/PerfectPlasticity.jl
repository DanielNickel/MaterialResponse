using Tensors
using Newton
import CairoMakie as CM
struct PerfectPlasticity
    G::Float64 # Shear modulus
    K::Float64 # Bulk modulus 
    Y::Float64 # Yield limit
end

elastic_stress(m::PerfectPlasticity, ϵe) = 2*m.G*dev(ϵe) + 3*m.K*vol(ϵe)

vonmises(σ) = sqrt((3/2)*dev(σ)⊡dev(σ))

function my_material_response(
    m::PerfectPlasticity, # material parameters
    ϵ::SymmetricTensor{2,3}, # current strain
    old_ϵp::SymmetricTensor{2,3}) # old value of plastic strain

    σ_trial = elastic_stress(m, ϵ-old_ϵp) # Stress assuming old_ϵp is constant
    Φ_trial = vonmises(σ_trial) - m.Y
    if Φ_trial <= 0 # Elastic response, plastic strain not changing
        return σ_trial, old_ϵp
    else # Plastic response
        # We need to find ϵp and Δλ such that 
        # (ϵp-old_ϵp) = Δλ ∂Φ/∂σ
        # Φ=0
        rf!(r, x) = residual!(r, x, m, ϵ, old_ϵp)
        # Make an initial guess
        x0 = zeros(7)
        tomandel!(x0, old_ϵp) # First guess is that ϵp = old_ϵp
        x, ∂r∂x, converged = newtonsolve(x0, rf!) # Use Newton.jl
        # Should check that converged == true 
        # If so, we have found an x such that r(x) = 0!
        ϵp = frommandel(SymmetricTensor{2,3}, x) # x[1:6]
        σ = elastic_stress(m, ϵ-ϵp)
        return σ, ϵp
    end
end

# Quick note about Mandel versus Voigt notation
# Mandel notation of tensor a: [a11, a22, a33, √2*a12, √2*a13, √2*a23]
# Voigt notation of tensor a: [a11, a22, a33, a12, a13, a23]

function residual!(r::Vector, x::Vector, m::PerfectPlasticity, ϵ, old_ϵp)
    # Translate x[1:6] in Mandel notation to SymmetricTensor
    ϵp = frommandel(SymmetricTensor{2,3}, x) # Only takes the first 6 components
    Δλ = x[7]
    
    # Do the calculations
    σ = elastic_stress(m, ϵ-ϵp)
    ∂Φ∂σ = gradient(vonmises, σ) # Using automatic differentiation
    r_ϵp = (ϵp - old_ϵp) - Δλ*∂Φ∂σ  # dϵp/dt = dλ/dt ∂Φ/∂σ
                                    # ϵp - old_ϵp = Δϵp = Δλ ∂Φ/∂σ
    r_Φ = vonmises(σ) - m.Y         # Φ = 0
    
    # Return the vector r_Φ
    tomandel!(r, r_ϵp)
    r[7] = r_Φ
end


m = PerfectPlasticity(80e3, 160e3, 300.0)
ϵ12_vector = collect(range(0, 0.01, 100))

function calculate_stress_history(m, ϵ12_vector)    
    ϵp = zero(SymmetricTensor{2,3})
    σ12_vector = Float64[]
    for ϵ12 in ϵ12_vector
        ϵ = SymmetricTensor{2,3}((i,j)->i==2  && j==1 ? ϵ12 : 0.0)
        old_ϵp = ϵp
        σ, ϵp = my_material_response(m, ϵ, old_ϵp)
        push!(σ12_vector, σ[1,2])
    end
    return σ12_vector
end

σ12_vector = calculate_stress_history(m, ϵ12_vector)

fig, ax, line = CM.lines(ϵ12_vector, σ12_vector)
CM.save("results.pdf", fig)
fig