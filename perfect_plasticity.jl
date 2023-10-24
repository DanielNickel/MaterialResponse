# Perfect plasticity code
struct PerfectPlasticity{T} <: AbstractMaterial
    G::T    # Shear modulus
    K::T    # Bulk modulus
    Y::T    # Yield limit
end
"""
    PerfectPlasticity(; G, K, Y)

Keyword constructor for `PerfectPlasticity`, makes 
`PerfectPlasticity(G=2, K=4, Y=1)` equivalent to 
`PerfectPlasticity(Y=1, G=4, K=2)`
"""
PerfectPlasticity(;G, K, Y) = PerfectPlasticity(G, K, Y)    

struct PerfectPlasticityState{T} <: AbstractMaterialState
    ϵp::SymmetricTensor{2,3,T,6}
end
MaterialModelsBase.initial_material_state(::PerfectPlasticity) = PerfectPlasticityState(zero(SymmetricTensor{2,3}))

## Elastic functions
"""
    elastic_stress(m::PerfectPlasticity, ϵe::SymmetricTensor{2,3})

Calculate the stress given the elastic strain `ϵe`
"""
function elastic_stress(m::PerfectPlasticity, ϵe::SymmetricTensor{2,3})
    return 2*m.G*dev(ϵe) + 3*m.K*vol(ϵe)
end

"""
    elastic_stiffness(m::PerfectPlasticity)

Return the elastic stiffness tensor considering the parameters in `m`
"""
function elastic_stiffness(m::PerfectPlasticity)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*m.G*(I4sym-I4vol) + 3*m.K*I4vol
end

"""
    extract_unknowns(::PerfectPlasticity, x)

Given the vector of unknowns, `x`, extract the tensor values of
those, specifically it returns `ϵp, Δλ`
"""
function extract_unknowns(::PerfectPlasticity, x)
    return frommandel(SymmetricTensor{2,3}, x), x[end]
end

"""
    sigma_from_x(m::PerfectPlasticity, ϵ, x)

Given the vector of unknowns, `x`, and the total strain `ϵ`,
return the stress.
"""
function sigma_from_x(m::PerfectPlasticity, ϵ, x)
    ϵp, _ = extract_unknowns(m, x)
    return elastic_stress(m, ϵ-ϵp)
end

# Main function
"""
    material_response(
        m::PerfectPlasticity, ϵ::SymmetricTensor{2,3}, 
        old_state::PerfectPlasticityState, Δt)

Calculate the updated stress, ATS-tensor and state variable 
given a new strain and the old state. For this specific model,
the time increment has no effect. 
"""
function MaterialModelsBase.material_response(
    m::PerfectPlasticity, ϵ::SymmetricTensor{2,3}, 
    old_state::PerfectPlasticityState, args...)
    σ_trial = elastic_stress(m, ϵ-old_state.ϵp)
    if vonmises(σ_trial) - m.Y < 0
        # Elastic
        return σ_trial, elastic_stiffness(m), old_state
    else
        # Plastic
        rf!(r_,x_) = residual!(r_, x_, m, ϵ, old_state)
        x0 = zeros(7)               # Allocate unknowns
        tomandel!(x0, old_state.ϵp) # Set starting guess to old plastic strain
        x, ∂r∂x, converged = newtonsolve(x0, rf!)
        if converged
            σ = sigma_from_x(m, ϵ, x)
            dσdϵ = calculate_ats(m, x, ϵ, old_state, ∂r∂x)
            ϵp, _ = extract_unknowns(m, x)
            new_state = PerfectPlasticityState(ϵp) 
            return σ, dσdϵ, new_state
        else
            throw(ErrorException("PerfectPlasticity did not converge"))
        end
    end
end

"""
    calculate_ats(
        m::PerfectPlasticity, x::Vector, 
        ϵ::SymmetricTensor{2,3}, old_state, ∂r∂x::Matrix)

Given the unknowns `x`, the total strain `ϵ`, the old state `old_state`,
and the derivative of the residual function wrt. to the unknowns, `∂r∂x`,
calculate the algorithmic tangent stiffness tensor. 
"""
function calculate_ats(m::PerfectPlasticity, x::Vector, ϵ::SymmetricTensor{2,3}, old_state, ∂r∂x::Matrix)
    # dσdϵ = ∂σ∂ϵ + ∂σ∂x ∂x∂ϵ
    # Problem: x is an implicit function of ϵ
    # drdϵ = 0 = ∂r∂ϵ + ∂r∂x ∂x∂ϵ => ∂x∂ϵ = -∂r∂x\∂r∂ϵ

    ∂σ∂ϵ = elastic_stiffness(m)
    rf_ϵ!(r_, ϵv_) = residual!(r_, x, m, frommandel(SymmetricTensor{2,3}, ϵv_), old_state)
    ∂r∂ϵ = ForwardDiff.jacobian(rf_ϵ!, zeros(7), tomandel(ϵ))
    ∂x∂ϵ = -∂r∂x\∂r∂ϵ
    
    ∂σ∂x = ForwardDiff.jacobian(x_ -> tomandel(sigma_from_x(m, ϵ, x_)), x)
    return ∂σ∂ϵ + frommandel(SymmetricTensor{4,3}, ∂σ∂x*∂x∂ϵ)
end

"""
    residual!(
        r::Vector, x::Vector, m::PerfectPlasticity, 
        ϵ::SymmetricTensor{2,3}, old_state::PerfectPlasticityState)

Calculate the residual `r` in-place, given a value for the unknowns `x`.
"""
function residual!(r::Vector, x::Vector, m::PerfectPlasticity, ϵ::SymmetricTensor{2,3}, old_state::PerfectPlasticityState)
    ϵp, Δλ = extract_unknowns(m, x)
    σ = elastic_stress(m, ϵ-ϵp)
    ν = gradient(vonmises, σ)
    rϵ = (ϵp - old_state.ϵp) - Δλ*ν
    rΦ = vonmises(σ) - m.Y
    tomandel!(r, rϵ)    # Equivalent to R[1:6] .= tomandel(Rϵ)
    r[7] = rΦ
end
