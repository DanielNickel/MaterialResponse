
using Tensors

# Define the Chaboche struct
struct Chaboche{T}<:AbstractMaterial
    G::T
    K::T
    Y::T
    Hiso::T
    κ∞::T
    Hkin::T
    β∞::T
end

struct ChabocheState{T} <: AbstractMaterialState
    β::SymmetricTensor{2,3,T,6}
    κ::T
    ϵp::SymmetricTensor{2,3,T,6}
end

"""
    Chaboche(; G, K, Y, Hiso, κ∞, Hkin, β∞)

Keyword constructor for the Chaboche struct
"""
Chaboche(; G, K, Y, Hiso, κ∞, Hkin, β∞) = Chaboche(G, K, Y, Hiso, κ∞, Hkin, β∞)
MaterialModelsBase.initial_material_state(::Chaboche) = ChabocheState(zero(SymmetricTensor{2,3}),0.0,zero(SymmetricTensor{2,3}))

function elastic_stress(m::Chaboche, ϵe::SymmetricTensor{2,3})
    return 2*m.G*dev(ϵe) + 3*m.K*vol(ϵe)
end

function elastic_stiffness(m::Chaboche)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*m.G*(I4sym-I4vol) + 3*m.K*I4vol
end

function extract_unknowns(::Chaboche, x)
    return frommandel(SymmetricTensor{2,3}, x),frommandel(SymmetricTensor{2,3},x;offset=6),x[end-1], x[end] #returns ϵp, Δλ
end

function sigma_from_x(m::Chaboche, ϵ, x)
    ϵp, _, _, _ = extract_unknowns(m, x)
    return elastic_stress(m, ϵ-ϵp) #Given the vector of unknowns, `x`, and the total strain `ϵ`, return the stress
end

function MaterialModelsBase.material_response(m::Chaboche, ϵ::SymmetricTensor{2,3}, old_state::ChabocheState, args...)
    σ_trial=elastic_stress(m,ϵ-old_state.ϵp)
    if vonmises(σ_trial-old_state.β)-(m.Y+old_state.κ)< 0 #elastic
        return σ_trial, elastic_stiffness(m), old_state
    else #plastic
        rf!(r_,x_)=residual!(r_, x_, m, ϵ, old_state)
        x0=zeros(14)
        tomandel!(x0,old_state.ϵp)
        tomandel!(x0, old_state.β; offset=6)
        x0[13]=old_state.κ
        x, ∂r∂x, converged = newtonsolve(x0, rf!)
            if converged
                σ=sigma_from_x(m, ϵ, x)
                dσdϵ=calculate_ats(m, x, ϵ, old_state, ∂r∂x)
                ϵp, β, κ, _ = extract_unknowns(m, x)
                new_state = ChabocheState(β, κ, ϵp) 
                return σ, dσdϵ, new_state
            else
                throw(ErrorException("Chaboche did not converge"))
            end
        end
    end
            

function calculate_ats(m::Chaboche, x::Vector, ϵ::SymmetricTensor{2,3}, old_state, ∂r∂x::Matrix)
    # dσdϵ = ∂σ∂ϵ + ∂σ∂x ∂x∂ϵ
    # Problem: x is an implicit function of ϵ
    # drdϵ = 0 = ∂r∂ϵ + ∂r∂x ∂x∂ϵ => ∂x∂ϵ = -∂r∂x\∂r∂ϵ

    ∂σ∂ϵ = elastic_stiffness(m)
    rf_ϵ!(r_, ϵv_) = residual!(r_, x, m, frommandel(SymmetricTensor{2,3}, ϵv_), old_state)
    ∂r∂ϵ = ForwardDiff.jacobian(rf_ϵ!, zeros(14), tomandel(ϵ))
    ∂x∂ϵ = -∂r∂x\∂r∂ϵ
    
    ∂σ∂x = ForwardDiff.jacobian(x_ -> tomandel(sigma_from_x(m, ϵ, x_)), x)
    return ∂σ∂ϵ + frommandel(SymmetricTensor{4,3}, ∂σ∂x*∂x∂ϵ)
end

function residual!(r::Vector, x::Vector, m::Chaboche, ϵ::SymmetricTensor{2,3}, old_state::ChabocheState)
    ϵp, β, κ, Δλ = extract_unknowns(m, x)
    σ = elastic_stress(m, ϵ-ϵp)
    σ_red=σ-β
    v=gradient(vonmises,σ_red)
    
    rϵ=(ϵp-old_state.ϵp)-(Δλ*v)
    rΦ=vonmises(σ_red) -(m.Y+κ)
    rβ=(β-old_state.β)- 2*(Δλ/3)*m.Hkin*(v-(3/2)*(β/m.β∞))
    rκ=(κ-old_state.κ)-Δλ*m.Hiso*(1-(κ/m.κ∞))

    tomandel!(r,rϵ)
    tomandel!(r,rβ;offset=6)
    r[13]=rκ
    r[14]=rΦ
end