
using Tensors

# Define the Chaboche struct
struct ChabocheNorton{T}<:AbstractMaterial
    G::T
    K::T
    Y::T
    Hiso::T
    κ∞::T
    Hkin::T
    β∞::T
    t_star::T
    n::T
end

struct ChabocheNortonState{T} <: AbstractMaterialState
    β::SymmetricTensor{2,3,T,6}
    κ::T
    ϵp::SymmetricTensor{2,3,T,6}
end

"""
    Chaboche(; G, K, Y, Hiso, κ∞, Hkin, β∞)

Keyword constructor for the Chaboche struct
"""
ChabocheNorton(; G, K, Y, Hiso, κ∞, Hkin, β∞, t_star, n) = ChabocheNorton(G, K, Y, Hiso, κ∞, Hkin, β∞, t_star, n)
MaterialModelsBase.initial_material_state(::ChabocheNorton) = ChabocheNortonState(zero(SymmetricTensor{2,3}),0.0,zero(SymmetricTensor{2,3}))

function elastic_stress(m::ChabocheNorton, ϵe::SymmetricTensor{2,3})
    return 2*m.G*dev(ϵe) + 3*m.K*vol(ϵe)
end

function elastic_stiffness(m::ChabocheNorton)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*m.G*(I4sym-I4vol) + 3*m.K*I4vol
end

function extract_unknowns(::ChabocheNorton, x)
    return frommandel(SymmetricTensor{2,3}, x),frommandel(SymmetricTensor{2,3},x;offset=6),x[end-1], x[end] #returns ϵp, Δλ
end

function sigma_from_x(m::ChabocheNorton, ϵ, x)
    ϵp, _, _, _ = extract_unknowns(m, x)
    return elastic_stress(m, ϵ-ϵp) #Given the vector of unknowns, `x`, and the total strain `ϵ`, return the stress
end

function MaterialModelsBase.material_response(m::ChabocheNorton, ϵ::SymmetricTensor{2,3}, old_state::ChabocheNortonState, Δt,displacement)
    σ_trial=elastic_stress(m,ϵ-old_state.ϵp)
    if vonmises(σ_trial-old_state.β)-(m.Y+old_state.κ)< 0 #elastic
        return σ_trial, elastic_stiffness(m), old_state
    else #plastic
        rf!(r_,x_)=residual!(r_, x_, m, ϵ, Δt, old_state)
        x0=zeros(14)
        tomandel!(x0,old_state.ϵp)
        tomandel!(x0, old_state.β; offset=6)
        x0[13]=old_state.κ
        x, ∂r∂x, converged = newtonsolve(x0, rf!,tol=1e-10)
            if converged
                σ=sigma_from_x(m, ϵ, x)
                dσdϵ=calculate_ats(m, x, ϵ, old_state, ∂r∂x, Δt)
                ϵp, β, κ, _ = extract_unknowns(m, x)
                new_state = ChabocheNortonState(β, κ, ϵp)
                return σ, dσdϵ, new_state
            else
                throw(ErrorException("Chaboche did not converge"))
            end
        end
    end 
            

function calculate_ats(m::ChabocheNorton, x::Vector, ϵ::SymmetricTensor{2,3}, old_state, ∂r∂x::Matrix, Δt)
    # dσdϵ = ∂σ∂ϵ + ∂σ∂x ∂x∂ϵ
    # Problem: x is an implicit function of ϵ
    # drdϵ = 0 = ∂r∂ϵ + ∂r∂x ∂x∂ϵ => ∂x∂ϵ = -∂r∂x\∂r∂ϵ

    ∂σ∂ϵ = elastic_stiffness(m)
    rf_ϵ!(r_, ϵv_) = residual!(r_, x, m, frommandel(SymmetricTensor{2,3}, ϵv_), Δt, old_state)
    ∂r∂ϵ = ForwardDiff.jacobian(rf_ϵ!, zeros(14), tomandel(ϵ))
    ∂x∂ϵ = -∂r∂x\∂r∂ϵ
    
    ∂σ∂x = ForwardDiff.jacobian(x_ -> tomandel(sigma_from_x(m, ϵ, x_)), x)
    return ∂σ∂ϵ + frommandel(SymmetricTensor{4,3}, ∂σ∂x*∂x∂ϵ)
end

macaulay(x)=max(zero(x),x)

function residual!(r::Vector, x::Vector, m::ChabocheNorton, ϵ::SymmetricTensor{2,3}, Δt, old_state::ChabocheNortonState)
    ϵp, β, κ, Δλ = extract_unknowns(m, x)
    σ = elastic_stress(m, ϵ-ϵp)
    σ_red=σ-β
    v=gradient(vonmises,σ_red)
    Φ=vonmises(σ_red) -(m.Y+κ)

    rϵ=(ϵp-old_state.ϵp)-(Δλ*v)
    rλ=Δλ-(Δt/m.t_star)*(macaulay(Φ)/(m.Y+κ))^m.n
    rβ=(β-old_state.β)- 2*(Δλ/3)*m.Hkin*(v-(3/2)*(β/m.β∞))
    rκ=(κ-old_state.κ)-Δλ*m.Hiso*(1-(κ/m.κ∞))
    

    tomandel!(r,rϵ)
    tomandel!(r,rβ;offset=6)
    r[13]=rκ
    r[14]=rλ
end

function twice_yield_limit_in_strain(m::AbstractMaterial)
    return (2*m.Y)/((9*m.K*m.G)/(3*m.K+m.G))
end