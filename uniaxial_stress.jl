"""
    function uniaxial_stress(m::AbstractMaterial, ϵ11::Vector, t::Vector)    
    
For a time history of the ``\\epsilon_{11}`` component, return 
history of the ``\\sigma_{11}`` component assuming a uniaxial stress
state. Requires that `initial_material_state(m)` and `material_response(m, ...)`
have been defined. 
"""
function uniaxial_stress(m::AbstractMaterial, ϵ11::Vector, t::Vector, displacement)
    state = initial_material_state(m)
    σ11 = zeros(length(ϵ11))
    ϵ = zero(SymmetricTensor{2,3})
    for i in 2:length(ϵ11)
        Δt = t[i]-t[i-1]
        ϵ_guess = SymmetricTensor{2,3}((k,l)->k==l==1 ? ϵ11[i] : ϵ[k,l])
        ϵ, σ, state = uniaxial_stress_iterations(m, state, Δt, ϵ_guess, displacement)
        σ11[i] = σ[1,1]
    end
    return σ11 
end

"""
    function uniaxial_stress_iterations(m::AbstractMaterial, old_state, Δt, ϵ_guess; maxiter=10, tol=1.e-6)

Iterate to find a uniaxial stress state (σ22=σ33=0). 
The guess `ϵ_guess` should contain the prescribed values for 
indicies `i=j=1` and `i≠j`, and the initial guess for `i=j∈{2,3}`.
Returns the resulting strain, `ϵ`, stress, `σ`, and state `new_state`.
"""
function uniaxial_stress_iterations(m::AbstractMaterial, old_state, Δt, ϵ_guess,displacement; maxiter=10, tol=1.e-6)    
    local ϵ, σ, new_state
    ϵ = ϵ_guess
    r = zeros(2)
    dσdϵᴹ = zeros(6,6)
    for iter_index in 1:maxiter
        σ, dσdϵ, new_state = material_response(m, ϵ, old_state, Δt, displacement)
        r[1] = σ[2,2]
        r[2] = σ[3,3]
        norm(r) < tol && break
        iter_index >= maxiter && throw(ErrorException("Did not converge in stress iterations"))
        tomandel!(dσdϵᴹ, dσdϵ)
        Δx = -dσdϵᴹ[2:3,2:3]\r
        Δϵ = ϵ_from_x(Δx, 0.0)
        ϵ += Δϵ
    end
    return ϵ, σ, new_state
end

"""
    function ϵ_from_x(x, ϵ11, ϵ12=0.0)

Given the vector of unknowns, `x=[ϵ22, ϵ33]` as well as 
`ϵ11` (and ϵ12, but not necessary as it defaults to zero), 
construct the full strain tensor assuming that 
all shear strain components are zero (except ϵ12 if given)
"""
function ϵ_from_x(x, ϵ11, ϵ12=0.0) 
    SymmetricTensor{2,3}((i,j)-> i==j==1 ? ϵ11 :
                                 i==j ? x[i-1] :
                                 (i==2 && j==1) ? ϵ12 : 0.0)
end



