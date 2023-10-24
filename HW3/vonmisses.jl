using Tensors
using ForwardDiff

function vonmises(σ::SymmetricTensor{2,3})
    σ_dev = dev(σ)
    return sqrt((3.0/2.0) * (σ_dev ⊡ σ_dev))
end

σ1 = SymmetricTensor{2,3}((i,j)->i==j==1 ? 1.0 : 0.0)
display(vonmises(σ1))

# Part 2
σ_dev = dev(σ1)
ϕ(x)= sqrt((3.0/2.0) * (dev(x) ⊡ dev(x)))

δϕ_δσ = Tensors.gradient(ϕ, σ1)
self_derivative = (3/2) * (σ_dev / sqrt((3/2) * (σ_dev ⊡ σ_dev)))

display(self_derivative)

if δϕ_δσ == self_derivative
    display("true")
else
    display("false")
end