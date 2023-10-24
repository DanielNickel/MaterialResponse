# Utility functions
"""
    vonmises(σ::SymmetricTensor{2,3})

Calculate the effective von Mises stress for the given 
stress tensor. 
"""
function vonmises(σ::SymmetricTensor{2,3})
    σ_dev = dev(σ)
    return sqrt((3.0/2.0) * (σ_dev ⊡ σ_dev))
end