using Tensors

x0=zeros(14)

ϵ=SymmetricTensor{2,3}([1 2 3; 2 1 2; 3 2 1])
tomandel!(x0,ϵ)
β=SymmetricTensor{2,3}([5 5 5; 5 5 5; 5 5 5])
tomandel!(x0,β;offset=6)

frommandel(SymmetricTensor{2,3},x0;offset=6)