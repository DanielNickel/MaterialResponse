# CA2 readme

## Getting started
With the included manifest, it should work to do run 
```julia
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```
in the REPL from inside the present folder. 
You can check where you are by using "print working directory": `pwd()`,
```julia
julia> pwd()
```

## References
See documentation for [MaterialModelsBase.jl](https://knutam.github.io/MaterialModelsBase.jl/dev/)
and [Newton.jl](https://knutam.github.io/Newton.jl/dev/)

## Check your implementation
The figure `demofigure.pdf` gives the results for the Chaboche model with the settings given in `demo.jl`.
It can be used to verify that your implementation is correct. 
Note that `demo.jl` will not run without implementing the `Chaboche` model, but if you comment out the 
`m_chaboche = ...` line and use `m_perfect` instead, you can plot the results for the supplied perfect
plasticity model instead. 

## Advanced usage
If you need to add the packages `MaterialModelsBase.jl` and 
`Newton.jl` to another environment (if you don't build on top
of the present environment) these can be added by running 
```julia
julia> using Pkg
julia> Pkg.add(;url="https://github.com/KnutAM/MaterialModelsBase.jl.git")
julia> Pkg.add(;url="https://github.com/KnutAM/Newton.jl.git")
```
