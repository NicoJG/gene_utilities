I made this a Julia package so that the dependencies can be activated easily. Just do the following:
```julia
using Pkg
Pkg.activate("/path/to/julia_gene_tools")
Pkg.instantiate()
```