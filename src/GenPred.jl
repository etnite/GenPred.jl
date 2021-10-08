module GenPred

# Packages included in Base install
using LinearAlgebra
using Statistics

# External packages
using Optim

export simple_imp, simple_imp!
export vanraden_grm
export emma_solve

include("impute.jl")
include("vanraden_grm.jl")
include("emma_solve.jl")

end ##module
