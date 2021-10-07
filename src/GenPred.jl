module GenPred

# Write your package code here.
using LinearAlgebra
using Statistics

using Optim

export simple_imp, simple_imp!
export vanraden_grm

include("impute.jl")
include("vanraden_grm.jl")

end ##module
