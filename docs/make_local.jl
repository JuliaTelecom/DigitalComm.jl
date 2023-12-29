# Assumes its being run from project root.

using Pkg
Pkg.activate("docs/")
Pkg.instantiate()

include("make.jl")