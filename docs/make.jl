push!(LOAD_PATH, "../src/")

using Documenter

# Running `julia --project docs/make.jl` can be very slow locally.
# To speed it up during development, one can use make_local.jl instead.
# The code below checks whether it's being called from make_local.jl or not.
const LOCAL = get(ENV, "LOCAL", "false") == "true"

if LOCAL
    include("../src/DigitalComm.jl")
    using .DigitalComm
else
    using DigitalComm
    ENV["GKSwstype"] = "100"	# Prevents warnings in the doc build on github actions.
end

DocMeta.setdocmeta!(DigitalComm, :DocTestSetup, :(using DigitalComm); recursive=true)

makedocs(
	modules = [DigitalComm],
	sitename="DigitalComm.jl",
	format = Documenter.HTML(),
	pages    = Any[
		"Introduction to DigitalComm"   => "index.md",
		"Function list"          		=> "base.md",
		"Examples"                		=> Any[ 
											"Examples/example_AWGN.md",
											"Examples/example_BER.md",
											"Examples/example_PSD.md",
										],
	],
	doctest  = true,
);

deploydocs(
    repo = "github.com/JuliaTelecom/DigitalComm.jl",
)
