push!(LOAD_PATH, "../src/")
using Documenter, DigitalComm

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
