push!(LOAD_PATH, "../src/")
using Documenter, DigitalComm

makedocs(sitename="DigitalComm.jl", 
		 format = Documenter.HTML(),
		 pages    = Any[
						"Introduction to DigitalComm"   => "index.md",
						"Function list"          => "base.md",
						"Examples"                => Any[ 
														 "Examples/example_AWGN.md",
														 "Examples/example_BER.md",
														 "Examples/example_PSD.md",
														 ],
						],
		 );

#makedocs(sitename="My Documentation", format = Documenter.HTML(prettyurls = false))

deploydocs(
    repo = "github.com/RGerzaguet/DigitalComm.jl",
)
