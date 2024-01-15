using Documenter, DigitalComm

# Runs the doctests (the jldoctest) in funciton documentation.
DocMeta.setdocmeta!(DigitalComm, :DocTestSetup, :(using DigitalComm); recursive=true)
Documenter.doctest(DigitalComm)