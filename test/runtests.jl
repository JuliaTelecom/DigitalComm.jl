# ---------------------------------------------------- 
# --- Tests  
# ---------------------------------------------------- 
import Printf, Test

# Run DocTests
include("test_doctest.jl")

# Binary sequence 
include("test_genBitSequence.jl");

# Bit Mapping 
include("test_bitMapping.jl");

# Bit Unmapping 
include("test_bitDemapping.jl");
include("test_hardConstellation.jl");

# NRZI Mapping 
include("test_nrzi.jl");

# Symbol demapper 
include("test_symbolDemapper.jl");

# AWGN 
include("test_addNoise.jl");

# Test Waveforms 
include("test_waveforms.jl");

# Test filters 
include("test_filters.jl")
