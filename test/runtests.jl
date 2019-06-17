# ---------------------------------------------------- 
# --- Tests  
# ---------------------------------------------------- 
import Printf
import Test

# Binary sequence 
include("test_genBitSequence.jl");


# Bit Mapping 
include("test_bitMapping.jl");

# Bit Unmapping 
include("test_bitDemapping.jl");
include("test_hardConstellation.jl");

# Symbol demapper 
include("test_symbolDemapper.jl");


# AWGN 
include("test_addNoise.jl");

# Test Waveforms 
include("test_waveforms.jl");
