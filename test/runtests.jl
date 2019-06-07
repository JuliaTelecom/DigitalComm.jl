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

# AWGN 
include("test_addNoise.jl");
