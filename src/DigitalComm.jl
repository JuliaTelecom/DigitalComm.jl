# ---------------------------------------------------- 
# --- Digital Communication module
# ---------------------------------------------------- 
__precompile__()


module DigitalComm
# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
# Only submodules are considered so don't need to overload packages here 
# --> Go to next section with submodules loading 


# ---------------------------------------------------- 
# --- Submodules inclusion  
# ---------------------------------------------------- 
# --- Binary managment 
include("genBitSequence.jl");
# Export 
export genBitSequence!, genBitSequence;


end
