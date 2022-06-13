using DigitalComm 


# ----------------------------------------------------
# --- Parameters 
# ---------------------------------------------------- 
nbBits  = 32768 
mcs     = 4 
nbUsers = 16

# ----------------------------------------------------
# --- CDMA Generation
# ---------------------------------------------------- 
# --- Binary sequence
bitSeq	      = genBitSequence(nbBits);
# Mapping
qamSeq		  = bitMappingQAM(mcs,bitSeq);
nbSymb        = length(qamSeq) ÷ nbUsers
# --- T/F matrix
qamMat		  = reshape(qamSeq,nbSymb,nbUsers);
# --- Signal
sigId		  = cdmaSigGen(qamMat,nbUsers,:ovsf)
