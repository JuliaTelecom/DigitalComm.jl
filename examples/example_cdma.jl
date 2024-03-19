using DigitalComm 


# ----------------------------------------------------
# --- Parameters 
# ---------------------------------------------------- 
nbBits  = 32768 
mcs     = 4 
nbUsers = 16
nbActiveUsers = 16

# ----------------------------------------------------
# --- CDMA Generation
# ---------------------------------------------------- 
# --- Binary sequence
bitSeq	      = genBitSequence(nbBits)
# Mapping
qamSeq		  = bitMappingQAM(mcs,bitSeq);
nbSymb        = length(qamSeq) ÷ nbActiveUsers
# --- T/F matrix
qamMat		  = reshape(qamSeq,nbActiveUsers,nbSymb);
# --- Signal
sigId		  = cdmaSigGen(qamMat,nbUsers,:ovsf)


# ----------------------------------------------------
# --- Ideal decoding
# ---------------------------------------------------- 
# --- Decoding 
qamRx = cdmaSigDecode(sigId,nbUsers,:ovsf,1:nbActiveUsers)
sir = getSIR(qamRx,qamMat)
println("SIR between Tx and Rx sequence (no noise) is $sir dB")


# ----------------------------------------------------
# --- Decoding with noise
# ---------------------------------------------------- 
snr = 30 
sigRx,_ = addNoise(sigId,snr)
qamRx = cdmaSigDecode(sigRx,nbUsers,:ovsf,1:nbActiveUsers)
sir = getSIR(qamRx,qamMat)
println("SIR between Tx and Rx sequence ($snr dB additive  noise) is $sir dB")
