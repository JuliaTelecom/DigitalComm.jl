# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
# Note --> The mapping system is described in bitMapping.jl

println("Tests for Hard symbol demapper");
for alphabet = [2;4;16;64;256]
	@testset  "$alphabet-QAM" begin 
		# Minimal Tx 
		mcs	  = alphabet;		  # Constellation size 
		n	  = Int(log2(mcs));		  # Bit per symbol 
		N	  = 4096;	  # Numbver of test symbols 
		bitSeq	= genBitSequence(N*n);
		qamSeq	= bitMappingQAM(mcs,bitSeq); 
		# Minimal Rx 
		bitDec	= bitDemappingQAM(mcs,qamSeq); 
		# Checking size of decoding vector 
		@test length(bitDec) == N*n; 
		# Checking that bang and classic methods are the same 
		buff  = zeros(UInt8,N*n);
		bitDemappingQAM!(buff,mcs,qamSeq);
		@test all(bitDec .== buff); 
		# Checking that decoded method is initial method 
		@test all( bitDec .≈ bitSeq)
	end
end

