# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
# Note --> The mapping system is described in bitMapping.jl

println("Tests for Max-Log soft  symbol demapper");
for alphabet = [4;16;64;256]
	@testset  "$alphabet-QAM" begin 
		# Minimal Tx 
		mcs	  = alphabet;		  # Constellation size 
		n	  = Int(log2(mcs));		  # Bit per symbol 
		N	  = 4096;	  # Numbver of test symbols 
		bitSeq	= genBitSequence(N*n);
		qamSeq	= bitMappingQAM(mcs,bitSeq); 
		# Assuming unitary channel 
		channel = ones(Complex{Float64},length(qamSeq));
		# Calculating LLR 
		llr		= symbolDemappingQAM(mcs,qamSeq,channel);
		hardB	= llrToHardBits(llr);
		# Checking size of decoding vector 
		@test length(llr)	== N*n; 
		@test length(hardB) == N*n; 
		# Checking that bang and classic methods are the same 
		buff  = zeros(Float64,N*n);
		symbolDemappingQAM!(buff,mcs,qamSeq,channel);
		buffB = zeros(UInt8,n*N);
		llrToHardBits!(buffB,buff);
		@test all(hardB .== buffB); 
		# Checking that decoded method is initial method 
		@test all( buffB .≈ bitSeq)
	end
end

