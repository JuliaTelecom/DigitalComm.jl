# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
using Statistics 

# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
println("Tests for noise addition");

@testset "Testing addNoise.jl..." begin 
	# Minimal Tx 
	mcs	  = 4;		  # Constellation size 
	n	  = Int(log2(mcs));		  # Bit per symbol 
	N	  = 4096;	  # Numbver of test symbols 
	snrVise	  = 10;
	bitSeq	= genBitSequence(N*n);
	qamSeq	= bitMappingQAM(mcs,bitSeq); 
	# Testing that SNR estimation corresponds 
	for snr = 0:5:30
		qamNoise,= addNoise(qamSeq,snrVise); 
		# Minimal Rx 
		bitDec	= bitDemappingQAM(mcs,qamSeq); 
		snr		= -10*log10(mean(abs2.(qamNoise .- qamSeq)));
		@test abs2(snr - snrVise) < 0.1
	end
	# Testing that bang method work
	snrVise	  = 10;
	mem	  = zeros(Complex{Float64},N);
	addNoise!(mem,qamSeq,snrVise);
	snr		= -10*log10(mean(abs2.(qamSeq .- mem)));
	@test abs2(snr - snrVise) < 0.1




end


