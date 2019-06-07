# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
println("Tests for Hard constellation");


for alphabet = [4;16;64;256]
	@testset  "$alphabet-QAM" begin 
		# Minimal Tx 
		mcs	  = alphabet;		  # Constellation size 
		n	  = Int(log2(mcs));		  # Bit per symbol 
		N	  = 4096;	  # Numbver of test symbols 
		bitSeq	= genBitSequence(N*n);
		qamSeq	= bitMappingQAM(mcs,bitSeq); 
		# Minimal Rx 
		qamDec	= hardConstellation(mcs,qamSeq);
		# Checking size of decoding vector 
		@test length(qamDec) == N;
		# Checking that bang and classic methods are the same 
		buff  = zeros(Complex{Float64},N);
		hardConstellation!(buff,mcs,qamSeq);
		@test all(qamDec .== buff); 
		# Checking that decoded method is initial method 
		@test all(qamDec .≈ qamSeq)
		# Checking that introduces a small biais does not change voronoi region 
		biais	= 1/sqrt(2/3*(mcs-1)) * 0.1 * (1+1im)
		@test all( hardConstellation(mcs, qamSeq .+ biais) .≈ qamSeq)
	end
end


