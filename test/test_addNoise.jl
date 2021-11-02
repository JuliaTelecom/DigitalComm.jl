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

@testset "Testing addCFO.jl..." begin 
    # We check with array of ones, to be sure rotation is Ok
    x = ones(2048)
    y = addCFO(x,0,1,0)
    @test length(x) == length(y)
    @test all(x .== y)
    addCFO!(y,x,0,1,0)
    @test all(x .== y)
    
    # Check it works with complex input 
    x = randn(ComplexF64,1024)
    y = addCFO(x,0,1,0)

    # Check that F32 is not promoted to F64 
    x2 = randn(ComplexF32,1024)
    y = addCFO(x2,100,1e3,0)
    @test y isa Vector{ComplexF32}
   
   
    # Check that CFO is applied 
    y = addCFO(x,100,1e3)
    @test y[1] == x[1]
    @test y[2] == x[2] * exp(2im*π*100/1e3)
    @test y[1024] == x[1024] * exp(2im*π*1023*100/1e3)
end

