# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
println("Tests for binary sequence generation");

@testset "genBitSequence" begin 
	# --- Check type and so on
	N		= 1000;
	rSeed	= 12;
	bitSeq	= genBitSequence(N,rSeed); 
	# --- Ensure that output is a UInt8 vector 
	@test typeof(bitSeq) == Array{UInt8,1}; 
	# --- Checking 0 and 1 populate the buffer 
	# Not elegant, commutative presence ? 
	@test unique(bitSeq) == [0x00,0x01] || unique(bitSeq) == [0x01,0x00] ;
	# --- Checking bang method 
	bitSeq2	= zeros(UInt8,N); 
	genBitSequence!(bitSeq2,N,rSeed); 
	# --- Ensure that output is a UInt8 vector 
	@test typeof(bitSeq2) == Array{UInt8,1}; 
	# --- Checking 0 and 1 populate the buffer 
	# Not elegant, commutative presence ? 
	@test unique(bitSeq2) == [0x00,0x01] || unique(bitSeq) == [0x01,0x00] ;
	# --- Seed should have generated same vector for both call 
	@test (bitSeq == bitSeq2)
	# --- Forcing a random seed and check taht vector is different 
	@test (genBitSequence(N) != bitSeq)
end


@testset "genByteSequence" begin 
	# --- Check type and so on
	N		= 1000;
	rSeed	= 12;
	byteSeq	= genByteSequence(N,rSeed); 
	# --- Ensure that output is a UInt8 vector 
	@test typeof(byteSeq) == Array{UInt8,1}; 
	# --- Checking bang method 
	byteSeq2	= zeros(UInt8,N); 
	genByteSequence!(byteSeq2,N,rSeed); 
	# --- Ensure that output is a UInt8 vector 
	@test typeof(byteSeq2) == Array{UInt8,1}; 
	# --- Seed should have generated same vector for both call 
	@test (byteSeq == byteSeq2)
	# --- Forcing a random seed and check taht vector is different 
	@test (genByteSequence(N) != byteSeq)
end
