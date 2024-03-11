# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
# Note --> The mapping system is described in bitMapping.jl
println("Tests for symbol mapper with NRZI sequences");

@testset  "NRZI" begin 
	# Create a bit squence (already tested) 
	nbBits	= 2 * 2048;
	bitSeq  = genBitSequence(nbBits);
	# Pass trough the function 
	buff	= zeros(Complex{Float64},nbBits)
	# Call 
	encodeNRZI!(buff,bitSeq); 
	buff2	= encodeNRZI(bitSeq);
	@test all( buff .== buff2)
        # Ensure Tx // Rx is Ok 
        @test all(bitSeq .== decodeNRZI(encodeNRZI(bitSeq)))
        @test all(bitSeq .== decodeNRZI(encodeNRZI(bitSeq,:high),:high))
	# Some manual check for both transitions
	buff  = [0x01;0x00;0x00;0x01;0x00;0x00;0x00;0x01 ];
        @test all( encodeNRZI(buff,:low) .== [0;1;0;0;1;0;1;1])  # Transitions on 0 
        @test all( encodeNRZI(buff,:high) .== [1;1;1;0;0;0;0;1]) # Transitions on 1
end
