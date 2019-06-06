# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
# Note --> The mapping system is described in bitMapping.jl
println("Tests for symbol mapper with QAM sequences");

@testset  "qpsk" begin 
	# Create a bit squence (already tested) 
	nbBits	= 2 * 2048;
	bitSeq  = genBitSequence(nbBits);
	# Pass trough the function 
	buff	= zeros(Complex{Float64},nbBits÷2)
	# Call 
	bitMappingQAM!(buff,4,bitSeq); 
	buff2	= bitMappingQAM(4,bitSeq);
	@test all( buff .== buff2)
	# --- Good const point
	# For QPSK, we should have |r|^2 = 1
	@test all(abs.(real(buff)) .== 1/sqrt(2)); 
	# Some manual check due to gray coding 
	buff  = [0x01 0x00 0x00 0x00 0x00 0x01 0x01 0x01 ];
	@test all(bitMappingQAM(4,buff) .≈ 1/sqrt(2)*[1-1im; -1-1im; -1+1im; 1+1im]); 
end


@testset "16-QAM" begin 
	# Create a bit squence (already tested) 
	mcs		= 16;
	n		= Int(log2(mcs));
	nbBits	= n * 2048;
	bitSeq  = genBitSequence(nbBits);
	# Pass trough the function 
	buff	= zeros(Complex{Float64},nbBits÷n)
	# Call 
	bitMappingQAM!(buff,mcs,bitSeq); 
	buff2	= bitMappingQAM(mcs,bitSeq);
	@test all( buff .== buff2)
	# Some manual check due to gray coding 
	buff  = [0x01 0x00 0x00 0x00 0x00 0x01 0x01 0x01 ];
	@test all(bitMappingQAM(mcs,buff) .≈ 1/sqrt(10)*[-3+3im;1-1im]);
end


@testset "64-QAM" begin 
	# Create a bit squence (already tested) 
	mcs		= 64;
	n		= Int(log2(mcs));
	nbBits	= n * 2048;
	bitSeq  = genBitSequence(nbBits);
	# Pass trough the function 
	buff	= zeros(Complex{Float64},nbBits÷n)
	# Call 
	bitMappingQAM!(buff,mcs,bitSeq); 
	buff2	= bitMappingQAM(mcs,bitSeq);
	@test all( buff .== buff2)
	# Some manual check due to gray coding 
	buff  = [0x01 0x00 0x00 0x00 0x00 0x01 0x01 0x01 0x01 0x00 0x01 0x00];
	@test all(bitMappingQAM(mcs,buff) .≈ 1/sqrt(42)*[-7+5im;-3-7im]);
end




@testset "256-QAM" begin 
	# Create a bit squence (already tested) 
	mcs		= 256;
	n		= Int(log2(mcs));
	nbBits	= n * 2048;
	bitSeq  = genBitSequence(nbBits);
	# Pass trough the function 
	buff	= zeros(Complex{Float64},nbBits÷n)
	# Call 
	bitMappingQAM!(buff,mcs,bitSeq); 
	buff2	= bitMappingQAM(mcs,bitSeq);
	@test all( buff .== buff2)
	# Some manual check due to gray coding 
	buff  = [0x01 0x00 0x00 0x00 0x00 0x01 0x01 0x01 0x01 0x00 0x01 0x00 0x01 0x01 0x00 0x00];
	@test all(bitMappingQAM(mcs,buff) .≈ 1/sqrt(170)*[3-5im;9-7im]);
end
