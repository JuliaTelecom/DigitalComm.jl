""" 
---  
Quadrature Amplitude Modulation (QAM) function
   Apply symbol mapping to a input binary sequence (of size 1xL) with constellation size M.
	Output is a vector (1xN) with N = L / log2(M)
	Conventional gray mapping is used. Output constellation is casted in float, with unitary average power
 Supported constellation
* QPSK
* 16-QAM
* 64-QAM
* 256-QAM
# --- Syntax 
  bitMappingQAM!(qamMat,M,bitSeq)
# --- Input parameters 
- qamMat	: Complex Vector to populate of size length(bitSeq) / log2(M) [Array{Complex{Float64}}]
- M			: Modulation size (i.e from 4 to 256) such as bit per symbol is log2(M) [Int]
- bitSeq	: Binary sequence to be transformed into QPSK symbols [Array{UInt8}]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function bitMappingQAM!(qamMat,M, bitSeq)
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------
	# --- Defining scaling factor
	# Constellation output should have an unitary average power
	# σ^2 = ∑ p(a_i)  ⃒ a_i ⃒ ^2 = 1
	scalingFactor	= sqrt(2/3*(M-1));
	nbSymb			= length(qamMat);
	# ----------------------------------------------------
	# --- Switch on modulation order
	# ----------------------------------------------------
	if M == 4
		## --- QPSK modulator
		# QPSK mapping
		# Odd bit are in I; Even bits are in Q = [I1 Q1 ...]
		# -------------------------#
		#   |      |
		#   0      1
		#  -1     -1
		# 0 -> -1
		# 1 -> 1
		#nbSymb	  = Int64(length(bitSeq)/2);
		#qamMat	  = zeros(Complex{Float64},nbSymb);
		alphabet  = [-1 1];
		for i = 1 : 1 : nbSymb
			## Binary to Symbol conversion
			# Conversion from 2 bit to one symbol (I path):
			classI	  = bitSeq[2(i-1)+1];
			# Conversion from 2 bit to one symbol (Q path)
			classQ	  = bitSeq[2(i-1)+2];
			## Gray encoding
			outI  = alphabet[1+classI];
			outQ  = alphabet[1+classQ];
			# Setting output
			qamMat[i] = (outI + 1im*outQ)/scalingFactor;
		end
	elseif M == 16
		# 16-QAM mapping
		# Odd bit are in I; Even bits are in Q = [I1 Q1 I2 Q2 ...]
		# I1 :> MSB || I2:> LSB
		# -------------------------#
		#   |      |     |     |
		#  10     11    01    00
		#  -3     -1     1     3
		# 0 -> 3
		# 1 -> 1
		# 2 -> -3
		# 3 -> -1
		nbSymb	  = Int64(length(bitSeq)/4);
		alphabet  = [3 1 -3 -1];
		for i = 1 : 1 : nbSymb
			## Binary to Symbol conversion
			# Conversion from 2 bit to one symbol (I path): {0,1,2,3}
			classI	  = 2bitSeq[4(i-1)+1] + bitSeq[4(i-1)+3];
			# Conversion from 2 bit to one symbol (Q path)
			classQ	  = 2bitSeq[4(i-1)+2] + bitSeq[4(i-1)+4];
			## Gray encoding
			outI  = alphabet[1+classI];
			outQ  = alphabet[1+classQ];
			# Setting output
			qamMat[i] = (outI + 1im*outQ)/scalingFactor;
		end
	elseif M == 64
		# 64-QAM mapping
		# Odd bit are in I; Even bits are in Q = [I1 Q1 I2 Q2  I3 Q3 ...]
		# I1 :> MSB || I4:> LSB
		# ------------------------------------------------#
		#   |      |     |     |     |     |      |      |
		#  100    101   111   110   010   011    001    000
		#   -7     -5    -3    -1    +1    +3     +5     +7
		# 0 -> 000 -> +7
		# 1 -> 001 -> +5
		# 2 -> 010 -> +1
		# 3 -> 011 -> +3
		# 4 -> 100 -> -7
		# 5 -> 101 -> -5
		# 6 -> 110 -> -1
		# 7 -> 111 -> -3
		#nbSymb	  = Int64(length(bitSeq)/6);
		#qamMat	  = zeros(Complex{Float64},nbSymb);
		alphabet  = [7 5 1 3 -7 -5 -1 -3];
		for i = 1 : 1 : nbSymb
			## Binary to Symbol conversion
			# Conversion from 2 bit to one symbol (I path): {0,1,2,3}
			classI	  = 4bitSeq[6(i-1)+1] + 2bitSeq[6(i-1)+3] + bitSeq[6(i-1)+5];
			# Conversion from 2 bit to one symbol (Q path)
			classQ	  = 4bitSeq[6(i-1)+2] + 2bitSeq[6(i-1)+4] + bitSeq[6(i-1)+6];
			## Gray encoding
			outI  = alphabet[1+classI];
			outQ  = alphabet[1+classQ];
			# Setting output
			qamMat[i] = (outI + 1im*outQ)/scalingFactor;
		end
	elseif M == 256
		# 256-QAM mapping
		# Inherited from 64-QAM with classic tree genereation
		# -> [0 Sequence(1:end)] [1 Sequence(end:-1:1)]
		# Odd bit are in I; Even bits are in Q = [I1 Q1 I2 Q2  I3 Q3 I4 Q4...]
		# I1 :> MSB || I4:> LSB
		# ------------------------------------------------#
		#	   |       |         |        |        |        |        |        |
		#	0100    0101 .    0111 .   0110 .   0010 .   0011 .   0001 .   0000
		#	 -15     -13 .     -11 .     -9 .     -7 .     -5 .     -3 .     -1
		#	   |       |         |        |        |        |        |        |
		#   1000 .   1001 .   1011 .   1010 .   1110 .   1111 .   1101 .   1100
		#     +1 .     +3 .     +5 .     +7 .     +9 .    +11 .    +13      +15.
		# 0 ->  0000 -> -1
		# 1 ->  0001 -> -3
		# 2 ->  0010 -> -7
		# 3 ->  0011 -> -5
		# 4 ->  0100 -> -15
		# 5 ->  0101 -> -13
		# 6 ->  0110 -> -9
		# 7 ->  0111 -> -11
		# 8 ->  1000 -> +1
		# 9 ->  1001 -> +3
		# 10 -> 1010 -> +7
		# 11 -> 1011 -> +5
		# 12 -> 1100 -> +15
		# 13 -> 1101 -> +13
		# 14 -> 1110 -> +9
		# 15 -> 1111 -> +11
		#nbSymb	  = Int64(length(bitSeq)/8);	#
		#qamMat	  = zeros(Complex{Float64},nbSymb);
		alphabet  = [-1 -3 -7 -5 -15 -13 -9 -11 +1 +3 +7 +5 +15 +13 +9 +11];
		for i = 1 : 1 : nbSymb
			## Binary to Symbol conversion
			# Conversion from 2 bit to one symbol (I path): {0,1,2,3}
			classI	  = 8bitSeq[8(i-1)+1] + 4bitSeq[8(i-1)+3] + 2bitSeq[8(i-1)+5] + bitSeq[8(i-1)+7];
			# Conversion from 2 bit to one symbol (Q path)
			classQ	  = 8bitSeq[8(i-1)+2] + 4bitSeq[8(i-1)+4] + 2bitSeq[8(i-1)+6] + bitSeq[8(i-1)+8];
			## Gray encoding
			outI  = alphabet[1+classI];
			outQ  = alphabet[1+classQ];
			# Setting output
			qamMat[i] = (outI + 1im*outQ)/scalingFactor;
		end
	end
end

# ---------------------------------------------------- 
# --- No allocation function  
# ---------------------------------------------------- 
function bitMappingQAM(M,bitSeq)
	nbBits	= length(bitSeq);
	nbSymb	= Int( nbBits ÷ log2(M)); 
	buffer	= zeros(Complex{Float64},nbSymb);
	bitMappingQAM!(buffer,M,bitSeq);
	return buffer;
end



# ----------------------------------------------------
# --- Multiple dispatch handling
# ----------------------------------------------------
# --- MD: String case for modulation order
function bitMappingQAM!(qamMat,M::String, bitSeq)
	# --- Casting modulation order to int8
	if M == "QPSK" || M == "4-QAM" || M == "QAM-4" || M=="4QAM" || M == "QAM4"
		bitMappingQAM!(qamMat,4,bitSeq);
	elseif M == "16-QAM" || M == "QAM-16" || M=="16QAM" || M == "QAM16"
		bitMappingQAM!(qamMat,16,bitSeq);
	elseif M == "64-QAM" || M == "QAM-64" || M=="64QAM" || M == "QAM64"
		bitMappingQAM!(qamMat,64,bitSeq);
	elseif M == "256-QAM" || M == "QAM-256" || M=="256QAM" || M == "QAM256"
		bitMappingQAM!(qamMat,256,bitSeq);
	end
end
function bitMappingQAM(M::String, bitSeq)
	# --- Casting modulation order to int8
	if M == "QPSK" || M == "4-QAM" || M == "QAM-4" || M=="4QAM" || M == "QAM4"
		bitMappingQAM(4,bitSeq);
	elseif M == "16-QAM" || M == "QAM-16" || M=="16QAM" || M == "QAM16"
		bitMappingQAM(16,bitSeq);
	elseif M == "64-QAM" || M == "QAM-64" || M=="64QAM" || M == "QAM64"
		bitMappingQAM(64,bitSeq);
	elseif M == "256-QAM" || M == "QAM-256" || M=="256QAM" || M == "QAM256"
		bitMappingQAM(256,bitSeq);
	end
end
