""" 
---  
Quadrature Amplitude Modulation (QAM) hard decoding function
Apply symbol hard demapping to a input symbol sequence (of size 1xN) with constellation size M.
Output is a binary (1xL) with N = L / log2(M)
Conventional gray demapping is used.
Input constellation is Array{Complex{Float64}}
Output bitsream is Array{Int8}
# --- Syntax 
	  bitDemappingQAM!(hardBits,M,qamVect)
# --- Input parameters 
- hardBits	: Vector of bits to populate [Array{UInt8}, length(qamVect)/log2(M)]
- M			: Constellation size (i.e from 4 to 256)
- qamVect	: Complex observation vector to decode.
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function bitDemappingQAM!(hardBits,M, qamVect)
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------
	# --- Bit sequence and initialisation
	nbBits	= Int( length(qamVect) * log2(M));
	# ----------------------------------------------------
	# --- Switch on modulation order
	# ----------------------------------------------------
	# Here all Voronoi region are defined with normalized constellation
	# Symbol demapper are recall from bitMappingQAM (without rescaling factor)
	if M == 2 
		## BPSK demodulator 
		hardBits .= Int.( real(qamVect) .> 0);
	elseif M == 4
		## --- QPSK demodulator
		# --- Gray encoding
		#  -1  1
		#   0  1
		# --- Voronoi region based on sign
		# < 0 --> -1
		# > 0 --> 1
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamVect)
			hardBits[(n-1)*2+1] = Int8(real(qamVect[n]) > 0);
			hardBits[(n-1)*2+2] = Int8(imag(qamVect[n]) > 0);
		end
	elseif M == 16
		## 16-QAM hard demapper
		# --- Gray encoding
		#  -3     -1     1     3
		#  10     11    01    00
		# --- Voronoi region definition
		bounds	= [-2 0 2]/(sqrt(2/3*15));
		#bounds	= [-sqrt(2)/2 0 sqrt(2)/2];
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamVect)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamVect[n]);
			# Getting area of interest: from area 1 to 4
			decI	= 1+ (realI > bounds[1]) + (realI>bounds[2]) + (realI>bounds[3]);
			# Convert into binary sequence based on Gray encoding scheme.
			if decI == 1
				hardBits[4(n-1)+1] = 1;
				hardBits[4(n-1)+3] = 0;
			elseif decI == 2
				hardBits[4(n-1)+1] = 1;
				hardBits[4(n-1)+3] = 1;
			elseif decI == 3
				hardBits[4(n-1)+1] = 0;
				hardBits[4(n-1)+3] = 1;
			elseif decI == 4
				hardBits[4(n-1)+1] = 0;
				hardBits[4(n-1)+3] = 0;
			end
			# --- Decision for imag part
			# Getting only imag part
			imagP	= imag(qamVect[n]);
			# Getting area of interest: from area 1 to 4
			decQ	= 1+ (imagP > bounds[1]) + (imagP>bounds[2]) + (imagP>bounds[3]);
			# Convert into binary sequence based on Gray encoding scheme.
			if decQ == 1
				hardBits[4(n-1)+2] = 1;
				hardBits[4(n-1)+4] = 0;
			elseif decQ == 2
				hardBits[4(n-1)+2] = 1;
				hardBits[4(n-1)+4] = 1;
			elseif decQ == 3
				hardBits[4(n-1)+2] = 0;
				hardBits[4(n-1)+4] = 1;
			elseif decQ == 4
				hardBits[4(n-1)+2] = 0;
				hardBits[4(n-1)+4] = 0;
			end
		end
	elseif M == 64
		## 64-QAM hard demapper
		# --- Gray encoding
		#   -7     -5    -3    -1    +1    +3     +5     +7
		#  100    101   111   110   010   011    001    000
		# --- Voronoi region definition
		bounds = [-6 -4 -2 0 2 4 6]/(sqrt(2/3*63));
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamVect)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamVect[n]);
			# Getting area of interest: from area 1 to 4
			decI	= 1+ (realI > bounds[1]) + (realI>bounds[2]) + (realI>bounds[3])+ (realI>bounds[4])+ (realI>bounds[5])+ (realI>bounds[6])+ (realI>bounds[7]);
			# Convert into binary sequence based on Gray encoding scheme.
			if decI == 1
				hardBits[6(n-1)+1] = 1;
				hardBits[6(n-1)+3] = 0;
				hardBits[6(n-1)+5] = 0;
			elseif decI == 2
				hardBits[6(n-1)+1] = 1;
				hardBits[6(n-1)+3] = 0;
				hardBits[6(n-1)+5] = 1;
			elseif decI == 3
				hardBits[6(n-1)+1] = 1;
				hardBits[6(n-1)+3] = 1;
				hardBits[6(n-1)+5] = 1;
			elseif decI == 4
				hardBits[6(n-1)+1] = 1;
				hardBits[6(n-1)+3] = 1;
				hardBits[6(n-1)+5] = 0;
			elseif decI == 5
				hardBits[6(n-1)+1] = 0;
				hardBits[6(n-1)+3] = 1;
				hardBits[6(n-1)+5] = 0;
			elseif decI == 6
				hardBits[6(n-1)+1] = 0;
				hardBits[6(n-1)+3] = 1;
				hardBits[6(n-1)+5] = 1;
			elseif decI == 7
				hardBits[6(n-1)+1] = 0;
				hardBits[6(n-1)+3] = 0;
				hardBits[6(n-1)+5] = 1;
			elseif decI == 8
				hardBits[6(n-1)+1] = 0;
				hardBits[6(n-1)+3] = 0;
				hardBits[6(n-1)+5] = 0;
			end
			# --- Decision for imag part
			# Getting only imag part
			realQ	= imag(qamVect[n]);
			# Getting area of interest: from area 1 to 8
			decQ	= 1+ (realQ > bounds[1]) + (realQ>bounds[2]) + (realQ>bounds[3])+ (realQ>bounds[4])+ (realQ>bounds[5])+ (realQ>bounds[6])+(realQ>bounds[7]);
			# Convert into binary sequence based on Gray encoding scheme.
			if decQ == 1
				hardBits[6(n-1)+2] = 1;
				hardBits[6(n-1)+4] = 0;
				hardBits[6(n-1)+6] = 0;
			elseif decQ == 2
				hardBits[6(n-1)+2] = 1;
				hardBits[6(n-1)+4] = 0;
				hardBits[6(n-1)+6] = 1;
			elseif decQ == 3
				hardBits[6(n-1)+2] = 1;
				hardBits[6(n-1)+4] = 1;
				hardBits[6(n-1)+6] = 1;
			elseif decQ == 4
				hardBits[6(n-1)+2] = 1;
				hardBits[6(n-1)+4] = 1;
				hardBits[6(n-1)+6] = 0;
			elseif decQ == 5
				hardBits[6(n-1)+2] = 0;
				hardBits[6(n-1)+4] = 1;
				hardBits[6(n-1)+6] = 0;
			elseif decQ == 6
				hardBits[6(n-1)+2] = 0;
				hardBits[6(n-1)+4] = 1;
				hardBits[6(n-1)+6] = 1;
			elseif decQ == 7
				hardBits[6(n-1)+2] = 0;
				hardBits[6(n-1)+4] = 0;
				hardBits[6(n-1)+6] = 1;
			elseif decQ == 8
				hardBits[6(n-1)+2] = 0;
				hardBits[6(n-1)+4] = 0;
				hardBits[6(n-1)+6] = 0;
			end
		end
	elseif M == 256
		## 64-QAM hard demapper
		# --- Gray encoding
		# ------------------------------------------------#
		#	   |       |         |        |        |        |        |        |
		#	 -15     -13 .     -11 .     -9 .     -7 .     -5 .     -3 .     -1
		#	0100    0101 .    0111 .   0110 .   0010 .   0011 .   0001 .   0000
		#	   |       |         |        |        |        |        |        |
		#     +1 .     +3 .     +5 .     +7 .     +9 .    +11 .    +13      +15.
		#   1000 .   1001 .   1011 .   1010 .   1110 .   1111 .   1101 .   1100
		# --- Voronoi region definition
		bounds = [-14 -12 -10 -8 -6 -4 -2 0 2 4 6 8 10 12 14]/(sqrt(2/3*255));
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamVect)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamVect[n]);
			# Getting area of interest: from area 1 to 15
			decI	= 1+ (realI > bounds[1]) + (realI>bounds[2]) + (realI>bounds[3])+ (realI>bounds[4])+ (realI>bounds[5])+ (realI>bounds[6])+ (realI>bounds[7])+ (realI>bounds[8])+(realI>bounds[9])+(realI>bounds[10])+(realI>bounds[11])+(realI>bounds[12])+(realI>bounds[13])+ (realI>bounds[14])+ (realI>bounds[15]) ;
			# Convert into binary sequence based on Gray encoding scheme.
			if decI == 1
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 2
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 3
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 4
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 5
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 6
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 7
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 8
				hardBits[8(n-1)+1] = 0;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 9
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 10
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 11
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 12
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 0;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 13
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 0;
			elseif decI == 14
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 1;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 15
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 1;
			elseif decI == 16
				hardBits[8(n-1)+1] = 1;
				hardBits[8(n-1)+3] = 1;
				hardBits[8(n-1)+5] = 0;
				hardBits[8(n-1)+7] = 0;
			end
			# --- Decision for imag part
			# Getting only imag part
			realQ	= imag(qamVect[n]);
			# Getting area of interest: from area 1 to 15
			decQ	= 1+ (realQ > bounds[1]) + (realQ>bounds[2]) + (realQ>bounds[3])+ (realQ>bounds[4])+ (realQ>bounds[5])+ (realQ>bounds[6])+ (realQ>bounds[7])+ (realQ>bounds[8])+(realQ>bounds[9])+(realQ>bounds[10])+(realQ>bounds[11])+(realQ>bounds[12])+(realQ>bounds[13])+ (realQ>bounds[14])+ (realQ>bounds[15]) ;
			# Convert into binary sequence based on Gray encoding scheme.
			if decQ == 1
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 2
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 3
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 4
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 5
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 6
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 7
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 8
				hardBits[8(n-1)+2] = 0;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 9
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 10
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 11
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 12
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 0;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 13
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 0;
			elseif decQ == 14
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 1;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 15
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 1;
			elseif decQ == 16
				hardBits[8(n-1)+2] = 1;
				hardBits[8(n-1)+4] = 1;
				hardBits[8(n-1)+6] = 0;
				hardBits[8(n-1)+8] = 0;
			end
		end
	end
	# ---  Output on decoded sequence
	return hardBits;
end



function bitDemappingQAM(M,qamVect);
	# --- Create receive bit vector 
	hardBits = zeros(UInt8, Int(length(qamVect)*log2(M))); 
	# --- Call bang method 
	bitDemappingQAM!(hardBits,M,qamVect);
	return hardBits;
end
