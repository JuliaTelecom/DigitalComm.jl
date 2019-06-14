"""
---  
Quadrature Amplitude Modulation (QAM) hard decoding function
Return the hard decoded constellation with voronoi baseds decision. The difference with bitDeMapping is that bitDeMapping returns the decoded bit sequence whereas hardConstellation returns the closest constellation point. This can be use to compute raw EVM estimation (assuming a sufficiently high SNR to avoid errors).
# --- Syntax 
	hardConstellation!(qamDec,M,qamMat)
# --- Input parameters 
- qamDec	: Vector to populate [Array{Complex{Float64},N}] with N = length(qamMat)
- M			: Constellation size (i.e 4 to 256) 
- qamMat	: Vector to decode
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function hardConstellation!(qamThres,M, qamMat)
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------
	scalingFactor	= sqrt(2/3*(M-1));
	#qamThres = zeros(Complex{Float64},length(qamMat));
	# ----------------------------------------------------
	# --- Switch on modulation order
	# ----------------------------------------------------
	# Here all Voronoi region are defined with normalized constellation
	# Symbol remapper are recall from bitMappingQAM (without rescaling factor)
	if M == 4
		## --- QPSK demodulator
		# --- Gray encoding
		#  -1  1
		#   0  1
		# --- Voronoi region based on sign
		# < 0 --> -1
		# > 0 --> 1
		# --- Iterative hard demapping
		alphabet  = [-1 1];
		bounds	  = 0;
		for n = 1 : 1 : length(qamMat)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamMat[n]);
			# Getting area of interest: from area 1 to 4
			decI	= 1+ (realI > bounds[1]);
			# --- Decision for imag part
			# Getting only imag part
			imagP	= imag(qamMat[n]);
			# Getting area of interest: from area 1 to 4
			decQ	= 1+ (imagP > bounds[1]) ;
			# --- TO hard constellation
			qamThres[n] = alphabet[decI] + 1im*alphabet[decQ];
		end
	elseif M == 16
		## 16-QAM hard demapper
		# --- Gray encoding
		#  -3     -1     1     3
		#  10     11    01    00
		# --- Voronoi region definition
		bounds	  = [-2 0 2]/(sqrt(2/3*15));
		alphabet  = [-3 -1 1 3];
		#bounds	= [-sqrt(2)/2 0 sqrt(2)/2];
		# FIXME Voronoi bug in decision
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamMat)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamMat[n]);
			# Getting area of interest: from area 1 to 4
			decI	= 1+ (realI > bounds[1]) + (realI>bounds[2]) + (realI>bounds[3]);
			# --- Decision for imag part
			# Getting only imag part
			imagP	= imag(qamMat[n]);
			# Getting area of interest: from area 1 to 4
			decQ	= 1+ (imagP > bounds[1]) + (imagP>bounds[2]) + (imagP>bounds[3]);
			# --- TO hard constellation
			qamThres[n] = alphabet[decI] + 1im*alphabet[decQ];
		end
	elseif M == 64
		## 64-QAM hard demapper
		# --- Gray encoding
		#   -7     -5    -3    -1    +1    +3     +5     +7
		#  100    101   111   110   010   011    001    000
		# --- Voronoi region definition
		bounds = [-6 -4 -2 0 2 4 6]/(sqrt(2/3*63));
		alphabet  = [-7 -5 -3 -1 1 3 5 7];
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamMat)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamMat[n]);
			# Getting area of interest: from area 1 to 4
			decI	= 1+ (realI > bounds[1]) + (realI>bounds[2]) + (realI>bounds[3])+ (realI>bounds[4])+ (realI>bounds[5])+ (realI>bounds[6])+ (realI>bounds[7]);
			# --- Decision for imag part
			# Getting only imag part
			realQ	= imag(qamMat[n]);
			# Getting area of interest: from area 1 to 8
			decQ	= 1+ (realQ > bounds[1]) + (realQ>bounds[2]) + (realQ>bounds[3])+ (realQ>bounds[4])+ (realQ>bounds[5])+ (realQ>bounds[6])+(realQ>bounds[7]);
			# --- TO hard constellation
			qamThres[n] = alphabet[decI] + 1im*alphabet[decQ];
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
		alphabet  = collect(-15:2:15);
		# --- Iterative hard demapping
		for n = 1 : 1 : length(qamMat)
			# --- Decision for real part
			# Getting only real part
			realI	= real(qamMat[n]);
			# Getting area of interest: from area 1 to 15
			decI	= 1+ (realI > bounds[1]) + (realI>bounds[2]) + (realI>bounds[3])+ (realI>bounds[4])+ (realI>bounds[5])+ (realI>bounds[6])+ (realI>bounds[7])+ (realI>bounds[8])+(realI>bounds[9])+(realI>bounds[10])+(realI>bounds[11])+(realI>bounds[12])+(realI>bounds[13])+ (realI>bounds[14])+ (realI>bounds[15]) ;
			# --- Decision for imag part
			# Getting only imag part
			realQ	= imag(qamMat[n]);
			# Getting area of interest: from area 1 to 15
			decQ	= 1+ (realQ > bounds[1]) + (realQ>bounds[2]) + (realQ>bounds[3])+ (realQ>bounds[4])+ (realQ>bounds[5])+ (realQ>bounds[6])+ (realQ>bounds[7])+ (realQ>bounds[8])+(realQ>bounds[9])+(realQ>bounds[10])+(realQ>bounds[11])+(realQ>bounds[12])+(realQ>bounds[13])+ (realQ>bounds[14])+ (realQ>bounds[15]) ;
			# --- TO hard constellation
			qamThres[n] = alphabet[decI] + 1im*alphabet[decQ];
		end
	end
	qamThres .= qamThres ./ scalingFactor;
	# ---  Output on decoded sequence
	return qamThres;
end

function hardConstellation(M,qamMat)
	qamThres  = zeros(Complex{Float64},length(qamMat));
	hardConstellation!(qamThres,M,qamMat);
	return qamThres;
end

# ----------------------------------------------------
# --- Multiple dispatch handling
# ----------------------------------------------------
# --- MD: String case for modulation order
function hardConstellation!(qamThres,M::String, qamMat)
	# --- Casting modulation order to int8
	if M == "QPSK" || M == "4-QAM" || M == "QAM-4" || M=="4QAM" || M == "QAM4"
		hardConstellation!(qamThres,4,qamMat);
	elseif M == "16-QAM" || M == "QAM-16" || M=="16QAM" || M == "QAM16"
		hardConstellation!(qamThres,16,qamMat);
	elseif M == "64-QAM" || M == "QAM-64" || M=="64QAM" || M == "QAM64"
		hardConstellation!(qamThres,64,qamMat);
	end
end
function hardConstellation!(M::String, qamMat)
	# --- Casting modulation order to int8
	if M == "QPSK" || M == "4-QAM" || M == "QAM-4" || M=="4QAM" || M == "QAM4"
		hardConstellation(4,qamMat);
	elseif M == "16-QAM" || M == "QAM-16" || M=="16QAM" || M == "QAM16"
		hardConstellation(16,qamMat);
	elseif M == "64-QAM" || M == "QAM-64" || M=="64QAM" || M == "QAM64"
		hardConstellation(64,qamMat);
	end

end
