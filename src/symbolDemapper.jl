
# 64-QAM
"""
---  
Calculate a euclidean distance based on input sample and alphabet distribution.
Used for symbolDemappingQAM function 
# --- Syntax 
	getDist!(out,sR,alphabet, a1...aN)
# --- Input parameters 
- out	  : Float64 array of size N 
- sR	  : Incoming real sequence 
- alphabet : Gray base for decision 
- a1...aN : N neigthbours for soft symbol demapper 
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
# 256-QAM
function getDist!(baseComp,sR,alphabet,a1,a2,a3,a4,a5,a6,a7,a8)
	baseComp[1] = abs2.(sR - alphabet[a1]);
	baseComp[2] = abs2.(sR - alphabet[a2]);
	baseComp[3] = abs2.(sR - alphabet[a3]);
	baseComp[4] = abs2.(sR - alphabet[a4]);
	baseComp[5] = abs2.(sR - alphabet[a5]);
	baseComp[6] = abs2.(sR - alphabet[a6]);
	baseComp[7] = abs2.(sR - alphabet[a7]);
	baseComp[8] = abs2.(sR - alphabet[a8]);
	return minimum(baseComp);
end

# 64-QAM
function getDist!(baseComp,sR,alphabet,a1,a2,a3,a4)
	baseComp[1] = abs2.(sR - alphabet[a1]);
	baseComp[2] = abs2.(sR - alphabet[a2]);
	baseComp[3] = abs2.(sR - alphabet[a3]);
	baseComp[4] = abs2.(sR - alphabet[a4]);
	return minimum(baseComp);
end
# 16-QAM
function getDist!(baseComp,sR,alphabet,a1,a2)
	baseComp[1] = abs2.(sR .- alphabet[a1]);
	baseComp[2] = abs2.(sR .- alphabet[a2]);
	return minimum(baseComp);
end

""" 
---  
Returns the final LLR value based on input distances 
# --- Syntax 
llr = calcLLR(e0,e1,c)
# --- Input parameters 
- e0  : Minimal distance 1 
- e1  : Minimal distance 2 
- c	  : Channel estimate 
# --- Output parameters 
- llr : Max likelihood estimate 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function calcLLR(e0,e1,c)
	return (e0-e1) * c;
end

"""
---
# --- Description
Returns the log likelihood ratio of the incoming sequence qamSeq based on the channel estimates channelIn. Max log approximation is considered
qamSeq is an input noisy QAM sequence with same size of channel estimate vector
Output is a vector of soft output binary sequence to be fed in a FEC 
# --- Syntax 
 output	  = :symbolDemappingQAM(mcs,qamSeq,channel)
# --- Input parameters 
- mcs	  : Constellation size (from 4 to 256) [Int]
- qamSeq  : Complex noisy received sequence (after equalization) [Array{Float64},N]
- channel : Complex channel estimate [Array{Float64},N]
# --- Output parameters 
- output  : Soft bits [Array{UInt8},N*log2(mcs)]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function symbolDemappingQAM(mcs,qamSeq,channel)
	# --- Creating output array for allocation 
	output = zeros(Float64,Int(length(qamSeq) * log2(mcs))); 
	# --- Call bang method 
	symbolDemappingQAM!(output,mcs,qamSeq,channel)
	# --- Return array
	return output;
end

""" 
---
# --- Description
Returns the log likelihood ratio of the incoming sequence qamSeq based on the channel estimates channelIn
qamSeq is an input noisy QAM sequence with same size of channel estimate vector
Output is populated by  soft output binary sequence to be fed in a FEC 
# --- Syntax 
 output	  = :symbolDemappingQAM(mcs,qamSeq,channel)
# --- Input parameters 
- output  : Soft bits [Array{UInt8},N*log2(mcs)]
- mcs	  : Constellation size (from 4 to 256) [Int]
- qamSeq  : Complex noisy received sequence (after equalization) [Array{Float64},N]
- channel : Complex channel estimate [Array{Float64},N]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function symbolDemappingQAM!(output,mcs,qamSeq,channel)
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------
	# --- Defining scaling factor
	# Constellation output should have an unitary average power
	# σ^2 = ∑ p(a_i)  ⃒ a_i ⃒ ^2 = 1
	scalingFactor	= sqrt(2/3*(mcs-1));
	# --- Getting size of input
	nbSymb	= length(qamSeq);
	# ----------------------------------------------------
	## --- Symbol demapper
	# ----------------------------------------------------
	if mcs == 4
		# ----------------------------------------------------
		# --- QPSK scheme
		# ----------------------------------------------------
		# Decision is based on euclidean distance with 1/sqrt(2) and -1/sqrt(2). 
		# We can do a direct form here.
		output[1:2:end] .= ((abs2.(real(qamSeq).+1/scalingFactor) .- abs2.(real(qamSeq).-1/scalingFactor)).*abs2.(channel))./2;
		output[2:2:end] .= ((abs2.(imag(qamSeq).+1/scalingFactor) .- abs2.(imag(qamSeq).-1/scalingFactor)).*abs2.(channel))./2;
	elseif mcs == 16
		# ----------------------------------------------------
		# --- 16-QAM scheme
		# ----------------------------------------------------
		#maximum log approximation
		# -------------------------%
		#   |      |     |     |
		#  10     11    01    00
		#  -3     -1     1     3
		alphabet  = [-3 -1 1 3]/ scalingFactor;
		# Trying to guess if the first bit  is a 0 or a 1
		# See bitMappingQAM: realSymb = 2bitSeq[4(i-1)+1] + bitSeq[4(i-1)+3]
		# * LSB is third bit for real part
		# * MSB is first bit for real part
		# --> 0 if -3 and 3
		# ---> 1 if -1 and 1
		# Compute here the maxlog approx --> minimum of distance
		# --- We create the container for comparison 
		baseComp = zeros(Float64,2);
		# --- Iterative decision 
		for iN = 1 : 1 : nbSymb 
			c	= abs2(channel[iN]);
			# MSB -> Re
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,3,4)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,2)
			output[4*(iN-1)+1]  = calcLLR(e_porteur0,e_porteur1,c);
			# LSB -> Re
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,4)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,2,3)
			output[4*(iN-1)+3]  =calcLLR(e_porteur0,e_porteur1,c);
			# MSB -> Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,3,4)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,2)
			output[4*(iN-1)+2]  =calcLLR(e_porteur0,e_porteur1,c);
			# LSB -> Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,4)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,2,3)
			output[4*(iN-1)+4]  =calcLLR(e_porteur0,e_porteur1,c);
		end
	elseif mcs == 64
		# ----------------------------------------------------
		# --- 64-QAM Soft demapper
		# ----------------------------------------------------
		#maximum log approximation
		# ------------------------------------------------#
		#   |      |     |     |     |     |      |      |
		#  100    101   111   110   010   011    001    000
		#   -7     -5    -3    -1    +1    +3     +5     +7
		alphabet  = [-7 -5 -3 -1 1 3 5 7]/ scalingFactor;
		# --- Creating comparison container 
		baseComp = zeros(Float64,4);
		# --- Iterative LLR computation 
		for iN = 1 : 1 : nbSymb 
			c	= abs2(channel[iN]);
			# --- Real part 
			# MSB0 - Re
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,5,6,7,8)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,2,3,4)
			output[6*(iN-1)+1]  = calcLLR(e_porteur0,e_porteur1,c);
			# MSB1 - Re 
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,2,7,8)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,3,4,5,6)
			output[6*(iN-1)+3]  =calcLLR(e_porteur0,e_porteur1,c);
			# LSB - Re 
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,4,5,8)
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,2,3,6,7)
			output[6*(iN-1)+5]  =calcLLR(e_porteur0,e_porteur1,c);
			# --- Imag part 
			# MSB0 - Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,5,6,7,8);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,2,3,4);
			output[6*(iN-1)+2]  = calcLLR(e_porteur0,e_porteur1,c);
			# MSB1 - Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,2,7,8);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,3,4,5,6);
			output[6*(iN-1)+4]  = calcLLR(e_porteur0,e_porteur1,c);
			# LSB - Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,4,5,8);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,2,3,6,7);
			output[6*(iN-1)+6]  =calcLLR(e_porteur0,e_porteur1,c);
		end
	elseif mcs == 256
		# ----------------------------------------------------
		# --- 256-QAM Soft Demapper
		# ----------------------------------------------------
		#maximum log approximation
		# ------------------------------------------------#
		#	   |       |         |        |        |        |        |        |
		#	0100    0101 .    0111 .   0110 .   0010 .   0011 .   0001 .   0000
		#	 -15     -13 .     -11 .     -9 .     -7 .     -5 .     -3 .     -1
		#	   |       |         |        |        |        |        |        |
		#   1000 .   1001 .   1011 .   1010 .   1110 .   1111 .   1101 .   1100
		#     +1 .     +3 .     +5 .     +7 .     +9 .    +11 .	   +13      +15
		#output    = zeros(Float64,nbSymb*8);
		alphabet  = (-15:2:15)./ scalingFactor;
		# --- Creating comparison container 
		baseComp = zeros(Float64,8);
		# --- Iterative LLR computation 
		for iN = 1 : 1 : nbSymb 
			c	= abs2(channel[iN]);
			# MSB0 - Re
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,2,3,4,5,6,7,8);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,9,10,11,12,13,14,15,16);
			output[8*(iN-1)+1]  = calcLLR(e_porteur0,e_porteur1,c);
			# MSB1 - Re
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,5,6,7,8,9,10,11,12);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,2,3,4,13,14,15,16);
			output[8*(iN-1)+3]  = calcLLR(e_porteur0,e_porteur1,c);
			# MSB2 - Re  
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,2,7,8,9,10,15,16);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,3,4,5,6,11,12,13,14);
			output[8*(iN-1)+5]  = calcLLR(e_porteur0,e_porteur1,c);
			# LSB - Re 
			e_porteur0 = getDist!(baseComp,qamSeq[iN].re,alphabet,1,4,5,8,9,12,13,16);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].re,alphabet,2,3,6,7,10,11,14,15);
			output[8*(iN-1)+7]  = calcLLR(e_porteur0,e_porteur1,c);
			# --- Imag part config
			# MSB0 - Im 
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,2,3,4,5,6,7,8);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,9,10,11,12,13,14,15,16);
			output[8*(iN-1)+2]  = calcLLR(e_porteur0,e_porteur1,c);
			# MSB1 - Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,5,6,7,8,9,10,11,12);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,2,3,4,13,14,15,16);
			output[8*(iN-1)+4]  = calcLLR(e_porteur0,e_porteur1,c);
			# MSB2 - Im 
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,2,7,8,9,10,15,16);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,3,4,5,6,11,12,13,14);
			output[8*(iN-1)+6]  = calcLLR(e_porteur0,e_porteur1,c);
			# LSB - Im
			e_porteur0 = getDist!(baseComp,qamSeq[iN].im,alphabet,1,4,5,8,9,12,13,16);
			e_porteur1 = getDist!(baseComp,qamSeq[iN].im,alphabet,2,3,6,7,10,11,14,15);
			output[8*(iN-1)+8]  = calcLLR(e_porteur0,e_porteur1,c);
		end
	end
	return output
end


# ----------------------------------------------------
# --- Convert LLR to Unsigned Char Array (soft decoder in C)
# ----------------------------------------------------
"""
---  
Convert a LLR floating array to a UInt LLR 
Some method (espcially C functions as in libfec) expect LLR to be UInt value from  0 (Likely a 0) to 255 (Likely a 1). 
symbolDemappingQAM create a Float64 array from -infty (likely a 0) to nfty (Likely a 1). This function take a Float64 LLR estimate and output a UInt8 vector
# --- Syntax 
	llrToUInt!(llrUInt,llr)
# --- Input parameters 
- llrUInt : Output LLR in UInt [Array{UInt},N]
- llr	  : Input LLR [Array{Float64},N]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function llrToUInt!(llrUInt,llr)
    # 0 ---> Likely a 0
    # 255 ---> Likely a 1
    # We have the llRin: >0 1 <0 -1, higher the better reliability.
    # We assume that 1  --> certainly a 1 (i.e 255)
    #                -1 --> certainly a 0 (i.e 0)
    # We need to convert LLR into UInt8: LLR is converted between 0 and 1 and map to 256 values
    # 0 is a 0 with high probability
    # 1 is a 1 with high probability
    # --- Minimal and maximal bounds (-1 +1)
    # Then, Switch to [0;2]
    # Then Convert to [0,1]
	# Finally extend to 255 
	llrUInt .= UInt8.( round.(255* (1 .+ map(x->min(x,1),map(x->max(x,-1),llr)))/2));
end

"""
---  
Convert a LLR floating array to a UInt LLR 
Some method (espcially C functions as in libfec) expect LLR to be UInt value from  0 (Likely a 0) to 255 (Likely a 1). 
symbolDemappingQAM create a Float64 array from -infty to infty. This function take a Float64 LLR estimate and output a UInt8 vector
# --- Syntax 
	llrUInt = llrToUInt(llr)
# --- Input parameters 
- llr	  : Input LLR [Array{Float64},N]
# --- Output parameters 
- llrUInt : Output LLR in UInt [Array{UInt},N]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function llrToUInt(llr)
	# --- Create output vector 
	llrUInt = zeros(UInt8,length(llr));
	# --- Call bang method 
	llrToUInt!(llrUInt,llr);
	# --- Return LLR 
    return llrUInt;
end


# ----------------------------------------------------
# --- Hard decoding from Soft LLRs
# ----------------------------------------------------
""" 
---  
Returns hard binary value from soft LLR estimate (no FEC decoder, only hard decision here!)
# --- Syntax 
	llrToHardBits(hardD,llr)
# --- Input parameters 
- hardD	  : Hard binary decision [Array{UInt},N] 
- llr	  : Input LLR (Float64 or UInt8 array) [Union{[Array{UInt}],[Array{Float64}]}]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function llrToHardBits!(hardD,llr::Array{Float64})
	# --- Direct LLR are given, sign is bit information
	hardD   .= map(x-> (x<0) ? 0 : 1,llr);
	return hardD
end

"""
---  
Returns hard binary value from soft LLR estimate (no FEC decoder, only hard decision here!)
# --- Syntax 
	llrToHardBits(hardD,llr)
# --- Input parameters 
- hardD	  : Hard binary decision [Array{UInt},N] 
- llr	  : Input LLR (Float64 or UInt8 array) [Union{[Array{UInt}],[Array{Float64}]}]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function llrToHardBits!(hardD,llr::Array{UInt8})
	# --- UInt LLR (input of soft decoder): rule is 127
	hardD   .= map(x-> (x<0x7f) ? 0 : 1,llr);
	return hardD
end

"""
---  
Returns hard binary value from soft LLR estimate (no FEC decoder, only hard decision here!)
# --- Syntax 
	hardD = llrToHardBits(llr)
# --- Input parameters 
- llr	  : Input LLR (Float64 or UInt8 array) [Union{[Array{UInt}],[Array{Float64}]}]
# --- Output parameters 
- hardD	  : Hard binary decision [Array{UInt},N] 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function llrToHardBits(llr::Array{Float64})
	# --- Direct LLR are given, sign is bit information
	hardD   = map(x-> (x<0) ? 0 : 1,llr);
	return hardD
end

""" 
---  
Returns hard binary value from soft LLR estimate (no FEC decoder, only hard decision here!)
# --- Syntax 
	hardD = llrToHardBits(llr)
# --- Input parameters 
- llr	  : Input LLR (Float64 or UInt8 array) [Union{[Array{UInt}],[Array{Float64}]}]
# --- Output parameters 
- hardD	  : Hard binary decision [Array{UInt},N] 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function llrToHardBits(llr::Array{UInt8})
	# --- UInt LLR (input of soft decoder): rule is 127
	hardD   = map(x-> (x<0x7f) ? 0 : 1,llr);
	return hardD
end

