""" 
---  
Apply OQAM demapping to incoming OQAM matrix qamMat of size nbSubcarriers x 2nbSymb
# --- Syntax 
qamMat = oqamDemapping(qamMat)
# --- Input parameters
- oqamMat : OQAM matrix (pure real and pure imag. alterns) [Array{Complex{Float64,nbSubcarriers,nbSymb*2}}]
# --- Output parameters
- qamMat : Output complex qam Matrix [Array{Complex{Float64,nbSubcarriers,nbSymb}}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function oqamDemapping(oQamMat)
	# --- Size of incoming matrix 
	(nbSubcarriers,nbOQAM)	  = size(oQamMat);
	nbSymb	  = nbOQAM ÷ 2;
	# Real and imag split:
	qamMat	  = zeros(Complex{Float64},nbSubcarriers,nbSymb);
	# --- Init matrix
	for iS = 1 : 1 : nbSubcarriers
		# Mask to keep real or imag part 
		pattern		= (-1im)^( mod( mod(iS,2) + 1,2));
		pattern2	= (-1im)^( mod( mod(iS,2) + 2,2));
		# 1 // 1i shift 
		scaleReIm	= (1im)^( mod( mod(iS,2) + 1,2));
		scaleReIm2	= (1im)^( mod( mod(iS,2) + 2,2));
		for iN = 1 : 1 : nbSymb
			# OQAM processing 
			qamMat[iS,iN]	= scaleReIm*real(pattern*oQamMat[iS,2*(iN-1)+1])+scaleReIm2*real(pattern2*oQamMat[iS,2*(iN-1)+2]);
		end
	end
	return qamMat;
end

""" 
---  
Demodulate FBMC waveform based on the dual operation of fbmcSigGen

FBMC is parametrized by its FFT size, its cyclic prefix length (in samples) and the allocSubcarriers vector
# Syntax
sigId	= fbmcSigDecode(sigRx,nFFT,K,allocSubcarriers)
# ---  Input parameters
- sigRx	  : Time domain FBMC signal  [Array{Complex{Float64},nbEch}]
- nFFT	  : FFT size [Int]
- K		  : Overlapping factor [Int]
allocSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
# ---  Output parameters
- qamMat  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function fbmcSigDecode(sigRx,nFFT,K,allocatedSubcarriers)
	# Set 4 core for FFT computation
	FFTW.set_num_threads(4)
	# --- Getting symbol size
	sizeSymb	   = nFFT * K;
	# --- Init FBMC filter 
	p			   = getFBMCFilter(K,nFFT);
	# --- Getting number of symbol to decode
	# Here, symbol are QAM symbols (need to demodulate twice of that)
	nbSymbDecode   = Int(floor(length(sigRx)/ nFFT - K + 1/2));
	sBeforeFFT	   = zeros(Complex{Float64},nFFT,2*nbSymbDecode);
	# --- Iterative decoding with overlap and sum undoing 
	for n = 1 : 1 : 2*nbSymbDecode
		# --- Received signal and Rx apodisation
		# Output is a KN vector multiplied by the filter 
		yTmp	= p  .* sigRx[(n-1)*nFFT÷2 .+ collect(1:sizeSymb)];
		# From this convert to N sequence. All samples delayed by nFFT are sumed
		sBeforeFFT[:,n] = sum(reshape(yTmp,nFFT,K),dims=2);
	end
	# FFT 
	oQamMat	= fft(sBeforeFFT,1);
	# --- OQAM to QAM conversion 
	qamMat	= oqamDemapping(oQamMat[allocatedSubcarriers,:]);
end

# MD for input structure 
function fbmcSigDecode(sigChan,fbmc::StrucFBMC)
	return fbmcSigDecode(sigChan,fbmc.nFFT,fbmc.K,fbmc.allocatedSubcarriers);
end
