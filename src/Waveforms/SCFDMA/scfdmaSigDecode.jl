""" 
---  
Apply SCFDMA demodulator to input signal and returns the T/F QAM matrix
# --- Syntax 
qamRx	= scfdmaSigDecode(sigId,nFFT,nCP,allocSubcarriers,sizeDFT);
# --- Input parameters 
# ---  Input parameters
- sigId	  : SCFDMA signal in time domain [Array{Complex{Float64},nbEch}]
- nFFT	  : FFT size [Int]
- nCp	  : Cyclic prefix size (in samples) [Int]
- allocSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
- sizeDFT	: Post processing size (DFT size)
- doPostProcessing : Do the post processing stage (IDFT): Default 1. In some case, we need the raw data (i.e data before postprocessing stage) such as in channel equalisation. In this case, the post processing should be done in a separate step (with the call of scfdmaPostProcessing function)
# ---  Output parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function scfdmaSigDecode(sigChan,nFFT,nCP,allocSubcarriers,sizeDFT,doPostProcessing=1)
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# --- Getting symbol size
	sizeSymb	   = nFFT + nCP;
	# --- Getting number of symbol to decode
	nbSymbDecode   = Int(floor(length(sigChan)/sizeSymb));
	# --- Initiate Rx T/F matrix
	qamRx		   = zeros(Complex{Float64},length(allocSubcarriers),nbSymbDecode);
	# --- SCFDMA demodulator
	for n = 1  : 1 : nbSymbDecode
		# --- Getting current symbol without CP
		currSymb  = sigChan[sizeSymb*(n-1).+collect(1+nCP:sizeSymb)];
		# --- TF
		fftEl	  = fft(currSymb,1);
		# --- Mapping in output matrix
		qamRx[:,n]		= fftEl[allocSubcarriers];
	end
	if Bool(doPostProcessing)
		qamRx	= scfdmaPostProcessing(qamRx,sizeDFT);
	end
	return qamRx;
end

"""
---  
Apply Post processing for SC-FDMA (i.e DFT post-processing stage). This function can be called if scfdmaSigDecode is called with post processing flag to 0 (if frequency egalisation is done for example)
# --- Syntax 
qamPost	  = scfdmaPostProcessing(qamMat,sizeDFT)
# --- Input parameters 
- qamMat  : T/F matrix after Rx FFT (and before Rx IDFT) [Array{Complex{Float64}},nbSubcarriers,nbSymb]
- sizeDFT	: Precoder bloc size (often 12)
# --- Output parameters 
- qamPost	: T/F matrix after IDFT processing [Array{Complex{Float64},nbSubcarriers,nbSymb}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function scfdmaPostProcessing(qamMat,sizeDFT)
	nbSymb		= size(qamMat,2);
	nbSubcarriers = size(qamMat,1);
	nbCarriers	= nbSubcarriers ÷ sizeDFT;
	qamIDFT		= zeros(Complex{Float64},nbSubcarriers,nbSymb);
	for iN = 1 : 1 : nbSymb 
		for iC = 1 : 1 : nbCarriers
			# --- Apply DFT to current carrier, at current symbol
			qamIDFT[(iC-1)*sizeDFT .+ (1:sizeDFT),iN] = ifft(qamMat[(iC-1)*sizeDFT  .+ (1:sizeDFT),iN]);
		end 
	end
	return qamIDFT
end


"""
---  
Apply SCFDMA demodulator to input signal and returns the T/F QAM matrix
# --- Syntax 
qamRx	= scfdmaSigDecode(sigId,nFFT,nCP,allocSubcarriers,sizeDFT);
# --- Input parameters 
# ---  Input parameters
- sigId	  : SCFDMA signal in time domain [Array{Complex{Float64},nbEch}]
- scfdma  : SCFDMA structure [StrucSCFDMA]
- doPostProcessing : Do the post processing stage (IDFT): Default 1. In some case, we need the raw data (i.e data before postprocessing stage) such as in channel equalisation. In this case, the post processing should be done in a separate step (with the call of scfdmaPostProcessing function)
# ---  Output parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function scfdmaSigDecode(sigChan,scfdma::StrucSCFDMA,doPostProcessing=1)
	return scfdmaSigDecode(sigChan,scfdma.nFFT,scfdma.nCP,scfdma.allocatedSubcarriers,scfdma.sizeDFT,doPostProcessing);
end
