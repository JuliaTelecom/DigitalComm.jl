""" 
---  
Structure for SCFDMA
# --- Syntax 
- nFFT		: FFT size [Int] 
- nCP		: Cyclic prefix size [Int] 
- allocatedSubcarriers	: Vector of allocated subbcarriers [Array{Int}]
- sizeDFT	: Size of DFT precoder (classic value is 12) [Int]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct StrucSCFDMA<: Waveform
	nFFT::Int;
	nCP::Int;
	allocatedSubcarriers::Array{Int};
	sizeDFT::Int;
end


"""
---  
Create SCFDMA structure 
# --- Syntax 
ofdm	= initSCFDMA(nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- nFFT	  : FFT size [Int]
- nCP	  : Cyclic prefix size [Int]
- allocatedSubcarrier	: Vector of allocated subbcarriers [Array{Int}]
- sizeDFT	: Size of DFT precoder (classic value is 12) [Int]
# --- Output parameters 
- scfdma: SCFDMA structure [StrucSCFDMA] 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function initSCFDMA(nFFT,nCP,allocatedSubcarriers,sizeDFT)
	# ---Checking FFT size
	if maximum(allocatedSubcarriers) > nFFT
		error("Subcarrier allocation is impossible");
	end
	# --- Create the SCFDMA structure
	return StrucSCFDMA(nFFT,nCP,allocatedSubcarriers,sizeDFT)
end

"""
---  
Create SCFDMA signal in time domain based on input T/F matrix and SCFDMA parameters. SCFDMA apply a precoder before the IFFT at the transmitter side in order to lower the signal fluctuation

qamMat is a complex symbol matrix (for instance QPSK) of size length(allocatedSubcarriers) x nbSymb 

The output signal in time domain is of size (nFFT+nCP)xnbSymb.
# --- Syntax 
  sigId = scfdmaSigGen(qamMat,nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- qamMat  : Complex T/F matrix to map [Array{Float64},length(allocatedSubcarriers),nbSymb]
- nFFT	  : FFT size [Int] 
- nCP	  : Cylic prefix size [Int] 
- allocatedSubcarrier	: Vector of allocated subbcarriers [Array{Int}]
- sizeDFT	: Size of DFT precoder (classic value is 12) [Int]
# --- Output parameters 
- sigId	  : Signal in time domain [Array{Complex{Float64}},(nFFT+nCP)xnbSymb]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function scfdmaSigGen(qamMat,nFFT,nCP,allocatedSubcarriers,sizeDFT)
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# --- Getting parameters
	nbSymb		= size(qamMat,2); 		# --- Applied on x symbols
	symbSize	= nFFT + nCP; 			# --- Symbol size
	nSize 		= symbSize*nbSymb; 		# --- Total burst length
	# --- DFT pre-processing  
	nbCarriers	= length(allocatedSubcarriers)÷sizeDFT;
	qamDFT		= zeros(Complex{Float64},length(allocatedSubcarriers),nbSymb);
	for iN = 1 : 1 : nbSymb 
		for iC = 1 : 1 : nbCarriers
			# --- Apply DFT to current carrier, at current symbol
			qamDFT[(iC-1)*sizeDFT .+ (1:sizeDFT),iN] = fft(qamMat[(iC-1)*sizeDFT  .+ (1:sizeDFT),iN]);
		end 
	end
	# --- Mapping to nFFT elements
	qamCurr		= zeros(Complex{Float64},nFFT,nbSymb);
	qamCurr[allocatedSubcarriers,:] = qamDFT;
	# --- Switch to time domain
	sTmp 	= ifft(qamCurr,1);
	# --- Inserting Cyclic prefix
	sig	  = [sTmp[end+1-nCP:end,:] ; sTmp][:];
end

""" 
---  
Create SCFDMA signal in time domain based on input T/F matrix and SCFDMA parameters. SCFDMA apply a precoder before the IFFT at the transmitter side in order to lower the signal fluctuation
qamMat is a complex symbol matrix (for instance QPSK) of size length(allocatedSubcarriers) x nbSymb 

The output signal in time domain is of size (nFFT+nCP)xnbSymb.
# --- Syntax 
sigId = scfdmaSigGen(qamMat,scfdma)
# --- Input parameters 
- qamMat  : Complex T/F matrix to map [Array{Float64},length(allocatedSubcarriers),nbSymb]
- scfdma  : SC-FDMA structure [StrucSCFDMA]
# --- Output parameters 
- sigId	  : Signal in time domain [Array{Complex{Float64}},(nFFT+nCP)xnbSymb]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
# --- MD is waveform structure is given
function scfdmaSigGen(qamMat,scfdma::StrucSCFDMA)
	return scfdmaSigGen(qamMat,scfdma.nFFT,scfdma.nCP,scfdma.allocatedSubcarriers,scfdma.sizeDFT);
end
