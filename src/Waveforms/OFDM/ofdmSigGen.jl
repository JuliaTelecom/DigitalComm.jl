"""
---  
Structure for OFDM
# --- Syntax 
- nFFT		: FFT size [Int] 
- nCP		: Cyclic prefix size [Int] 
- allocatedSubcarriers	: Vector of allocated subbcarriers [Array{Int}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct StrucOFDM <: Waveform
	nFFT::Int;
	nCP::Int;
	allocatedSubcarriers::Array{Int};
end

""" 
---  
Create OFDM structure 
# --- Syntax 
ofdm	= initOFDM(nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- nFFT	  : FFT size [Int]
- nCP	  : Cyclic prefix size [Int]
- allocatedSubcarrier	: Vector of allocated subbcarriers [Array{Int}]
# --- Output parameters 
- ofdm	  : OFDM structure [StrucOFDM] 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function initOFDM(nFFT,nCP,allocatedSubcarriers)
	# ---Checking FFT size
	if maximum(allocatedSubcarriers) > nFFT || length(allocatedSubcarriers) > nFFT 
		error("Subcarrier allocation is impossible");
	end
	# --- Create the OFDM structure
	return StrucOFDM(nFFT,nCP,allocatedSubcarriers)
end

""" 
---  
Create OFDM signal in time domain based on input T/F matrix and OFDM parameters

qamMat is a complex symbol matrix (for instance QPSK) of size length(allocatedSubcarriers) x nbSymb 

The output signal in time domain is of size (nFFT+nCP)xnbSymb.
# --- Syntax 
  sigId = ofdmSigGen(qamMat,nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- qamMat  : Complex T/F matrix to map [Array{Float64},length(allocatedSubcarriers),nbSymb]
- nFFT	  : FFT size [Int] 
- nCP	  : Cylic prefix size [Int] 
- allocatedSubcarrier	: Vector of allocated subbcarriers [Array{Int}]
# --- Output parameters 
- sigId	  : Signal in time domain [Array{Complex{Float64}},(nFFT+nCP)xnbSymb]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigGen(qamMat,nFFT,nCP,allocatedSubcarriers)
	nbSymb = size(qamMat,2);
	sigId  = zeros(Complex{Float64},(nFFT+nCP)*nbSymb);
	ofdmSigGen!(sigId,qamMat,nFFT,nCP,allocatedSubcarriers);
	return sigId;
end 


""" 
---  
Populate a OFDM signal in time domain based on input T/F matrix and OFDM parameters

qamMat is a complex symbol matrix (for instance QPSK) of size length(allocatedSubcarriers) x nbSymb 

The output signal in time domain is of size (nFFT+nCP)xnbSymb.
# --- Syntax 
  ofdmSigGen!(sigId,qamMat,nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- sigId	  : Signal in time domain [Array{Complex{Float64}},(nFFT+nCP)xnbSymb]
- qamMat  : Complex T/F matrix to map [Array{Float64},length(allocatedSubcarriers),nbSymb]
- nFFT	  : FFT size [Int] 
- nCP	  : Cylic prefix size [Int] 
- allocatedSubcarrier	: Vector of allocated subbcarriers [Array{Int}]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigGen!(sigId,qamMat,nFFT,nCP,allocatedSubcarriers)
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# --- Getting parameters
	nbSymb		= size(qamMat,2); 		# --- Applied on x symbols
	symbSize	= nFFT + nCP; 			# --- Symbol size
	nSize 		= symbSize*nbSymb; 		# --- Total burst length
	# --- Mapping to nFFT elements
	qamCurr		= zeros(Complex{Float64},nFFT,nbSymb);
	qamCurr[allocatedSubcarriers,:] .= qamMat;
	# --- Switch to time domain
	ifft!(qamCurr,1);
	# --- Inserting Cyclic prefix
	#sigId	 .= reshape([qamCurr[end+1-nCP:end,:] ; qamCurr], symbSize*nbSymb,1);
	for iN = 1 : 1 : nbSymb 
		# --- Classic cyclic 
		sigId[ (iN-1)*symbSize +  nCP .+ (1:nFFT)] .= qamCurr[:,iN];
		# --- CP insertion 
		sigId[ (iN-1)*symbSize .+ (1:nCP)] .= qamCurr[end-nCP+1:end,iN];
	end
end

# --- MD is waveform structure is given
"""
---  
Create OFDM signal in time domain based on input T/F matrix and OFDM structure  

qamMat is a complex symbol matrix (for instance QPSK) of size length(allocatedSubcarriers) x nbSymb 

The output signal in time domain is of size (nFFT+nCP)xnbSymb.
# --- Syntax 
  sigId = ofdmSigGen(qamMat,nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- qamMat  : Complex T/F matrix to map [Array{Float64},length(allocatedSubcarriers),nbSymb]
- ofdm	  : OFDM structure [StrucOFDM]
# --- Output parameters 
- sigId	  : Signal in time domain [Array{Complex{Float64}},(nFFT+nCP)xnbSymb]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigGen(qamMat,ofdm::StrucOFDM)
	return ofdmSigGen(qamMat,ofdm.nFFT,ofdm.nCP,ofdm.allocatedSubcarriers);
end


# --- MD is waveform structure is given
""" 
---  
Populate a  OFDM signal in time domain based on input T/F matrix and OFDM structure  

qamMat is a complex symbol matrix (for instance QPSK) of size length(allocatedSubcarriers) x nbSymb 

The output signal in time domain is of size (nFFT+nCP)xnbSymb.
# --- Syntax 
  ofdmSigGen!(sigId,qamMat,nFFT,nCP,allocatedSubcarriers)
# --- Input parameters 
- sigId	  : Signal in time domain [Array{Complex{Float64}},(nFFT+nCP)xnbSymb]
- qamMat  : Complex T/F matrix to map [Array{Float64},length(allocatedSubcarriers),nbSymb]
- ofdm	  : OFDM structure [StrucOFDM]
# --- Output parameters 
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigGen!(sigId,qamMat,ofdm::StrucOFDM)
	return ofdmSigGen!(sigId,qamMat,ofdm.nFFT,ofdm.nCP,ofdm.allocatedSubcarriers);
end
