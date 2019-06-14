""" ofdmSigDecode
---  
Decodes a time domain OFDM signal and returns a T/F matrix with decoded QAM symbols
# --- Syntax 
qamRx	= ofdmSigDecode(sigId,nFFT,nCP,allocSubcarriers);
# ---  Input parameters
- sigId	  : OFDM signal in time domain [Array{Complex{Float64},nbEch}]  ,(nbEch	: Number of samples: nbSymb*(nFFT+nCp))
- nFFT	  : FFT size [Int]
- nCp	  : Cyclic prefix size (in samples) [Int]
allocSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
# ---  Output parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigDecode(sigChan,nFFT,nCP,allocSubcarriers)
	# --- Getting number of symbol to decode
	nbSymbDecode   = Int(floor(length(sigChan)/(nFFT+nCP)));
	# --- Allocation 
	qamRx	= zeros(Complex{Float64},length(allocSubcarriers),nbSymbDecode);
	# --- OFDM decoder
	return ofdmSigDecode!(qamRx,sigChan,nFFT,nCP,allocSubcarriers); 
end

""" ofdmSigDecode! 
---  
Decodes a time domain OFDM signal and populate the  T/F matrix with decoded QAM symbols
# --- Syntax 
ofdmSigDecode(qamRx,sigId,nFFT,nCP,allocSubcarriers);
# ---  Input parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
- sigId	  : OFDM signal in time domain [Array{Complex{Float64},nbEch}]  ,(nbEch	: Number of samples: nbSymb*(nFFT+nCp))
- nFFT	  : FFT size [Int]
- nCp	  : Cyclic prefix size (in samples) [Int]
allocSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
# ---  Output parameters
- []
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigDecode!(qamRx,sigChan,nFFT,nCP,allocSubcarriers)
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	currSymb  = zeros(Complex{Float64},nFFT);
	# --- OFDM demodulator
	for n = 1  : 1 : size(qamRx,2)
		# --- Getting current symbol without CP
		currSymb  .= sigChan[(nFFT+nCP)*(n-1) .+ nCP .+ (1:nFFT)];
		# --- TF
		fft!(currSymb,1);
		# --- Mapping in output matrix
		qamRx[:,n]		= currSymb[allocSubcarriers];
	end
	return qamRx;
end


""" ofdmSigDecode
---  
Decodes a time domain OFDM signal and returns a T/F matrix with decoded QAM symbols
# --- Syntax 
qamRx	= ofdmSigDecode(sigId,ofdm);
# ---  Input parameters
- sigId	  : OFDM signal in time domain [Array{Complex{Float64},nbEch}]  ,(nbEch	: Number of samples: nbSymb*(nFFT+nCp))
- ofdm	  : OFDM structure [StrucOFDM]
# ---  Output parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigDecode(sigChan,ofdm::StrucOFDM)
	return ofdmSigDecode(sigChan,ofdm.nFFT,ofdm.nCP,ofdm.allocatedSubcarriers);
end


""" ofdmSigDecode!
---  
Decodes a time domain OFDM signal and returns a T/F matrix with decoded QAM symbols
# --- Syntax 
qamRx	= ofdmSigDecode!(qamRx,sigId,ofdm);
# ---  Input parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
- sigId	  : OFDM signal in time domain [Array{Complex{Float64},nbEch}]  ,(nbEch	: Number of samples: nbSymb*(nFFT+nCp))
- ofdm	  : OFDM structure [StrucOFDM]
# ---  Output parameters
- []
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function ofdmSigDecode!(qamRx,sigChan,ofdm::StrucOFDM)
	return ofdmSigDecode!(qamRx,sigChan,ofdm.nFFT,ofdm.nCP,ofdm.allocatedSubcarriers);
end
