""" 
---  
WOLA waveform structure 
# --- Syntax 
- nFFT	  : FFT size [Int] 
- nCP	  : CP size [Int] 
- allocatedSubcarriers	: Vector of allocated subcarriers [Array{Int}] 
- windowTx	: Window used @ transmitter side 
- windowRx	: Window used @ recevier side 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct StrucWOLA<: Waveform
	nFFT::Int;
	nCP::Int;
	allocatedSubcarriers::Array{Int};
	windowTx::Window;
	windowRx::Window;
end


""" 
---  
Create and initiate a WOLA waveform structure 
# --- Syntax 
 wola = initWOLA(nFFT,nCP,allocatedSubcarriers,winFuncTx,winLengthTx,winFuncRx,winLengthRx;windowTx=[],windowRx=[])
# --- Input parameters 
- nFFT	  : FFT size [Int] 
- nCP	  : CP size [Int] 
- allocatedSubcarriers	: Vector of allocated subcarriers [Array{Int}] 
- winFuncTx	  : Name of window used @Tx 
- winLengthTx : Size of window @Tx 
- winFuncRx	  : Name of window used @Rx 
- winLengthRx : Size of window @Rx 
- windowTx	  : Coefficient of Tx window. By default it is empty. To force a given window, populate this vector. If let empty, the window will be created based on winLengthTx and winFuncTx
- windowRx	  : Coefficient of Rx window. By default it is empty. To force a given window, populate this vector. If let empty, the window will be created based on winLengthRx and winFuncRx
# --- Output parameters 
- wola		  : Waveform structure
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function initWOLA(nFFT,nCP,allocatedSubcarriers,winFuncTx,winLengthTx,winFuncRx,winLengthRx;windowTx=[],windowRx=[])
	# ----------------------------------------------------
	# --- Tx window
	# ----------------------------------------------------
	# --- Creating Wola Window
	if isempty(windowTx)
		windowTx	= getWolaWindow(winFuncTx,nFFT,nCP,winLengthTx);
	else
		if !isa(windowTx,Window)
			# --- Proposed window is not Window
			# Convert taps to a Window structure
			windowTx = Window("custom",length(windowTx),windowTx);
		end
	end
	# ----------------------------------------------------
	# --- Rx window
	# ----------------------------------------------------
	# --- Creating Wola Window
	# For receiver, no CP is considered
	if isempty(windowRx)
		windowRx	= getWolaWindow(winFuncRx,nFFT,0,winLengthRx);
	else
		if !isa(windowRx,Window)
			# --- Proposed window is not Window
			# Convert taps to a Window structure
			windowRx = Window("custom",length(windowRx),windowRx);
		end
	end
	# ----------------------------------------------------
	# --- Create objects
	# ----------------------------------------------------
	# --- Create structure
	return StrucWOLA(nFFT,nCP,allocatedSubcarriers,windowTx,windowRx);
end

""" 
---  
Apply Weighted Overlap and Add Orthogonal Frequency Division Multiplexing (WOLA-OFDM) to the time frequency matrix qamMat and returns a time domain OFDM signal
OFDM is parametrized by its FFT size, its cyclic prefix length (in samples) and the allocSubcarriers vector.
The WOLA part is parametrized by the window applied at each beginning and ending of symbols. The window size can be (and is likely to be) higher than the length of the CP and pure OFDM compatibility is ensured by overlapping the symbols. The interested reader can refer to [1] [2] and [3] for WOLA principle description.
# ---
# Syntax
sigId	= wolaSigGen(qamMat,nFFT,nCP,allocSubcarriers,winFunc,winLength=0)
# ---  Input parameters
- qamMat  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
- nbSymb			: Number of OFDM symbol tro be transmitted
- nbSubcarriers	: Number of allocated subcarriers (shall be < nFFT)
- nFFT	  : FFT size [Int]
- nCP	  : Cyclic prefix size (in samples) [Int]
- allocSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
- winFunc : Type of window. Can be a string of supported window or directly the window taps.
- winLength : Length of window. This parameter is not used if winFunc is an array of the window tap.
# ---  Output parameters
- sigId	  : WOLA-OFDM signal in time domain [Array{Complex{Float64},nbEch}]. nbEch	: Number of samples: nbSymb*(nFFT+nCP)
# ---
# Supported window
-	"Triangle"		: Triangle window
-	"srrc"			: Square Root Raised Cosine
-	"Meyer"			: Meyer window (See [1])
# ---
# References
- [1] R. Zayani, Y. Medjahdi, H. Shaiek and D. Roviras, "WOLA-OFDM: A Potential Candidate for Asynchronous 5G," 2016.
- [2] Y. Medjahdi and al, "On the road to 5G: Comparative study of Physical layer in MTC context", 2017.
- [3] R. Gerzaguet and al, "Comparison of Promising Candidate Waveforms for 5G: WOLA-OFDM Versus BF-OFDM", 2017.
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function wolaSigGen(qamMat,nFFT::Union{Int,Float64},nCP::Union{Int,Float64},allocatedSubcarriers,winFunc,winLength=0)
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------
	# --- Getting parameters
	sizeSymb		= nFFT + nCP;
	# --- Extended size
	sizeExt			= sizeSymb +  winLength;
	# --- WOLA signal parameters
	nbSymb			= size(qamMat,2);
	sigWOLA			= zeros(Complex{Float64},nbSymb*sizeSymb+winLength);
	# ----------------------------------------------------
	# --- Getting window
	# ----------------------------------------------------
	if isa(winFunc,String)
		# --- We create a window with a dedicated function
		# No control on window and parameters, it will be done in the function
		window = getWolaWindow(winFunc,nFFT,nCP,winLength);
	else
		# --- We have a array with the filter taps, no processing here
		window = winFunc;
	end
	# ----------------------------------------------------
	# --- Modulator stage
	# ----------------------------------------------------
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# --- Getting parameters
	nbSymb		= size(qamMat,2); 		# --- Applied on x symbols
	sizeSymb	= nFFT + nCP; 			# --- Symbol size
	nSize 		= sizeSymb*nbSymb; 		# --- Total burst length
	# --- Mapping to nFFT elements
	qamCurr		= zeros(Complex{Float64},nFFT,nbSymb);
	qamCurr[allocatedSubcarriers,:] = qamMat;
	# --- Switch to time domain
	sTmp 	= ifft(qamCurr,1);
	# --- Inserting Cyclic prefix and cyclic suffix
	sigCP	  = [sTmp[end+1-nCP-winLength÷2:end,:]; sTmp; sTmp[1:winLength÷2,:] ];
	# --- Apodisation
	sigApod	  = mapslices(x -> x .* window, sigCP,dims=1);
	# --- Overlap and add
	y = zeros(Complex{Float64},sizeSymb*nbSymb + winLength)
	for iB = 1 : 1 : nbSymb
		y[ 1 + (iB-1)*sizeSymb : iB*sizeSymb+ winLength ] = y[ 1 + (iB-1)*sizeSymb : iB*sizeSymb+ winLength ] + sigApod[:,iB];
	end
	return y;
end


# MD Tx
function wolaSigGen(qamMat::Array{Complex{Float64}},wola::StrucWOLA)
	return wolaSigGen(qamMat,wola.nFFT,wola.nCP,wola.allocatedSubcarriers,wola.windowTx.window,wola.windowTx.winLength);
end
