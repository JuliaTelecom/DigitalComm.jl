""" 
---  
Apply WOLA-OFDM demodulator to input signal and returns the T/F QAM matrix
# ---
Syntax
qamRx	= wolaSigDecode(sigId,nFFT,nCP,allocSubcarriers);
# ---  Input parameters
- sigId	  : OFDM signal in time domain [Array{Complex{Float64},nbEch}]
- nbEch	: Number of samples: nbSymb*(nFFT+nCP)
- nFFT	  : FFT size [Int]
- nCP	  : Cyclic prefix size (in samples) [Int]
- allocSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
- winLengthTx	: Window size @Tx side (FFT rotation)
- doTailBiting	: Do tail biting approach (by default 1) [Int]
- winFunc : Type of window. Can be a string of supported window or directly the window taps.
- winLength : Length of window. This parameter is not used if winFunc is an array of the window tap.
# ---  Output parameters
- qamMat  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
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
function wolaSigDecode(sigChan,nFFT,nCP,allocatedSubcarriers,winLengthTx,doTailBiting,winFunc::String,winLength=0,window=[])
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# --- WOLA signal parameters
	sizeSymb		= nFFT + nCP;
	# ----------------------------------------------------
	# --- Getting window
	# ----------------------------------------------------
	if isempty(window)
		# --- We create a window with a dedicated function
		# No control on window and parameters, it will be done in the function
		objWindow = getWolaWindow(winFunc,nFFT,0,winLength);
		window	  = objWindow.window;
	end
	# --- Rotation applied after FFT
	# As we have overlapping between symbol, FFT position is set @ middle (and not at the end
	# as for CP-OFDM)
	# --> We need to calculate this position (function of Tx side)
	# --> We need to apply a rotation factor
	#delta           = (nFFT+nCP)÷2 - (nFFT + winLength)÷2;
	delta            = winLengthTx÷2 + nCP + nFFT -nFFT - winLength;	# End of apod.
	rotFactor       = -winLengthTx÷2 - nCP + delta +winLength÷2 +1;
	# --- Getting number of symbol to decode
	nbSymbDecode   = Int(floor((length(sigChan)-winLength)/sizeSymb));
	# --- Initiate Rx T/F matrix
	qamRx		   = zeros(Complex{Float64},length(allocatedSubcarriers),nbSymbDecode);
	sigW		   = zeros(Complex{Float64},nFFT + winLength);
	# --- OFDM demodulator
	for n = 1  : 1 : nbSymbDecode
		# --- Getting current extnded symbol
		currSymb  = sigChan[ 1+delta + (n-1)*sizeSymb .+ (1:nFFT+winLength)];
		# --- Rx windowing
		currSymb = currSymb .* window;
		# --- Tail biting
		if Bool(doTailBiting)
			# --- Append end to begining
			sigW[ 1:winLength ] = currSymb[ end+1-winLength:end ];
			# --- Append beginning to end
			sigW[ end+1-winLength:end ] = currSymb[ 1:winLength ];
			# --- Add
			sigTB	  = currSymb + sigW;
		else
			sigTB     = currSymb;
		end
		# -- Removing CP
		currSymb  = sigTB[ 1+winLength÷2:end-winLength÷2 ]
		# --- FFT and rotation
		fftEl	  = fft(currSymb,1) .* exp.(-2im*pi*rotFactor/nFFT.*(0:nFFT-1));
		# --- Mapping in output matrix
		qamRx[:,n]		= fftEl[allocatedSubcarriers];
	end
	return qamRx;
end


function wolaSigDecode(sigChan,wola::StrucWOLA,doTailBiting=1);
	return wolaSigDecode(sigChan,wola.nFFT,wola.nCP,wola.allocatedSubcarriers,wola.windowTx.winLength,doTailBiting,wola.windowRx.winFunc,wola.windowRx.winLength,wola.windowRx.window);
end
