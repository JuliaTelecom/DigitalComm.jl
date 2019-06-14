""" 
---  
Apply Universal Filtered Orthogonal Frequency Division Multiplexing (UF-OFDM) to the time frequency matrix qamMat and returns a time domain ufofdm signal [1]

ufofdm is parametrized by its FFT size, the filter length (in samples) and the allocatedSubcarriers vector. Optional parameters are carrier size in subcarrier (by default RB size which is 12)
# --- Syntax 
qamRx	= genereSignalufofdm(sigRx,nFFT,nCp,allocatedSubcarriers;sizeRB=12)
# ---  Input parameters
- sigRx	  : Time UF-OFDM recevied signal [Array{Complex{Float64},nbEch}]
- nFFT	  : FFT size [Int]
- L		  : Dolph Chebyshev filter length [Int]
- allocatedSubcarriers : Vector of index of allocated subcarriers [Array{Int,nbSubcarriers}]
- sizeRB  : Carrier size in subcarriers (default : LTE RB size: 12) [Int]
# ---  Output parameters
- qamRx  : Time frequency matrix : [Array{Complex{Float64},nbSubcarriers,nbSymb}]
- nbSymb			: Number of OFDM symbol tro be transmitted
- nbSubcarriers	: Number of allocated subcarriers (shall be < nFFT)
# --- 
# v 1.0 - Robin Gerzaguet.
"""

function ufofdmSigDecode(sigRx,nFFT,L,allocatedSubcarriers;sizeRB=12,window=0)
	# ----------------------------------------------------
	# --- UF-OFDM parameters
	# ----------------------------------------------------
	# --- Symbol size (classic convolution)
	sizeSymb	  = nFFT + L - 1;
	nbSymbDecode		  = Int( floor( length(sigRx)/sizeSymb  ));
	# --- Number of physical RB (i.e carriers)
	nbRB		  = Int(floor(length(allocatedSubcarriers) / sizeRB));
	nbDataSubcarrier = length(allocatedSubcarriers);
	if window == 1
		# TODO Create windowing on UFOFDM Rx
	end
	qamRx	= zeros(Complex{Float64},nbDataSubcarrier,nbSymbDecode);
	# --- Iterative decoding
	for iB = 1 : 1 :nbSymbDecode
		# --- Extract subsignal
		subSig	  = sigRx[ 1+(iB-1)*sizeSymb:iB*sizeSymb ];
		# --- Apodisation
		if window != 0
			subSig	= subSig .* filterTimeDomain;
		end
		# --- Zero Padding to get a 2nFFT
		sigZP	  = [subSig;zeros(Complex{Float64},nFFT-L+1)];
		# --- Apply FFT
		sigFFT	  = fft(sigZP);
		# --- Getting only odd subcarrier
		sigFFT	  = sigFFT[1:2:end];
		# --- Getting received QPSK symbols
		qamRx[ :,iB ]  = sigFFT[allocatedSubcarriers];
	end
	return qamRx;
end
function ufofdmSigDecode(sigChan,ufofdm::StrucUFOFDM;window=0)
	return ufofdmSigDecode(sigChan,ufofdm.nFFT,ufofdm.L,ufofdm.allocatedSubcarriers,sizeRB=ufofdm.sizeRB,window=window);
end
