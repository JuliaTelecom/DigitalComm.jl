""" 
---  
Demodulate a BF-OFDM signal and returns the obtained T/F matrix composed of QAM symbols (without equalisation).
Based on bfofdmSigGen.jl
# --- Syntax 
qamRx	= bfofdmSigDecode(sigRx,nFBMC,nOFDM,K,GI,δ,allocatedSubcarriers;posWindow="end")
# --- Input parameters 
- sigRx		: Complex baseband signal to decode [Array{Complex{Float64}}]
- nFBMC		: PPN size [Int]
- nOFDM		: FFT precoder size [Int]
- K			: Overlapping factor [Int]
- GI		: CP size of precoder [Int]
- δ			: Compression rate  [Float64]
- allocatedSubcarriers : Vector of allocated subcarriers [Array{Int}]
- posWindow	  : Receiver window position. By default window is at same place as OFDM (i.e drop CP). For BF-OFDM, middle window position can also be considered with small delay spread channel (reduce ISI incuded by PPN). In that case, a phase rotation must be applied (see [1])
# --- Output parameters 
- qamRx		: Decoded constellation [Array{Complex{Float64}}]
# References 
- [1]	: Demmer, D.; Zakaria, R.; Gerzaguet, R.; Doré, J. & Le Ruyet, D. Study of OFDM Precoded Filter-Bank Waveforms, IEEE Transactions on Wireless Communications, 2019.
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function bfofdmSigDecode(sigRx,nFBMC,nOFDM,K,GI,δ,allocatedSubcarriers;posWindow="end")
	# ----------------------------------------------------
	# --- Getting number of decoded symbol and create output matrix
	# ----------------------------------------------------
	# --- PPN output symbols
	nbSymb	  = Int( length(sigRx)/nFBMC/δ +1 -K/δ );
	# --- QAM output symbols
	sizeSymb		= nOFDM + GI;
	nbSymbQamDecode = Int( nbSymb ÷ sizeSymb);
	qamRx			= zeros(Complex{Float64},length(allocatedSubcarriers),nbSymbQamDecode);
	# --- Decoding parameters
	sizeFFT			= Int( nFBMC * nOFDM * δ );
	sizeL			= Int( nFBMC * (nOFDM+GI) * δ );
	# --- Window selector
	# Defining CP as pure CP (and not half CP+CS)
	gil				= GI;				  # --- Left part appended (this is a CP)
	gir				= 0;				  # --- Right part appended (this is a CS -> 0)
	perf_GI			= Int(GI*nFBMC*δ);	  # --- Offset for Rx FFT window
	# Calculating Rx position
	giLeft		    = Int(floor(-(sizeFFT-((nOFDM+gil+gir-1)*nFBMC*δ+K*nFBMC))*δ));

	# ----------------------------------------------------
	# --- Rw Windowing
	# ----------------------------------------------------
	#TODO To be continued
	ww				  = 0;
	windowBFOFDM	  = 0;
	tmpVect			  = zeros(Complex{Float64}, sizeFFT + ww);

	# ----------------------------------------------------
	# --- Decoding
	# ----------------------------------------------------
	if posWindow == "middle"
		# ----------------------------------------------------
		# --- Best position for Rx window
		# ----------------------------------------------------
		# PPN acts as a channel FIR. But, high power taps is not the first h = . : | : .
		# As a consequence, ISI occurs (see [1] for mathematical explanation)
		# Contrary to OFDM where FFT window is get at the end of CP, for us as PPN filter raises and fall off
		# It is better to get the middle position for the Rx FFT.
		for n = 1 : 1 : nbSymbQamDecode
			# --- Getting current signal
			currSig = sigRx[ (n-1)*sizeL + giLeft - ww÷2 .+ collect((1:sizeFFT+ww)) ];
			# --- Windowing
			if windowBFOFDM == 1
				# --- Apply Windowing
				currSig	  = currSig .* filterRx;
				# --- Apply tail biting
				tmpVect[1:winLenght]		  = currSig[ end+1-winLenght:end ];
				tmpVect[ end+1-winLenght:end ]  = currSig[ 1:winLenght ];
				currSig						  = currSig + tmpVect;
				currSig						  = currSig[ 1+winLenght÷2:end-winLenght÷2 ];
			end
			# --- FFT
			# Adding a rotor term as the FFT window is not set at the end of CP but in the middle of the extended symbol
			xFFT	= fft(currSig,1) .* exp.(2im*pi.*collect(0:length(currSig)-1) .* (perf_GI-giLeft ) ./length(currSig));
			# --- OQAM decoding
			qamRx[:,n]	  = xFFT[allocatedSubcarriers];
		end
	else
		# ----------------------------------------------------
		# --- Position at the end of CP
		# ----------------------------------------------------
		# Not a clever idea for interference managment but it works (let's be honest, it is really good except in strong multipath channels).
		# Pure CP-OFDM receiver here.
		for n = 1 : 1 : nbSymbQamDecode
			# --- Getting current signal
			currSig = sigRx[ (n-1)*sizeL+perf_GI.+ collect((1:sizeFFT)) ];
			# --- FFT
			# Note that no need to rotor the signal as window is set on the pure symbol.
			xFFT    = fft(currSig,1);
			# --- Restrict to allocated subcarriers
			qamRx[:,n]	  = xFFT[allocatedSubcarriers];
		end
	end
	# --- Output on demodulated constellation
	return qamRx;
end


# ----------------------------------------------------
# --- Multiple dispatch with structure
# ----------------------------------------------------
function bfofdmSigDecode(sigId::Array{Complex{Float64}},bfofdm::StrucBFOFDM;posWindow="end")
	  return bfofdmSigDecode(sigId,bfofdm.nFBMC,bfofdm.nOFDM,bfofdm.K,bfofdm.GI,bfofdm.δ,bfofdm.allocatedSubcarriers,posWindow=posWindow);
end
