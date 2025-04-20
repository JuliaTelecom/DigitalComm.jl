module example_BER_cdma



# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
using DigitalComm 
# --- External Modules
using Plots
using Printf
using FFTW
using Statistics



# ---------------------------------------------------- 
# --- Core functions  
# ---------------------------------------------------- 
function  getBER(waveform,mcs,snrVect,nbSymb,nbIt)
    if waveform isa DigitalComm.StrucCDMA 
        # --- Get code allocation 
        nbSubcarriers       = length(waveform.userMask)
    else 
        # --- Getting frequency allocation 
        allocatedSubcarriers  = waveform.allocatedSubcarriers;
        nbSubcarriers	      = length(allocatedSubcarriers);
    end
	# Deduce number of required bits 
	nbBits			      = nbSymb * nbSubcarriers * Int(log2(mcs));
	# --- Init BER vector
	nbSNR			= length(snrVect);
	ber				= zeros(Float64,nbSNR);
	# --- Iterative PSD calculation
	for k = 1 : 1 : nbSNR
		# --- Update counters
		nbC		= 0;
		nbE		= 0;
		powSig	= 0;
		for iN = 1 : 1 : nbIt
			# --- Binary sequence
			bitSeq	      = genBitSequence(nbBits);
			# Mapping
			qamSeq		  = bitMappingQAM(mcs,bitSeq);
			# --- T/F matrix
			qamMat		  = reshape(qamSeq,nbSubcarriers,nbSymb);
			# --- Signal
			sigId		  = genSig(qamMat,waveform);
			# ---------------------------------------------------- 
			# --- Channel  
			# ---------------------------------------------------- 
			#  --- AWGN
			if iN == 1 
				# We compute the power based on generated signal 
				# We troncate the beginning and end of signal to avoid 
				# estimation biais induced by (potential) filter tails
				powSig = mean(abs2.( @views sigId[1+ end÷4 : end - end÷4]));
			end
			sigNoise,  = addNoise(sigId,snrVect[k],powSig);
			# ----------------------------------------------------
			# --- Rx Stage
			# ----------------------------------------------------
			# --- Waveform demodulator 
			qamDec	  = decodeSig(sigNoise,waveform);
			# --- Binary demapper
			bitDec	  = bitDemappingQAM(mcs,qamDec[:]);
			# --- BER measure
			nbE	 += sum(xor.(bitDec,bitSeq));
			nbC	 += length(bitSeq);
		end
		# --- BER measure
		# Adding 1e-10 to avoid log plot of zero errors
		ber[k]		= 1e-10 .+ nbE / nbC;
	end
	return ber
end

# ---------------------------------------------------- 
# --- Main routine  
# ---------------------------------------------------- 
function main()
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------
	# --- Overall PHY parameters
	nbIt			= 20;			  # --- Iteration number
	nbSymb 			= 1400;			  # --- Number of symbols (one frame)
	nFFT 			= 1024;			  # --- Base FFT size
	samplingFreq	= 3.84e6;		  # --- Frequency value (MHz)
	mcs				= 16;			  # --- 16-QAM 
	snrVect			= (-10:30);
	# --- Frequency allocation
	allocatedSubcarriers= getLTEAlloc(nFFT);

	# ----------------------------------------------------
	# --- Waveform contender
	# ----------------------------------------------------
	# --- Init OFDM structure
    cdma4 = initCDMA(
                    4,
                    :ovsf
                   )
    cdma16 = initCDMA(
                    16,
                    :ovsf
                   )
    cdma32 = initCDMA(
                    32,
                    :ovsf
                   )
    cdma64 = initCDMA(
                    64,
                    :ovsf
                   )
	# ----------------------------------------------------
	# --- Merging structures
	# ----------------------------------------------------
	# Create  a dictionnary to rule them all 
	waveforms 	= initWaveforms(
                                cdma4,
                                cdma16,
                                cdma32,
                                cdma64,
								);


	# ---------------------------------------------------- 
	# --- BER main calculation  
	# ---------------------------------------------------- 
	# --- Init plot container 
	plt	  = plot(reuse=false,yscale=:log10,legend=:bottomleft);
	# --- Iterative PSD generation
	for (name,struc) in waveforms 
		# --- Calculate PSD for the configuration 
		ber	= getBER(struc,mcs,snrVect,nbSymb,nbIt);
		# --- Plot stuff
		plot!(plt,snrVect,ber,label=struc.nbUsers);
	end
	# --- Adding metata do plot curve
	xlabel!("SNR [dB]");
	ylabel!("Bit Error Rate");
	ylims!(1e-6,1);
	display(plt);
end



end



