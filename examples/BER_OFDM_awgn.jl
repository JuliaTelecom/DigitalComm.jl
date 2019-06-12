module BER_OFDM_awgn 

# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
using DigitalComm
using PGFPlotsX


function doPlot(snrVect,ber,qamVect,ofdm)
	a = 0;
	@pgf a = Axis({
				   ymode	  = "log",
				   height      ="3in",
				   width       ="4in",
				   grid,
				   xlabel      = "SNR [dB]",
				   ylabel      = "Bit Error Rate ",
				   ymax 	   = 1,
				   ymin 	   = 10.0^(-5),
				   title       = "AWGN BER for OFDM",
				   legend_style="{at={(0,0)},anchor=south west,legend cell align=left,align=left,draw=white!15!black}"
				   },
				  Plot({color="red",mark="square*"},Table([snrVect,ber[1,:]])),
				  LegendEntry("QPSK"),
				  Plot({color="green",mark="*"},Table([snrVect,ber[2,:]])),
				  LegendEntry("16-QAM"),

				  Plot({color="purple",mark="triangle*"},Table([snrVect,ber[3,:]])),
				  LegendEntry("64-QAM"),
				  Plot({color="blue",mark="diamond*"},Table([snrVect,ber[4,:]])),
				  LegendEntry("256-QAM"),
				  );
	# ---  Adding theoretical curve
	snrLin  = (10.0).^(snrVect/10)
	for qamScheme = qamVect
		# Taking spectral efficiency into account
		ebNo 	= snrLin / length(ofdm.allocatedSubcarriers) * ofdm.nFFT / log2(qamScheme);
		# This approximation is only valid for high SNR (one symbol error is converted to one bit error with Gray coding).
		berTheo	  = 4 * ( 1 - 1 / sqrt(qamScheme)) / log2(qamScheme) * qFunc.(sqrt.( 2*ebNo * 3 * log2(qamScheme) / (2*(qamScheme-1)  )));
		@pgf push!(a,Plot({color="black"},Table([snrVect,berTheo])));
	end
	display(a);
end 


function main()
	# ----------------------------------------------------
	# --- Overall parameters
	# ----------------------------------------------------

	# --- Defining key parameters
	qamVect	        = [4,16,64,256];		# --- Constellation size
	nbSymb			= 14;					# --- Number of symbols per MC run
	nbIt	        = 10; 					# --- MC runs
	snrVect	        = (-10:30);				# --- SNR range
	# --- OFDM parameters 
	nFFT			= 1024;					# --- LTE-10MHz 
	nCP				= 72; 
	allocatedSubcarriers  = getLTEAlloc(nFFT); 
	nbSubcarriers	= length(allocatedSubcarriers);
	# --- Init vector
	nbSNR			= length(snrVect);
	ber				= zeros(Float64,length(qamVect),nbSNR);
	# --- Create OFDM structure 
	ofdm			= initOFDM(nFFT,nCP,allocatedSubcarriers); 
	# --- Preallocate buffers
	sigId			= zeros(Complex{Float64},(nFFT+nCP)*nbSymb);
	sigNoise		= zeros(Complex{Float64},(nFFT+nCP)*nbSymb);
	qamSeq			= zeros(Complex{Float64},nbSubcarriers*nbSymb);
	qamDec			= zeros(Complex{Float64},nbSubcarriers,nbSymb);
	# --- MC run
	for iN = 1 : 1 : length(qamVect)
		# ---  Setting MCS
		mcs = qamVect[iN];
		n	= Int(log2(mcs));
		# ----------------------------------------------------
		# --- Generating data
		# ----------------------------------------------------
		# --- Calculate number of bits
		nbBits		= nbSymb * nbSubcarriers * n;
		# --- Init MC buffers 
		bitSeq	    = zeros(UInt8,nbBits);
		bitDec	    = zeros(UInt8,nbBits);
		# ----------------------------------------------------
		# --- Iterative BER measure
		# ----------------------------------------------------
		for k = 1 : 1 : nbSNR
			# --- Update counters
			nbC		= 0;
			nbE		= 0;
			for iN = 1 : 1 : nbIt
				# ---------------------------------------------------- 
				# --- Tx stage  
				# ---------------------------------------------------- 
				# --- Create random sequence 
				# Forcing seed
				genBitSequence!(bitSeq,nbBits);
				# --- QPSK mapping
				bitMappingQAM!(qamSeq,mcs,bitSeq);
				# --- Moove to T/F matrix 
				qamSeq	= reshape(qamSeq,nbSubcarriers,nbSymb);
				# --- To OFDM modulator 
				ofdmSigGen!(sigId,qamSeq,nFFT,nCP,allocatedSubcarriers);
				# ---------------------------------------------------- 
				# --- Channel  
				# ---------------------------------------------------- 
				#  --- AWGN
				addNoise!(sigNoise,sigId,snrVect[k]);
				# ----------------------------------------------------
				# --- Rx Stage: SRRC
				# ----------------------------------------------------
				# --- OFDM demodulator 
				ofdmSigDecode!(qamDec,sigNoise,ofdm);
				# --- Binary demapper
				bitDemappingQAM!(bitDec,mcs,@view qamDec[:]);
				# --- BER measure
				nbE	 += sum(xor.(bitDec,bitSeq));
				nbC	 += length(bitSeq);
			end
			# --- BER measure
			ber[iN,k]		= nbE / nbC;
		end
	end
	# --- Plotting routine
	doPlot(snrVect,ber,qamVect,ofdm);
end

end

