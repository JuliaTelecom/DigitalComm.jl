module example_BER_awgn 

# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
using DigitalComm
using PGFPlotsX


function doPlot(snrVect,ber,qamVect)
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
				   title       = "AWGN BER for QAM",
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
		ebNo 	= snrLin / log2(qamScheme);
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
	nbSymb			= 1000;					# --- Number of symbols per MC run
	nbIt	        = 10; 					# --- MC runs
	snrVect	        = (-10:30);				# --- SNR range
	# --- Init vector
	nbSNR			= length(snrVect);
	ber				= zeros(Float64,length(qamVect),nbSNR);
	qamSeq			= zeros(Complex{Float64},nbSymb);
	qamNoise		= zeros(Complex{Float64},nbSymb);
	# --- MC run
	for iN = 1 : 1 : length(qamVect)
		# ---  Setting MCS
		mcs = qamVect[iN];
		n	= Int(log2(mcs));
		# ----------------------------------------------------
		# --- Generating data
		# ----------------------------------------------------
		# --- Calculate number of bits
		nbBits		= nbSymb *n;
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
				genBitSequence!(bitSeq,nbBits,124+iN);
				# --- QPSK mapping
				bitMappingQAM!(qamSeq,mcs,bitSeq);
				# ---------------------------------------------------- 
				# --- Channel  
				# ---------------------------------------------------- 
				#  --- AWGN
				# Theoretical power is 1 (normalized constellation)
				addNoise!(qamNoise,qamSeq,snrVect[k],1);
				# ----------------------------------------------------
				# --- Rx Stage: SRRC
				# ----------------------------------------------------
				# --- Binary demapper
				bitDemappingQAM!(bitDec,mcs,qamNoise);
				# --- BER measure
				nbE	 += sum(xor.(bitDec,bitSeq));
				nbC	 += length(bitSeq);
			end
			# --- BER measure
			ber[iN,k]		= nbE / nbC;
		end
	end
	# --- Plotting routine
	doPlot(snrVect,ber,qamVect);
end

end

