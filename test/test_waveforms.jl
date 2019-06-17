# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
println("Tests for multicarrier waveforms");

# ----------------------------------------------------
# --- Overall parameters
# ----------------------------------------------------
# --- Overall PHY parameters
nbSymb 			= 14;			  # --- Number of symbols (one frame)
nFFT 			= 1024;			  # --- Base FFT size
samplingFreq	= 15.36;		  # --- Frequency value (MHz)
snrVect			= (-10:30);
# --- Frequency allocation
allocatedSubcarriers= getLTEAlloc(nFFT);
nbSubcarriers	      = length(allocatedSubcarriers);

# ----------------------------------------------------
# --- Waveform contender
# ----------------------------------------------------
# --- Init OFDM structure
ofdm  = initOFDM(
				 nFFT,						        # --- nFFT                 : FFT size
				 72,						        # --- nCP                  : CP size
				 allocatedSubcarriers		        # --- allocatedSubcarriers : Subcarrier allocation
				 );
# --- Init SCFDMA structure
scfdma  = initSCFDMA(
					 nFFT,						        # --- nFFT                 : FFT size
					 72,						        # --- nCP                  : CP size
					 allocatedSubcarriers,		        # --- allocatedSubcarriers : Subcarrier allocation
					 12;								# --- sizeDFT			   : DFT preprocessing size
					 );
# --- Init UF-OFDM structure
ufofdm  = initUFOFDM(
					 nFFT,					        # --- nFFT                 : FFT size
					 73,					        # --- L                    : Filter length (same size +1 due to conv)
					 allocatedSubcarriers,	        # --- allocatedSubcarriers : Subcarrier allocation
					 applyPD=1,				        # --- applyPD              : Do predistortion at Tx stage
					 attenuation=40,		        # --- attenuation          : Filter attenuation in dB
					 );
# --- Init BF-OFDM structure
bfofdm	= initBFOFDM(
					 32,				            # --- nFBMC                : PPN size (max number of carriers)
					 64,				            # --- nOFDM                : Precoder size (OFDM sizer)
					 3,				            	# --- K                    : Overlapping factor
					 9,					            # --- GI                   : CP size of precoder
					 0.5,				            # --- δ                    : compression factor
					 allocatedSubcarriers,          # --- allocatedSubcarriers : Subcarrier allocation
					 "gaussian",	            # --- filterName           : Pulse shape name
					 BT=0.36,				        # --- BT                   : Potential BT value for Gaussian
					 filterStopBand = 110,			# --- filterStopBand       : DC stopband value
					 fS=[],				            # --- fS                   : Potential frequency coefficient for FS filter
					 nFFT= 1024,		            # --- nFFT                 : associated FFT value in Rx
					 nCP= 72,			            # --- nCP                  : extended CP size
					 );
# --- Init WOLA-OFDM structure
wola  = initWOLA(
				 nFFT,						        # --- nFFT                 : FFT size
				 72,						        # --- nCP                  : CP size
				 allocatedSubcarriers,		        # --- allocatedSubcarriers : Subcarrier allocation
				 "triangle",						# --- Window type @Tx side
				 20,								# --- Window size @Tx side
				 "triangle",						# --- Window type @Rx side
				 20,								# --- Window size @Rx side
				 );
fbmc  = initFBMC(
				 nFFT,						        # --- nFFT                 : FFT size
				 4,									# --- K					   : Overlapping factor
				 allocatedSubcarriers		        # --- allocatedSubcarriers : Subcarrier allocation
				 );

# ----------------------------------------------------
# --- Merging structures
# ----------------------------------------------------
# Create  a dictionnary to rule them all 
waveforms 	= initWaveforms(ofdm,
							scfdma,
							ufofdm,
							bfofdm,
							wola,
							fbmc,
							);


for  (name,struc) in waveforms
	@testset  "$name" begin 
		for mcs=[4,16,64,256]
			nbBits			      = nbSymb * nbSubcarriers * Int(log2(mcs));
			# --- Binary sequence
			bitSeq	      = genBitSequence(nbBits);
			# Mapping
			qamSeq		  = bitMappingQAM(mcs,bitSeq);
			# --- T/F matrix
			qamMat		  = reshape(qamSeq,nbSubcarriers,nbSymb);
			# --- Signal
			sigId		  = genSig(qamMat,struc);
			# --- Waveform demodulator 
			qamDec	  = decodeSig(sigId,struc);
			# --- Binary demapper
			bitDec	  = bitDemappingQAM(mcs,qamDec[:]);
			# --- BER measure
			nbE	 = sum(xor.(bitDec,bitSeq));
			@test nbE == 0
		end
	end
end

