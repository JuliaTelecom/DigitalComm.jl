# --- testingChannelApplication.jl
# ---
# Description
#    Testing how wo apply time varying channel model
# ---
# Syntax
#    include("testingChannelApplication.jl");
# ---
# v 1.0 - Robin Gerzaguet.


module example_multipathChannel 

# ----------------------------------------------------
# --- Modules
# ----------------------------------------------------
# --- External
using Plots
using Printf
using DSP
using FFTW
# --- Custom
using DigitalComm

# ----------------------------------------------------
## --- Functions
# ----------------------------------------------------
# --- Tx scheme
function minimalPhy(nbSymb,mcs,ofdm)
	# Minimal transmitter scheme
	nbSubcarriers = length(ofdm.allocatedSubcarriers);		  # Subcarrier length
	nbBits		  = Int(nbSymb * nbSubcarriers * log2(mcs));		  # nbBits
	bitSeq	      = genBitSequence(nbBits);					  # Binary sequence
	qamSeq		  = bitMappingQAM(mcs,bitSeq);				  # QAM generation
	qamMat		  = reshape(qamSeq,nbSubcarriers,nbSymb);	  # QAM matrix
	sigId		  = genSig(qamMat,ofdm);					  # OFDM signal
	return sigId;
end
function zeroPad(x,n)
	if n > length(x) 
		out = [x;zeros(n-length(x))];
	else 
		out = x;
	end
	return out;
end
# --- Receiver equalizer
function minimalEqualizer(qamRx,ofdm,strucChan)
	# Simple ZF equalisation scheme:
	qamEq	= zeros(Complex{Float64},size(qamRx));
	firAll	= zeros(Complex{Float64},ofdm.nFFT,size(qamRx,2));
	for iN = 1 : 1 : size(qamRx,2)
		# --- Getting all channel samples in current symbol
		# We get a nFFT x nbTaps matrix (all channel taps in current symbol)
		# We will consider that channel is constant, and that channel of interest
		# is at the middle of the OFDM symbol
		if Bool(strucChan.timeVarying)
			currChan  = strucChan.cir[(iN-1)*(ofdm.nFFT+ofdm.nCP) + ofdm.nCP + ofdm.nFFT÷2,:];
		else
			# COnstant FIR
			currChan = strucChan.cir;
		end
		# Switch to frequency domain
		firFreq   = fft(zeroPad(currChan,ofdm.nFFT));
		qamEq[:,iN] = qamRx[:,iN] ./ firFreq[ofdm.allocatedSubcarriers];
		firAll[:,iN] = firFreq;
	end
	return qamEq,firAll;
end

function plotFIR(fir,ofdm,samplingFreq) 
	# --- Calculate freq axis 
	freqAx	= (((0:1:ofdm.nFFT-1) ./ ofdm.nFFT) .- 0.5) .* samplingFreq;
	firC	= circshift(10*log10.(abs2.(fir)), ofdm.nFFT ÷2);
	pObj	= plot(freqAx,firC[:,1],label="First Channel");
	plot!(pObj,freqAx,firC[:,end],label="Last Channel");
	xlabel!("Frequency [Hz]");
	ylabel!("Channel magnitude");
	return pObj;
end


function main()
	# ----------------------------------------------------
	# --- Parameters
	# ----------------------------------------------------
	sizeFFT	  		= 1024;
	carrierFreq		= 5.2e9;
	speed	  		= 1300;
	samplingFreq	= 15.36e6;
	allocatedSubcarriers= getLTEAlloc(sizeFFT);

	# ----------------------------------------------------
	# --- Minimal Phy
	# ----------------------------------------------------
	# --- Init OFDM structure
	ofdm  = initOFDM(
								 sizeFFT,						    # --- nFFT                 : FFT size
								 72,						        # --- nCP                  : CP size
								 allocatedSubcarriers,		        # --- allocatedSubcarriers : Subcarrier allocation
								 );
	# --- Init parameters
	sigId	  = minimalPhy(14,4,ofdm);
	# ----------------------------------------------------
	# --- Multipath Channel
	# ----------------------------------------------------
	# --- profile 
	profile         = "etu";
	# --- Channel object
	channelObj      = initChannel(profile,carrierFreq,samplingFreq,speed);
	# --- Channel instance
	strucChan0       = getChannel(length(sigId),channelObj,-1);
	# ---------------------------------------------------- 
	# --- Channel application  
	# ---------------------------------------------------- 
	sigChan	  = applyChannel(sigId,strucChan0);
	qamRx	  = ofdmSigDecode(sigChan,ofdm);
	scatter(qamRx[:])
	# ----------------------------------------------------
	# --- Setting a minimal equalizer
	# ----------------------------------------------------
	# --- Assuming a perfect CSI
	qamEq,firFreq	= minimalEqualizer(qamRx,ofdm,strucChan0)
	#
	pObj = scatter(qamRx[:],xlims=(-1,1),ylims=(-1,1),label="before eq.");
	scatter!(pObj,qamEq[:],xlims=(-1,1),ylims=(-1,1),label="after eq.");
	display(pObj);
	# 
	pObj2 = plotFIR(firFreq,ofdm,samplingFreq);
	display(pObj2);


end

end
