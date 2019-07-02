# ----------------------------------------------------
## --- Channel definition
# ----------------------------------------------------
# Object with all channel attribute (but no realisations)
""" 
---  
Channel model object
# --- Syntax 
- profile		  : Name of profile [String]
- carrierFreq	  : Carrier frequency [Union{Int,Float64}]
- samplingFreq	  : Sampling frequency [Union{Int,Float64}]
- speed			  : Desired speed (km/h) [Union{Int,Float64}]
- powerProfile	  : Distribution of power values in dB [Array{Float64}]
- delayProfile	  : Distribution of delay values in s [Array{Float64}] (same size as powerProfile) 
- dopplerFreq	  : Doppler frequency (inherited from samplingFreq  and speed) 
- delaySpread	  : Max support of CIR in samples [Int]
# --- 
"""
struct ChannelModel
	profile::String;					  # Name of profile
	carrierFreq::Union{Int,Float64};	  # Carrier frequency (Hz)
	samplingFreq::Union{Int,Float64};	  # Sampling frequency (Hz)
	speed::Union{Int,Float64};			  # Velocity (km/h)
	powerProfile::Union{Array{Int},Array{Float64}};		  # Vector of power value (dB)
	delayProfile::Union{Array{Int},Array{Float64}};		  # Vector of delay value (s or samples)
	dopplerFreq::Union{Int,Float64};	  # Doppler frequency in Hz
	#randSeed::Int;					  	  # Seed for random variables %
	delaySpread::Int;					  # Delay spread value (samples)
end


# 
""" 
--- 
Object with channel realisations
# --- Syntax 
- timeVarying	: Flag for constant vs time varying channel [Int] 
- cir			: CIR matrix (nbTap x nbChannel) [Union{Array{Int},Array{Float64},Array{Complex{Float64}}}]
- channelModel  : Model use for generation [DigitalComm.ChannelModel]
- powerLin		: Power distribution with the rayleigh distribution [Array{Complex{Float64}}] 
- randSeed		: Seed use for rayleigh generation 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct ChannelImpl
	timeVarying::Int;					  # Time varing or constant channel
	cir::Union{Array{Int},Array{Float64},Array{Complex{Float64}}}; # CIR matrix (nbTap x nbReal)
	channelModel::ChannelModel;			  # Channel model baseline
	powerLin::Array{Complex{Float64}};	  # Power distribution (with random rayleigh application)
	randSeed::Int;						  # Used seed
end

# ----------------------------------------------------
## ---  Channel initialisation
# ----------------------------------------------------
""" 
---  
 Create a channel object associated to physical parameters and propagation profile.
# --- Syntax 
channObj  = initChannel(profile,carrierFreq,samplingFreq,dopplerFreq,powerProfile,delayProfile,randSeed=-1)
# --- Input parameters 
# --- Input parameters
- profile : Multipath profile (see below) [String]
- carrierFreq : Carrier frequency in Hz [Union{Int,Float64}]
- samplingFreq: baseband sampling frequency [Union{Int,Float64}]
- speed : Velocity (km/h) [Float64]
- randSeed : Seed if necessary (default: -1) [Float64]
# --- Output parameters
- channel : Julia Channel Object [Channel]
# --- Channel models 
-	ETU	  : Extended Typical Urban
-	EVA	  : Extended Vehicular area
-	TDLC  :
-	Rayleigh : Single tap rayleigh channel model
# --- 
"""
function initChannel(profile::String,carrierFreq::Union{Int,Float64},samplingFreq::Union{Int,Float64},speed::Union{Int,Float64};powerProfile=[],delayProfile=[],tau=1)
	# ----------------------------------------------------
	# --- Init parameters
	# ----------------------------------------------------
	# --- Creating profile without case sensitive
	profile = lowercase(profile);
	# ----------------------------------------------------
	# --- Channel model switch
	# ----------------------------------------------------
	if profile == "multipath"
		# ----------------------------------------------------
		# --- Custom multipath channel profile
		# ----------------------------------------------------
		# In this case, the profile is direclty inherited from powerProfile and delayProfile array
		# --- Calculating doppler frequency
		dopplerFreq	  = speed * 1e3 / 3600 * carrierFreq / 3e8;
		# --- Creating object
		return ChannelModel("multipath",carrierFreq,samplingFreq,speed,powerProfile,delayProfile,dopplerFreq,max(delayProfile));
	elseif profile == "none" || profile == "awgn"
		# ----------------------------------------------------
		# --- AWGN channel model
		# ----------------------------------------------------
		# AWGN channel model, constant with no doppler frequency
		return ChannelModel("none",carrierFreq,samplingFreq,0,[1.0],[1.0],0,1);
	elseif profile == "exponential"
		# ----------------------------------------------------
		# --- Exponential rayleigh profile
		# ----------------------------------------------------
		# --- Setting doppler frequency
		dopplerFreq	  = speed * 1e3 / 3600 * carrierFreq / 3e8;
		# --- Creating random profile
		# We use tau as a specific parameter for exponential decay profile
		# This is a constant time --> Channel spread on 5tau
		# tau can be given in s or in samples. We assume that if tau > 1, it is in samples.
		if tau < 1
			# --- second to samples conversion
			tau	  = tau * samplingFreq;
		end
		# --- Channel delay spread
		delaySpread	  = 5 * tau;
		# --- Enveloppe
		enveloppe	  = 10*log10.(exp.(-(0:1:delaySpread-1) ./ tau));
		# --- Creating object
		return ChannelModel("exponential",carrierFreq,samplingFreq,speed,enveloppe,collect(1:delaySpread),dopplerFreq,delaySpread);
	elseif profile == "flatrayleigh"
		# ----------------------------------------------------
		# --- Flat fading channel
		# ----------------------------------------------------
		# --- Calculating doppler frequency
		dopplerFreq	  = speed * 1e3 / 3600 * carrierFreq / 3e8;
		# --- Creating object
		return ChannelModel("flatRayleigh",carrierFreq,samplingFreq,speed,[1],[1],dopplerFreq,1);
	elseif profile == "eva"
		# ----------------------------------------------------
		# ---EVA channel
		# ----------------------------------------------------
		# --- Calculating doppler frequency
		dopplerFreq	    = speed * 1e3 / 3600 * carrierFreq / 3e8;
		# --- EVA power and delay profile
		# second
		powerProfile			= [0.0; -1.5; -1.4; -3.6; -0.6; -9.1; -7.0; -12.0; -16.9];
		tapDelaySecond		= [0; 30; 150; 310; 370; 710; 1090; 1730; 2510] * 1e-9;
		delaySpread		    = ceil(maximum(tapDelaySecond) * samplingFreq);
		return ChannelModel("EVA",carrierFreq,samplingFreq,speed,powerProfile,tapDelaySecond,dopplerFreq,delaySpread);
	elseif profile == "consteva"
		# ----------------------------------------------------
		# ---EVA channel
		# ----------------------------------------------------
		# --- EVA power and delay profile
		# second
		powerProfile			= [0.0; -1.5; -1.4; -3.6; -0.6; -9.1; -7.0; -12.0; -16.9];
		tapDelaySecond		= [0; 30; 150; 310; 370; 710; 1090; 1730; 2510] * 1e-9;
		delaySpread		    = ceil(maximum(tapDelaySecond) * samplingFreq);
		return ChannelModel("EVA",carrierFreq,samplingFreq,0,powerProfile,tapDelaySecond,0,delaySpread);

	elseif profile == "etu"
		# ----------------------------------------------------
		# --- ETU fading channel
		# ----------------------------------------------------
		# --- Calculating doppler frequency
		dopplerFreq	    = speed * 1e3 / 3600 * carrierFreq / 3e8;
		# --- EVA power and delay profile
		# second
		powerProfile		= [-1.0; -1.0; -1.0; +0.0; +0.0; +0.0; -3.0; -5.0; -7.0];
		tapDelaySecond		= [0;   50;  120;  200;  230;  500; 1600;  2300; 5000] * 1e-9;
		delaySpread		    = ceil(maximum(tapDelaySecond) * samplingFreq);
		return ChannelModel("ETU",carrierFreq,samplingFreq,speed,powerProfile,tapDelaySecond,dopplerFreq,delaySpread);
	elseif profile == "constetu"
		# ----------------------------------------------------
		# --- ETU fading channel: w/o modification even with speed
		# ----------------------------------------------------
		# --- EVA power and delay profile
		# second
		powerProfile			= [-1.0;-1.0; -1.0; +0.0; +0.0; +0.0; -3.0; -5.0; -7.0];
		tapDelaySecond		= [0;   50;  120;  200 ; 230;  500; 1600;  2300; 5000] * 1e-9;
		delaySpread		    = ceil(maximum(tapDelaySecond) * samplingFreq);
		return ChannelModel("ETU",carrierFreq,samplingFreq,0,powerProfile,tapDelaySecond,0,delaySpread);

	end
end

""" 
---  
Returns a complete FIR response based on the power profile (in linear scale), the delay profile (in s) and the sampling frequency.
# --- Syntax 
cir = getFIRResponse(delayProfile:,powerProfile,freq,sincSupport=5,interpSystem=0)
# --- Input parameters 
- delayProfile	: Vector of delay in second [Union{Array{Int},Array{Float64}}]
- powerProfile	: Vector of attenuation with same size of delayProfile in linear scale with rayleigh distribution
- freq			: Sampling frequency
- sincSupport	: Size of interpolator (default 5)
- interpSystem	: Forcing extra delay for all samples to have proper interpolation of the FIR beginning (default 0)
# --- Ouput parameters
- cir			: Finite impulse response [Array{Complex{Float64}},max(delayProfile)+sincSupport]
# --- 
"""
function getFIRResponse(delayProfile::Union{Array{Int},Array{Float64}},powerProfile,freq,sincSupport=5,interpSystem=0)
	# --- Instantiate grid
	# FIXME: A signal processing question here.
	# % We should do fractional interpolation as the delay provided in s
	# % will not match the time grid (and depends on the sampling frequency).
	# % I propose Shannon interpolation (and not Lagragian here) in order to have something
	# % slightly different with the farrow interpolation that can be done in practical Rx
	# % However, it leads to the definition of a sinc support as the spread of the interpolator is infinite.
	# % 5 additional taps for sinc should be sufficient
	# % Two ways to to this:
	# %	  (0) - Forcing 5 taps before and after: proper interpolator but all delays will be shifted by 5Te
	# %	  (1) - Keep all the same delays values: it means that the interpolation will not be perfect @ the first samples
	# % interpSystem = 0 or 1
	if Bool(interpSystem)
		# No additional delay
		sS		= sincSupport;			# Half support (only append @ end)
		fdGap	= 0;					# We keep the delay grid
	else
		# Grid is shifted to enable FD (fractional delay) interpolation
		sS		= 2*sincSupport;		# Support is doubled (before and after delay grid)
		fdGap	= sincSupport / freq;	# We shift all the delay value to be sure first samples are FD.
	end
	# Final CIR with size of nbTaps + overhead due to sinc. interpolation
	nbTaps	  = length(delayProfile);
	delaySpread = Int(ceil(delayProfile[end]*freq));
	cir		  = zeros(Complex{Float64}, delaySpread + sS);
	# --- Ensuring additional delay for proper interpolation
	delay	  = delayProfile .+ fdGap;
	# sinc. interpolation due to delay in s (not in grid)
	for k = 1 : 1 : delaySpread + sS
		# --- Each index is superimposition of all timed compondent
		for ∂ = 1 : 1 : nbTaps
			# --- Shannon interpolation with sinc. on sincSupport
			cir[k] += powerProfile[∂] * sinc.(pi*(k/freq - delay[∂]) * freq);
		end
	end
	return cir;
end






# ----------------------------------------------------
# --- Create and instanciate channel model
# ----------------------------------------------------
""" 
---  
Create a channelImpl based on desired channel model and number of realisations 
# --- Syntax 
channelImpl	  = getChannel(nbSamples,channelModel,randSeed=-1)
# --- Input parameters 
- nbSamples		: Number of samples on which channel will be applied [Int]
- channelModel	: Channel object [ChannelModel]
- randSeed		: Desired seed (default -1)
# --- Output parameters 
- channelImpl	: Channel implementation 
# --- 
"""
function getChannel(nbSamples::Int,channelModel::ChannelModel,randSeed=-1)
	# ----------------------------------------------------
	# --- Constant or varying channel
	# ----------------------------------------------------
	if channelModel.speed == 0 || channelModel.dopplerFreq == 0
		timeVarying = 0;
	else timeVarying = 1;
	end
	# ----------------------------------------------------
	# --- Other parameters
	# ----------------------------------------------------
	sincSupport	  = 5;	  # Extra support for sinc pulse shape (in samples)
	# --- Setting the random seed if given
	if randSeed != -1
		# --- Seed as the second parameter
		Random.seed!(randSeed);
	end
	if channelModel.profile == "none"
		# ----------------------------------------------------
		# --- Special case of none channel
		# ----------------------------------------------------
		# Channel is unitary, no power loss (no rayleigh distrib)
		cir = [1.0];
		return ChannelImpl(0,cir,channelModel,[1.0],randSeed);
	else
		if Bool(timeVarying)
			# ----------------------------------------------------
			# --- Time varying channel
			# ----------------------------------------------------
			# FIXME: A signal processing question here.
			# % We should do fractional interpolation as the delay provided in s
			# % will not match the time grid (and depends on the sampling frequency).
			# % I propose Shannon interpolation (and not Lagragian here) in order to have something
			# % slightly different with the farrow interpolation that can be done in practical Rx
			# % However, it leads to the definition of a sinc support as the spread of the interpolator is infinite.
			# % 5 additional taps for sinc should be sufficient
			# % Two ways to to this:
			# %	  (0) - Forcing 5 taps before and after: proper interpolator but all delays will be shifted by 5Te
			# %	  (1) - Keep all the same delays values: it means that the interpolation will not be perfect in the first samples
			# % interpSystem = 0 or 1
			interpSystem  = 1;
			if Bool(interpSystem)
				# No additional delay
				sS		= sincSupport;			# Half support (only append @ end)
				fdGap	= 0;					# We keep the delay grid
			else
				# Grid is shifted to enable FD (fractional delay) interpolation
				sS		= 2*sincSupport;		# Support is doubled (before and after delay grid)
				fdGap	= sincSupport / freq;	# We shift all the delay value to be sure first samples are FD.
			end
			# --- Generate realisation of rayleigh distribution
			nbTaps	  = length(channelModel.powerProfile);
			delaySpread = Int(ceil(channelModel.delayProfile[end]*channelModel.samplingFreq));
			# --- Creating amplitudes
			powerLin  = 10 .^(channelModel.powerProfile / 10);
			# --- Correlated rayleigh amplitude
			raylAmpl  = rayleighChan(nbTaps,nbSamples,channelModel.samplingFreq,channelModel.dopplerFreq,randSeed);
			# --- Getting sampling frequency
			freq	  = channelModel.samplingFreq;
			# --- Ensuring additional delay for proper interpolation
			delay	  = channelModel.delayProfile .+ fdGap;
			cir		  = zeros(Complex{Float64}, nbSamples, delaySpread + sS);
			for indexTime = 1 : 1: nbSamples
				for k = 1 : 1 : delaySpread + sS
					# --- Each index is superimposition of all timed compondent
					for ∂ = 1 : 1 : nbTaps
						# --- Shannon interpolation with sinc. on sincSupport
						cir[indexTime,k] += powerLin[∂] * raylAmpl[indexTime,∂] * sinc.(pi*(k/freq - delay[∂]) * freq);
					end
				end
			end
			return ChannelImpl(1,cir,channelModel,powerLin,randSeed);
		else
			# ----------------------------------------------------
			# --- Constant FIR channel
			# ----------------------------------------------------
			# --- Generate realisation of rayleigh distribution
			nbTaps	  = length(channelModel.powerProfile);
			raylComp  = 1/sqrt(2)*(randn(nbTaps) + 1im .* randn(nbTaps));
			# --- Creating amplitudes
			powerLin  = 10 .^(channelModel.powerProfile / 10) .* raylComp;
			# --- Create FIR
			cir	  = getFIRResponse(channelModel.delayProfile,powerLin,channelModel.samplingFreq);
			# --- Create object
			return ChannelImpl(0,cir,channelModel,powerLin,randSeed);
		end
	end
end



# --- Channel model
""" 
---  
Apply a channel implementation  to an input signal
# --- Syntax 
sigChan	  : applyChannel(sigId,channelImpl)
# --- Input parameters 
- sigId	  : Input signal [Array{Any}]
- channelImpl : Channel implementation [ChannelImpl]
# --- Output parameters 
- sigChan	: Output signal [Array{Any}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function applyChannel(sigId,channelImpl)
	if !Bool(channelImpl.timeVarying)
		# --- Classic convolution
		return sigChan = conv(sigId,channelImpl.cir);
	else
		# --- Time vraying convolution
		# Check size
		if length(sigId) != size(channelImpl.cir,1)
			error("Time varying channel should have same size as incoming signal");
		end
		# Init outut vector: Output is a convolution
		nbTaps	=size(channelImpl.cir,2);
		nbOut	=length(sigId)+nbTaps
		sigChan = zeros(Complex{Float64},nbOut)
		# n is time index of output, δ is lag index
		for n = 1 : 1 : length(sigId)
			# --- At beginning, not possible to have all index
			nbMaxTap  = min(n,nbTaps);
			# --- Classic time domain convolution
			for δ = 1 : 1 : nbMaxTap
				sigChan[n] += channelImpl.cir[n,δ] * sigId[n-δ+1];
			end
		end
		# At the end we have zero
		buffEnd = zeros(Complex{Float64},nbTaps *2);
		buffEnd[1:nbTaps] = sigChan[1:nbTaps];
		for n = 1:nbTaps
			# --- Classic time domain convolution
			for δ = 1 : 1 : nbTaps
				sigChan[length(sigId)+n] += channelImpl.cir[end,δ] * buffEnd[n+nbTaps-δ+1];
			end
		end
		return sigChan;
	end
end


