# ---------------------------------------------------- 
# --- Digital Communication module
# ---------------------------------------------------- 
__precompile__()


module DigitalComm
# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
# Only submodules are considered so don't need to overload too much packages here 
# --> Go to next section with submodules loading 
using FFTW 


# ---------------------------------------------------- 
# --- Submodules inclusion  
# ---------------------------------------------------- 
# --- Binary managment 
include("genBitSequence.jl");
# Export 
export genBitSequence!  , genBitSequence;
export genByteSequence! , genByteSequence;

# --- QAM bit mapping 
include("bitMapping.jl");
export bitMappingQAM! , bitMappingQAM;

# --- QAM Hard demapper 
include("bitDeMapping.jl");
export bitDemappingQAM! ,  bitDemappingQAM;
include("hardConstellation.jl");
export hardConstellation! , hardConstellation;

# --- QAM soft demapper 
include("symbolDemapper.jl");
export symbolDemappingQAM! , symbolDemappingQAM; 
export llrToHardBits! , llrToHardBits; 
export llrToUInt! , llrToUInt; 

# --- AWGN channel model 
include("Channel/addNoise.jl");
export addNoise, addNoise!
include("Channel/rayleighChan.jl");
export rayleighChan; 
include("Channel/getChannel.jl");
export initChannel 
export getChannel 
export getFIRResponse
export applyChannel 



# ---------------------------------------------------- 
# --- Function definition  
# ---------------------------------------------------- 
# --- Q function definition 
import SpecialFunctions.erfc
"""
---  
Returns the Q function used for Bit error rate computation in digital system 
 Q(x)= 1/2 erfc(x/sqrt(2))
erfc is the Complexmplementary error function, i.e. the accurate version of 1-erf(x) for large x
erfc is inherited from DSP
# --- Syntax 
y = qfunc(x)
# --- Input parameters 
- x: Input [Float64]
# --- Output parameters 
- y: Q(x)[Float64] 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function qFunc(x) 
	return 1/2 * erfc.( x / sqrt(2));
end
export qFunc;



"""
--- 
Returns the Signal to interference ratio expressed in dB (or in linear) between a obersvation signal d(n) and a reference signal u(n)
The ratio is expressed as 10*log10( E[ || d(n) - u(n) || / E[||u(n)||^2]  )
 with E the expectation wrt to time
 The 2 vectors d and u should have the same length L
# --- Syntax 
  sir = getSIR( d, u , type="dB")
# ---  Input parameter 
- d	: Observation signal [Array{Any}]
- u	: Reference signal [Array{Any}]	
- type: Output unit [String]: "dB" or "Linear" (default, "dB")
# --- Output parameters
- sir	: Signal to interference ratio in unit `type`
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getSIR(d,u,type="dB")
	# --- Setting type
	if lowercase(type) == "linear"
		# --- flag for casting
		castype = 1;
	else
		# --- By default SIR is in dB
		castype = 0;
	end
	# --- Moving vector to column system
	ur	= u[:];
	dr  = d[:];
	# --- Checking length
	if length(ur) != length(dr)
		# --- Raise an error
		error("getSIR.jl: Observation and reference signals should have the same length")
	end
	# --- Calculating SIR
	sirLinear	= mean( (abs2.( ur - dr))) / custMean( (abs2.(dr)));
	# --- Setting output
	if castype == 1
		# --- Return linear SIR
		return sirLinear
	else
		# --- Return SIR in dB scale
		return 10*log10(sirLinear)
	end
end
export getSIR;


# --- Complex Convolution definition 
# In DSP conv is only defined as Re * Re -> Extend to C
import DSP.conv
function conv(x::Array{Complex{T}},h::Array{Complex{T}}) where T
	 y = conv(real(x),real(h)).-conv(imag(x),imag(h)).+ 1im*(conv(imag(x),real(h)).+conv(real(x),imag(h)));
end
function conv(x::Array{Complex{T}},h) where T
	y = conv(real(x),h).+1im*conv(imag(x),h);
end
function conv(x,h::Array{Complex{T}}) where T
	y = conv(x,real(h)).+1im*conv(x,imag(h));
end


# ---------------------------------------------------- 
# --- Waveform type definition  
# ---------------------------------------------------- 
""" 
---  
Abstract type gathering all waveform configuration 
# v 1.0 - Robin Gerzaguet.
"""
abstract type Waveform end
export Waveform 

# --- Waveforms functions
# ZC sequence
include("./Waveforms/genZCSequence.jl");
export genZCSequence;
# LTE Allocation
include("./Waveforms/getLTEAlloc.jl");
export getLTEAlloc;

# OFDM generation
include("./Waveforms/OFDM/ofdmSigGen.jl");
include("./Waveforms//OFDM/ofdmSigDecode.jl");
export ofdmSigGen, initOFDM;
export ofdmSigGen!;
export ofdmSigDecode;
export ofdmSigDecode!;

# SC-FDMA generation
include("./Waveforms//SCFDMA/scfdmaSigGen.jl");
include("./Waveforms//SCFDMA/scfdmaSigDecode.jl");
export scfdmaSigGen, initSCFDMA;
export scfdmaSigDecode, scfdmaPostProcessing;

# UF-OFDM generation
include("./Waveforms//UFOFDM/filterUFOFDM.jl")
include("./Waveforms//UFOFDM/ufofdmSigGen.jl");
include("./Waveforms//UFOFDM/ufofdmSigDecode.jl");
export ufofdmSigGen, initUFOFDM;
export ufofdmSigDecode;

# BF-OFDM generation
include("./Waveforms/BFOFDM/BFOFDM_filter.jl")
include("./Waveforms/BFOFDM/carrierManipulation.jl");
include("./Waveforms/BFOFDM/bfofdmSigGen.jl");
include("./Waveforms/BFOFDM/bfofdmSigDecode.jl");
export getBFOFDMFilter;
export getBFOFDM_subCarrierFromFBMCCarriers, getBFOFDM_carrierFromSubcarriers,getBFOFDM_oversampledGridSubcarriers;
export bfofdmSigGen, initBFOFDM;
export bfofdmSigDecode;

# WOLA generation 
include("./Waveforms/WOLA/getWolaWindow.jl");
export getWolaWindow; 
include("./Waveforms/WOLA/wolaSigGen.jl");
export initWOLA,wolaSigGen
include("./Waveforms/WOLA/wolaSigDecode.jl");
export wolaSigDecode

# FBMC Generation 
include("./Waveforms/FBMC/fbmcSigGen.jl");
export fbmcSigGen 
export getFBMCFilter 
export oqamMapping 
export initFBMC 
include("./Waveforms/FBMC/fbmcSigDecode.jl"); 
export oqamDemapping 
export fbmcSigDecode



# --- Global waveform alias
include("./Waveforms/genSig.jl");
export genSig
export getWaveformName;
export initWaveforms 
export decodeSig

end
