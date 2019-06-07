# ---------------------------------------------------- 
# --- Digital Communication module
# ---------------------------------------------------- 
__precompile__()


module DigitalComm
# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
# Only submodules are considered so don't need to overload packages here 
# --> Go to next section with submodules loading 
using SpecialFunctions

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


# --- AWGN channel model 
include("addNoise.jl");
export addNoise, addNoise!



# ---------------------------------------------------- 
# --- Function definition  
# ---------------------------------------------------- 
# --- Q function definition 
"""
# --- qfunc.jl 
# --- 
# Description 
#    Returns the Q function used for Bit error rate computation in digital system 
# Q(x)= 1/2 erfc(x/sqrt(2))
#erfc is the Complexmplementary error function, i.e. the accurate version of 1-erf(x) for large x
#erfc is inherited from 
# --- 
# Syntax 
#y = qfunc(x)
## --- Input parameters 
#x: Input [Float64]
#y: Q(x)[Float64] 
# --- 
# v 1.0 - Robinn Gerzaguet.
"""
function qFunc(x) 
	return 1/2 * erfc.( x / sqrt(2));
end
export qFunc;


# --- SIR computation 
# --- getSIR
# ---
# Description
#    Returns the Signal to interference ratio expressed in dB (or in linear) between a obersvation signal d(n) and a reference signal u(n)
#	 The ratio is expressed as 10*log10( E[ || d(n) - u(n) || / E[||u(n)||^2]  )
#	  with E the expectation wrt to time
#	  The 2 vectors d and u should have the same length L
# ---
# Syntax
#		  sir	= getSIR( d, u , type="dB")
#				# ---  Input parameters
#					d	: Observation signal [Array{Any}]
#					u	: Reference signal [Array{Any}]
#					type: Output unit [String]: "dB" or "Linear"
#				# --- Output parameters
#					sir	: Signal to interference ratio in unit `type`
# ---
# v 1.0 - Robin Gerzaguet.
""" getSIR 
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
function conv(x::Array{Complex{Float64}},h::Array{Complex{Float64}})
	y = conv(real(x),real(h)).-conv(imag(x),imag(h)).+ 1im*(conv(imag(x),real(h)).+conv(real(x),imag(h)));
end
function conv(x::Array{Complex{Float64}},h)
	y = conv(real(x),h).+1im*conv(imag(x),h);
end
function conv(x,h::Array{Complex{Float64}})
	y = conv(x,real(h)).+1im*conv(x,imag(h));
end



end
