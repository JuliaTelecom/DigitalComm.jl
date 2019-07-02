""" 
---  
 Add a white additive circular gaussian noise to input signal
Added noise is real if input signal is real and computedlex is input signal is complex
Noise level is controled by the second input parameter which is the signal to noise ratio (SNR)
Signal power powSig is computed based in the input sequence power (average power in time domain)
although a third parameter (theoretical power of input signal is given). In that case, snr is computed
based on the value powSig
output parameters are the signal with noise, and the noise samples
See the bang methods for non-buffer alloc.
# --- Syntax 
[y,n]   = addNoise(x,snr,powSig);
# --- Input parameters 
- x= Input signal [Array{Real{Float64}}, Array{Complex{Float64}}] of size N
- snr= Desired signal to noise ratio [Float64]
- powSig= Power of input signal. If not given, power is evaluated based on input signal x
# --- Output parameters
- y= Signal with noise [Array{Real{Float64}}, Array{Complex{Float6464}}] of size N
- n= noise samples [Array{Real{Float64}}, Array{Complex{Float64}}] of size N
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function addNoise(x::Array{T},snr::Any,powSig=0) where T
	# --- If powSig is equal to 0, power is computed with x input
	# Expectation is taken with respect to time
	if powSig == 0
		powSig  = sum(abs.(x).^2) / length(x);
	end
	# --- Evaluation of noise power
	powNoise  = sqrt(powSig) *  10^(-snr/20);
	# --- Create a complex random sequence as additive noise
	n		  = randn(T,length(x))* powNoise;
	# --- Create output signal
	y		  = x + n;
	return y,n
end

function addNoise!(y::Array{T},x::Array{T},snr::Any,powSig=0.0) where T
	# --- If powSig is equal to 0, power is computed with x input
	# Expectation is taken with respect to time
	if powSig == 0.0
		powSig  = sum(abs2.(x)) / length(x);
	end
	# --- Evaluation of noise power
	# Additional /sqrt(2) to divide power between I and Q paths
	# No sqrt(2)--> Complex randn is unitary variance in complex plan
	powNoise  = sqrt(powSig) * 10^(-snr/20);
	# --- Create a complex random sequence as additive noise
	y  .= x .+ randn!(y) .* powNoise;
end
