""" 
---  
Constant Modulus Zero Autocorrelation (CAZAC) sequence generation. Sequence can be used as preamble sequences for OFDM systems. The function generate the ZC sequence in frequency domain.
Generates a Zadoff-Chu  sequence on allocated subcarrier mapped on nFFT vector with generated kernel muPSS and power boost zcBoost
# --- Syntax 
zcSeq	  = genZCSequence(nFFT,allocatedSubcarrier,muPSS=0,zcBoost=0)
# --- Input parameters 
- nFFT				  : FFT size for output [Int]
- allocatedSubcarrier : vector of allocated subcarrier [Array{Int,L}]. 
- muPSS				  : Kernel for ZC sequence [Int] -- default 0
- zcBoost			  : Power boost (in dB) applied to sequence [Float32] -- default 0.
# --- Output parameters
- zcSeq				  : ZC sequence [Array{Complex{Float64}},nFFT]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genZCSequence(nFFT::Int,allocatedSubcarrier::Array{Int},muPss=0::Int,zcBoost=0)
	# --- Getting ZC length
	nbSubcarriers	= length(allocatedSubcarrier);
	# --- Create a subcarrier index vector
	nI				= (1:nbSubcarriers);
	# --- Create frequency container 
	seqZC	= zeros(Complex{Float64},nFFT);
	# --- Fill frequency location with ZC sequence 
	seqZC[allocatedSubcarrier] = (10).^(zcBoost./20) .* exp.(im*pi*muPss*nI.*(nI.+mod(nbSubcarriers,2))/nbSubcarriers); 
	# --- Return generated sequence 
	return seqZC;
end
