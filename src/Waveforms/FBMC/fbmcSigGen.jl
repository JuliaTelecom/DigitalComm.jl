""" 
---  
FBMC waveform structure
# --- Syntax 
- nFFT	  : FFT size 
- K		  : Overlapping factor 
- allocatedSubcarriers	: Vector of allocated subcarriers [Array{Int}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct StrucFBMC<: Waveform
	nFFT::Int;
	K::Int;
	allocatedSubcarriers::Array{Int};
end

""" 
---  
Initiate FBMC structure 
# --- Syntax 
  fbmc = initFBMC(nFFT,K,allocatedSubcarriers)
# --- Input parameters 
- nFFT	  : FFT size 
- K		  : Overlapping factor 
- allocatedSubcarriers	: Vector of allocated subcarriers [Array{Int}]
# --- Output parameters 
- fbmc	  : FBMC structure [StrucFBMC]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function initFBMC(nFFT,K,allocatedSubcarriers)
	# ---Checking FFT size
	if maximum(allocatedSubcarriers) > nFFT
		error("Subcarrier allocation is impossible");
	end
	# --- Create the FBMC structure
	return StrucFBMC(nFFT,K,allocatedSubcarriers)
end

"""
---  
Apply OQAM pre-processing to incoming matrix qamMat of size nbSubcarriers x nbSymb
# --- Syntax 
oqamMat = oqamConversion(qamMat)
# --- Input parameters
- qamMat : Input complex qam Matrix [Array{Complex{Float64,nbSubcarriers,nbSymb}}]
# --- Output parameters
- oqamMat : OQAM matrix (pure real and pure imag. alterns) [Array{Complex{Float64,nbSubcarriers,nbSymb*2}}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function oqamMapping(qamMat)
	# --- Size of incoming matrix :w
	(nbSubcarriers,nbSymb)	  = size(qamMat);
	# Real and imag split:
	oqamMat	  = zeros(Complex{Float64},nbSubcarriers,2*nbSymb);
	# --- Init matrix
	for iS = 1 : 1 : nbSubcarriers
		# Mask to keep real or imag part 
		pattern		= (-1im)^( mod( mod(iS,2) + 1,2));
		pattern2	= (-1im)^( mod( mod(iS,2) + 2,2));
		# 1 // 1i shift 
		scaleReIm	= (1im)^( mod( mod(iS,2) + 1,2));
		scaleReIm2	= (1im)^( mod( mod(iS,2) + 2,2));
		for iN = 1 : 1 : nbSymb
			# OQAM processing 
			oqamMat[iS,2*(iN-1)+1]	= scaleReIm*real(pattern * qamMat[iS,iN]);
			oqamMat[iS,2*(iN-1)+2]	= scaleReIm2*real(pattern2* qamMat[iS,iN]);
		end
	end
	return oqamMat;
end

"""
---  
Returns the PHYDIAS time domain impulse response of desired FBMC filter parametrized by its overlapping factor and FFT size
# --- Syntax 
p = getFBMCFilter(K,nFFT,type);
# --- Input parameters 
- K	  : Overlapping factor
- nFFT: FFT size 
# --- Output parameters 
- p	  : FBMC time impulse response 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getFBMCFilter(K,nFFT)
	if K == 2
		p1			= sqrt(2)/2;
		P			= p1;
	elseif K== 3
		p1			= 0.91143783;
		P			= [p1;sqrt(1-p1^2)];
	elseif K== 4
		p1			= 0.97195983;
		p2			= 1/sqrt(2);
		p3          = sqrt(1-p1^2);
		P			= [p1;p2;p3];
	elseif K== 5
		# Parameters values are arbitrarily chosen
		p0          = 1;
		P1          = 0.998;
		P2          = 0.844;
		P           = [P1;P2;sqrt(1-P2^2);sqrt(1-P1^2)];
	elseif K== 6
		p0          = 1;
		p1          = 0.99722723;
		p2          = 0.94136732;
		P           = [ p1;p2;sqrt(2)/2;sqrt(1-p2^2);sqrt(1-p1^2)];
	elseif K== 8
		p0          = 1;
		p1          = 0.99988389;
		p2          = 0.99315513;
		p3          = 0.92708081;
		P          = [ p1;p2;p3;sqrt(2)/2;sqrt(1-p3^2);sqrt(1-p2^2);sqrt(1-p1^2)];
	else 
		error("Overlapping factor should be 2, 3, 4, 6 or 8");
	end
	mSize = K*nFFT;
	t	  = collect(1:mSize);
	p	  = zeros(Float64,mSize);
	kV	  = 1:K-1;
	for n = t 
		p[n] = 1 + 2 * sum( (-1).^(kV) .* P .* cos.( 2*pi*kV ./ (K*nFFT) * (n-1) )   );
	end 
	return p./K;
end


"""
---  
Generate a FBMC-OQAM signal in time domain, based on input complex matrix (before OQAM processing) and FBMC parameters.

Transmitter is based on Polyphase Network implementation, with PHYDIAS filter (parametrized by overlapping factor).
# --- Syntax 
fbmcSigGen(qamMat,nFFT,K,allocatedSubcarriers)
# --- Input parameters 
- qamMat	: Complex QAM Time-Frequency matrix [Array{Complex{Float64}}]
- nFFT		: FFT size  [Int]
- K			: Overlapping factor [Int] 
- allocatedSubcarriers	: Vector of allocated subcarriers [Array{Int}]
# --- Output parameters 
- sigId		: Complex baseband signal in time domain [Array{Complex{Float64}},nbEch]; nbEch = (2*nbSymb-1)*nFFT*/2+nFFT*K)
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function fbmcSigGen(qamMat,nFFT,K,allocatedSubcarriers)
	# Set 4 core for FFT computation
	#FFTW.set_num_threads(4)
	# --- Getting parameters
	nbSymb		= size(qamMat,2); 		# --- Applied on x symbols
	# --- Mapping to nFFT elements
	oqamMat		= oqamMapping(qamMat);
	oqamCurr	= zeros(Complex{Float64},nFFT,2*nbSymb);
	oqamCurr[allocatedSubcarriers,:] = oqamMat;
	# --- Switch to time domain
	sTmp 	= ifft(oqamCurr,1);
	# --- Init FBMC filter 
	p = getFBMCFilter(K,nFFT);
	δ = 1/2;					# Space between 2 symbols
	y = zeros(Complex{Float64},Int((2*nbSymb-1)*nFFT*δ +nFFT*K));
	for iN = 1 : 1 : 2*nbSymb
		# Repeat IFFT and Phydias apodisation
		sigApod	  = p .* repeat(sTmp[:,iN],K);
		# Overlap and add 
		y[ (iN-1)*Int(nFFT*δ)  .+ collect(1:( nFFT*K)) ] = y[ (iN-1)*Int(nFFT*δ)  .+ collect(1:( nFFT*K)) ] .+ sigApod;
	end
	return y
end

# --- MD is waveform structure is given
function fbmcSigGen(qamMat,fbmc::StrucFBMC)
	return fbmcSigGen(qamMat,fbmc.nFFT,fbmc.K,fbmc.allocatedSubcarriers);
end
