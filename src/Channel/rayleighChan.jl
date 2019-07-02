""" 
---  
Generates a correlated random rayleigh sequence of size N, parametrized by the doppler frequency fd and the sampling frqeuency fs.
The method uses IFFT approach described in [1]
# --- Syntax 
α		= rayleighChan(nbPath,N,fs,fd,seed=-1)
# --- Input parameters 
- nbPath : Number of path to be generated (>0) [Int]
- N	: Size of desired output (size of FFT) [Union{Int,Float64}]
- fs	: Sampling frequency [Union{Int,Float64}]
- fd	: doppler frequency [Union{Int,Float64}]
- seed: Seed for random generation ([Int], default =-1)
# --- Output parameters
- α	: Complex rayleigh coefficient [Array{Complex{Float64},N}]
# --- References 
[1] D. J. Young and N. C. Beaulieu, "The generation of correlated Rayleigh random variates by inverse discrete Fourier transform," in IEEE Transactions on Communications

# --- 
"""
function rayleighChan(nbPath,sizeFFT,fs,fd,randSeed=-1)
	# ----------------------------------------------------
	# --- Parameters init
	# ----------------------------------------------------
	# --- FFT size and FFT parameters
	# Core number
	#FFTW.set_num_threads(4);
	# FFT size
	# Forcing power of 2 for FFT computation
	N	= nextpow(2,sizeFFT);
	# --- Compute normalized doppler frequency
	normFd	= fd / fs;
	# --- Deduce grid-based index associated to normFd
	# This is denoted as km in [1]
	km		= Int(floor(normFd * N));
	# Defined as F_M[k] in [1]
	fM	= zeros(Complex{Float64},N);
	if km == 0
		# --- Target doppler frequency cannot be
		# obtained with frequency grid
		# Residual doppler frequency
		# We set the spectrum to have a constant output
		fM[1] = 1;
	else
		# ----------------------------------------------------
		# --- Filter definition
		# ----------------------------------------------------
		# --- Init filter
		# First area
		for k = 0:km-1
			fM[1+km] = sqrt( 1/ (2*sqrt(1 - (k/(normFd*N)^2))));
		end
		# @ km +1
		fM[km+1]       = sqrt( km/2 *( pi/2 - atan( (km-1)/sqrt(2km-1))));
		# Second area
		fM[km+2:N-km]  .= 0;
		# @ N - km +1
		fM[N - km + 1] = fM[km+1];
		# Third area
		fM[ N-km+2:N]  = fM[km:-1:2];
	end
	# ----------------------------------------------------
	# --- Casting output
	# ----------------------------------------------------
	# --- Seed
	if randSeed != -1
		# ---
		Random.seed!(randSeed);
	end
	# --- Creating random sequence
	# We need to have unitary variance @output of system model
	# /2 for I/Q equivalent power
	# Σ | F[k] | ^2 -> Normalisation wrt doppler filter
	variance  = sqrt(1/2/sum(abs2.(fM)));
	seqPhase  = variance* randn(N,nbPath) .* fM;	  # --- In phase
	seqQuad	  = variance* randn(N,nbPath) .* fM;	  # --- Quadrature
	α 		  = sizeFFT * ifft(seqPhase + 1im* seqQuad,1);
	# --- Creating random sequence
	return α
end
