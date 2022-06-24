""" 
Returns the vector of allocated subcarriers associated to Long Term evolution frequency mapping.
In LTE, depending on FFT size, only few subcarriers are allocated (45-55%). This function takes a FFT size as input 
and returns an array of size nbSubcarriers  \\
allocatedSubcarrier = getLTEAlloc(nFFT) \\
Input parameters 
- nFFT	  : Desired FFT size (128, 256, 512, 1024, 1536, 2048, 4096) [Int]
Output parameters  \\
- allocatedSubcarrier : Vector of subcarriers index [Array{Int}]
"""
function getLTEAlloc(nFFT)
	# --- Define RB size 
	# Subcarrier alloc is by multiple of 12 (rbSize)
	rbSize	= 12; 
	# --- RB allocation 
	# Each alloc. has a specific RB alloc.
	rbAlloc	  = Dict( 128	=> 6,256	=> 15, 512	=> 25, 1024	=> 50, 1536 	=> 75,2048 	=> 100, 4096 	=> 200 );
	# --- Getting current allocation 
	currRB		= rbAlloc[nFFT];
	# Handle left/right alloc if odd allocation
	# Negative frequencies have lower subcarrier 
	currRBL		= currRB ÷ 2;
	currRBR		= currRB - currRBL; 
	# Deduce allocation 
	# First half is positive, start at 2 (DC offset removal)
	allocatedSubcarrier = [(1 .+ (1:rbSize*currRBR));(nFFT - currRBL*rbSize + 1:nFFT)];
    return allocatedSubcarrier
end




""" 
Returns the vector of allocated subcarriers associated to 5G New Radio frequency mapping.
In 5G-NR, depending on FFT size, only few subcarriers are allocated (55-65%). This function takes a FFT size as input 
and returns an array of size nbSubcarriers. The output has more subcarriers thanj getLTEAlloc (for 4G-LTE)  \\
allocatedSubcarrier = getLTEAlloc(nFFT) \\
Input parameters 
- nFFT	  : Desired FFT size (128, 256, 512, 1024, 1536, 2048, 4096) [Int]
Output parameters  \\
- allocatedSubcarrier : Vector of subcarriers index [Array{Int}]
"""
function get5GNRAlloc(nFFT)
	# --- Define RB size
	# Subcarrier alloc is by multiple of 12 (rbSize)
	rbSize	= 12;
	# --- RB allocation
	# Each alloc. has a specific RB alloc.
	rbAlloc	  = Dict( 128	=> 6,256	=> 15, 512	=> 25, 1024	=> 52, 1536 	=> 79,2048 	=> 106, 4096 	=> 216 );
	# --- Getting current allocation
	currRB		= rbAlloc[nFFT];
	# Handle left/right alloc if odd allocation
	# Negative frequencies have lower subcarrier
	currRBL		= currRB ÷ 2;
	currRBR		= currRB - currRBL;
	# Deduce allocation
	# First half is positive, start at 2 (DC offset removal)
	allocatedSubcarrier = [(1 .+ (1:rbSize*currRBR));(nFFT - currRBL*rbSize + 1:nFFT)];
    return allocatedSubcarrier
end
