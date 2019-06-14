""" 
---  
Structure describing window parameters 
# --- Syntax 
- winFunc	: Name of window  
- winLength : Size of window 
- window	: Window coefficients 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct Window
	winFunc::String;
	winLength::Int;
	window::Array{Float64};
end

""" 
---  
Returns the desired window used in transmitter or receiver weigthed overlap and add methods (WOLA). See [1] for window reference design.
# --- Syntax 
window = getWolaWindow(winFunc,nFFT,nCp,winlength)
# --- Input parameters
- winFunc	: Window type (see supported format below) [String]
- nFFT		: FFT size [Int]
- nCP		: CP size [Int]
- winlength : length of window [Int]
# --- Output parameters
window	: Window [Window]
# ---
# Supported window
- "Triangle"		: Triangle window
- "srrc"			: Square Root Raised Cosine
- "Meyer"			: Meyer window (See [1])
# ---
# References
- [1] R. Zayani, Y. Medjahdi, H. Shaiek and D. Roviras, "WOLA-OFDM: A Potential Candidate for Asynchronous 5G," 2016.
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getWolaWindow(winFunc,nFFT,nCP,winlength)
	# --- Overall parameters
	sizeSymb		= nFFT + nCP;
	# --- Switch to lowercase to avoid case mistakes.
	winName = lowercase(winFunc);
	if winFunc == "triangle"
		# ----------------------------------------------------
		# --- Triangle window
		# ----------------------------------------------------
		# Classic ramp up // ramp down window
		filterTimeUp  = 1/(winlength-1) * collect(0:1:winlength-1);
		filterTimeDown= 1/(winlength-1) * collect(winlength-1:-1:0);
		window      = [filterTimeUp ; ones(Float64,sizeSymb-winlength) ; filterTimeDown];
	elseif 	winFunc == "meyer"
		# ----------------------------------------------------
		# --- Meyer window
		# ----------------------------------------------------
		ww      = Int.(floor(winlength)/2);
		# --- Up
		n       = collect(0:ww-1) ./(ww-1);
		V       = n.^4 .*(35 .- 84n .+ (70 .*n.^2) .- (20 .*n.^3));
		wUp     = 1/2 .- 1/2 .* cos.( pi .* V);
		# --- Down
		wDown     = collect(wUp[ end:-1:1 ]);
		# --- Complete filter
		window                = [wUp;ones(Float64,sizeSymb) ;wDown];
	else
		error("Unsupported apoodisation window");
	end
	return Window(winFunc,winlength,window)
end
