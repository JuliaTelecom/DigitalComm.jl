# ----------------------------------------------------
# --- BF-OFDM structure
# ----------------------------------------------------
"""
---  
BF-OFDM structure 
# --- Syntax 
- nFBMC		: Number of carriers (PPN size)
- nOFDM		: Number of subbcarriers per carrier (OFDM precoder size)
- K			: Overlapping factor of the PPN 
- GI		: CP size of the precoder 
- δ			: Rate factor (compression factor)  
- allocatedSubcarriers	: Vector of allocated subcarriers   
- filterName  : Type of filter used (name)  
- filterTaps  : FIR coefficient values 
- BT		: Gaussian parameters (if gaussian filter is used)
- filterStopBand : DolphChebyshev attenuation factor (if DC is used)
- fS		: Frqeuency spreading coefficients (if FS filter is used)
- nFFT		: Equivalent OFDM FFT size 
- nCP		: Equivalent OFDM CP size 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
struct StrucBFOFDM<: Waveform
	nFBMC::Int;
	nOFDM::Int;
	K::Int;
	GI::Int;
	δ::Float64;
	allocatedSubcarriers::Array{Int};
	filterName::String;
	filterTaps::Array{Float64};
	BT::Float64;
	filterStopBand::Float64;
	fS::Array{Float64};
	nFFT::Int;
	nCP::Int;
end

""" 
---  
BF-OFDM initialisation
# --- Syntax 
bfofdm = initBFOFDM(nFBMC::Int,nOFDM::Int,K::Int,GI::Int,δ::Float64,allocatedSubcarriers::Array{Int},filterType::String;BT=-1,filterStopBand=-1,fS=[],nFFT=-1,nCP=-1)
# --- Input parameters 
- nFBMC		: Number of carriers (PPN size)
- nOFDM		: Number of subbcarriers per carrier (OFDM precoder size)
- K			: Overlapping factor of the PPN 
- GI		: CP size of the precoder 
- δ			: Rate factor (compression factor)  
- allocatedSubcarriers	: Vector of allocated subcarriers   
- filterName  : Type of filter used (name)  
- filterTaps  : FIR coefficient values 
- BT		: Gaussian parameters (if gaussian filter is used)
- filterStopBand : DolphChebyshev attenuation factor (if DC is used)
- fS		: Frqeuency spreading coefficients (if FS filter is used)
- nFFT		: Equivalent OFDM FFT size 
- nCP		: Equivalent OFDM CP size 
# --- Output parameters 
bfofdm		: BF-OFDM structure 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function initBFOFDM(nFBMC,nOFDM,K,GI,δ,allocatedSubcarriers,filterType;BT=-1,filterStopBand=-1,fS=[],nFFT=-1,nCP=-1)
	# --- Create direclty the filter taps + filter name.
	(filterTaps,fS) = getBFOFDMFilter(nFBMC,K,filterType,GI=GI,BT=BT,filterStopBand=filterStopBand,fS=fS);
	# --- Create the BF-OFDM structure
	return StrucBFOFDM(nFBMC,nOFDM,K,GI,δ,allocatedSubcarriers,filterType,filterTaps,BT,filterStopBand,fS,nFFT,nCP)
end



# ----------------------------------------------------
# --- Main BF-OFDM signal generation method
# ----------------------------------------------------
"""  
--- 
Creates a Block Filtered - OFDM (BF-OFDM) signal parametrized by its numerlogy and its waveform parameter. Generate a time domain signal based on the input matrix qamMat. The input matrix is a T/F matrix of size (nRe x nbSymb) with nRe the number of allocated subcarrier and nbSymb the number of symbols.

The waveform is parametrized by the PPN size (number of carriers) nFBMC, the size of the precoding stage nOFDM, the CP size of the precoding stage GI and the compression factor parameter δ. See [1], [2] for waveform description and [3,5] for importance and design of compression rate (set to 0.5 in former works on BF-OFDM).

BF-OFDM is characterized by its pulse shape filter. Pulse shape can be automatically tuned to optimal (in terms of  innterference management) [4]. In such a case, filterType should be a string (can be "gaussian_opt", "dc_opt" and "fs_opt"). For a specific pulse shape filter, an array of the nFBMC*K filter taps should be given (and generated for example  with getBFOFDMFilter).
# --- Syntax
sigId	  =bfofdmSigGen(qamMat,nFBMC,nOFDM,K,GI,δ,allocatedSubcarriers,filterType;typeTx="PPN")
# --- Input parameters
- qamMat		: T/F symbols to transmit (QAM symbols) [Array{Complex{Float64},nRe,nbSymb] with nbSymb number of BF-OFDM symbols and nRe number of allocated subcarriers.
- nFBMC		: PPN size
- nOFDM		: OFDM precoding size
- GI			: CP size of precoding stage
- K			: Overlapping factor
- δ			: Compression rate
- allocatedSubcarriers : Vector or allocated subcarriers Array{Int,nRe} with maximum(nRe) < nOFDM*nFBMC*δ
- filterType	: Pulse shape type. Can be a string (see getBFOFDMFilter) or a vector of size K x nFBMC with filter taps.
- typeTx		: Filterbank architecture. Can be "PPN" (polyphase network) or "FS" (frequency spreading) if filterType is of type "fs" or "fs_opt".
# --- Output parameters
- sigId		: Time domainn BF-OFDM signal. [Array{Complex{Float64},nChip}.
# ---
# References
- [1] Gerzaguet, R and Demmer, D and Doré, J-B and Le Ruyet, D. and Kténas, D, "Block-Filtered OFDM: A new Promising Waveform for Multi-service Scenarios", 2017
- [2] Demmer, D and Gerzaguet, R and Doré, J-B and Le Ruyet, D. and Kténas, D, "Block-filtered OFDM: A novel waveform for future wireless technologies", 2017.
- [3] Demmer, D and Rostom, Z and Gerzaguet, R and Doré, J-B and Le Ruyet, D. "Study of OFDM Precoded Filter-Bank Waveforms", 2018
- [4] Demmer, D and Gerzaguet, R and Doré, J-B and Le Ruyet, D. and Kténas, D, "Filter Design for 5G BF-OFDM Waveform", 2017.
- [5] Demmer, D.; Zakaria, R.; Gerzaguet, R.; Doré, J. & Le Ruyet, D. Study of OFDM Precoded Filter-Bank Waveforms, IEEE Transactions on Wireless Communications, 2019.
# v 1.0 - Robin Gerzaguet.
"""
 function bfofdmSigGen(qamMat::Array{Complex{Float64}},nFBMC::Int,nOFDM::Int,K::Int,GI::Int,δ::Float64,allocatedSubcarriers::Array{Int},filterType;typeTx="PPN")
	 # ----------------------------------------------------
	 # --- Overall parameters
	 # ----------------------------------------------------
	 # Set 4 core for FFT computation
	 #FFTW.set_num_threads(4)
	 # --- Define compression alloc
	 nOPB				= Int(nOFDM * δ );
	 nbDataSubcarriers	= length(allocatedSubcarriers);
	 nbSymb				= size(qamMat,2);

	 # ----------------------------------------------------
	 # --- Subcarrier and carrier management
	 # ----------------------------------------------------
	 # We have the subcarrier indexes, as with CP-OFDM
	 # --> We should deduce which carriers is enable (i.e with PPN stage holds data)
	 # --> We also should deduce how the oversampled grid (nOFDM*nFBMC) bears data.
	 # --- Getting PPN alloc.
	 fbmcCarriers	= getBFOFDM_carrierFromSubcarriers(nFBMC, nOFDM, δ , allocatedSubcarriers);
	 nbDataCarriers	= length(fbmcCarriers);
	 # --- Getting oversampling grid index
	 # Indead in BF-OFDM the precoding stage are IFFT blocks of size nOFDM on which nOFDM*δ are allocated.
	 # BF-OFDM can be seen as a global grid of nOFDM*nFBMC, with nRe allocated subcarriers.
	 # Purpose is to get the nRe vector that points which elemnts are allocated on the nOFDM*nFBMC grid
	 allSubcarriers	= getBFOFDM_oversampledGridSubcarriers(nFBMC,nOFDM,δ,fbmcCarriers);
	 # In subcarrier mode, we can have only part of a carrier that is allocated
	 # In this case, size(qamMap,2) is not equal to nOFDM*δ*nbDataCarriers
	 # The given matrix have always lower number of subcarriers, so zero padding approach here
	 # We should remap with 0, and extract only submatrixe of appropriate size
	 if size(qamMat,1) != nOFDM*δ*nbDataCarriers
		 qTmp		= zeros(Complex{Float64},Int(nOFDM*δ*nFBMC),nbSymb);
		 # --- Get all the index associadted to fbmcCarriers alloc.
		 carrierAllmapped = getBFOFDM_subCarrierFromFBMCCarriers(nFBMC, nOFDM, δ, fbmcCarriers)
		 # --- Mapp in the total T/F matrix
		 qTmp[ allocatedSubcarriers,: ]	= qamMat;
		 # --- Getting only full carrier allocated part
		 qamMat = qTmp[ carrierAllmapped,: ];
	 end
	 # At this stage qamMat should be [Complex{Float64},nbDataCarriers*nOFDM*δ,nbSymb]

	 # ----------------------------------------------------
	 # --- Getting filter  and Pre-distortion coefficients
	 # ----------------------------------------------------
	 # --- Time impulse response
	 if isa(filterType,String)
		 # --- Filter is defined as optimal, we launch getBFOFDMFilter without additional parameters
		 # This can only be done with optimal filters, but this is checked in BFOFDM_filters.
		 # Lazy approach here, launch directly without proper check.
		 pulseShape,	  = getBFOFDMFilter(nFBMC,K,filterType,GI=GI);
		 pulseShape	  = 1/(K*nFBMC) .*pulseShape;
	 else
		 # --- Filter is given as filter taps
		 # We check the length
		 if length(filterType) != K * nFBMC
			 error("Unapproriate filter for BF-OFDM: It seems that  filter taps are given, but the length of the filter does not match the PPN expected size (getting %d instead of %d)",length(filterType),K*nFBMC);
		 end
		 # Map to pulseShape
		 pulseShape	  = 1/(K*nFBMC) .* filterType;
	 end
	 # --- Frequency response
	 # Frequency response is necessary to compute the pre-distortion coefficient (based on inverse of frequency response)
	 # Calculate the frequency response of the filter
	 # FFT is done on nOFDM*n FBMC *δ , should be doing zero padding
	 # We assume that nOB > K as otherwise Near Perfect Reconstruction (NPR) cannnot be achieved.
	 pZP	= [pulseShape;zeros(Complex{Float64},nFBMC*nOPB - K * nFBMC)];
	 G		= fft(pZP,1);
	 # --- Calculate pre-distortion coefficient
	 # Get frequency location where subcarrier are effectively mapped and calculate inverse. See [2].
	 zoneInt	= [G[ end-nOPB÷2+1:end ] ;G[ 1:nOPB÷2 ]];
	 Gcorr		= 1 ./abs.(zoneInt).^2 .*conj.(zoneInt);
	 # Same PD coefficients for each FBMC carrier
	 GCorrFreq  = repeat(Gcorr,nbDataCarriers);
	 # --- Apply Predistortion
	 qamMap		= mapslices( x-> x.* GCorrFreq,qamMat,dims=1);

	#----------------------------------------------------
	## --- OFDM preprocessing
	#----------------------------------------------------
	#TODO Adding window for BF-OFDM
	winLenght   = 0;
	# --- OFDM modulation + GI insertion
	output_ifft_full			  = zeros(Complex{Float64},nbDataCarriers,(GI+nOFDM)*nbSymb+winLenght);
	# --- Init a signal with Extended size for windowing
    sizeSymb         = nOFDM+GI;
	sizeExt			 = sizeSymb + winLenght;
	sigExtend		 = zeros(Complex{Float64},nbDataCarriers,sizeExt);
	filterTx         = ones(Int, sizeExt);
	wLD2             = winLenght÷2;
	apodWindow       = transpose(repeat(filterTx,1,nbDataCarriers));
	for idx_j=1:nbSymb
		# --- Getting data of poreprocessing blocks
		input_ifft                                = reshape(qamMap[ :,idx_j ],nOPB ,nbDataCarriers);
		# Mapping carriers depending on delta value
		input_fft_padded                          = zeros(Complex{Float64},nOFDM,nbDataCarriers);
		# --- Map data to appropriate location for each carrier
		input_fft_padded[ allSubcarriers ]        = input_ifft;
		# IFFT preprocessing
		output_ifft			                      = transpose(ifft(input_fft_padded,1));

		#FIXME Speed up with CP insertion with matrix approach.
		# --- Adding CP
		sigExtend[ :,1:GI+wLD2 ]	              = output_ifft[ :,1+end-(GI+wLD2):end ];
		# --- Adding signal
		sigExtend[ :,GI+wLD2 .+ collect(1:nOFDM) ]= output_ifft;
		# --- Adding CS
		sigExtend[ :,end+1-wLD2:end ]             = output_ifft[ 1:wLD2 ];

		# --- Apodisation
		sigApod                                   = sigExtend .* apodWindow;
		# --- Folding in output vector
		output_ifft_full[ :,1 + (idx_j-1)*sizeSymb: idx_j*sizeSymb+ winLenght ] = output_ifft_full[ :,1 + (idx_j-1)*sizeSymb : idx_j*sizeSymb+ winLenght ] .+ sigApod;
	end
	# --- Map preprocessing block to synthesis block
	synthesisMat					= zeros(Complex{Float64},nFBMC, (nOFDM+GI)*nbSymb);
	synthesisMat[fbmcCarriers,:]  	= output_ifft_full;

	# ----------------------------------------------------
	## --- Filter bank stage
	# ----------------------------------------------------
	# --- Init vector
	#TODO Frequency spreading architecture
	# --- Windowing is applied for each symbol
	pKron	= repeat(nFBMC*pulseShape,1,(nOFDM+GI)*nbSymb+winLenght);
	G       =  [];
	# --- Repeat and copy in polyphase filter
	x		= repeat(ifft(synthesisMat,1),K,1);
	# --- Apply filter
	x		= x .* pKron;
	# --- Overlap and Add
	nbFBMCSymbTotal     = (nOFDM+GI) * nbSymb + winLenght;
	y = zeros(Complex{Float64},Int((nbFBMCSymbTotal-1)*nFBMC*δ +nFBMC*K));
	for k=1: nbFBMCSymbTotal
		#FIXME Heavy line in term of computation -> switch to C ?
		y[ (k-1)*Int(nFBMC*δ)  .+ collect(1:( nFBMC*K)) ] = y[ (k-1)*Int(nFBMC*δ)  .+ collect(1:( nFBMC*K)) ] .+ (x[:,k]);
	end

	# ----------------------------------------------------
	# --- Output
	# ----------------------------------------------------
	return y;
 end


 # ----------------------------------------------------
 # ---  Multiple dispatch for structure managment
 # ----------------------------------------------------
 function bfofdmSigGen(qamMat,bfofdm::StrucBFOFDM)
	 return bfofdmSigGen(qamMat,bfofdm.nFBMC,bfofdm.nOFDM,bfofdm.K,bfofdm.GI,bfofdm.δ,bfofdm.allocatedSubcarriers,bfofdm.filterTaps);
 end
