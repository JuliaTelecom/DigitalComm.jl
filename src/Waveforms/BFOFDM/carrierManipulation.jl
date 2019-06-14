# Usefull functions for subcarrier and carrier mappign in BF-OFDM

""" 
---  
Get the index of the oversampled carrier (i.e subcarriers index in the precoding field) based on the index of the allocated carriers. It is assumed that for each carriers, all subcarriers are enable (with orthogonality, i.e Nδ allocated subcarriers).
# --- Syntax 
subcarriers	= getBFOFDM_subCarrierFromFBMCCarriers(nFBMC,nOFDM,δ,fbmcCarriers)
# --- Input parameters
- nFBMC	  : FBMC carrier size (number of carriers) [Int]
- nOFDM	  : OFDM subcarrier size [Int]
- δ		  : Rate factor [Float64]
- fbmcCarriers : Vector of allocated carrier [Array{Int,L}], L < nFBMC
# --- Output parameters
- subcarriers : Vector of allocated subcarriers [Array{Int, L δ ]
# v 1.0 - Robin Gerzaguet.
"""
function getBFOFDM_subCarrierFromFBMCCarriers(nFBMC, nOFDM, δ, fbmcCarriers)
	# ----------------------------------------------------
	# --- Getting parameters
	# ----------------------------------------------------
	nbDataCarrier	= length(fbmcCarriers);
	subcarrierPerCarrier  = Int(nOFDM*δ);
	# --- Index of occupied carriers in the oversampled plan
	carrierOversamp	      = zeros(Int,nbDataCarrier*subcarrierPerCarrier);
	for idx = 1:nbDataCarrier
		carrierOversamp[ collect(1+(idx-1)*Int(nOFDM*δ):idx*Int(nOFDM*δ)) ]  = 1 .+ mod.((fbmcCarriers[idx]-1)*Int(nOFDM*δ)-Int(nOFDM*δ/2)  .+  collect(1:Int(nOFDM*δ)) .-1,Int(nOFDM*nFBMC*δ));
		end
	return carrierOversamp
end

""" 
---  
Get the index of the allocated carriers based on the index of the subcarriers. We consider here that a FBMC carriers bear data (i.e is allocated) is at least one of its subcarrier is enable.
# --- Syntax 
fbmcCarriers	= getBFOFDM_carrierFromSubcarriers(nFBMC,nOFDM,δ,subcarrier)
# --- Input parameters 
- nFBMC	  : FBMC carrier size (number of carriers) [Int]
- nOFDM	  : OFDM subcarrier size [Int]
- δ		  : Rate factor [Float64]
- fbmcCarriers : Vector of allocated carrier [Array{Int,L}], L < nFBMC*nOFDM*δ 
# --- Output parameters 
- fbmcCarriers : Vector of allocated carriers [Array{Int, P ] P < nFBMC
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getBFOFDM_carrierFromSubcarriers(nFBMC,nOFDM,δ,subcarriers)
	# ----------------------------------------------------
	# --- Getting parameters
	# ----------------------------------------------------
	# --- Allocate worst case scenario
	subcarrierPerCarrier =  nOFDM * δ;
	allCarriers = zeros(Int,nFBMC);		# --- Carrier buffer
	nbC         = 1;					# --- Carrier counter
	# --- Iterative loop
	for iN = 1 : 1 : nFBMC
		# --- For each carrier, we calculate associated "subcarrier index length
		# If at least 1 index of the subcarrier index is in subcarrier, we decide to allocate the index.
		# --- Calculation of subcarrier index for current FBMC carrier
		tmp	 = 1 .+ mod.( (iN-1)*Int(nOFDM*δ) .- Int(nOFDM*δ/2) .+ collect(1:nOFDM*δ) .-1,nOFDM*nFBMC*δ );
		# --- If at least one index belongs to subcarrier and tmp, store as an allocated carrier
		tmpInt	= intersect(tmp,subcarriers);
		if length(tmpInt) != 0
			# --- Savign currennt carrier index
			allCarriers[nbC] = iN;
			# --- Setting next element for saving next carrier
			nbC = nbC + 1;
		end
	end
	# --- Returnns onnly allocated vector
	return allCarriers[1:nbC-1];
end

""" 
---  
Returns the oversampled grid at the subcarrier level for BF-OFDM. In BF-OFDM as we use a nOFDM*nFBMC grid with nOFDM*δ allocated subcarrier per carriers (for a precodinng stage of size nOFDM), it can be usefull to get all the index of the allocated subcarrier per carrier. The output is a vector of the allocated subcarrier index in the oversampled frequency grid.
# --- Syntax 
subcarrierTx	= getCarrierFromSubcarriers(nOFDM,δ,fbmcCarriers)
# --- Input parameters
- nOFDM	  : OFDM subcarrier size [Int]
- δ		  : Rate factor [Float64]
- fbmcCarriers : Vector of allocated carrier [Array{Int,L}], L < nFBMC*nOFDM*δ
# --- Output parameters
- subcarrierTx : Vector of allocated carriers [Array{Int, L*nOFDM*δ  ]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getBFOFDM_oversampledGridSubcarriers(nFBMC,nOFDM,δ,fbmcCarriers)
	# --- Getting alloc. parameters
	nbDataCarrier	= length(fbmcCarriers);
	subcarrierPerCarrier  = Int(nOFDM*δ);
	# --- Index of occupied carriers in the oversampled plan
	subcarrierTx          = zeros(Int,nbDataCarrier*subcarrierPerCarrier);
	for idx = 1:nbDataCarrier
		# --- Calculating index of subcarriers index into carriers
		# It depends on delta (delta mapping pattern exists)
		# We also have a shift which is unknown
		# --- Specific pattern index for given carrier
		tmp  = mod.((-Int(nOFDM*δ/2)+(fbmcCarriers[idx]-1)*subcarrierPerCarrier: Int(nOFDM*δ/2) +(fbmcCarriers[idx]-1)*subcarrierPerCarrier-1), nOFDM) .+1;
		# --- Store this pattern as function of carrier index
		# No need to have a absolute positionning as after that we will use
		# only matrix with data
		subcarrierTx[collect(1 + (idx-1)*Int(nOFDM*δ):idx*Int(nOFDM*δ)) ] = tmp .+ nOFDM*(idx-1);
	end
	return subcarrierTx
end
