""" 

Returns prototype fitlering matrix for BF-OFDM  with specified by overlapping factor K, PPN size nFBMC, filter type. Additional parameters are the CP size of the precoding stage (set to optimal value by default) and filter structure for specific parameter overset.
# --- Syntax 
(p,pF) = getBFOFDMFilter(K,nFBMC,filter;GI,BT,filterStopBand,fS)>
# --- Input parameters
- K	  : Overlapping factor [Int]
- nFBMC : PPN size (Number of FBMC carriers) [Int]
- filter: Filter type: "Gaussian","phydyas","FS","DC": [String]
- BT	  : Bandwidth time product for custom gaussian filter (used with "gaussian") [Float64]
- filterStopBand : Attenuation with custom DC window (used with "DC") [Float64]
- fS	  : Frequency sampling custom  coefficients (Array{Float64,K}) (used with "fs")
# --- Output parameters
- p		: Filter impulse response [Array{Float64,K*nFBMC}].
- pF		: Filter coefficient in frequency domain (For Frequency sampling implementation)
---
Filter types
All filter type can be called with 2 ways. "filter" and "filter_opt". When "filter_opt" is used, the filter coefficient are obtained through numerical optimisation (SIR optimisation). Otherwise, it uses its associated parameter to extract appropriate coefficients. See [1] and [2] for filter optimisation principles. See [4] for special cases about 5G-NR compatibilty. 
The supported filter/windows are:
- Gaussian	: Gaussian filter shape, specified by BT. Optimized Gaussian shape can be used [1]
- phydyas		: Classic FBMC pulse shape defined in frequency domain [3]
- FS			: Coefficient defined in frequency domain. Filter coefficient are based on intrinsic SIR optimisation [1]
- DC			: Dolph Chebyshev window. Can be based on window optimisation (SIR optimisation) or defined by filterStopBand parameter.
# ---
# References
- [1]	Demmer, D and Gerzaguet, R and Doré, J-B and Le Ruyet, D. and Kténas, D, "Filter Design for 5G BF-OFDM Waveform", 2017.
- [2] A. Sahin, I. Guvenc, and H. Arslan, “A Survey on Multicarrier Communications: Prototype Filters, Lattice Structures, and Implementation Aspects", 2014
- [3] Phydyas project, "Deliverable D5.1 : Prototype filter and structure optimization", 2009
- [4]	Demmer, D and Rostom, Z and Gerzaguet, R and Doré, J-B and Le Ruyet, D. "Study of OFDM Precoded Filter-Bank Waveforms", 2018
# Examples
Get a classic FBMC filter with phydyas and an overlapping factor of 4, with a input PPN of size 64. The precodinng stage is a OFDM of size 64 and a GI of size 4.

	# --- phydyas filter
	( hphydyas,hFphydyas )	 = getBFOFDMFilter(4,64,"phydyas");

	# --- Gaussian shape
	( hGaussian,hFGaussian )	 = getBFOFDMFilter(4,64,"Gaussian_opt",GI=4);

	# --- Custom Gaussian window
	( hGaussianC,hFGaussianC ) = getBFOFDMFilter(4,64,"Gaussian",GI=4,BT=0.5)

# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getBFOFDMFilter(nFBMC,K,filter;GI=-1,BT=-1,filterStopBand =-1,fS=[])
	# --- Avoid type issue with lowercase casting
	filterSwitch = lowercase(filter);
	# ----------------------------------------------------
	# --- Filter switch
	# ----------------------------------------------------
	if filterSwitch == "phydyas"
		# ----------------------------------------------------
		# --- phydyas prototype filter [3]
		# ----------------------------------------------------
		# --- Getting frequency coefficients
		if K == 2
			p1			= sqrt(2)/2;
			P			= p1;
		elseif K == 3
			p1			= 0.91143783;
			P			= [p1;sqrt(1-p1^2)];
		elseif K == 4
			p1			= 0.97195983;
			p2			= 1/sqrt(2);
			p3          = sqrt(1-p1^2);
			P			= [p1;p2;p3];
		elseif K == 6
			p1          = 0.99722723;
			p2          = 0.94136732;
			P           = [ p1;p2 sqrt(2)/2;sqrt(1-p2^2);sqrt(1-p1^2)];
		elseif  K== 8
			p1          = 0.99988389;
			p2          = 0.99315513;
			p3          = 0.92708081;
			P          = [ p1;p2;p3;sqrt(2)/2;sqrt(1-p3^2);sqrt(1-p2^2);sqrt(1-p1^2)];
		else
			error("For Phydyas prototype filter, supported overlapping factor are 2, 3, 4, 6 or 8");
		end
		# --- Getting time impulse response
		# Time indexes
		t	      = collect(1:K*nFBMC);
		# Frequency vector
		P_K	      = [P[end:-1:1];1;P];
		# Init filter time impulse reponse
		p         = zeros(Float64,K*nFBMC);
		for k     = K+1 : 2*length(P)+1
			# --- Iterative generation (see [3])
			p[ t ]  = p[ t ] .+  2*P_K[ k ]*(-1)^(k-K) .* cos.(2*pi/(K*nFBMC)*(k-K) .*(t .- 1));
		end
		p         =p .+ P_K[ K ];
		# --- Final impulse response
		imp_resp  = (-1)^K * p;
		# --- Final frequency response
		freq_coef = P_K;
		# --- Output
		return (imp_resp,freq_coef)
	elseif filterSwitch == "gaussian_opt"
		#FIXME: Optimize value for any (K,GI) tuple
		# ----------------------------------------------------
		# --- Gaussian optimized window [1]
		# ----------------------------------------------------
		# To use an optimized Gaussian function as it depends on GI size, we have to use a double switch (on K and on GI)
		if K == 2
			if GI ==  2
				a = sqrt(log(2)/2)/(0.40);
			elseif GI ==  1
				a = sqrt(log(2)/2)/(0.45);
			elseif GI == 9 
				# 5G-NR value [3]
				a = sqrt(log(2)/2)/(1.08);
			else
				# --- Of other value (>), it corresponds to the optimized Gaussian shape without ISI (2K-1)
				a = sqrt(log(2)/2)/(0.38);
			end
		elseif K == 3
			if GI ==  9 || GI == 18 || GI == 27 GI == 36
				# --- 5G NR with 0.25 compression rate
				a = sqrt(log(2)/2)/(0.3);
			elseif GI ==  5
				a = sqrt(log(2)/2)/(0.32);
			elseif GI ==  4
				a = sqrt(log(2)/2)/(0.33);
			elseif GI ==  2
				a = sqrt(log(2)/2)/(0.40);
			elseif GI ==  1
				a = sqrt(log(2)/2)/(0.40);
			end
		elseif K == 4
			if GI ==  7
				a = sqrt(log(2)/2)/(0.28);
			elseif GI ==  4
				a = sqrt(log(2)/2)/(0.33);
			elseif GI ==  3
				a = sqrt(log(2)/2)/(0.357);
			elseif GI ==  2
				a = sqrt(log(2)/2)/(0.40);
			else
				# --- Raise an error for other values
				error("unknown optimized filter for GI (gaussian)");
			end
		elseif K == 6
			a = sqrt(log(2)/2)/(0.327);
		else
			error("Unknown optimized filter for K (gaussian)");
		end
		# --- Time index
		t = range(-K/2,stop=K/2,length=K*nFBMC);
		# --- Gaussian filter shape
		p = sqrt(pi)/a*exp.(-(pi/a .* t).^2);
		# --- Output
		imp_resp  = p;
		freq_coef = [];
		return(imp_resp,freq_coef)
	elseif filterSwitch == "gaussian"
		# ----------------------------------------------------
		# --- Manual Gaussian filter
		# ----------------------------------------------------
		if BT == -1
			# --- Corresponds to the default value for BT.
			# Raise an error, as BT should be defined for a specific Gaussian filter
			error("For custom Gaussian filter, BT should be defined. For an optimized Gaussian shape, use gaussian_opt instead");
		end
		a = sqrt(log(2)/2)/(BT);
		# --- Time index
		t = range(-K/2,stop=K/2,length=K*nFBMC);
		# --- Gaussian filter shape
		p = sqrt(pi)/a*exp.(-(pi/a .* t).^2);
		# --- Output
		imp_resp  = p;
		freq_coef = [];
		return(imp_resp,freq_coef)
	elseif filterSwitch == "fs_opt"
		# ----------------------------------------------------
		# --- Frequency sampling optimized filter
		# ----------------------------------------------------
		# --- Get frequency coefficients
		if K == 4
			if GI ==  7
				P = [0.72071;0.25452;0.03380];
			elseif GI ==  4
				P = [0.792;0.375;0.082];
			elseif GI ==  2
				P = [0.870;0.581;0.211];
			else
				error("unknown optimized filter for GI (fs_opt)");
			end
		elseif K == 3
			if GI == 5
				P = [0.633;0.133];
			elseif GI ==  4
				P = [0.656;0.154];
			elseif GI ==  2
				P = [0.806;0.301];
			else
				error("unknown optimized filter for GI (fs_opt)");
			end
		elseif K == 2
			if GI == 3
				P = 0.475;
			elseif GI ==  2
				P = 0.485;
			elseif GI ==  1
				P = 0.524;
			else
				error("unknown optimized filter for GI (fs_opt)");
			end
		else
			error("unknown optimized filter for K (fs_opt)");
		end
		# --- Deduce filter taps
		t	  = collect(1:K*nFBMC);
		P_K	  = [P[ end:-1:1 ];1;P];
		p = zeros(Float64,K*nFBMC);
		for k = K+1 : 2*length(P)+1
			p[ t ] = p[ t ] +2*P_K[ k ]*(-1)^(k-K) .* cos.(2*pi/(K*nFBMC)*(k-K) .*(t .- 1));
		end
		p =p .+ P_K(K);
		# --- Final impulse response
		imp_resp  = (-1)^K * p;
		# --- Final frequency response
		freq_coef = P_K;
		# --- Output
		return (imp_resp,freq_coef)
	elseif filterSwitch == "fs"
		# ----------------------------------------------------
		# --- Customized Frequency sampling filter
		# ----------------------------------------------------
		if isempty(fS)
			# --- Corresponds to the default value for freuency points.
			# Raise an error, as the vector should be defined properly.
			error("For customFrequency sampling filter, fS should be defined. For an optimized FS system, use fs_opt instead");
		end
		t	  = collect(1:K*nFBMC);
		P_K	  = fS;
		p = zeros(Float64,K*nFBMC);
		for k = K+1 : 2*length(P)+1
			p[ t ] = p[ t ] +2*P_K[ k ]*(-1)^(k-K) .* cos.(2*pi/(K*nFBMC)*(k-K) .*(t .- 1));
		end
		p =p .+ P_K(K);
		# --- Final impulse response
		imp_resp  = (-1)^K * p;
		# --- Final frequency response
		freq_coef = P_K;
		# --- Output
		return (imp_resp,freq_coef)
	elseif filterSwitch == "dc_opt"
		# ----------------------------------------------------
		# --- Dolph Chebyshev optimized window
		# ----------------------------------------------------
		# --- Getting optimized stop band attenuation
		if K == 6
			if G == 7
				filterStopBand = 159;
			elseif GI ==  4
				filterStopBand = 340;
			elseif GI ==  2
				filterStopBand = 317;
			else
				error("unknown optimized filter for GI (DC)");
			end
		elseif K == 4
			if GI ==  7
				filterStopBand = 159;
			elseif GI ==  4
				filterStopBand = 220;
			elseif GI ==  2
				filterStopBand = 317;
			else
				error("unknown optimized filter for GI (DC)");
			end
		elseif K == 3
			if GI == 5
				filterStopBand = 119.30;
			elseif GI ==  4
				filterStopBand = 129.19;
			elseif GI ==  2
				filterStopBand = 182.00;
			else
				error("unknown optimized filter for GI (DC)");
			end
		elseif K ==2
			if GI ==  3
				filterStopBand = 80.04;
			elseif GI ==  2
				filterStopBand = 86.31;
			elseif GI ==  1
				filterStopBand = 106.18;
			else
				error("unknown optimized filter for GI (DC)");
			end
		else
			error("unknown optimized filter for K (DC)");
		end
		# --- Final impulse response
		imp_resp  = dolphChebyshev(K*nFBMC,Float64(filterStopBand));
		# --- Final frequency response
		freq_coef = [];
		# --- Output
		return (imp_resp,freq_coef)
	elseif filterSwitch == "dc"
		# ----------------------------------------------------
		# --- Customized Dolph Chebyshev window
		# ----------------------------------------------------
		if filterStopBand == -1
			# --- Corresponds to the default value for BT.
			# Raise an error, as BT should be defined for a specific Gaussian filter
			error("For custom Dolph Chebyshev window, filterStopBand should be defined. For an optimized DC window, use dc_opt instead");
		end
		# --- Final impulse response
		imp_resp  = dolphChebyshev(K*nFBMC,Float64(filterStopBand));
		# --- Final frequency response
		freq_coef = [];
		# --- Output
		return (imp_resp,freq_coef)
	end
end



