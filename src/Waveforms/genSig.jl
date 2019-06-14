""" 
---  
Generate a time domain signal for a given multicarrier waveform. Input are the Time-Frequency matrix qamMat and the Waveform structure.
# --- Syntax 
sigId	= genSig(qamMat,struc::Waveform)
# --- Input parameters 
- qamMat  : Time frequency complex QAM matrix [Array{Complex{Float64}}]
- struc	  : Waveform structure 
# --- Output parameters 
- sigId	  : Time domain signal 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genSig(qamMat,struc::StrucOFDM)
	return ofdmSigGen(qamMat,struc);
end
function genSig(qamMat,struc::StrucSCFDMA)
	return scfdmaSigGen(qamMat,struc);
end
function genSig(qamMat,struc::StrucUFOFDM)
	return ufofdmSigGen(qamMat,struc);
end
function genSig(qamMat,struc::StrucBFOFDM)
	return bfofdmSigGen(qamMat,struc);
end
function genSig(qamMat,struc::StrucWOLA)
	return wolaSigGen(qamMat,struc);
end
function genSig(qamMat,struc::StrucFBMC)
	return fbmcSigGen(qamMat,struc);
end


""" 

Decode a time domain signal for a given multicarrier waveform and returns a T/F QAM constellation. Input are the compelx baseband signal  and the Waveform structure.
# --- Syntax 
	qamDec	= decodeSig(signal,struc);
# --- Input parameters 
- signal  : Time domain signal 
- struc	  : Waveform structure 
# --- Output parameters 
- qamDec  : Time frequency complex QAM matrix [Array{Complex{Float64}}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function decodeSig(signal,struc::StrucOFDM)
	return ofdmSigDecode(signal,struc);
end
function decodeSig(signal,struc::StrucSCFDMA)
	return scfdmaSigDecode(signal,struc);
end
function decodeSig(signal,struc::StrucUFOFDM)
	return ufofdmSigDecode(signal,struc);
end
function decodeSig(signal,struc::StrucBFOFDM)
	return bfofdmSigDecode(signal,struc);
end
function decodeSig(signal,struc::StrucWOLA)
	return wolaSigDecode(signal,struc);
end
function decodeSig(signal,struc::StrucFBMC)
	return fbmcSigDecode(signal,struc);
end


""" 
---  
Returns the waveform name based on input type structure
# --- Syntax 
 name = getWaveformName(struc)
# --- Input parameters 
- struc	  : Waveform structure [Waveform]
# --- Output parameters 
- name	  : String associated to waveform name
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function getWaveformName(s::Waveform)
	if typeof(s) == StrucOFDM
		name = "OFDM";
	elseif typeof(s) == StrucSCFDMA
		name = "SCFDMA";
	elseif typeof(s) == StrucUFOFDM
		name = "UFOFDM";
	elseif typeof(s) == StrucWOLA
		name = "WOLA";
	elseif typeof(s) == StrucFBMC 
		name = "FBMC";
	elseif typeof(s) == StrucBFOFDM 
		name = "BFOFDM";
	end
end

""" 
---  
Create a dictionnary of waveform configurations. To compare and use the same script for different waveform configuration, we propose to add a dictionnary to have a container that contains all waveform configuraton. The function is called with every desired waveform structure. The dictionnary as a key associated to the waveform name, and a field associated to the waveform structure. If the waveform is present several times (several configuration with same waveform type, for instance FBMC with different overlapping factor values), a counter index is added to the waveform key.
# --- Syntax 
	dWav  = initWaveforms(x1,x2,...)
# --- Input parameters 
- x		: Waveform structure [Waveform]
# --- Output parameters 
- dWav	: Dictionnary of waveforms.
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function initWaveforms(x...)
	# --- Init dictionnary 
	waveformDict = Dict{String,Waveform}();
	nbElem		 = nfields(x);
	#@show nbElem
	# --- Iterative extraction of waveform configuration
	for i = 1 : 1 : nbElem
		cnt = 1;
		# --- Getting name of configuration
		charW	= getWaveformName(x[i]);
		# --- Iterative earch of other config of same waveform 
		tC	= charW;
		al  = collect(keys(waveformDict));
		for iN = 1 : 1 : length(waveformDict) 
			# --- 
			if al[iN] == tC 
				tC = "$charW-$cnt";
				cnt+=1;
			end
		end
		# --- Pushing the configuration in the dictionnary
		waveformDict[tC] = x[i];
	end
	return waveformDict;
end

""" 
---  
Create a signal based on a waveform dictionnary and a desired configuration (i.e the key)
# --- Syntax 
sigId	= genSig(qamMat,dWav,key)
# --- Input parameters 
- qamMat  : Time frequency complex QAM matrix [Array{Complex{Float64}}]
- dWav	  : Waveform dictionnary (see initWaveforms)
- key	  : Desired waveform configuration
# --- Output parameters 
- sigId	  : Time domain signal 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genSig(qamMat::Array{Complex{Float64}},dict::Dict{String,Waveform},mod::String)
	return genSig(qamMat,dict[mod]);
end

""" 
---  
Create a signal based on a waveform dictionnary and a desired configuration (i.e the key)
# --- Syntax 
sigId	= decodeSig(qamMat,dWav,key)
# --- Input parameters 
- signal	  : Time domain signal 
- qamMat  : Time frequency complex QAM matrix [Array{Complex{Float64}}]
- dWav	  : Waveform dictionnary (see initWaveforms)
- key	  : Desired waveform configuration
# --- Output parameters 
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function decodeSig(signal,dict::Dict{String,Waveform},mod::String)
	return decodeSig(signal,dict[mod]);
end


