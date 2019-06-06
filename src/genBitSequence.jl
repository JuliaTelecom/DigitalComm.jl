# --- Binary managment  

# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
using Random


# ---------------------------------------------------- 
# --- Source code 
# ---------------------------------------------------- 

""" genbBitsequence!
---  
Create a binary sequence and populate input buffer with nbBits bits
The array is of type UInt8 with x00 or x01)
If stated, randSeed controls the seed of the random generator
# --- Syntax 
	  genbBitsequence!(buffer,nbBits,randSeed=-1);
# --- Input parameters 
- buffer	: Buffer to populate [Array{UInt8,nbBits}]
- nbBits	: Number of bits to generate [Int]
- randSeed	: Seed of random process (default -> -1)
# --- Output parameters 
- buffer	: Populated buffer [Array{UInt8}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genBitSequence!(buffer::Array{UInt8},nbBits,randSeed=-1)
	# --- Setting the random seed if given
	if randSeed != -1
		# --- Seed as the second parameter
		Random.seed!(randSeed);
	end
	# --- Generating binary sequence
	buffer .= Int8.(rand([0,1],Int(nbBits)));
	# --- Switch to "pure" random system
	RandomDevice();
	# --- Export buffer 
	return buffer;
end

""" genbBitsequence
---  
Create a binary sequence and return a buffer with nbBits bits
The array is of type UInt8 with x00 or x01)
If stated, randSeed controls the seed of the random generator
# --- Syntax 
	  genbBitsequence!(nbBits,randSeed=-1);
# --- Input parameters 
- nbBits	: Number of bits to generate [Int]
- randSeed	: Seed of random process (default -> -1) [Int]
# --- Output parameters 
- buffer	: Populated buffer [Array{UInt8}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genBitSequence(nbBits,randSeed=-1)
	# --- Create the input buffer 
	buffer	= zeros(UInt8,nbBits);
	# --- Call the bang method 
	genBitSequence!(buffer,nbBits,randSeed);
	# --- Return the buffer 
	return buffer;
end




#function genByteSequence(nbBytes::Any,randSeed=-1)
	## --- Setting the random seed if given
	#if randSeed != -1
		## --- Seed as the second parameter
		#Random.seed!(randSeed);
	#end
	## --- Generating binary sequence
	#return bitSeq  = convert(Array{UInt8}, rand(collect(0x00:0x01:0xff),Int(nbBytes)));
	## --- Switch to "pure" random system
	#RandomDevice();
#end
