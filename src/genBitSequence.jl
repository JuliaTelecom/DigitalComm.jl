# --- Binary managment  

# ---------------------------------------------------- 
# --- Modules  
# ---------------------------------------------------- 
using Random


# ---------------------------------------------------- 
# --- Source code 
# ---------------------------------------------------- 

"""
---  
Create a binary sequence and populate input buffer with nbBits bits
The array is of type UInt8 with x00 or x01)
If stated, randSeed controls the seed of the random generator
# --- Syntax 
	  genBitsequence!(buffer,nbBits,randSeed=-1);
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
	#buffer .= Int8.(rand([0,1],Int(nbBits)));
	Random.rand!(buffer,[0,1]);
	# --- Switch to "pure" random system
	if randSeed != -1
		RandomDevice();
	end
end

"""
---  
Create a binary sequence and return a buffer with nbBits bits
The array is of type UInt8 with x00 or x01)
If stated, randSeed controls the seed of the random generator
# --- Syntax 
	  genBitsequence!(nbBits,randSeed=-1);
# --- Input parameters 
- nbBits	: Number of bits to generate [Int]
- randSeed	: Seed of random process (default -> -1) [Int]
# --- Output parameters 
- buffer	: Populated buffer [Array{UInt8}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genBitSequence(nbBits,randSeed=-1)
	# --- Setting the random seed if given
	if randSeed != -1
		# --- Seed as the second parameter
		Random.seed!(randSeed);
	end
	# --- Create the input buffer 
	buffer	= zeros(UInt8,nbBits);
	# --- Call the bang method 
	genBitSequence!(buffer,nbBits,randSeed);
	# --- Return the buffer 
	return buffer;
end



""" 
---  
Create a byte sequence and populate input buffer with nbytes bytes
The array is of type UInt8 with x00 to xff
If stated, randSeed controls the seed of the random generator
# --- Syntax 
	  genByteSequence!(buffer,nbBytes,randSeed=-1);
# --- Input parameters 
- buffer	: Buffer to populate [Array{UInt8,nbByte}]
- nbBytes	: Number of byte to generate [Int]
- randSeed	: Seed of random process (default -> -1)
# --- Output parameters 
- buffer	: Populated buffer [Array{UInt8}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genByteSequence!(buffer::Array{UInt8},nbBytes,randSeed=-1)
	# --- Setting the random seed if given
	if randSeed != -1
		# --- Seed as the second parameter
		Random.seed!(randSeed);
	end
	# --- Generating byte sequence
	Random.rand!(buffer, 0x00:0x01:0xff);
	# --- Switch to "pure" random system
	RandomDevice();
	if randSeed != -1
		RandomDevice();
	end
	# --- Export buffer 
	return buffer;
end

""" 
---  
Create a byte sequence and return a populated  buffer with nbytes bytes
The array 
# --- Syntax 
	  genByteSequence(nbBytes,randSeed=-1);
# --- Input parameters 
- nbBytes	: Number of byte to generate [Int]
- randSeed	: Seed of random process (default -> -1)
# --- Output parameters 
- buffer	: Populated buffer [Array{UInt8}]
# --- 
# v 1.0 - Robin Gerzaguet.
"""
function genByteSequence(nbBytes,randSeed=-1)
	# --- Create buffer 
	buffer = zeros(UInt8,nbBytes); 
	# --- Call ! 
	genByteSequence!(buffer,nbBytes,randSeed);
end


