""" Orthogonal Variable Spreading Factor (OVSF) code generation
It generates N codes that are orthogonal and follows the OVSF code pattern 
codes = ovsf(N,initState::UInt8=1) \\
Input 
- N = Number of generated codes (should be a power of 2)
- initState : Initial state value (default 1)
Output
- codes : Matrix of size code_size x N with type UInt8
[1] Yuh-Shyan Chen and Ting-Lung Lin, "Code placement and replacement schemes for WCDMA Rotated-OVSF code tree management," in IEEE Transactions on Mobile Computing, vol. 5, no. 3, pp. 224-239, March 2006, doi: 10.1109/TMC.2006.30.
"""
function ovsf(N,initState::Int8=Int8(1))
    n = Int(log2(N))
    C = [initState]
    for _ in 0 : n - 1
        newC = Int8[]
        for c in eachrow(C)
            c_0 = vcat(c,c)
            c_1 = vcat(c,-c)
            newC = _hcat_empty(newC,c_0)
            newC = _hcat_empty(newC,c_1)
        end 
        C = newC
    end
    return C
end

# A custom concat operator that handle empty input vector
function _hcat_empty(a,b)
    isempty(a) ? b : hcat(a,b)
end


""" Generates a Code Division Multiplexing Access signal from a Symbol matrix `qamMat` and using code family `code`. \\
cdmaSigGen(qamMat::AbstractMatrix{T},nbUsers::Number,code::Symbol,userMask=nothing) \\
Input parameters 
- qamMat : Matrix of symbols of size nbActiveUser x nbSymb. Each line corresponds to a user stream and there is at least nbActiveUser colunm. nbActiveUser should be ≤ than the number of code we have 
- nbUsers : Maximal number of user that can be multiplexed. Used for code generation 
- code : Type of code used as a Symbol. We support `:ovsf`
- userMask: Optional parameter, a vector of the code used. Should be of size nbActiveUser. For example, if we have 16 max users but we only want to encode 1000 symbols for  user 2 and user 12, we should call `cdmaSigGen` with a matrix `qamMat` of size (1000,2), `nbUsers` = 16 and `userMask` = [2;12]
"""
function cdmaSigGen!(sigOut::AbstractVector,qamMat::AbstractMatrix{T},nbUsers::Number,code::Symbol,userMask=nothing) where T
    # ----------------------------------------------------
    # --- Code generation 
    # ---------------------------------------------------- 
    if code == :ovsf
        # --- Generate the code base 
        c = ovsf(nbUsers)
    else 
        @error "Code $code is unsupported. We currently inly support OVSF"
    end
    nbActiveUser = size(qamMat,1)
    nbSymb       = size(qamMat,2)
    # Check we can encode users (more codes than uers)
    @assert nbActiveUser ≤ nbUsers "Number of active users ($nbActiveUser) should be lower than the code size ($nbUsers)."
    # Check output is of enough size 
    @assert length(sigOut) == nbUsers * nbSymb "Size of output vector ($(length(sigOut)) does not match input size x code size ($(nbSymb)x$(nbUsers)= $(nbUsers*nbSymb))"
    # Code selection
    if isnothing(userMask) 
        # We use the first codes as codeset
        userMask = 1 : nbActiveUser 
    end
    # ----------------------------------------------------
    # --- Scramling
    # ---------------------------------------------------- 
    for n = 1 : 1 : nbActiveUser
        # Apply the spreading and accumulate
        _spread_accum!(sigOut,qamMat[n,:],c[:,userMask[n]],nbUsers)
    end
    return sigOut
end
function cdmaSigGen(qamMat::AbstractMatrix{T},nbUsers::Number,code::Symbol,userMask=nothing) where T
    sigOut  = zeros(T,size(qamMat,2)*nbUsers)
    cdmaSigGen!(sigOut,qamMat,nbUsers,code,userMask)
return sigOut
end
    
""" Apply a spreading of input signal `seq` using `code` and accumulate the result in `tmp`
"""
function _spread_accum!(tmp,seq,code,sF)
    @inbounds @simd for n ∈ eachindex(seq)
        for k ∈ 1 : sF 
            tmp[(n-1)*sF + k] += seq[n] * code[k]
        end
    end
end

