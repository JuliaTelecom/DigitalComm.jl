

function cdmaSigDecode(sig::AbstractVector{T},nbUsers,code::Symbol,userMask=nothing) where T
    # --- Define constant 
    nbChips = length(sig)
    nbSymb = nbChips  ÷ nbUsers
    # --- Define code 
    if code == :ovsf
        c = ovsf(nbUsers)
    else 
        @error "Code $code is not supported"
    end
    if isnothing(userMask)
        userMask = 1:nbUsers
    end
    nbActiveUsers = length(userMask)
    qamRx = zeros(T,nbActiveUsers,nbSymb)
    tmp = zeros(T,nbSymb) # To be optimized as a view of qamRx
    for n ∈ eachindex(userMask)
        despreading!(tmp,sig,nbUsers,c[:,userMask[n]])
        qamRx[n,:] .= tmp 
    end
    return qamRx 
end

function despreading!(tmp,sig,nbUsers,code)
    for n ∈ eachindex(tmp)
        tmp[n] = 0
        for k ∈ 1 : nbUsers
            tmp[n] += sig[(n-1)*nbUsers + k] * code[k]
        end
        tmp[n] = tmp[n] / nbUsers
    end
end


# Dispatch with StrucCDMA 
cdmaSigDecode(sig::AbstractVector,cdma::StrucCDMA) = cdmaSigDecode(sig,cdma.nbUsers,cdma.code,cdma.userMask)
