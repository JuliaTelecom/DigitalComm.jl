

function cdmaSigDecode(sig::AbstractVector{T},nbUser,code::Symbol,userMask=nothing) where T
    # --- Define constant 
    nbChips = length(sig)
    nbSymb = nbChips  ÷ nbUser
    # --- Define code 
    if code == :ovsf
        c = ovsf(nbUser)
    else 
        @error "Code $code is not supported"
    end
    if isnothing(userMask)
        userMask = 1:nbUser
    end
    nbActiveUsers = length(userMask)
    qamRx = zeros(T,nbActiveUsers,nbSymb)
    tmp = zeros(T,nbSymb) # To be optimized as a view of qamRx
    for n ∈ eachindex(userMask)
        despreading!(tmp,sig,nbUser,c[:,userMask[n]])
        qamRx[n,:] .= tmp 
    end
    return qamRx 
end

function despreading!(tmp,sig,nbUser,code)
    for n ∈ eachindex(tmp)
        tmp[n] = 0
        for k ∈ 1 : nbUser
            tmp[n] += sig[(n-1)*nbUser + k] * code[k]
        end
        tmp[n] = tmp[n] / nbUser
    end
end

