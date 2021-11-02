# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using DSP 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
println("Tests for filters");
@testset "Raised Cosine filter test" begin 
    # Create the vector 
    h = raisedCosine(12,0.5,16)
    # Check type 
    @test h isa Vector{Float64}
    # Check size 
    sH = 12*16*2+1
    @test length(h)== sH
    # Check it is 1 in the middle 
    @test h[1+ (sH-1)÷2] == 1
    # Check Nyquist criterion
    for n ∈ 0 : 11 
        @test ≈(h[1 + n*16],0,atol=1e-8)
    end
end

@testset "Raised Cosine filter test" begin 
    # Create the vector 
    h = sqrtRaisedCosine(12,0.5,16)
    # Check type 
    @test h isa Vector{Float64}
    # Check size 
    sH = 12*16*2+1
    @test length(h)== sH
    # Check it is 1 in the middle 
    @test h[1+ (sH-1)÷2] == 1
    # Check Nyquist criterion
    # sqrt is not a Nyquist filter but h*h is !
    p = conv(h,h)
    p = p /maximum(p)
    for n ∈ 0 : 2*11 
        @test ≈(p[1 + n*16],0,atol=1e-4)
    end
end
