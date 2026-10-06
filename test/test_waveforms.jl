# ---------------------------------------------------- 
# --- Import modules  
# ---------------------------------------------------- 
using DigitalComm 
using Test
# ---------------------------------------------------- 
# --- Tests 
# ---------------------------------------------------- 
println("Tests for multicarrier waveforms");

# ----------------------------------------------------
# --- Overall parameters
# ----------------------------------------------------
# --- Overall PHY parameters
nbSymb 			= 14;			  # --- Number of symbols (one frame)
nFFT 			= 1024;			  # --- Base FFT size
samplingFreq	= 15.36;		  # --- Frequency value (MHz)
snrVect			= (-10:30);
# --- Frequency allocation
allocatedSubcarriers= getLTEAlloc(nFFT);
nbSubcarriers	      = length(allocatedSubcarriers);

# ----------------------------------------------------
# --- Waveform contender
# ----------------------------------------------------
# --- Init OFDM structure
ofdm  = initOFDM(
                 nFFT,						        # --- nFFT                 : FFT size
                 72,						        # --- nCP                  : CP size
                 allocatedSubcarriers		        # --- allocatedSubcarriers : Subcarrier allocation
                );
# --- Init SCFDMA structure
scfdma  = initSCFDMA(
                     nFFT,						        # --- nFFT                 : FFT size
                     72,						        # --- nCP                  : CP size
                     allocatedSubcarriers,		        # --- allocatedSubcarriers : Subcarrier allocation
                     12;								# --- sizeDFT			   : DFT preprocessing size
                    );
# --- Init UF-OFDM structure
ufofdm  = initUFOFDM(
                     nFFT,					        # --- nFFT                 : FFT size
                     73,					        # --- L                    : Filter length (same size +1 due to conv)
                     allocatedSubcarriers,	        # --- allocatedSubcarriers : Subcarrier allocation
                     applyPD=1,				        # --- applyPD              : Do predistortion at Tx stage
                     attenuation=40,		        # --- attenuation          : Filter attenuation in dB
                    );
# --- Init BF-OFDM structure
bfofdm	= initBFOFDM(
                     32,				            # --- nFBMC                : PPN size (max number of carriers)
                     64,				            # --- nOFDM                : Precoder size (OFDM sizer)
                     3,				            	# --- K                    : Overlapping factor
                     9,					            # --- GI                   : CP size of precoder
                     0.5,				            # --- δ                    : compression factor
                     allocatedSubcarriers,          # --- allocatedSubcarriers : Subcarrier allocation
                     "gaussian",	            # --- filterName           : Pulse shape name
                     BT=0.36,				        # --- BT                   : Potential BT value for Gaussian
                     filterStopBand = 110,			# --- filterStopBand       : DC stopband value
                     fS=[],				            # --- fS                   : Potential frequency coefficient for FS filter
                     nFFT= 1024,		            # --- nFFT                 : associated FFT value in Rx
                     nCP= 72,			            # --- nCP                  : extended CP size
                    );
# --- Init WOLA-OFDM structure
wola  = initWOLA(
                 nFFT,						        # --- nFFT                 : FFT size
                 72,						        # --- nCP                  : CP size
                 allocatedSubcarriers,		        # --- allocatedSubcarriers : Subcarrier allocation
                 "triangle",						# --- Window type @Tx side
                 20,								# --- Window size @Tx side
                 "triangle",						# --- Window type @Rx side
                 20,								# --- Window size @Rx side
                );
fbmc  = initFBMC(
                 nFFT,						        # --- nFFT                 : FFT size
                 4,									# --- K					   : Overlapping factor
                 allocatedSubcarriers		        # --- allocatedSubcarriers : Subcarrier allocation
                );

# ----------------------------------------------------
# --- Merging structures
# ----------------------------------------------------
# Create  a dictionnary to rule them all 
waveforms 	= initWaveforms(ofdm,
                            scfdma,
                            ufofdm,
                            bfofdm,
                            wola,
                            fbmc,
                           );


for  (name,struc) in waveforms
    @testset  "$name" begin 
        for mcs=[4,16,64,256]
            nbBits			      = nbSymb * nbSubcarriers * Int(log2(mcs));
            # --- Binary sequence
            bitSeq	      = genBitSequence(nbBits);
            # Mapping
            qamSeq		  = bitMappingQAM(mcs,bitSeq);
            # --- T/F matrix
            qamMat		  = reshape(qamSeq,nbSubcarriers,nbSymb);
            # --- Signal
            sigId		  = genSig(qamMat,struc);
            # --- Waveform demodulator 
            qamDec	  = decodeSig(sigId,struc);
            # --- Binary demapper
            bitDec	  = bitDemappingQAM(mcs,qamDec[:]);
            # --- BER measure
            nbE	 = sum(xor.(bitDec,bitSeq));
            @test nbE == 0
        end
    end
end



@testset "OVSF CDMA " begin 
    # Testing OVSF code 
    for N = [4,8,16,64]
        c = ovsf(N)
        for cn ∈ 1 : N
            for ck ∈ 1 : N 
                # --- Code correlation 
                γ = sum(c[cn,:] .* conj(c[ck,:]) )
                if cn == ck 
                    # Autocorrelation 
                    @test γ == N 
                else 
                    # Zero correlation 
                    @test γ == 0
                end
            end
        end
    end
end


@testset "CDMA generator" begin 
    N = 16
    C = ovsf(N)
    # --- Spreading
    for n ∈ 1 : N 
        c = C[:,n]
        # Check spreading with one as input, we should obtain the code 
        in = [1]
        out = zeros(16)
        DigitalComm._spread_accum!(out,in,c,N)
        @test all(out .== c)
        # Check that accumulation is Ok 
        DigitalComm._spread_accum!(out,in,c,N)
        @test all(out .== 2 .*c)
    end
    # --- True generator 
    nS = 100
    qamMat = ones(N,nS)
    sigId = cdmaSigGen(qamMat,N,:ovsf)
    @test sigId isa Vector 
    @test length(sigId) == nS * N
    for n = 1 : 1 : N
        @test sigId[n] == sum(C[n,:])
    end
    # --- Generate only user 8
    qamMat = ones(1,nS) # Equivalent to a one subcarrier OFDM system
    sigId = cdmaSigGen(qamMat,N,:ovsf,[8])
    @test sigId isa Vector 
    @test length(sigId) == nS * N
    for n = 1 : 1 : N
        @test sigId[n] == C[n,8]
    end
    # Testing with struct 
    cdma = initCDMA(N,:ovsf,[8])
    @test cdma isa Waveform 
    @test cdma isa DigitalComm.StrucCDMA
    sigId2 = cdmaSigGen(qamMat,cdma)
    @test all(sigId2 .≈ sigId)
end


@testset "CDMA receiver" begin 
    nS = 32768
    N  = 16 
    mcs = 4
    nbSymb = nS ÷ N ÷ Int(log2(mcs))
    # Mapping
    bitSeq	      = genBitSequence(nS);
    qamSeq		  = bitMappingQAM(mcs,bitSeq);
    # --- T/F matrix
    qamMat		  = reshape(qamSeq,N,nbSymb);
    sigId = cdmaSigGen(qamMat,N,:ovsf)
    # Testing global decoding method
    qamRx = cdmaSigDecode(sigId,N,:ovsf,1:N)
    @test all(qamRx .≈ qamMat)
    # Testing one user decoding 
    for n = 1 : 1 : N 
        qamRx = cdmaSigDecode(sigId,N,:ovsf,[n])
        @test all(qamRx[1,:] .≈ qamMat[n,:])
    end
    # Testing with struct 
    cdma = initCDMA(N,:ovsf,[8])
    @test cdma isa Waveform 
    @test cdma isa DigitalComm.StrucCDMA
    qamRx0 = cdmaSigDecode(sigId,N,:ovsf,[8])
    qamRx1 = cdmaSigDecode(sigId,cdma)
    @test all(qamRx0 .≈ qamRx1)
end


@testset "Gold code generator" begin
    for n = [5,6,7]
        N = 2^n - 1
        C = gold(n)
        @test size(C) == (N,N)
        @test all(abs.(C) .== 1) # bipolar codes, no zero entries
        # --- Autocorrelation: peak at zero shift should equal N (perfect match)
        for cn ∈ 1 : N
            @test sum(C[:,cn] .* C[:,cn]) == N
        end
        # --- Correlation under cyclic shift stays inside the Gold bound
        # {-1, -t(n), t(n)-2} ∪ {N} -- unlike OVSF, Gold codes are NOT
        # perfectly orthogonal, so we check boundedness rather than exact
        # zero. Exhaustive for n=5 (cheap); randomly sampled for n=6,7 to
        # keep runtime reasonable (full search is O(N^4)).
        t     = 2^(div(n+2,2)) + 1
        bound = Set([-1, -t, t-2, N])
        pairs = n == 5 ? [(cn,ck,s) for cn in 1:N, ck in 1:N, s in 0:N-1] :
                          [(rand(1:N),rand(1:N),rand(0:N-1)) for _ in 1:2000]
        for (cn,ck,s) in pairs
            γ = sum(C[:,cn] .* circshift(C[:,ck],s))
            @test γ ∈ bound
        end
    end
    # --- Unsupported register length without explicit taps should error
    @test_throws ErrorException gold(4)
end

@testset "CDMA generator with Gold codes" begin
    n  = 5
    N  = 2^n - 1   # 31
    C  = gold(n)
    nS = 50
    # --- Single active user (DSSS: one user, one spreading code)
    qamMat = ones(1,nS)
    sigId  = cdmaSigGen(qamMat,N,:gold,[3])
    @test sigId isa Vector
    @test length(sigId) == nS * N
    for k = 1 : 1 : N
        @test sigId[k] == C[k,3]
    end
    # --- Full multi-user load, noiseless synchronous channel: every user's
    # bits should be recovered exactly, same as with :ovsf
    qamMat = Float64.(rand((-1,1),N,nS))
    sigId  = cdmaSigGen(qamMat,N,:gold)
    qamRx  = cdmaSigDecode(sigId,N,:gold,1:N)
    @test all(sign.(qamRx) .≈ qamMat)
    # --- Dispatch through the StrucCDMA / genSig / decodeSig path
    cdma  = initCDMA(N,:gold)
    sigId2 = genSig(qamMat,cdma)
    @test all(sigId2 .== sigId)
    qamDec = decodeSig(sigId2,cdma)
    @test all(sign.(qamDec) .≈ qamMat)
end

@testset "CDMA through genSig/decodeSig" begin
    nbSymb = 14
    nbUsers = 16
    cdma = initCDMA(nbUsers,:ovsf)
    for mcs=[4,16,64,256]
        nbBits			      = nbSymb * nbUsers * Int(log2(mcs));
        # --- Binary sequence
        bitSeq	      = genBitSequence(nbBits);
        # Mapping
        qamSeq		  = bitMappingQAM(mcs,bitSeq);
        # --- T/F matrix
        qamMat		  = reshape(qamSeq,nbUsers,nbSymb);
        # --- Signal
        sigId		  = genSig(qamMat,cdma);
        # --- Waveform demodulator
        qamDec	  = decodeSig(sigId,cdma);
        # --- Binary demapper
        bitDec	  = bitDemappingQAM(mcs,qamDec[:]);
        # --- BER measure
        nbE	 = sum(xor.(bitDec,bitSeq));
        @test nbE == 0
    end
end
