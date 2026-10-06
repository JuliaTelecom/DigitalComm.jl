using DigitalComm

# ----------------------------------------------------
# --- DSSS vs CDMA : same spreading machinery, different use
# ----------------------------------------------------
# DSSS (Direct Sequence Spread Spectrum) spreads ONE user's data with ONE
# chip sequence - the goal is processing gain / LPI / anti-jam, not sharing
# the channel with anyone else.
# CDMA (Code Division Multiple Access) is DSSS used as a MULTIPLE ACCESS
# scheme: every user gets its own code from the same code family, and all
# users share the same chip-rate channel simultaneously. Decoding a user is
# just DSSS despreading with that user's code; the other users' spread
# signals act as (controlled) interference on top of the noise.
#
# OVSF codes (already in this PR) are only orthogonal when every user's
# transmission is chip-aligned (synchronous). Gold codes trade perfect
# orthogonality for a correlation that stays *bounded* under arbitrary
# relative shifts - the classical choice for asynchronous / multipath CDMA,
# and also for ranging (good autocorrelation).

n = 5                   # LFSR register length -> code length N = 2^n - 1 = 31
N = 2^n - 1

# ----------------------------------------------------
# --- 1) DSSS: a single user, single spreading code
# ----------------------------------------------------
# This is the direct analogue of the CommunicationsSequences.jl LFSR/DSSS
# snippet, but using DigitalComm's own Gold code generator instead of a bare
# m-sequence - any one column of gold(n) is a valid chip sequence on its own.
code = gold(n)[:, 1]         # pick one code out of the Gold family
nbSymb = 10
txData = rand((0,1), nbSymb)
qamMat = reshape(2 .* txData .- 1, 1, nbSymb)   # 0/1 -> -1/+1, 1 user x nbSymb

sigId = cdmaSigGen(qamMat, N, :gold, [1])        # spread with code #1 only
qamRx = cdmaSigDecode(sigId, N, :gold, [1])      # despread with the same code
rxData = Int.(sign.(qamRx[1,:]) .> 0)

println("DSSS (1 user, Gold code #1): all bits recovered = ", all(rxData .== txData))

# ----------------------------------------------------
# --- 2) CDMA: every code slot used by a different simultaneous user
# ----------------------------------------------------
nbUsers = N                                      # full load: N users, N codes
qamMatMU = Int8.(rand((-1,1), nbUsers, nbSymb))  # one data stream per user

sigIdMU = cdmaSigGen(qamMatMU, N, :gold)          # all users summed in the channel
qamRxMU = cdmaSigDecode(sigIdMU, N, :gold, 1:N)   # despread user-by-user

println("CDMA (N=$nbUsers users, Gold codes, synchronous, noiseless): ",
        "all users recovered exactly = ", all(sign.(qamRxMU) .≈ qamMatMU))

# Residual multi-user interference, if you want to see it (will be exactly
# 0 here because synchronous + noiseless + Gold's -1 zero-shift
# cross-correlation keeps |interference| < 1, see test_waveforms.jl):
err = qamRxMU .- qamMatMU
println("max |decode error| per user = ", maximum(abs.(err)))
