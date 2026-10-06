using DigitalComm
using Plots
using DSP: welch_pgram, power, freq
gr()

outdir = joinpath(@__DIR__, "cdma_plots")
mkpath(outdir)

# ----------------------------------------------------
# --- 1) What the codes themselves look like
# ----------------------------------------------------
No = 16            # OVSF code size
Co = ovsf(No)       # No x No, columns are codes

ng = 5
Ng = 2^ng - 1       # 31, Gold code size
Cg = gold(ng)        # Ng x Ng, columns are codes

p1 = heatmap(Co, c=:grays, colorbar=false, title="OVSF codes (N=$No)",
             xlabel="code #", ylabel="chip #", yflip=true)
p2 = heatmap(Cg, c=:grays, colorbar=false, title="Gold codes (N=$Ng)",
             xlabel="code #", ylabel="chip #", yflip=true)
plot(p1, p2, layout=(1,2), size=(900,400))
savefig(joinpath(outdir, "1_code_matrices.png"))

# A single chip sequence, zoomed in, as a stairs plot (what one code "looks like")
p3 = plot(Cg[:,1], seriestype=:steppost, lw=2, legend=false,
          title="One Gold chip sequence (code #1, N=$Ng)",
          xlabel="chip index", ylabel="amplitude")
savefig(joinpath(outdir, "2_single_chip_sequence.png"))

# ----------------------------------------------------
# --- 2) DSSS: spreading a single user's bits over time
# ----------------------------------------------------
nbSymb = 8
bits   = genBitSequence(nbSymb)
qam    = bitMappingQAM(2, bits)             # BPSK: +-1 symbols
qamMat = reshape(qam, 1, nbSymb)

sigDSSS = cdmaSigGen(qamMat, Ng, :gold, [1])  # spread with Gold code #1

p4 = plot(real.(sigDSSS), lw=1, legend=false,
          title="DSSS chip-level signal (single user, $nbSymb symbols)",
          xlabel="chip index (time)", ylabel="amplitude")
# mark symbol boundaries
vline!(p4, collect(0:nbSymb) .* Ng, ls=:dash, lc=:red, alpha=0.4)
savefig(joinpath(outdir, "3_dsss_timedomain.png"))

# ----------------------------------------------------
# --- 3) CDMA: several users summed on the same channel
# ----------------------------------------------------
nbActiveUsers = 4
qamMU  = reshape(bitMappingQAM(2, genBitSequence(nbSymb*nbActiveUsers)), nbActiveUsers, nbSymb)
sigCDMA = cdmaSigGen(qamMU, Ng, :gold, 1:nbActiveUsers)

p5 = plot(real.(sigDSSS), lw=1, label="1 user", alpha=0.8)
plot!(p5, real.(sigCDMA), lw=1, label="$nbActiveUsers users summed",
      title="Single user vs multi-user CDMA (chip domain)",
      xlabel="chip index", ylabel="amplitude")
savefig(joinpath(outdir, "4_single_vs_multiuser.png"))

# ----------------------------------------------------
# --- 4) Spectrum: spreading flattens/widens the occupied band
# ----------------------------------------------------
nbSymbLong = 2000
qamLong    = bitMappingQAM(2, genBitSequence(nbSymbLong))
sigUnspread = real.(qamLong)                                       # un-spread BPSK, 1 sample/symbol
sigSpread   = real.(cdmaSigGen(reshape(qamLong,1,nbSymbLong), Ng, :gold, [1]))  # spread, Ng samples/symbol

pxxU = welch_pgram(sigUnspread; fs=1.0)
pxxS = welch_pgram(sigSpread;   fs=Float64(Ng))  # chip rate = Ng x symbol rate

p6 = plot(freq(pxxU), 10 .* log10.(power(pxxU)), label="Unspread BPSK",
          title="PSD: unspread vs DSSS-spread", xlabel="Normalized frequency", ylabel="PSD [dB]")
plot!(p6, freq(pxxS), 10 .* log10.(power(pxxS)), label="DSSS spread (gain N=$Ng)")
savefig(joinpath(outdir, "5_psd_comparison.png"))

# ----------------------------------------------------
# --- 5) Despreading: correlate the multi-user signal against every
# code, in a noisy channel - active users should stand out
# ----------------------------------------------------
snr = 5  # dB, deliberately low to make the contrast visible
sigNoisy, _ = addNoise(sigCDMA, snr)
corrPerCode = [sum(abs2, cdmaSigDecode(sigNoisy, Ng, :gold, [c])) for c in 1:Ng]

p7 = bar(1:Ng, corrPerCode, legend=false,
         title="Despread energy per code (SNR=$(snr) dB, users 1:$nbActiveUsers active)",
         xlabel="code index", ylabel="despread energy")
vline!(p7, [nbActiveUsers+0.5], ls=:dash, lc=:red, label="active | inactive")
savefig(joinpath(outdir, "6_despread_energy_per_code.png"))

println("Saved plots to: $outdir")
for f in readdir(outdir)
    println("  - ", f)
end
