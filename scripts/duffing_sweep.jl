# Előfeltétel: includet("src/runner.jl")  (betölti a run_duffing-ot)
@isdefined(run_duffing) || error("Előbb töltsd be: includet(\"src/runner.jl\")")
using PyPlot, DelimitedFiles, Printf
PyPlot.ioff()

out = joinpath(@__DIR__, "..", "data", "Sweeps", "Duffing");
mkpath(out);
lg(a, b, n) = 10.0 .^ range(a, b, length=n)
const NAMES = ("rms", "maxe", "maxu")

# ---------- Rács ----------
rows = Any[]
function sweep(p1, v1, p2, v2)
    R = fill(NaN, length(v2), length(v1), 3)             # rms, maxe, maxu
    for (j, b) in enumerate(v2), (i, a) in enumerate(v1)
        m = run_duffing(; Dict(p1 => a, p2 => b)..., metrics=true)
        m.ok && (R[j, i, :] = [m.rms, m.maxe, m.maxu])
        push!(rows, (p1, a, p2, b, m.ok, m.rms, m.maxe, m.maxu))
    end
    R
end

# ---------- Heatmap és mátrix-CSV ----------
function heat(R, p1, v1, p2, v2)
    cm = get_cmap("viridis")
    cm.set_bad("lightgray")      # szürke = divergált
    for (k, nm) in enumerate(NAMES)
        figure()
        imshow(log10.(R[:, :, k]), origin="lower", aspect="auto", cmap=cm)
        colorbar(label="log10($nm)")
        xticks(0:length(v1)-1, [@sprintf("%.2g", x) for x in v1], rotation=45)
        yticks(0:length(v2)-1, [@sprintf("%.2g", x) for x in v2])
        xlabel(string(p1))
        ylabel(string(p2))
        title("Duffing – $nm")
        tight_layout()
        savefig(joinpath(out, "heat_$(p1)_$(p2)_$nm.png"))
        close()
    end
end
mat_csv(io, M, v1, v2) =
    writedlm(io, vcat(permutedims(vcat("", v1)), hcat(v2, M)), ',')

# ---------- Futtatás ----------
b = run_duffing(metrics=true)
@printf("Alapeset: RMS=%.3g  max|e|=%.3g  max|u|=%.3g\n", b.rms, b.maxe, b.maxu)

pairs = [(:K_VSSM, lg(0, 4, 9), :w, lg(-2, 1, 7)),
    (:A, lg(-7, -3, 9), :K, lg(3, 7, 9)),
    (:Λ, lg(-1, 1, 7), :K_VSSM, lg(0, 4, 9))]
mats = Any[]
for (p1, v1, p2, v2) in pairs
    R = sweep(p1, v1, p2, v2)
    heat(R, p1, v1, p2, v2)
    for (k, nm) in enumerate(NAMES)
        mat_csv(joinpath(out, "matrix_$(p1)_$(p2)_$nm.csv"), R[:, :, k], v1, v2)
    end
    push!(mats, (p1, v1, p2, v2, R))
end

# ---------- results.csv és summary.txt (minden egy helyen) ----------
hdr = ["p1" "v1" "p2" "v2" "ok" "rms" "maxe" "maxu"]
writedlm(joinpath(out, "results.csv"),
    vcat(hdr, reduce(vcat, [permutedims(collect(r)) for r in rows])), ',')
good = sort(filter(r -> r[5], rows), by=r -> r[6])
open(joinpath(out, "summary.txt"), "w") do f
    println(f, "Futások: $(length(rows)), divergált: $(count(r -> !r[5], rows))")
    println(f, "Alapeset: RMS=$(b.rms) max|e|=$(b.maxe) max|u|=$(b.maxu)")
    println(f, "Legjobb 10 (RMS szerint):")
    for r in good[1:min(10, end)]
        @printf(f, "%s=%.3g, %s=%.3g | RMS=%.3g max|e|=%.3g max|u|=%.3g\n",
            r[1], r[2], r[3], r[4], r[6], r[7], r[8])
    end
    for (p1, v1, p2, v2, R) in mats, (k, nm) in enumerate(NAMES)
        println(f, "\n# $nm  (oszlop: $p1, sor: $p2; NaN = divergált)")
        mat_csv(f, R[:, :, k], v1, v2)
    end
end
PyPlot.ion()
println("Kész: $out")