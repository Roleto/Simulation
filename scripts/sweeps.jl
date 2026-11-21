# Előfeltétel: includet("src/runner.jl") 
# Használat:   includet("scripts/sweeps.jl")
# using Pkg
# Pkg.status()                  # mi van a projektben
# Pkg.add(["PyPlot", "Revise"]) # ha hiányzik valamelyik
# Pkg.instantiate()             # a Manifest alapján mindent letölt
# using Revise
# includet("src/runner.jl")
# includet("scripts/sweeps.jl")
# sw_duffing_L_KVSSM()
# scripts/sweeps.jl
# Előfeltétel: includet("src/runner.jl")  (betölti a run_duffing-ot)
# Használat:   includet("scripts/sweeps.jl");  sw_duffing_L_KVSSM()
using PyPlot, DelimitedFiles, Printf, Dates
PyPlot.ioff()

const NAMES = ("rms_ss", "maxe_ss", "rms_all", "maxe_all", "maxu")
lg(a, b, n) = 10.0 .^ range(a, b, length=n)
const SAVED = String[]     # a sweep ide gyűjti a képek útvonalát; csak a logged írja ki, a konzolra (a logba nem)

# ---------- Motor: bármely run_fn(; kwargs..., metrics=true) -> (ok,rms,maxe,maxu) ----------
function sweep(run_fn, name, p1, v1, p2, v2; tail=0.2)
    out = joinpath(@__DIR__, "..", "data", "Sweeps", name)
    mkpath(out)
    tg = "tail$(tail)"
    R = fill(NaN, length(v2), length(v1), length(NAMES))
    rows = Any[]
    for (j, b) in enumerate(v2), (i, a) in enumerate(v1)
        m = run_fn(; Dict(p1 => a, p2 => b)..., metrics=true, tail=tail)
        vals = [m.rms_ss, m.maxe_ss, m.rms_all, m.maxe_all, m.maxu]
        m.ok && (R[j, i, :] = vals)                     # divergált -> NaN (szürke)
        push!(rows, (a, b, m.ok, vals...))
    end
    # heatmap + mátrix CSV
    cm = get_cmap("viridis")
    cm.set_bad("lightgray")
    for (k, nm) in enumerate(NAMES)
        figure()
        imshow(log10.(R[:, :, k]), origin="lower", aspect="auto", cmap=cm)
        colorbar(label="log10($nm)")
        xticks(0:length(v1)-1, [@sprintf("%.2g", x) for x in v1], rotation=45)
        yticks(0:length(v2)-1, [@sprintf("%.2g", x) for x in v2])
        xlabel(string(p1))
        ylabel(string(p2))
        title("$name – $nm (tail=$tail)")
        tight_layout()
        png = joinpath(out, "heat_$(p1)_$(p2)_$(nm)_$tg.png")
        savefig(png)
        push!(SAVED, normpath(png))
        close()
        writedlm(joinpath(out, "matrix_$(p1)_$(p2)_$(nm)_$tg.csv"),
            vcat(permutedims(vcat("", v1)), hcat(v2, R[:, :, k])), ',')
    end
    # összefoglaló CSV (minden futás egy sor)
    writedlm(joinpath(out, "summary_$(p1)_$(p2)_$tg.csv"),
        vcat(permutedims(Any[string(p1), string(p2), "ok", NAMES...]),
            permutedims(hcat(map(collect, rows)...))), ',')
    println("CSV → ", normpath(out))
    nok = count(r -> r[3], rows)
    @printf("%s %s×%s: %d/%d stabil\n", name, p1, p2, nok, length(rows))
    return R
end

# ---------- Duffing: egysoros sweepek ----------
sw_duffing_KVSSM_w(; kw...) = sweep(run_duffing, "Duffing", :K_VSSM, lg(0, 4, 9), :w, lg(-2, 1, 7); kw...)
sw_duffing_L_KVSSM(; kw...) = sweep(run_duffing, "Duffing", :Λ, lg(0, 2, 9), :K_VSSM, lg(0, 4, 9); kw...)  # Λ: 1–100

function sw_duffing_all(; kw...)
    sw_duffing_KVSSM_w(; kw...)
    sw_duffing_L_KVSSM(; kw...)
end

# A·K (és A×K) sweep minden modellre: lásd sweep_AK_1d / sweep_params a runner.jl-ben.
# Másik rendszer Λ×K_VSSM sweepje egy sor, pl.:
# sw_lorenz_L_KVSSM(; kw...) = sweep(run_lorenz, "Lorenz", :Λ, lg(0,2,9), :K_VSSM, lg(0,3,9); kw...)

# ---------- Logolás: a kimenet + a visszatérési érték fájlba megy, a konzolra csak egy sor ----------
# Használat:  logged(() -> sw_duffing_L_KVSSM())
# logged(() -> sw_duffing_L_KVSSM(); name="duffing_L_KVSSM")
function logged(f; name="sweeps")
    path = joinpath(@__DIR__, "..", "data", "Sweeps", "$name.log")
    mkpath(dirname(path))
    res = nothing
    empty!(SAVED)
    open(path, "w") do io
        println(io, "=== ", Dates.now(), " ===")
        redirect_stdio(() -> (res = f()); stdout=io, stderr=io)
        if res !== nothing                       # teljes tömb, csonkolás nélkül
            show(IOContext(io, :limit => false), MIME("text/plain"), res)
            println(io)
        end
    end
    foreach(f -> println("kép → ", f), SAVED)
    println("kész → ", normpath(path))
    return nothing                               # a REPL ne írja ki újra
end