# scripts/run_all_sweeps.jl
# Mind a 4 modell összes sweepje egy paranccsal; a kimenet a data/Sweeps/<modell>/ mappába kerül.
#
# Előfeltétel:  using Revise; includet("src/runner.jl")
# Használat:    includet("scripts/run_all_sweeps.jl")
#               run_all_sweeps()                 # minden modell
#               run_all_sweeps(models=["Lorenz"]) # csak egy modell
#               run_all_sweeps(cap=1.0)          # szigorúbb divergálás-küszöb
#
# A hiba a modellek skalár visszatérési értéke (max követési hiba). Divergált / NaN / Inf /
# cap fölötti futás -> NaN (szürke a képen). Nincs korai kilépés.

@isdefined(sweep_params) || error("Előbb töltsd be: includet(\"src/runner.jl\")")

lg(a, b, n) = 10.0 .^ range(a, b, length=n)

# (név, run függvény, K alapérték, alapértelmezett cap)
# Lorenz/Rössler: referencia-amplitúdó norm = sqrt(2²+3²+1²) ≈ 3.7 → cap=10 választja el
#   a jól követő (hiba << 1) és a rosszul követő (hiba > 10) tartományokat.
# Duffing/VdP: SISO, cap=1e3 megőrizve (ott a my_duffing ok-feltétel szűr).
const SWEEP_MODELS = [
    ("Duffing", run_duffing, 1e5, 1e3),
    ("VanDerPol", run_vanderpol, 1e5, 1e3),
    ("Lorenz", run_lorenz, 1e2, 10.0),
    ("Rossler", run_rossler, 1e2, 10.0),
]

# Rácsok (minden modellre ugyanazok)
const GRID_AK = lg(-3, 0.5, 25)      # A·K = 1e-3 … ~3,16
const GRID_L = lg(0, 2, 9)         # Λ = 1 … 100
const GRID_KV = lg(0, 4, 9)         # K_VSSM = 1 … 1e4
const GRID_W = lg(-2, 1, 7)         # w = 0,01 … 10

# cap=nothing → a SWEEP_MODELS-ből vett modell-specifikus érték használatos
function sweeps_one_model(name, fn, K, default_cap; cap=nothing)
    c = cap === nothing ? default_cap : cap
    println("=== $name  (cap=$c) ===")
    # 1) A·K (1D)
    P, E = sweep_AK_1d(fn, GRID_AK, K; cap=c)
    save_AK_1d(name, P, E)
    # 2) Λ × K_VSSM
    a, b, M = sweep_params(fn, :Λ, GRID_L, :K_VSSM, GRID_KV; cap=c, A=0.2 / K)
    save_sweep_2d(name, a, b, M, :Λ, :K_VSSM)
    # 3) K_VSSM × w
    a, b, M = sweep_params(fn, :K_VSSM, GRID_KV, :w, GRID_W; cap=c, A=0.2 / K)
    save_sweep_2d(name, a, b, M, :K_VSSM, :w)
    return nothing
end

# cap=nothing → minden modell a saját SWEEP_MODELS-beli cap-jét használja
# cap=<érték> → felülírja az összeset (pl. run_all_sweeps(cap=5.0) minden modellnél)
function run_all_sweeps(; models=nothing, cap=nothing)
    for (name, fn, K, default_cap) in SWEEP_MODELS
        (models === nothing || name in models) || continue
        sweeps_one_model(name, fn, K, default_cap; cap=cap)
    end
    println("kész.")
    return nothing
end
