# ═══════════════════════════════════════════════════════════════════════
# initial_condition_sweeps.jl
#
# Kezdeti feltétel sweep (sweep_2d) különböző paraméter-konfigurációkkal.
# Futtatás: julia --project=. (a Simulation mappából)
#   include("scripts/initial_condition_sweeps.jl")
#
# PARAMÉTER REFERENCIA:
#   DEFAULT paraméterek (my_xxx.jl fájlokban):
#     Duffing   : Λ=1,  K_VSSM=500, w=1
#     VanDerPol : Λ=1,  K_VSSM=500, w=1
#     Lorenz    : Λ=1,  K_VSSM=50,  w=1
#     Rössler   : Λ=1,  K_VSSM=50,  w=1
#
#   OPTIMÁLIS paraméterek (sweep elemzésből):
#     Duffing   : Λ=100, K_VSSM=1000, w=1  → max error ≈ 4.887e-7
#     VanDerPol : Λ=56,  K_VSSM=1000, w=1  → max error ≈ 1.257e-6
#     Lorenz    : Λ=10,  K_VSSM=100,  w=1  → max error ≈ 2.532e-3 (attraktor korlát)
#     Rössler   : Λ=10,  K_VSSM=1000, w=1  → max error ≈ 1.002e-4
# ═══════════════════════════════════════════════════════════════════════

include(joinpath(@__DIR__, "..", "src", "runner.jl"))  # betölti a sweep_2d, plot_heatmap_2d stb. függvényeket
# és a my_xxx.jl modulokat is

# ───────────────────────────────────────────────────────────────────────
# 0. Segédfüggvények – egy futás / egy kezdeti feltétel sweep
#
#   ic_run(:duffing)                                     # default kezdőérték + default paraméterek
#   ic_run(:duffing, 5.0, -3.0)                          # saját kezdőérték
#   ic_run(:duffing, 5.0, -3.0; Λ=100, K_VSSM=1000)      # + paraméter felülírás
#   ic_sweep(:duffing, -20:2:20, -20:2:20)               # állandósult hiba (maxe_ss), log színskála
#   ic_sweep(:duffing, -20:2:20, -20:2:20; metric=:maxe_all, Λ=100, K_VSSM=1000)
#
#   Összehasonlítás egy ábrán (tag = felirat):
#   ic_run(:duffing, 5.0, 0.0; Adaptive=0, tag="csak VSSM")
#   ic_run(:duffing, 5.0, 0.0; Robust=0,   tag="csak RFPT")
#   ic_run(:duffing, 5.0, 0.0;             tag="hibrid")
#   close("all")   # ábrák törlése két összehasonlítás között
#
#   metric lehet: :maxe_ss, :rms_ss, :maxe_all, :rms_all, :maxu
# ───────────────────────────────────────────────────────────────────────
# paraméterek olvasható formában a címhez / kiíráshoz: "default" vagy "Λ=100, K_VSSM=1000"
pstr(kw) = isempty(kw) ? "default" : join(["$k=$v" for (k, v) in pairs(kw)], ", ")

const MODELS = Dict(
    :duffing => run_duffing,
    :vanderpol => run_vanderpol,
    :lorenz => run_lorenz,        # ic_run(:lorenz, x0, y0, z0); ic_sweep: x0–y0 sík, z0 = 0
    :rossler => run_rossler,      # ic_run(:rossler, x0, y0, z0); ic_sweep: x0–y0 sík, z0 = 0
)

function ic_run(model::Symbol, ic...; kw...)
    m = MODELS[model](ic...; metrics=true, do_plot=true, kw...)
    println("$model  ic=$(ic)  ($(pstr(kw)))")
    println("  maxe_ss=$(m.maxe_ss)  rms_ss=$(m.rms_ss)  maxe_all=$(m.maxe_all)  maxu=$(m.maxu)")
    show()
    return m
end

# ───────────────────────────────────────────────────────────────────────
# Robusztusság: egy mutató a modellhiba (vagy zaj) függvényében, több esetre (abláció)
#
#   robust_sweep(:duffing)                                        # modellhiba ε = -0.9 … 2.0
#   robust_sweep(:duffing, 0:0.25:3; Λ=100, K_VSSM=1000)           # saját szintek + paraméterek
#   robust_sweep(:duffing, [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]; param=:noise)   # mérési zaj
#   robust_sweep(:duffing; ic=(3.0, -1.0), metric=:rms_ss)        # más kezdőérték / mutató
#
#   ε: a becsült paraméterek pₐ = (1+ε)·pₑ  →  ε=0 pontos modell, ε=-0.5 fele, ε=1 kétszerese
# ───────────────────────────────────────────────────────────────────────
const ABLATION = [("hibrid", (;)), ("csak RFPT", (Robust=0,)), ("csak VSSM", (Adaptive=0,))]

function robust_sweep(model::Symbol, levels=-0.9:0.1:2.0; param=:mismatch, metric=:maxe_ss,
                      ic=(), cases=ABLATION, kw...)
    res = Dict{String,Vector{Float64}}()
    for (name, ckw) in cases
        v = Float64[]
        for lv in levels
            m = MODELS[model](ic...; metrics=true, do_plot=false, ckw..., kw...,
                              NamedTuple{(param,)}((lv,))...)
            e = getproperty(m, metric)
            push!(v, isfinite(e) ? e : NaN)
        end
        res[name] = v
    end

    # táblázat a konzolra
    println("\n$model – $metric vs $param  ($(pstr(kw)))")
    println(rpad(string(param), 10), join([rpad(n, 14) for (n, _) in cases]))
    for (j, lv) in enumerate(levels)
        println(rpad(string(lv), 10), join([rpad(string(round(res[n][j], sigdigits=3)), 14) for (n, _) in cases]))
    end

    # ábra
    figure("$model - Robustness ($param)")
    grid(true, which="both")
    for (n, _) in cases
        semilogy(collect(levels), res[n], marker="o", markersize=3, label=n)
    end
    param == :noise && xscale("symlog", linthresh=1e-6)
    xlabel(param == :mismatch ? "modellhiba ε  (pₐ = (1+ε)·pₑ)" : "mérési zaj szórása")
    ylabel(string(metric))
    title("$(uppercasefirst(string(model))) – robusztusság  ($(pstr(kw)))")
    legend()
    show()
    return res
end

function ic_sweep(model::Symbol, r1, r2; metric=:maxe_ss, kw...)
    f = (a, b; do_plot=false) -> getproperty(MODELS[model](a, b; metrics=true, do_plot=false, kw...), metric)
    q, qp, E = sweep_2d(f, r1, r2)
    E[.!isfinite.(E)] .= NaN          # divergált futás → szürke cella
    Ef = filter(isfinite, E)
    println("$model  $metric:  max=$(isempty(Ef) ? NaN : maximum(Ef))  min=$(isempty(Ef) ? NaN : minimum(Ef))  nem véges: $(count(!isfinite, E))/$(length(E))")
    is3d = model in (:lorenz, :rossler)          # 3D modell: x0–y0 sík (z0 ≈ 0)
    plot_heatmap_2d(q, qp, E; logscale=true, label=string(metric), title_str="$(uppercasefirst(string(model))) – $metric  ($(pstr(kw)))",
        xlab=is3d ? "y₀  (initial y)" : "q̇₀  (initial velocity)", ylab=is3d ? "x₀  (initial x)" : "q₀  (initial position)")
    show()
    return q, qp, E
end

# ───────────────────────────────────────────────────────────────────────
# 1. DEFAULT paraméterekkel – kezdeti feltétel sweep
# ───────────────────────────────────────────────────────────────────────
function rundefault()

    # ── Duffing (default: Λ=1, K_VSSM=500) ──
    ic_sweep(:duffing, -20:2:20, -20:2:20)   # állandósult hiba (maxe_ss), log skála

    # ── Van der Pol (default: Λ=1, K_VSSM=500) ──
    ic_sweep(:vanderpol, -5:0.5:5, -5:0.5:5)   # állandósult hiba (maxe_ss), log skála

    # ── Lorenz (default: Λ=1, K_VSSM=50) ──
    ic_sweep(:lorenz, -5:0.5:5, -5:0.5:5)   # x0–y0 sík, z0 = 0; állandósult hiba, log skála

    # # ── Rössler (default: Λ=1, K_VSSM=50) ──
    ic_sweep(:rossler, -5:0.5:5, -5:0.5:5)   # x0–y0 sík, z0 ≈ 0; állandósult hiba, log skála

end
# ───────────────────────────────────────────────────────────────────────
# 2. OPTIMÁLIS paraméterekkel – kezdeti feltétel sweep (kommentezd ki, ha kell)
# ───────────────────────────────────────────────────────────────────────

# ── Duffing (optimális: Λ=100, K_VSSM=1000) ──
# q_vec, qp_vec, E = sweep_2d(
#     (q0, qp0; kw...) -> run_duffing(q0, qp0; Λ=100, K_VSSM=1000, kw...),
#     -20:5:20, -20:5:20)
# plot_heatmap_2d(q_vec, qp_vec, E; title_str="Duffing – OPTIMAL (Λ=100, K_VSSM=1000)")
# show()

# ── Van der Pol (optimális: Λ=56, K_VSSM=1000) ──
# q_vec, qp_vec, E = sweep_2d(
#     (q0, qp0; kw...) -> run_vanderpol(q0, qp0; Λ=56, K_VSSM=1000, kw...),
#     -5:1:5, -5:1:5)
# plot_heatmap_2d(q_vec, qp_vec, E; title_str="VanDerPol – OPTIMAL (Λ=56, K_VSSM=1000)")
# show()

# ── Lorenz (optimális: Λ=10, K_VSSM=100) ──
# q_vec, qp_vec, E = sweep_2d(
#     (x0, y0; kw...) -> run_lorenz(x0, y0, 0.0; Λ=10, K_VSSM=100, kw...),
#     -5:1:5, -5:1:5)
# plot_heatmap_2d(q_vec, qp_vec, E; title_str="Lorenz – OPTIMAL (Λ=10, K_VSSM=100)")
# show()

# ── Rössler (optimális: Λ=10, K_VSSM=1000) ──
# q_vec, qp_vec, E = sweep_2d(
#     (x0, y0; kw...) -> run_rossler(x0, y0, 0.0; Λ=10, K_VSSM=1000, kw...),
#     -5:1:5, -5:1:5)
# plot_heatmap_2d(q_vec, qp_vec, E; title_str="Rössler – OPTIMAL (Λ=10, K_VSSM=1000)")
# show()

# ───────────────────────────────────────────────────────────────────────
# 3. Egyedi futtatás – single run vizsgálathoz
# ───────────────────────────────────────────────────────────────────────

# run_duffing(10.0, -10.0; Λ=100, K_VSSM=1000, do_plot=true); show()   # optimális, messze induló
# run_duffing(10.0, -10.0; do_plot=true); show()                         # default, messze induló
# run_lorenz(1.0, 0.5, 0.2; Λ=10, K_VSSM=100, do_plot=true); show()    # Lorenz optimális
