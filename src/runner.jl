using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("TkAgg")
PyPlot.ion() # interactive mode: allows plots to update without blocking the REPL
# ─────────────────────────────────────────
# Load simulation modules
# ─────────────────────────────────────────
include("my_duffing.jl")
include("my_lorenz.jl")
include("my_rossler.jl")
include("my_vanderpool.jl")

using .DuffingModule
using .LorenzModule
using .RosslerModule
using .VanDerPolModule

# ─────────────────────────────────────────
# sweep_2d
#
# Runs any 2-argument simulation function over a grid of initial conditions.
# Returns (q_vec, qp_vec, errors_matrix).
#
# Usage for SISO systems:
#   q_vec, qp_vec, E = sweep_2d(run_duffing,   -20:5:20, -20:5:20)
#   q_vec, qp_vec, E = sweep_2d(run_vanderpol,  -5:1:5,   -5:1:5)
#
# Usage for MIMO systems (fix third coordinate via a closure):
#   q_vec, qp_vec, E = sweep_2d((x0, y0) -> run_lorenz(x0, y0, 0.0), -5:1:5, -5:1:5)
# ─────────────────────────────────────────
function sweep_2d(run_fn, q_range, qp_range)
    q_vec = collect(q_range)
    qp_vec = collect(qp_range)
    n_q = length(q_vec)
    n_qp = length(qp_vec)
    errors = zeros(n_q, n_qp)
    total = n_q * n_qp
    done = 0
    for (i, q0) in enumerate(q_vec)
        for (j, qp0) in enumerate(qp_vec)
            errors[i, j] = run_fn(q0, qp0; do_plot=false)
            done += 1
            # print("\rProgress: $done / $total  (q0=$q0, qp0=$qp0)")
        end
    end
    println()
    return q_vec, qp_vec, errors
end

# ─────────────────────────────────────────
# plot_scatter_3d
#
# 3-D scatter: x = initial position, y = initial velocity, z = maxhiba
# ─────────────────────────────────────────
function plot_scatter_3d(q_vec, qp_vec, errors; title_str="Error Scatter")
    X = vec([q for q in q_vec, _ in qp_vec])
    Y = vec([qp for _ in q_vec, qp in qp_vec])
    Z = vec(errors)

    fig = figure()
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(X, Y, Z, c=Z, cmap="viridis", s=25)
    fig.colorbar(sc, ax=ax, label="maxhiba")
    ax.set_xlabel("q₀  (initial position)")
    ax.set_ylabel("q̇₀  (initial velocity)")
    ax.set_zlabel("maxhiba")
    ax.set_title(title_str)
    return fig
end

# ─────────────────────────────────────────
# plot_surface_3d
#
# 3-D surface: x = initial position, y = initial velocity, z = maxhiba
# ─────────────────────────────────────────
function plot_surface_3d(q_vec, qp_vec, errors; title_str="Error Surface")
    X = [q for q in q_vec, _ in qp_vec]
    Y = [qp for _ in q_vec, qp in qp_vec]

    fig = figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, errors, cmap="viridis", alpha=0.85)
    ax.set_xlabel("q₀  (initial position)")
    ax.set_ylabel("q̇₀  (initial velocity)")
    ax.set_zlabel("maxhiba")
    ax.set_title(title_str)
    return fig
end

# ─────────────────────────────────────────
# plot_heatmap_2d
#
# 2-D heatmap: x axis = initial velocity, y axis = initial position, colour = maxhiba
# ─────────────────────────────────────────
function plot_heatmap_2d(q_vec, qp_vec, errors; title_str="Error Heatmap")
    fig, ax = subplots()
    pc = ax.pcolormesh(qp_vec, q_vec, errors, cmap="viridis", shading="auto")
    fig.colorbar(pc, ax=ax, label="maxhiba")
    ax.set_xlabel("q̇₀  (initial velocity)")
    ax.set_ylabel("q₀  (initial position)")
    ax.set_title(title_str)
    return fig
end

# ─────────────────────────────────────────
# PARAMÉTER-SWEEP  (A, K, Λ, K_VSSM, w ...)  – minden modellre
#
# Bármely run_xxx-szel megy, ami kulcsszavakat fogad és skalár hibát ad vissza
# (run_duffing, run_vanderpol, run_lorenz, run_rossler). A modellfájlokhoz nem kell nyúlni.
# Divergált / NaN / Inf / cap fölötti hiba -> NaN (szürke a képen). Nincs korai kilépés.
#
# Usage:
#   a, b, E = sweep_params(run_duffing, :Λ, 10 .^ range(0,2,length=9), :K_VSSM, 10 .^ range(0,4,length=9))
#   plot_heatmap_params(a, b, E, :Λ, :K_VSSM; title_str="Duffing")
# ─────────────────────────────────────────
function sweep_params(run_fn, p1::Symbol, v1, p2::Symbol, v2; cap=1e3, kw...)
    a_vec = collect(v1)
    b_vec = collect(v2)
    E = fill(NaN, length(b_vec), length(a_vec))          # E[j,i]: j = p2, i = p1
    for (j, b) in enumerate(b_vec), (i, a) in enumerate(a_vec)
        e = run_fn(; Dict(p1 => a, p2 => b)..., kw..., do_plot=false)
        E[j, i] = (isfinite(e) && e < cap) ? e : NaN
    end
    println("$p1×$p2: $(count(!isnan, E))/$(length(E)) stabil")
    return a_vec, b_vec, E
end

# 2-D hőtérkép (log10 hiba, log tengelyek). A×K párnál A·K = 0.1 / 1 / 3 vonalakat is rajzol.
function plot_heatmap_params(a_vec, b_vec, E, p1::Symbol, p2::Symbol;
    title_str="Parameter sweep", ak_lines=true, savepath=nothing)
    cm = get_cmap("viridis")
    cm.set_bad("lightgray")
    fig, ax = subplots()
    pc = ax.pcolormesh(a_vec, b_vec, log10.(E), cmap=cm, shading="nearest")
    fig.colorbar(pc, ax=ax, label="log10(maxhiba)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    if ak_lines && p1 == :A && p2 == :K
        for c in (0.1, 1.0, 3.0)
            ax.plot(a_vec, c ./ a_vec, color="white", linestyle="--", linewidth=1)
            ax.annotate("A·K=$c", xy=(a_vec[end], c / a_vec[end]), color="white", fontsize=8,
                ha="right", va="bottom")
        end
        ax.set_xlim(minimum(a_vec), maximum(a_vec))
        ax.set_ylim(minimum(b_vec), maximum(b_vec))
    end
    ax.set_xlabel(string(p1))
    ax.set_ylabel(string(p2))
    ax.set_title(title_str)
    fig.tight_layout()
    savepath === nothing || fig.savefig(savepath)
    return fig
end

# 1-D sweep a P = A·K szorzaton (K fix, A = P/K) – az A×K rács helyett.
function sweep_AK_1d(run_fn, P_vec, K; cap=1e3, kw...)
    P = collect(P_vec)
    E = fill(NaN, length(P))
    for (i, p) in enumerate(P)
        e = run_fn(; A=p / K, K=K, kw..., do_plot=false)
        E[i] = (isfinite(e) && e < cap) ? e : NaN
    end
    println("A·K sweep: $(count(!isnan, E))/$(length(E)) stabil")
    return P, E
end

function plot_AK_1d(P, E; title_str="A·K sweep", savepath=nothing)
    fig, ax = subplots()
    ax.loglog(P, E, "o-")
    ax.set_xlabel("A·K")
    ax.set_ylabel("maxhiba")
    ax.set_title(title_str)
    ax.grid(true, which="both", alpha=0.3)
    fig.tight_layout()
    savepath === nothing || fig.savefig(savepath)
    return fig
end

using DelimitedFiles

# A·K sweep mentése: CSV + PNG a data/Sweeps/<name>/ mappába.
function save_AK_1d(name, P, E; title_str="$name – A·K")
    out = joinpath(@__DIR__, "..", "data", "Sweeps", name)
    mkpath(out)
    csv = joinpath(out, "AK_1d.csv")
    writedlm(csv, vcat(["A*K" "maxhiba"], hcat(P, E)), ',')
    png = joinpath(out, "AK_1d.png")
    plot_AK_1d(P, E; title_str=title_str, savepath=png)
    println("CSV → ", normpath(csv))
    println("kép → ", normpath(png))
end

# 2D sweep mentése: CSV (mátrix) + hőtérkép PNG a data/Sweeps/<name>/ mappába.
function save_sweep_2d(name, a_vec, b_vec, E, p1::Symbol, p2::Symbol)
    out = joinpath(@__DIR__, "..", "data", "Sweeps", name)
    mkpath(out)
    csv = joinpath(out, "sweep_$(p1)_$(p2).csv")
    writedlm(csv, vcat(permutedims(vcat("", a_vec)), hcat(b_vec, E)), ',')
    png = joinpath(out, "sweep_$(p1)_$(p2).png")
    plot_heatmap_params(a_vec, b_vec, E, p1, p2; title_str="$name – $p1 × $p2", savepath=png)
    println("CSV → ", normpath(csv))
    println("kép → ", normpath(png))
end

# ═══════════════════════════════════════════════════════════
# USAGE EXAMPLES  –  uncomment the block you want to run
# ═══════════════════════════════════════════════════════════

# ── 1. Single run – default initial conditions ─────────────
# err = run_duffing()
# println("Duffing  default maxhiba: ", err)
# show()

# err = run_vanderpol()
# println("VanDerPol default maxhiba: ", err)
# show()

# err = run_lorenz()
# println("Lorenz  default maxhiba: ", err)
# show()

# err = run_rossler()
# println("Rössler default maxhiba: ", err)
# show()

# ── 2. Single run – custom initial conditions ───────────────
# err = run_duffing(5.0, -3.0; do_plot=true)
# println("Duffing  custom maxhiba: ", err)
# show()

# err = run_vanderpol(1.0, 0.5; do_plot=true)
# println("VanDerPol custom maxhiba: ", err)
# show()

# err = run_lorenz(1.0, 0.5, 0.2; do_plot=true)
# println("Lorenz  custom maxhiba: ", err)
# show()

# err = run_rossler(0.5, 1.0, 0.0; do_plot=true)
# println("Rössler custom maxhiba: ", err)
# show()

# ── 3. Grid sweep – Duffing ────────────────────────────────
# q_range  = -20:5:20
# qp_range = -20:5:20
# q_vec, qp_vec, E = sweep_2d(run_duffing, q_range, qp_range)
# plot_scatter_3d(q_vec, qp_vec, E; title_str="Duffing – Error Scatter")
# plot_surface_3d(q_vec, qp_vec, E; title_str="Duffing – Error Surface")
# show()

# ── 4. Grid sweep – Van der Pol ───────────────────────────
# q_range  = -5:1:5
# qp_range = -5:1:5
# q_vec, qp_vec, E = sweep_2d(run_vanderpol, q_range, qp_range)
# plot_scatter_3d(q_vec, qp_vec, E; title_str="VanDerPol – Error Scatter")
# plot_surface_3d(q_vec, qp_vec, E; title_str="VanDerPol – Error Surface")
# show()

# ── 5. Grid sweep – Lorenz (sweep x0, y0; fix z0=0) ──────
# x_range = -5:1:5
# y_range = -5:1:5
# q_vec, qp_vec, E = sweep_2d((x0, y0) -> run_lorenz(x0, y0, 0.0), x_range, y_range)
# plot_scatter_3d(q_vec, qp_vec, E; title_str="Lorenz – Error Scatter")
# plot_surface_3d(q_vec, qp_vec, E; title_str="Lorenz – Error Surface")
# show()

# ── 6. Grid sweep – Rössler (sweep x0, y0; fix z0=0) ────
# x_range = -5:1:5
# y_range = -5:1:5
# q_vec, qp_vec, E = sweep_2d((x0, y0) -> run_rossler(x0, y0, 0.0), x_range, y_range)
# plot_scatter_3d(q_vec, qp_vec, E; title_str="Rössler – Error Scatter")
# plot_surface_3d(q_vec, qp_vec, E; title_str="Rössler – Error Surface")
# show()

# ── 7. Paraméter-sweep (A·K, Λ, K_VSSM, w ...) – bármely modell ──
# lg(a, b, n) = 10.0 .^ range(a, b, length=n)
# P, E = sweep_AK_1d(run_duffing, 10 .^ range(-3, 0.5, length=25), 1e5)      # Duffing/VdP: K=1e5
# plot_AK_1d(P, E; title_str="Duffing – A·K")
# P, E = sweep_AK_1d(run_lorenz, 10 .^ range(-3, 0.5, length=25), 1e2)       # Lorenz/Rössler: K=1e2
# a, b, E = sweep_params(run_duffing, :Λ, 10 .^ range(0, 2, length=9), :K_VSSM, 10 .^ range(0, 4, length=9))
# plot_heatmap_params(a, b, E, :Λ, :K_VSSM; title_str="Duffing – Λ × K_VSSM")
# a, b, E = sweep_params(run_vanderpol, :A, 10 .^ range(-7, -3, length=9), :K, 10 .^ range(3, 7, length=9))
# plot_heatmap_params(a, b, E, :A, :K; title_str="VdP – A × K")   # A·K vonalakkal
# show()
