using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("Qt5Agg")

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
    q_vec  = collect(q_range)
    qp_vec = collect(qp_range)
    n_q    = length(q_vec)
    n_qp   = length(qp_vec)
    errors = zeros(n_q, n_qp)
    total  = n_q * n_qp
    done   = 0
    for (i, q0) in enumerate(q_vec)
        for (j, qp0) in enumerate(qp_vec)
            errors[i, j] = run_fn(q0, qp0; do_plot=false)
            done += 1
            print("\rProgress: $done / $total  (q0=$q0, qp0=$qp0)")
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
    X = vec([q  for q  in q_vec,  _  in qp_vec])
    Y = vec([qp for _  in q_vec,  qp in qp_vec])
    Z = vec(errors)

    fig = figure()
    ax  = fig.add_subplot(111, projection="3d")
    sc  = ax.scatter(X, Y, Z, c=Z, cmap="viridis", s=25)
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
    X = [q  for q  in q_vec,  _  in qp_vec]
    Y = [qp for _  in q_vec,  qp in qp_vec]

    fig = figure()
    ax  = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, errors, cmap="viridis", alpha=0.85)
    ax.set_xlabel("q₀  (initial position)")
    ax.set_ylabel("q̇₀  (initial velocity)")
    ax.set_zlabel("maxhiba")
    ax.set_title(title_str)
    return fig
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
