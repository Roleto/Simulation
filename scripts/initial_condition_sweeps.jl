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
# 1. DEFAULT paraméterekkel – kezdeti feltétel sweep
# ───────────────────────────────────────────────────────────────────────
function rundefault()

    # ── Duffing (default: Λ=1, K_VSSM=500) ──
    q_vec, qp_vec, E = sweep_2d((q0, qp0; kw...) -> run_duffing(q0, qp0; kw...), -20:5:20, -20:5:20)
    plot_heatmap_2d(q_vec, qp_vec, E; title_str="Duffing – DEFAULT (Λ=1, K_VSSM=500)")
    show()

    # ── Van der Pol (default: Λ=1, K_VSSM=500) ──
    q_vec, qp_vec, E = sweep_2d((q0, qp0; kw...) -> run_vanderpol(q0, qp0; kw...), -5:1:5, -5:1:5)
    plot_heatmap_2d(q_vec, qp_vec, E; title_str="VanDerPol – DEFAULT (Λ=1, K_VSSM=500)")
    show()

    # ── Lorenz (default: Λ=1, K_VSSM=50) ──
    q_vec, qp_vec, E = sweep_2d(
        (x0, y0; kw...) -> run_lorenz(x0, y0, 0.0; kw...),
        -5:1:5, -5:1:5)
    plot_heatmap_2d(q_vec, qp_vec, E; title_str="Lorenz – DEFAULT (Λ=1, K_VSSM=50)")
    show()

    # # ── Rössler (default: Λ=1, K_VSSM=50) ──
    q_vec, qp_vec, E = sweep_2d(
        (x0, y0; kw...) -> run_rossler(x0, y0, 0.0; kw...),
        -5:1:5, -5:1:5)
    plot_heatmap_2d(q_vec, qp_vec, E; title_str="Rössler – DEFAULT (Λ=1, K_VSSM=50)")
    show()

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
