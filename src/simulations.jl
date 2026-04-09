# julia --project=.
PyPlot.ioff()

const SWEEP_RANGE = -60:2:60
const SAVE_PLOTS = false   # set to false to skip saving

function log_run(model_name, ic_str, maxhiba)
    println("$model_name, $ic_str,   maxhiba = $maxhiba")
    show()
end

if SAVE_PLOTS
    mkpath(joinpath(@__DIR__, "..", "data", "Errors", "Duffing"))
    mkpath(joinpath(@__DIR__, "..", "data", "Errors", "VanDerPol"))
    mkpath(joinpath(@__DIR__, "..", "data", "Errors", "Lorenz"))
    mkpath(joinpath(@__DIR__, "..", "data", "Errors", "Rossler"))
end

# ── Duffing (10 runs) ────────────────────────────────────────
# err = run_duffing(; do_plot=true);
# log_run("Duffing", "q0=default, qp0=default", err);
# err = run_duffing(0.0, 0.0; do_plot=true);
# log_run("Duffing", "q0=0.0, qp0=0.0", err);
# err = run_duffing(5.0, 45.0; do_plot=true);
# log_run("Duffing", "q0=5.0, qp0=45.0", err);
# err = run_duffing(15.0, 30.0; do_plot=true);
# log_run("Duffing", "q0=15.0, qp0=30.0", err);
# err = run_duffing(23.0, 67.0; do_plot=true);
# log_run("Duffing", "q0=23.0, qp0=67.0", err);
# err = run_duffing(42.0, 7.0; do_plot=true);
# log_run("Duffing", "q0=42.0, qp0=7.0", err);
# err = run_duffing(0.0, 125.0; do_plot=true);
# log_run("Duffing", "q0=0.0, qp0=125.0", err);
# err = run_duffing(0.0, 150.0; do_plot=true);
# log_run("Duffing", "q0=0.0, qp0=150.0", err);
# err = run_duffing(50.0, 90.0; do_plot=true);
# log_run("Duffing", "q0=50.0, qp0=90.0", err);
# err = run_duffing(61.0, 38.0; do_plot=true);
# log_run("Duffing", "q0=61.0, qp0=38.0", err);
# err = run_duffing(73.0, 12.0; do_plot=true);
# log_run("Duffing", "q0=73.0, qp0=12.0", err);
# err = run_duffing(88.0, 55.0; do_plot=true);
# log_run("Duffing", "q0=88.0, qp0=55.0", err);

# # err = run_duffing(0.0, 1000.0; do_plot=true);
# # show()

# q_vec, qp_vec, E = sweep_2d(run_duffing, SWEEP_RANGE, SWEEP_RANGE)
# println("errors: ", E)
# fig = plot_heatmap_2d(q_vec, qp_vec, E; title_str="Duffing – Error Heatmap")
if SAVE_PLOTS
    fig.savefig(joinpath(@__DIR__, "..", "data", "Errors", "Duffing", "heatmap_duffing.png"), dpi=150, bbox_inches="tight")
end
show()

# ── Van der Pol (10 runs) ────────────────────────────────────
# err = run_vanderpol(; do_plot=true);
# log_run("VanDerPol", "q0=default, qp0=default", err);
# err = run_vanderpol(1.0, 0.5; do_plot=true);
# log_run("VanDerPol", "q0=1.0, qp0=0.5", err);
# err = run_vanderpol(5.0, 70.0; do_plot=true);
# log_run("VanDerPol", "q0=5.0, qp0=70.0", err);
# err = run_vanderpol(10.0, 25.0; do_plot=true);
# log_run("VanDerPol", "q0=10.0, qp0=25.0", err);
# err = run_vanderpol(20.0, 95.0; do_plot=true);
# log_run("VanDerPol", "q0=20.0, qp0=95.0", err);
err = run_vanderpol(35.0, 15.0; do_plot=true);
log_run("VanDerPol", "q0=35.0, qp0=15.0", err);
# # run_vanderpol(35.0, 15.0; K_VSSM=2000, Λ=2, do_plot=true)
# # show()
# err = run_vanderpol(48.0, 33.0; do_plot=true);
# log_run("VanDerPol", "q0=48.0, qp0=33.0", err);
# err = run_vanderpol(60.0, 80.0; do_plot=true);
# log_run("VanDerPol", "q0=60.0, qp0=80.0", err);
# err = run_vanderpol(77.0, 62.0; do_plot=true);
# log_run("VanDerPol", "q0=77.0, qp0=62.0", err);
# err = run_vanderpol(90.0, 45.0; do_plot=true);
# log_run("VanDerPol", "q0=90.0, qp0=45.0", err);

q_vec, qp_vec, E = sweep_2d(run_vanderpol, SWEEP_RANGE, SWEEP_RANGE)
fig = plot_heatmap_2d(q_vec, qp_vec, E; title_str="VanDerPol – Error Heatmap")
if SAVE_PLOTS
    fig.savefig(joinpath(@__DIR__, "..", "data", "Errors", "VanDerPol", "heatmap_vanderpol.png"), dpi=150, bbox_inches="tight")
end
show()

# # ── Lorenz (10 runs, sweep fixes z0=0) ──────────────────────
# err = run_lorenz(; do_plot=true);
# log_run("Lorenz", "x0=default, y0=default, z0=default", err);
# err = run_lorenz(1.0, 0.5, 0.2; do_plot=true);
# log_run("Lorenz", "x0=1.0, y0=0.5, z0=0.2", err);
# err = run_lorenz(5.0, 70.0, 85.0; do_plot=true);
# log_run("Lorenz", "x0=5.0, y0=70.0, z0=85.0", err);
# err = run_lorenz(10.0, 25.0, 50.0; do_plot=true);
# log_run("Lorenz", "x0=10.0, y0=25.0, z0=50.0", err);
# err = run_lorenz(20.0, 95.0, 65.0; do_plot=true);
# log_run("Lorenz", "x0=20.0, y0=95.0, z0=65.0", err);
# err = run_lorenz(35.0, 15.0, 70.0; do_plot=true);
# log_run("Lorenz", "x0=35.0, y0=15.0, z0=70.0", err);
# err = run_lorenz(48.0, 33.0, 12.0; do_plot=true);
# log_run("Lorenz", "x0=48.0, y0=33.0, z0=12.0", err);
# err = run_lorenz(60.0, 80.0, 30.0; do_plot=true);
# log_run("Lorenz", "x0=60.0, y0=80.0, z0=30.0", err);
# err = run_lorenz(77.0, 62.0, 40.0; do_plot=true);
# log_run("Lorenz", "x0=77.0, y0=62.0, z0=40.0", err);
# err = run_lorenz(90.0, 45.0, 20.0; do_plot=true);
# log_run("Lorenz", "x0=90.0, y0=45.0, z0=20.0", err);
# q_vec, qp_vec, E = sweep_2d((x0, y0) -> run_lorenz(x0, y0, 0.0), SWEEP_RANGE, SWEEP_RANGE)
# fig = plot_heatmap_2d(q_vec, qp_vec, E; title_str="Lorenz – Error Heatmap")
# if SAVE_PLOTS
#     fig.savefig(joinpath(@__DIR__, "..", "data", "Errors", "Lorenz", "heatmap_lorenz.png"), dpi=150, bbox_inches="tight")
# end
# show()

# # ── Rössler (10 runs, sweep fixes z0=0) ─────────────────────
# err = run_rossler(; do_plot=true);
# log_run("Rossler", "x0=default, y0=default, z0=default", err);
# err = run_rossler(0.5, 1.0, 0.0; do_plot=true);
# log_run("Rossler", "x0=0.5, y0=1.0, z0=0.0", err);
# err = run_rossler(8.0, 22.0, 55.0; do_plot=true);
# log_run("Rossler", "x0=8.0, y0=22.0, z0=55.0", err);
# err = run_rossler(12.0, 68.0, 90.0; do_plot=true);
# log_run("Rossler", "x0=12.0, y0=68.0, z0=90.0", err);
# err = run_rossler(25.0, 92.0, 60.0; do_plot=true);
# log_run("Rossler", "x0=25.0, y0=92.0, z0=60.0", err);
# err = run_rossler(40.0, 18.0, 75.0; do_plot=true);
# log_run("Rossler", "x0=40.0, y0=18.0, z0=75.0", err);
# err = run_rossler(52.0, 30.0, 10.0; do_plot=true);
# log_run("Rossler", "x0=52.0, y0=30.0, z0=10.0", err);
# err = run_rossler(65.0, 78.0, 32.0; do_plot=true);
# log_run("Rossler", "x0=65.0, y0=78.0, z0=32.0", err);
# err = run_rossler(72.0, 60.0, 45.0; do_plot=true);
# log_run("Rossler", "x0=72.0, y0=60.0, z0=45.0", err);
# err = run_rossler(85.0, 50.0, 18.0; do_plot=true);
# log_run("Rossler", "x0=85.0, y0=50.0, z0=18.0", err);
# q_vec, qp_vec, E = sweep_2d((x0, y0) -> run_rossler(x0, y0, 0.0), SWEEP_RANGE, SWEEP_RANGE)
# fig = plot_heatmap_2d(q_vec, qp_vec, E; title_str="Rössler – Error Heatmap")
# if SAVE_PLOTS
#     fig.savefig(joinpath(@__DIR__, "..", "data", "Errors", "Rossler", "heatmap_rossler.png"), dpi=150, bbox_inches="tight")
# end
# show()

PyPlot.ion()

