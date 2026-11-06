# 1) Duffing
using DelimitedFiles

q_range = -100:10:100
qp_range = -100:10:100

q_vec, qp_vec, E = sweep_2d(run_duffing, q_range, qp_range)
writedlm("../data/Errors/Duffing/duffing_heatmap.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Duffing – Error Heatmap")
savefig("../data/Errors/Duffing/duffing_heatmap.png", dpi=150, bbox_inches="tight")

# 2) Van der Pol
q_vec, qp_vec, E = sweep_2d(run_vanderpol, q_range, qp_range)
writedlm("../data/Errors/VanDerPol/vanderpol_heatmap.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="VanDerPol – Error Heatmap")
savefig("../data/Errors/VanDerPol/vanderpol_heatmap.png", dpi=150, bbox_inches="tight")

# 3) Lorenz
q_vec, qp_vec, E = sweep_2d((x0, y0; do_plot=false) -> run_lorenz(x0, y0, 0.0; do_plot=do_plot), q_range, qp_range)
writedlm("../data/Errors/Lorenz/lorenz_heatmap.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Lorenz – Error Heatmap (XY plane, z0=0)")
savefig("../data/Errors/Lorenz/lorenz_heatmap.png", dpi=150, bbox_inches="tight")

q_vec, qp_vec, E = sweep_2d((x0, y0; do_plot=false) -> run_lorenz(0.0, x0, y0; do_plot=do_plot), q_range, qp_range)
writedlm("../data/Errors/Lorenz/lorenz_heatmap_yz.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Lorenz – Error Heatmap (YZ plane, x0=0)")
savefig("../data/Errors/Lorenz/lorenz_heatmap_yz.png", dpi=150, bbox_inches="tight")

q_vec, qp_vec, E = sweep_2d((x0, y0; do_plot=false) -> run_lorenz(x0, 0.0, y0; do_plot=do_plot), q_range, qp_range)
writedlm("../data/Errors/Lorenz/lorenz_heatmap_xz.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Lorenz – Error Heatmap (XZ plane, y0=0)")
savefig("../data/Errors/Lorenz/lorenz_heatmap_xz.png", dpi=150, bbox_inches="tight")

# 6) Rossler
q_vec, qp_vec, E = sweep_2d((x0, y0; do_plot=false) -> run_rossler(x0, y0, 0.0; do_plot=do_plot), q_range, qp_range)
writedlm("../data/Errors/Rossler/rossler_heatmap.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Rössler – Error Heatmap (XY plane, z0=0)")
savefig("../data/Errors/Rossler/rossler_heatmap.png", dpi=150, bbox_inches="tight")

q_vec, qp_vec, E = sweep_2d((x0, y0; do_plot=false) -> run_rossler(0.0, x0, y0; do_plot=do_plot), q_range, qp_range)
writedlm("../data/Errors/Rossler/rossler_heatmap_yz.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Rössler – Error Heatmap (YZ plane, x0=0)")
savefig("../data/Errors/Rossler/rossler_heatmap_yz.png", dpi=150, bbox_inches="tight")

q_vec, qp_vec, E = sweep_2d((x0, y0; do_plot=false) -> run_rossler(x0, 0.0, y0; do_plot=do_plot), q_range, qp_range)
writedlm("../data/Errors/Rossler/rossler_heatmap_xz.csv", E, ',')
plot_heatmap_2d(q_vec, qp_vec, E; title_str="Rössler – Error Heatmap (XZ plane, y0=0)")
savefig("../data/Errors/Rossler/rossler_heatmap_xz.png", dpi=150, bbox_inches="tight")