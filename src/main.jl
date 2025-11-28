# Main entry: include model files and run selected simulations
using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("Qt5Agg")

# Toggle which simulations to run (1 = run, 0 = skip)
run_duffing = 1
run_lorenz = 0
run_rossler = 0
run_vanderpol = 0

include("my_duffing.jl")
using .DuffingModule
# Include model source files
# include("my_lorenz.jl")
# include("my_rossler.jl")
# include("my_vanderpool.jl")

# Execute simulations (each file currently auto-runs on include; optional explicit calls)
if run_duffing == 1 && isdefined(@__MODULE__, :duffing_single_run)
	println("Running Duffing single run...")

	# egyszerű hívás
	err = duffing_single_run(6.5, 1.0; do_plot = true)

	# saját paraméter
	# p = DuffingParams()
	# p.δt = 1e-3
	# p.N = Int(2e4)
	# p.save_pdf = false
	# p.pdf_dir = "data/Duffing/"
	# q_range = -120:1:120.0
	# n = length(q_range)

	# errors = zeros(n, n)

	# for (ix, i) in enumerate(q_range)
	# 	for (iy, y) in enumerate(q_range)

	# 		err = duffing_single_run(p, i, y; do_plot = false)

	# 		errors[ix, iy] = err
	# 	end
	# end

	# using PyPlot
	# using PyPlot: plot_surface

	# figure("Image_plot")
	# imshow(errors, extent = (first(q_range), last(q_range), first(q_range), last(q_range)), origin = "lower")
	# colorbar()
	# title("Hibatér (Grid search)")
	# xlabel("q₀ (kezdő pozíció)")
	# xlabel("q0")
	# ylabel("q̇₀ (kezdő sebesség)")
	# ylabel("q_p0")

	# X = repeat(collect(q_range)', n, 1)
	# Y = repeat(collect(q_range), 1, n)

	# fig = figure("Surface_plot")
	# ax = fig.add_subplot(111, projection = "3d")
	# ax.plot_surface(X, Y, errors, cmap = "viridis")
	# xlabel("q0")
	# ylabel("q_p0")
	# zlabel("hiba")
	# title("Duffing hibafelület")


	# q_vec = collect(q_range)
	# X = repeat(q_vec, inner = n)
	# Y = repeat(q_vec', outer = n)
	# Z = vec(errors)

	# figure("Scatter_plot")
	# scatter(X, Y, c = Z, cmap = "viridis")
	# colorbar()
	# xlabel("q₀ (kezdő pozíció)")
	# ylabel("q̇₀ (kezdő sebesség)")
	# title("Duffing kezdőállapot → maximális hiba")


	# show()
	# ax = fig.add_subplot(111, projection = "3d")

	# ax.scatter(X, Y, Z, c = Z, cmap = "viridis")
	# xlabel("q₀ (kezdő pozíció)")
	# ylabel("q̇₀ (kezdő sebesség)")
	# zlabel("hiba")
	# title("Duffing hibafelület – 3D scatter")

	# figure("Test4")



	# err2 = duffing_single_run(p, 0.0, 0.0; do_plot = true)

	# # grid search
	# res, csvfile = grid_search(p, -1.0:0.2:1.0, -1.0:0.2:1.0; plot = false)

	# duffing_single_run(0, 0; plot = true)
end

if run_lorenz == 1 && isdefined(@__MODULE__, :simulate_lorenz)
	println("Running Lorenz simulation...")
	simulate_lorenz() # already auto-run
end

if run_rossler == 1 && isdefined(@__MODULE__, :simulate_rossler)
	println("Running Rössler simulation...")
	simulate_rossler() # already auto-run
end

if run_vanderpol == 1 && isdefined(@__MODULE__, :vanderpol_single_run)
	println("Running VanDerPol simulation...")
	vanderpol_single_run() # underlying call in file through simulate_vanderpol
end

# println("Minden szimuláció lefutott (vagy ki volt kapcsolva).")
