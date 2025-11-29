# Main entry: include model files and run selected simulations
using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("Qt5Agg")

# Toggle which simulations to run (1 = run, 0 = skip)
run_duffing = 1
run_lorenz = 1
run_rossler = 1
run_vanderpol = 1


# Execute simulations (each file currently auto-runs on include; optional explicit calls)
if run_duffing == 1
	include("my_duffing.jl")
	using .DuffingModule
	println("Running Duffing...")

	# egyszerű hívás
	err = duffing_single_run(-60, 9.0; do_plot = true)
	println("Max követési hiba: ", err)

	# saját paraméter
	# p = DuffingParams()
	# p.save_pdf = true
	# p.Robust = false
	# p.Adaptive = false

	# err = duffing_single_run(p, 15.0, 0.0; do_plot = true)

	# q_range = -60:3:60.0
	# n = length(q_range)

	# errors = zeros(n, n)

	# for (ix, i) in enumerate(q_range)
	# 	for (iy, y) in enumerate(q_range)

	# 		err = duffing_single_run(p, i, y; do_plot = false)
	# 		println("q0=", i, ", q_p0=", y, " → hiba=", err)
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
end

if run_lorenz == 1
	include("my_lorenz.jl")
	using .LorenzModule

	println("Running Lorenz simulation...")

	p = LorenzParams()
	p.save_pdf = true
	q0 = (2.0, 3.0, 1.0)
	q_p0 = (0.5, 0.0, -0.3)
	# q_pp0 = (0.0, 0.0, 0.0)

	err = simulate_lorenz(p, q0, q_p0; do_plot = true)
	println("Max követési hiba: ", err)
end

if run_rossler == 1
	include("my_rossler.jl")
	using .RosslerModule
	println("Running Rössler simulation...")
	p = RosslerParams()
	p.save_pdf = true
	q0 = (0.0, 0.0, 0.0)
	q_p0 = (1.0, 1.0, 1.0)

	err = simulate_rossler(p, q0, q_p0; do_plot = true)
	println("Max követési hiba: ", err)
end

if run_vanderpol == 1
	include("my_vanderpool.jl")
	using .VanDerPolModule

	println("Running VanDerPol simulation...")

	p = VanDerPolParams()
	p.save_pdf = true

	err = simulate_vanderpol(p, 1.0, 1.0; do_plot = true)
	println("Max követési hiba: ", err)
end

println("Minden szimuláció lefutott (vagy ki volt kapcsolva).")
