# Main entry: include model files and run selected simulations
using PyPlot
using Pkg
PyPlot.matplotlib.use("TkAgg")
Pkg.activate(".")
# Toggle which simulations to run (1 = run, 0 = skip)
run_duffing = 1
run_lorenz = 1
run_rossler = 1
run_vanderpol = 1

# Include model source files
include("my_duffing.jl")
include("my_lorenz.jl")
include("my_rossler.jl")
include("my_vanderpool.jl")

# Execute simulations (each file currently auto-runs on include; optional explicit calls)
if run_duffing == 1 && isdefined(@__MODULE__, :simulate_duffing)
	println("Running Duffing simulation...")
	# Already executed inside file; could re-run if desired
	simulate_duffing()
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

println("Minden szimuláció lefutott (vagy ki volt kapcsolva).")
