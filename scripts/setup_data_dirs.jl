# Script to create required data and plotting directory structure
# Run with: julia scripts/setup_data_dirs.jl

using Dates

# Base directories
models = [
	"Duffing",
	"Lorenz",
	"Rossler",
	"VanDerPol",
]

base_data = joinpath(@__DIR__, "..", "data")

mkpath(base_data)

# Per-model folders
for m in models
	model_dir = joinpath(base_data, m)
	mkpath(model_dir)
	log_path = joinpath(model_dir, "log.txt")
	if !isfile(log_path)
		open(log_path, "w") do io
			println(io, "# Log file for $m created at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
		end
	end
end
