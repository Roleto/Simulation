# julia --project=.
cd("src")
include("my_vanderpool.jl")
using .VanDerPolModule

err = run_vanderpol(0.0, 1.0; do_plot=true)
println("VanDerPol maxhiba: ", err)