"""
DuffingModule.jl

Refactored, self-contained Duffing oscillator simulator.
- No globals
- Clear parameter struct
- Fast, type-stable arrays
- Optional plotting and PDF export
- Convenience wrappers for single-run and grid-run

Usage (from your main.jl or REPL):

	# save this file somewhere, e.g. /mnt/data/DuffingModule.jl
	include("/mnt/data/DuffingModule.jl")
	using .DuffingModule

	# simple single run with default params
	err = duffing_single_run(0.0, 0.0; plot=true)

	# or pass custom params
	p = DuffingParams(); p.δt = 1e-3; p.N = Int(2e4)
	err = duffing_single_run(p, 0.0, 0.0; plot=false)

	# grid search on initial conditions (returns array of errors and saves CSV)
	results = grid_search(p, -1.0:0.1:1.0, -1.0:0.1:1.0; plot=false)

"""

module DuffingModule

using Dates
using PyPlot
using DelimitedFiles

export DuffingParams, duffing_single_run, grid_search, simulate_duffing

# -----------------------------
# Parameter container
# -----------------------------
mutable struct DuffingParams
	# Control flags
	Adaptive::Bool
	Robust::Bool
	# Time
	δt::Float64
	N::Int
	# Control params
	K::Float64
	B::Float64
	A::Float64
	# Kinematic block
	Λ::Float64
	w::Float64
	K_VSSM::Float64
	# Nominal trajectory
	ω::Float64
	Amp::Float64
	# Exact/Approx params
	αe::Float64
	δe::Float64
	βe::Float64
	αa::Float64
	δa::Float64
	βa::Float64
	# Output options
	save_pdf::Bool
	pdf_dir::String
end

# default constructor
function DuffingParams(; Adaptive = true, Robust = true, δt = 1e-3, N = Int(2e4), K = 1e5, B = -1.0, A = 1e-5,
	Λ = 1.0, w = 1.0, K_VSSM = 500.0, ω = 0.5, Amp = 2.0,
	αe = 1.0, δe = 0.2, βe = 1.0, αa = 0.8, δa = 0.1, βa = 0.9,
	save_pdf = false, pdf_dir = "data/Duffing")
	return DuffingParams(Adaptive, Robust, δt, N, K, B, A, Λ, w, K_VSSM, ω, Amp,
		αe, δe, βe, αa, δa, βa, save_pdf, pdf_dir)
end

# -----------------------------
# Helper model functions
# -----------------------------
Exact(p::DuffingParams, q, q_p, u) = p.αe * q + p.δe * q_p + p.βe * q^3 + u
Approx(p::DuffingParams, q, q_p, q_Np) = q_Np - p.αa * q - p.δa * q_p - p.βa * q^3
ErrorMetrics(p::DuffingParams, h_int, h, h_p) = p.Λ^2 * h_int + 2 * p.Λ * h + h_p
KinBlock(p::DuffingParams, S, h, h_p, qN_pp) = p.K_VSSM * tanh(S / p.w) + p.Λ^2 * h + 2 * p.Λ * h_p + qN_pp
G(p::DuffingParams, past_input, past_response, xDnow) = (p.K + past_input) * (1 + p.B * tanh(p.A * (past_response - xDnow))) - p.K

# nominal trajectory at time t
function nominalTraj(p::DuffingParams, t)
	qN = p.Amp * sin(p.ω * t)
	q_pN = p.ω * p.Amp * cos(p.ω * t)
	q_ppN = -p.Amp * p.ω^2 * sin(p.ω * t)
	return qN, q_pN, q_ppN
end

# -----------------------------
# Core simulator (type-stable, no globals)
# -----------------------------
function simulate_duffing(p::DuffingParams, q0::Float64, q_p0::Float64; do_plot::Bool = true)
	δt = p.δt
	N = p.N

	LONG = N - 1

	# preallocate arrays (Float64)
	time_mem = zeros(Float64, N)
	q_mem = zeros(Float64, N)
	q_p_mem = zeros(Float64, N)
	q_pp_mem = zeros(Float64, N)
	qN_mem = zeros(Float64, N)
	qN_p_mem = zeros(Float64, N)
	qN_pp_mem = zeros(Float64, N)

	u_mem = zeros(Float64, N)
	qDes_pp_mem = zeros(Float64, N)
	qDef_pp_mem = zeros(Float64, N)

	# initial conditions
	q_mem[1] = q0
	q_p_mem[1] = q_p0
	q_pp_mem[1] = -p.Amp * p.ω^2 * sin(p.ω * δt)
	max_index = 1
	h_int = 0.0
	l = LONG - 1
	for t in 1:l
		time_mem[t] = (t - 1) * δt

		qN_mem[t], qN_p_mem[t], qN_pp_mem[t] = nominalTraj(p, time_mem[t])

		h = qN_mem[t] - q_mem[t]
		h_p = qN_p_mem[t] - q_p_mem[t]

		h_max = qN_mem[max_index] - q_mem[max_index]
		if (h*h) < (h_max*h_max)
			max_index = t
		end

		if p.Robust
			S = ErrorMetrics(p, h_int, h, h_p)
			qDes_pp_mem[t] = KinBlock(p, S, h, h_p, qN_pp_mem[t])
		else
			qDes_pp_mem[t] = qN_pp_mem[t] + p.Λ^3 * h_int + 3 * p.Λ^2 * h + 3 * p.Λ * h_p
		end

		if p.Adaptive && t > 3
			qDef_pp_mem[t] = G(p, qDef_pp_mem[t-1], q_pp_mem[t-1], qDes_pp_mem[t])
		else
			qDef_pp_mem[t] = qDes_pp_mem[t]
		end

		u_mem[t] = Approx(p, q_mem[t], q_p_mem[t], qDef_pp_mem[t])

		q_pp_mem[t] = Exact(p, q_mem[t], q_p_mem[t], u_mem[t])

		# Integrate with Euler
		q_p_mem[t+1] = q_p_mem[t] + δt * q_pp_mem[t]
		q_mem[t+1] = q_mem[t] + δt * q_p_mem[t]

		h_int += δt * h
	end

	# # compute l for plotting bounds safely
	# l = LONG - 1

	# compute final max error value
	max_error = qN_mem[max_index] - q_mem[max_index]

	if do_plot
		# Trajectory
		figure("Trajectory_tracking")
		clf()
		grid(true)
		title("Pályakövetés az idő függvényében")
		xlabel("Idő [s]")
		ylabel("Pozíció [m]")
		plot(time_mem[1:l], qN_mem[1:l], label = "Nominális")
		plot(time_mem[1:l], q_mem[1:l], linestyle = "--", label = "Megvalósult")
		legend(loc = 1)

		# Acceleration
		figure("Acceleration")
		clf()
		grid(true)
		title("Accelerations vs Time")
		xlabel("Time [s]")
		ylabel("Acceleration [m/s²]")
		plot(time_mem[1:l], qN_pp_mem[1:l], label = "Névleges")
		plot(time_mem[1:l], q_pp_mem[1:l], linestyle = "--", label = "Realizált")
		plot(time_mem[1:l], qDes_pp_mem[1:l], linestyle = "-.", label = "Desired")
		legend(loc = 1)

		# Tracking error
		figure("Tracking_Error")
		clf()
		grid(true)
		title("Követési hiba az idő függvényében")
		xlabel("Idő [s]")
		ylabel("Követési hiba [m]")
		plot(time_mem[1:l], qN_mem[1:l] .- q_mem[1:l], label = "Hiba")
		legend(loc = 1)

		# Phase space
		figure("Phase_Space")
		clf()
		grid(true)
		title("Fázistér")
		xlabel("Pozíció [m]")
		ylabel("Sebesség [m/s]")
		plot(qN_mem[1:l], qN_p_mem[1:l], label = "Nominális")
		plot(q_mem[1:l], q_p_mem[1:l], linestyle = "--", label = "Megvalósult")
		legend(loc = 1)

		if p.save_pdf
			# combine and save a timestamped PDF
			try
				ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
				fname = joinpath(p.pdf_dir, "duffing_$(ts).pdf")
				# ensure directory
				isdir(p.pdf_dir) || mkpath(p.pdf_dir)
				savefig(fname)
			catch e
				@warn "Could not save PDF: $e"
			end
		end

		# show plots
		show()
	end

	return max_error
end

# Convenience wrappers to preserve original API
function duffing_single_run(q0::Real, q_p0::Real; do_plot::Bool = true)
	p = DuffingParams()
	return simulate_duffing(p, float(q0), float(q_p0); do_plot = do_plot)
end

function duffing_single_run(p::DuffingParams, q0::Real, q_p0::Real; do_plot::Bool = true)
	return simulate_duffing(p, float(q0), float(q_p0); do_plot = do_plot)
end

# Grid search helper
function grid_search(p::DuffingParams, q0_range::AbstractRange, q_p0_range::AbstractRange; do_plot::Bool = false)
	nx = length(q0_range)
	ny = length(q_p0_range)
	errors = zeros(Float64, nx, ny)
	for (i, q0) in enumerate(q0_range)
		for (j, qp0) in enumerate(q_p0_range)
			errors[i, j] = simulate_duffing(p, float(q0), float(qp0); do_plot = do_plot)
		end
	end
	# save CSV for later inspection
	ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
	fname = "duffing_grid_$(ts).csv"
	writedlm(fname, errors, ',')
	return (errors, fname)
end

end # module
