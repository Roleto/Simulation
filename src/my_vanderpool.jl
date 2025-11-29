module VanDerPolModule

using PyPlot
using PDFmerger
using Dates

export VanDerPolParams, simulate_vanderpol, vanderpol_single_run

#######################
# Paraméter struktúra #
#######################

mutable struct VanDerPolParams
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
	K_VSSM::Float64  # bent hagyom, ha később használni szeretnéd

	# Nominal trajectory
	ω::Float64
	Amp::Float64

	# Exact model params
	mₑ::Float64
	μₑ::Float64
	ωₑ::Float64
	αₑ::Float64
	λₑ::Float64

	# Approx model params
	mₐ::Float64
	μₐ::Float64
	ωₐ::Float64
	αₐ::Float64
	λₐ::Float64
	ll::Float64

	# Flags
	Adaptive::Bool
	Robust::Bool

	# Output
	save_pdf::Bool
	pdf_dir::String
end

function VanDerPolParams(;
	δt = 1e-3,
	N   = Int(2e4),

	K = 1e5,
	B = -1.0,
	A = 1e-5,

	Λ = 1.0,
	w = 1.0,
	K_VSSM = 500.0,

	ω = 0.5,
	Amp = 2.0,

	mₑ = 1.0,
	μₑ = 0.4,
	ωₑ = 0.46,
	αₑ = 1.0,
	λₑ = 0.1,

	mₐ = 0.8,
	μₐ = 0.5,
	ωₐ = 0.42,
	αₐ = 0.9,
	λₐ = 0.09,
	ll = 0.015,

	Adaptive = true,
	Robust = true,

	save_pdf = false,
	pdf_dir = "data/VanDerPol",
)
	return VanDerPolParams(
		δt, N,
		K, B, A,
		Λ, w, K_VSSM,
		ω, Amp,
		mₑ, μₑ, ωₑ, αₑ, λₑ,
		mₐ, μₐ, ωₐ, αₐ, λₐ, ll,
		Adaptive, Robust,
		save_pdf, pdf_dir,
	)
end

###################
# Segédfüggvények #
###################

ErrorMetrics(p::VanDerPolParams, h_int, h, h_p) =
	p.Λ^2 * h_int + 2 * p.Λ * h + h_p

KinBlock(p::VanDerPolParams, S, h, h_p, qN_pp) =
	p.K * tanh(S / p.w) + p.Λ^2 * h + 2 * p.Λ * h_p + qN_pp

function Exact(p::VanDerPolParams, u, q, q_p)
	return (u + p.μₑ * (1 - q^2) * q_p - p.ωₑ^2 * q - p.αₑ * q^3 - p.λₑ * q^5) / p.mₑ
end

function Approx(p::VanDerPolParams, q_pp)
	u = p.mₐ * q_pp
	return u * p.ll
end

function G(p::VanDerPolParams, past_input, past_response, xDnow)
	(p.K + past_input) * (1 + p.B * tanh(p.A * (past_response - xDnow))) - p.K
end

function nominalTraj(p::VanDerPolParams, t_idx::Int)
	τ = t_idx * p.δt
	qN = p.Amp * sin(p.ω * τ)
	q_pN = p.ω * p.Amp * cos(p.ω * τ)
	q_ppN = -p.Amp * p.ω^2 * sin(p.ω * τ)
	return qN, q_pN, q_ppN
end

###########################
# Fő szimulációs függvény #
###########################

"""
simulate_vanderpol(p, q0, q_p0; q_pp0=nothing, do_plot=true)

- p     : VanDerPolParams
- q0    : kezdő pozíció
- q_p0  : kezdő sebesség
- q_pp0 : opcionális kezdő gyorsulás (ha nothing, akkor nominálisból számoljuk)
- do_plot: ha true, kirajzol + (ha save_pdf) PDF-et is ment

Visszatérés: egyetlen skalár hiba (Duffing-stílus, 1D eset)
"""
function simulate_vanderpol(p::VanDerPolParams,
	q0::Real,
	q_p0::Real;
	q_pp0::Union{Nothing, Real} = nothing,
	do_plot::Bool = true)

	δt = p.δt
	N = p.N
	l = N - 1

	# Memóriák
	time_mem = zeros(Float64, N)
	q_mem    = zeros(Float64, N)
	q_p_mem  = zeros(Float64, N)
	q_pp_mem = zeros(Float64, N)

	qN_mem    = zeros(Float64, N)
	qN_p_mem  = zeros(Float64, N)
	qN_pp_mem = zeros(Float64, N)

	u_mem = zeros(Float64, N)
	qDes_pp_mem = zeros(Float64, N)
	qDef_pp = zeros(Float64, N)

	# Kezdeti feltételek
	q_mem[1]   = float(q0)
	q_p_mem[1] = float(q_p0)

	if q_pp0 === nothing
		q_pp_mem[1] = -p.Amp * p.ω^2 * sin(p.ω * p.δt)
	else
		q_pp_mem[1] = float(q_pp0)
	end

	# hiba integrál
	h_int = 0.0

	# max_index inicializálás
	max_index = 1

	for t in 1:l
		time_mem[t] = t * δt

		# nominális pálya
		qN_mem[t], qN_p_mem[t], qN_pp_mem[t] = nominalTraj(p, t)

		# hiba
		h   = qN_mem[t] - q_mem[t]
		h_p = qN_p_mem[t] - q_p_mem[t]

		# kinematikus blokk
		if p.Robust
			S = ErrorMetrics(p, h_int, h, h_p)
			qDes_pp_mem[t] = KinBlock(p, S, h, h_p, qN_pp_mem[t])
		else
			qDes_pp_mem[t] = qN_pp_mem[t] + p.Λ^3 * h_int + 3 * p.Λ^2 * h + 3 * p.Λ * h_p
		end

		# deformáció (adaptív G blokk)
		if p.Adaptive && t > 10
			qDef_pp[t] = G(p, qDef_pp[t-1], q_pp_mem[t-1], qDes_pp_mem[t])
		else
			qDef_pp[t] = qDes_pp_mem[t]
		end

		# irányítójel
		u_mem[t] = Approx(p, qDef_pp[t])

		# rendszer válasza
		q_pp_mem[t] = Exact(p, u_mem[t], q_mem[t], q_p_mem[t])

		# Euler integrálás
		q_p_mem[t+1] = q_p_mem[t] + δt * q_pp_mem[t]
		q_mem[t+1]   = q_mem[t] + δt * q_p_mem[t]

		# hiba integrál frissítése
		h_int += δt * h

		h_now = qN_mem[t] - q_mem[t]
		h_max_now = qN_mem[max_index] - q_mem[max_index]
		if (h_now * h_now) < (h_max_now * h_max_now)
			max_index = t
		end
	end

	#Plotting
	if do_plot
		# Trajectory tracking
		figure("Trajectory_tracking")
		grid(true)
		title("Pályakövetés az idő függvényében")
		xlabel("Idő [s]")
		ylabel("Pozíció [m]")
		plot(time_mem[1:l], qN_mem[1:l], color = "red", label = "Nominális")
		plot(time_mem[1:l], q_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
		legend(loc = 1, borderaxespad = 0)

		if p.save_pdf
			savefig("temp_vander.pdf")
			append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)
		end

		# Velocities
		figure("Velocities")
		grid(true)
		title("Sebesség az idő függvényében")
		xlabel("Idő [s]")
		ylabel("Sebesség [m/s]")
		plot(time_mem[1:l], qN_p_mem[1:l], color = "red", label = "Nominális")
		plot(time_mem[1:l], q_p_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
		legend(loc = 1, borderaxespad = 0)

		if p.save_pdf
			savefig("temp_vander.pdf")
			append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)
		end

		# Accelerations
		figure("Accelerations")
		grid(true)
		title("Gyorsulás az idő függvényében")
		xlabel("Idő [s]")
		ylabel("Gyorsulás [m/s²]")
		plot(time_mem[1:l], qN_pp_mem[1:l], color = "red", label = "Nominális")
		plot(time_mem[1:l], q_pp_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
		legend(loc = 1, borderaxespad = 0)

		if p.save_pdf
			savefig("temp_vander.pdf")
			append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)
		end

		# Tracking error
		figure("Tracking_Error")
		title("Követési hiba az idő függvényében")
		grid(true)
		xlabel("Idő [s]")
		ylabel("Követési hiba [m]")
		plot(time_mem[1:l], qN_mem[1:l] .- q_mem[1:l], color = "red")

		if p.save_pdf
			savefig("temp_vander.pdf")
			append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)
		end

		# Control signal
		figure("Control_Signal")
		title("Irányítójel az idő függvényében")
		grid(true)
		xlabel("Idő [s]")
		ylabel("Irányítójel [N]")
		plot(time_mem[1:l], u_mem[1:l], color = "red")

		if p.save_pdf
			savefig("temp_vander.pdf")
			append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)
		end

		# Phase space
		figure("Phase_Space")
		title("Fázistér")
		xlabel("Pozíció [m]")
		ylabel("Sebesség [m/s]")
		grid(true)
		plot(qN_p_mem[1:l], qN_mem[1:l], color = "red", label = "Nominális")
		plot(q_p_mem[1:l], q_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
		legend(loc = 1, borderaxespad = 0)

		if p.save_pdf
			savefig("temp_vander.pdf")
			append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)
			# timestampelt mentés
			try
				# timestamp
				ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

				# Ensure output dir exists
				isdir(p.pdf_dir) || mkpath(p.pdf_dir)

				# Path of the merged PDF in project root
				project_root = normpath(joinpath(@__DIR__, ".."))
				merged_pdf = joinpath(project_root, "allplots_vander.pdf")

				# Expected final file
				final_pdf = joinpath(p.pdf_dir, "vanderpol_$(ts).pdf")

				# MOVE the merged file → correct location!
				mv(merged_pdf, final_pdf; force = true)
			catch e
				@warn "Could not save PDF: $e"
			end
		end

		show()
	end

	return abs(qN_mem[max_index] - q_mem[max_index])
end

end # module
