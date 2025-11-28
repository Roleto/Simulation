module RosslerModule

using LinearAlgebra
using PyPlot
using PDFmerger
using Dates

export RosslerParams, simulate_rossler, rossler_single_run

#######################
# Paraméter struktúra #
#######################

mutable struct RosslerParams
	δt::Float64
	N::Int

	# Control parameters
	K::Float64
	B::Float64
	A::Float64

	# Kinematic block
	Λ::Float64
	K_VSSM::Float64
	w::Float64

	# Nominal trajectory parameters
	A1::Float64
	ω1::Float64
	A2::Float64
	ω2::Float64
	A3::Float64
	ω3::Float64

	# Exact model params
	ae::Float64
	be::Float64
	ce::Float64

	# Approx model params
	aa::Float64
	ba::Float64
	ca::Float64

	Adaptive::Bool
	Robust::Bool

	save_pdf::Bool
	pdf_dir::String
end


function RosslerParams(;
	δt = 1e-3,
	N = Int(2e4),

	K = 1e2,
	B = -1.0,
	A = 1.97e-3,

	Λ = 1.0,
	K_VSSM = 50.0,
	w = 1.0,

	A1 = 2.0,
	ω1 = 0.5,
	A2 = 3.0,
	ω2 = 0.7,
	A3 = 1.0,
	ω3 = 1.0,

	ae = 0.01,
	be = 0.2,
	ce = 5.7,

	aa = 0.1,
	ba = 0.3,
	ca = 5.5,

	Adaptive = true,
	Robust = true,

	save_pdf = false,
	pdf_dir = "data/Rossler",
)
	return RosslerParams(
		δt, N,
		K, B, A,
		Λ, K_VSSM, w,
		A1, ω1, A2, ω2, A3, ω3,
		ae, be, ce,
		aa, ba, ca,
		Adaptive, Robust,
		save_pdf, pdf_dir,
	)
end


###################
# Segédfüggvények #
###################

ErrorMetric(hint, h, Λ) = Λ .* hint .+ h

function G_MIMO(past_input, past_response, desired, err_limit, p::RosslerParams)
	Amatr_h = past_response - desired
	error_norm = norm(Amatr_h)

	if error_norm > err_limit
		e_direction = Amatr_h / error_norm
		B_factor = p.B * tanh(p.A * error_norm)
		G = (1 + B_factor) * past_input + B_factor * p.K * e_direction
	else
		G = past_input
	end
	return G
end

###########################
# Fő szimulációs függvény #
###########################

function simulate_rossler(p::RosslerParams,
	q0::NTuple{3, Real},
	q_p0::NTuple{3, Real};
	q_pp0 = nothing,
	do_plot = true)

	δt = p.δt
	N = p.N
	l = N - 1

	# Allocation
	time_mem = zeros(N)

	x = zeros(N);
	y = zeros(N);
	z = zeros(N)
	x_p = zeros(N);
	y_p = zeros(N);
	z_p = zeros(N)

	xN = zeros(N);
	xN_p = zeros(N)
	yN = zeros(N);
	yN_p = zeros(N)
	zN = zeros(N);
	zN_p = zeros(N)

	xDes_p = zeros(N);
	yDes_p = zeros(N);
	zDes_p = zeros(N)
	xDef_p = zeros(N);
	yDef_p = zeros(N);
	zDef_p = zeros(N)

	S_x = zeros(N);
	S_y = zeros(N);
	S_z = zeros(N)

	past_input = zeros(3)
	past_response = zeros(3)
	past_responses = zeros(N, 3)

	hint_x = 0.0
	hint_y = 0.0
	hint_z = 0.0
	error_limit = 1e-3

	# Initial conditions
	x[1]   = q0[1];
	y[1]   = q0[2];
	z[1]   = q0[3]
	x_p[1] = q_p0[1];
	y_p[1] = q_p0[2];
	z_p[1] = q_p0[3]

	max_index = 1

	for t in 1:l
		time_mem[t] = (t - 1) * δt

		# Nominal traj
		xN[t] = p.A1 * sin(p.ω1 * time_mem[t])
		xN_p[t] = p.A1 * p.ω1 * cos(p.ω1 * time_mem[t])

		yN[t] = p.A2 * sin(p.ω2 * time_mem[t])
		yN_p[t] = p.A2 * p.ω2 * cos(p.ω2 * time_mem[t])

		zN[t] = p.A3 * sin(p.ω3 * time_mem[t])
		zN_p[t] = p.A3 * p.ω3 * cos(p.ω3 * time_mem[t])

		# Errors
		h_x = xN[t] - x[t]
		h_y = yN[t] - y[t]
		h_z = zN[t] - z[t]

		# Robust
		if p.Robust
			S = ErrorMetric([hint_x, hint_y, hint_z], [h_x, h_y, h_z], p.Λ)
			S_x[t] = S[1];
			S_y[t] = S[2];
			S_z[t] = S[3]

			xDes_p[t] = xN_p[t] + p.Λ*h_x + p.K_VSSM*tanh(S_x[t]/p.w)
			yDes_p[t] = yN_p[t] + p.Λ*h_y + p.K_VSSM*tanh(S_y[t]/p.w)
			zDes_p[t] = zN_p[t] + p.Λ*h_z + p.K_VSSM*tanh(S_z[t]/p.w)

			desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
		else
			xDes_p[t] = p.Λ^2*hint_x + 2*p.Λ*h_x + xN_p[t]
			yDes_p[t] = p.Λ^2*hint_y + 2*p.Λ*h_y + yN_p[t]
			zDes_p[t] = p.Λ^2*hint_z + 2*p.Λ*h_z + zN_p[t]

			desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
		end

		# Adaptive deformation
		if p.Adaptive && t > 3
			past_input = G_MIMO(past_input, past_response, desired, error_limit, p)
		else
			past_input = desired
		end

		xDef_p[t] = past_input[1]
		yDef_p[t] = past_input[2]
		zDef_p[t] = past_input[3]

		# Control signal
		u_x = xDef_p[t] + y[t] + z[t]
		u_y = yDef_p[t] - x[t] - p.aa*y[t]
		u_z = zDef_p[t] - p.ba - z[t]*(x[t]-p.ca)

		# Model dynamics
		x_p[t] = -y[t] - z[t] + u_x
		y_p[t] = x[t] + p.ae*y[t] + u_y
		z_p[t] = p.be + z[t]*(x[t]-p.ce) + u_z

		past_response .= [x_p[t], y_p[t], z_p[t]]
		past_responses[t, :] .= past_response

		# Integrate
		x[t+1] = x[t] + δt*x_p[t]
		y[t+1] = y[t] + δt*y_p[t]
		z[t+1] = z[t] + δt*z_p[t]

		hint_x += δt*h_x
		hint_y += δt*h_y
		hint_z += δt*h_z

		h2 = h_x*h_x + h_y*h_y + h_z*h_z

		h_max_x = xN[max_index] - x[max_index]
		h_max_y = yN[max_index] - y[max_index]
		h_max_z = zN[max_index] - z[max_index]
		h2_max = h_max_x*h_max_x + h_max_y*h_max_y + h_max_z*h_max_z

		if h2 > h2_max
			max_index = t
		end
	end

	hx = xN[max_index] - x[max_index]
	hy = yN[max_index] - y[max_index]
	hz = zN[max_index] - z[max_index]

	# Plot
	if do_plot
		# Nominal–Realized trajectories
		fig = figure("nominal_realized_rossler")
		title("Nominal vs Realized — Rössler")
		subplot(311);
		grid(true)
		ylabel("X")
		plot(time_mem[1:l], xN[1:l], "r")
		plot(time_mem[1:l], x[1:l], "g--")

		subplot(312);
		grid(true)
		ylabel("Y")
		plot(time_mem[1:l], yN[1:l], "r")
		plot(time_mem[1:l], y[1:l], "g--")

		subplot(313);
		grid(true)
		ylabel("Z");
		xlabel("t")
		plot(time_mem[1:l], zN[1:l], "r")
		plot(time_mem[1:l], z[1:l], "g--")

		tight_layout()
		if p.save_pdf
			savefig("allplots_rossler.pdf")
		end


		# TRACKING ERROR PLOT
		fig = figure("tracking_error_rossler")
		subplot(311);
		grid(true)
		title("Tracking Errors")
		ylabel("Xerr")
		plot(time_mem[1:l], xN[1:l]-x[1:l], "r")

		subplot(312);
		grid(true)
		ylabel("Yerr")
		plot(time_mem[1:l], yN[1:l]-y[1:l], "r")

		subplot(313);
		grid(true)
		ylabel("Zerr")
		plot(time_mem[1:l], zN[1:l]-z[1:l], "r")
		xlabel("t")

		tight_layout()
		if p.save_pdf
			savefig("temp_rossler.pdf")
			append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)
		end


		# CONTROL SIGNALS
		fig = figure("control_signal_rossler")
		grid(true)
		title("Control Signals")
		xlabel("t");
		ylabel("u")
		plot(time_mem[1:l], xDef_p[1:l], "r", label = "u_x")
		plot(time_mem[1:l], yDef_p[1:l], "g--", label = "u_y")
		plot(time_mem[1:l], zDef_p[1:l], "b", label = "u_z")
		legend()
		tight_layout()

		if p.save_pdf
			savefig("temp_rossler.pdf")
			append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)
		end


		# PHASE PLOTS
		fig = figure("phase_x_rossler")
		grid(true)
		title("Phase X")
		xlabel("x");
		ylabel("x_p")
		plot(xN[1:l], xN_p[1:l], "r")
		plot(x[1:l], x_p[1:l], "g--")
		if p.save_pdf
			savefig("temp_rossler.pdf")
			append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)
		end

		fig = figure("phase_y_rossler")
		grid(true)
		title("Phase Y")
		xlabel("y");
		ylabel("y_p")
		plot(yN[1:l], yN_p[1:l], "r")
		plot(y[1:l], y_p[1:l], "g--")
		if p.save_pdf
			savefig("temp_rossler.pdf")
			append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)
		end

		fig = figure("phase_z_rossler")
		grid(true)
		title("Phase Z")
		xlabel("z");
		ylabel("z_p")
		plot(zN[1:l], zN_p[1:l], "r")
		plot(z[1:l], z_p[1:l], "g--")
		if p.save_pdf
			savefig("temp_rossler.pdf")
			append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)
		end


		if p.save_pdf
			try
				ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
				isdir(p.pdf_dir) || mkpath(p.pdf_dir)
				fname = joinpath(p.pdf_dir, "rossler_$(ts).pdf")
				savefig(fname)
			catch e
				@warn "Could not save PDF: $e"
			end
		end

		show()
	end

	return sqrt(hx*hx + hy*hy + hz*hz)
end

end # module
