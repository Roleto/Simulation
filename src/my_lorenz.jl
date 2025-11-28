module LorenzModule

using LinearAlgebra
using PyPlot
using PDFmerger
using Dates

export LorenzParams, simulate_lorenz, lorenz_single_run

#######################
# Paraméter struktúra #
#######################
mutable struct LorenzParams
	# Idő
	δt::Float64
	N::Int

	# Kontroll paraméterek
	K::Float64
	B::Float64
	A::Float64

	# Kinematikus blokk
	Λ::Float64
	K_VSSM::Float64
	w::Float64

	# Nominális pálya paraméterek
	A1::Float64
	A2::Float64
	A3::Float64
	ω1::Float64
	ω2::Float64
	ω3::Float64

	# Pontos modell paraméterei
	βe::Float64
	σe::Float64
	ρe::Float64

	# Közelítő modell paraméterei
	βa::Float64
	σa::Float64
	ρa::Float64

	# Logikai flag-ek
	Adaptive::Bool
	Robust::Bool

	# Output options
	save_pdf::Bool
	pdf_dir::String
end

function LorenzParams(;
	δt::Float64 = 1e-2,
	N::Int = Int(2e4),

	K::Float64 = 1e2,
	B::Float64 = -1.0,
	A::Float64 = 1.97e-3,

	Λ::Float64 = 1.0,
	K_VSSM::Float64 = 50.0,
	w::Float64 = 1.0,

	A1::Float64 = 2.0,
	ω1::Float64 = 0.5,
	A2::Float64 = 3.0,
	ω2::Float64 = 0.7,
	A3::Float64 = 1.0,
	ω3::Float64 = 1.0,

	βe::Float64 = 8/3,
	σe::Float64 = 5.0,
	ρe::Float64 = 40.0,

	βa::Float64 = 7/3,
	σa::Float64 = 4.0,
	ρa::Float64 = 36.0,

	Adaptive::Bool = true,
	Robust::Bool = true,

	save_pdf::Bool = false,
	pdf_dir::String = "data/Lorenz",
)
	return LorenzParams(
		δt, N,
		K, B, A,
		Λ, K_VSSM, w,
		A1, A2, A3, ω1, ω2, ω3,
		βe, σe, ρe,
		βa, σa, ρa,
		Adaptive, Robust,
		save_pdf, pdf_dir,
	)
end

###################
# Segédfüggvények #
###################

"ErrorMetric az eredeti formulával: S, Ṡ"
function ErrorMetric(p::LorenzParams, hint::NTuple{3, Float64},
	h::NTuple{3, Float64}, h_p::NTuple{3, Float64})
	Λ = p.Λ
	S = (
		Λ * hint[1] + h[1],
		Λ * hint[2] + h[2],
		Λ * hint[3] + h[3],
	)
	S_p = (
		Λ * h[1] + h_p[1],
		Λ * h[2] + h_p[2],
		Λ * h[3] + h_p[3],
	)
	return S, S_p
end

"""
Multiváltozós adaptív blokk (G_MIMO) – az eredeti kódból átemelve.

r_prev       – előző bemenet (v)
r_prev_last  – r_prev(t) - r_prev(t-1) szerepben van
f_prev       – előző válasz (sebességek)
des_now      – aktuális kívánt sebesség vektor
err_limit    – hiba küszöb
"""
function G_MIMO(p::LorenzParams,
	r_prev::Vector{Float64},
	r_prev_last::Vector{Float64},
	f_prev::Vector{Float64},
	des_now::Vector{Float64},
	err_limit::Float64)

	K = p.K
	B = p.B
	A = p.A

	Amatr_h = (f_prev .- des_now)
	error_norm = norm(Amatr_h)

	if error_norm > err_limit
		e_direction = Amatr_h / error_norm
		B_factor = B * tanh(A * error_norm)
		G = (1 + B_factor) .* r_prev .+ B_factor * K .* e_direction
	else
		G = copy(r_prev)
	end

	G_p = G .- r_prev_last
	return G, G_p
end

###########################
# Fő szimulációs függvény #
###########################

"""
simulate_lorenz(p, q0, q_p0; q_pp0 = nothing, do_plot = true)

- p      : LorenzParams
- q0     : (x0, y0, z0)
- q_p0   : (ẋ0, ẏ0, ż0)
- q_pp0  : opcionális (ẍ0, ÿ0, z̈0), ha nem adod meg, az eredeti nominális gyorsulást használjuk
- do_plot: ha true, az összes eredeti ábrát kirajzoljuk és PDF-be mentjük

Visszatérés: max tracking error (pozíciós hiba max. normája az időben)
"""
function simulate_lorenz(p::LorenzParams,
	q0::NTuple{3, Real},
	q_p0::NTuple{3, Real};
	q_pp0::Union{Nothing, NTuple{3, Real}} = nothing,
	do_plot::Bool = true)

	δt  = p.δt
	LONG = p.N
	l    = LONG - 1

	time_mem = zeros(Float64, LONG)

	xN    = zeros(Float64, LONG);
	xN_p  = zeros(Float64, LONG);
	xN_pp = zeros(Float64, LONG)
	yN    = zeros(Float64, LONG);
	yN_p  = zeros(Float64, LONG);
	yN_pp = zeros(Float64, LONG)
	zN    = zeros(Float64, LONG);
	zN_p  = zeros(Float64, LONG);
	zN_pp = zeros(Float64, LONG)

	xDes_p = zeros(Float64, LONG)
	yDes_p = zeros(Float64, LONG)
	zDes_p = zeros(Float64, LONG)

	xDef_p = zeros(Float64, LONG)
	yDef_p = zeros(Float64, LONG)
	zDef_p = zeros(Float64, LONG)

	xDef_pp = zeros(Float64, LONG)
	yDef_pp = zeros(Float64, LONG)
	zDef_pp = zeros(Float64, LONG)

	u_x = zeros(Float64, LONG)
	u_y = zeros(Float64, LONG)
	u_z = zeros(Float64, LONG)

	x = zeros(Float64, LONG)
	y = zeros(Float64, LONG)
	z = zeros(Float64, LONG)

	x_p = zeros(Float64, LONG)
	y_p = zeros(Float64, LONG)
	z_p = zeros(Float64, LONG)

	x_pp = zeros(Float64, LONG)
	y_pp = zeros(Float64, LONG)
	z_pp = zeros(Float64, LONG)

	S_x   = zeros(Float64, LONG)
	S_p_x = zeros(Float64, LONG)
	S_y   = zeros(Float64, LONG)
	S_p_y = zeros(Float64, LONG)
	S_z   = zeros(Float64, LONG)
	S_p_z = zeros(Float64, LONG)

	past_input = zeros(Float64, 3)
	past_input_p = zeros(Float64, 3)
	past_inputs = zeros(Float64, LONG, 3)
	past_inputs_p = zeros(Float64, LONG, 3)
	past_response = zeros(Float64, 3)
	past_responses = zeros(Float64, LONG, 3)
	accelerations = zeros(Float64, LONG, 3)

	x[1]   = float(q0[1])
	y[1]   = float(q0[2])
	z[1]   = float(q0[3])
	x_p[1] = float(q_p0[1])
	y_p[1] = float(q_p0[2])
	z_p[1] = float(q_p0[3])

	if q_pp0 === nothing
		x_pp[1] = -p.A1 * p.ω1^2 * sin(p.ω1 * δt)
		y_pp[1] = -p.A2 * p.ω2^2 * sin(p.ω2 * δt)
		z_pp[1] = -p.A3 * p.ω3^2 * sin(p.ω3 * δt)
	else
		x_pp[1] = float(q_pp0[1])
		y_pp[1] = float(q_pp0[2])
		z_pp[1] = float(q_pp0[3])
	end

	# Integrálok a hibametrikához
	hint_x = 0.0
	hint_y = 0.0
	hint_z = 0.0

	error_limit = 1e-3  # Hiba limit az adaptív blokkhoz
	max_index = 1

	for t ∈ 1:l
		time_mem[t] = δt * t

		# Nominális pálya 
		xN[t] = p.A1 * sin(p.ω1 * time_mem[t])
		xN_p[t] = p.A1 * p.ω1 * cos(p.ω1 * time_mem[t])
		xN_pp[t] = -p.A1 * p.ω1^2 * sin(p.ω1 * time_mem[t])

		yN[t] = p.A2 * sin(p.ω2 * time_mem[t])
		yN_p[t] = p.A2 * p.ω2 * cos(p.ω2 * time_mem[t])
		yN_pp[t] = -p.A2 * p.ω2^2 * sin(p.ω2 * time_mem[t])

		zN[t] = p.A3 * sin(p.ω3 * time_mem[t])
		zN_p[t] = p.A3 * p.ω3 * cos(p.ω3 * time_mem[t])
		zN_pp[t] = -p.A3 * p.ω3^2 * sin(p.ω3 * time_mem[t])

		# Hibák
		h_x = xN[t] - x[t]
		h_y = yN[t] - y[t]
		h_z = zN[t] - z[t]

		h2 = h_x*h_x + h_y*h_y + h_z*h_z

		h_max_x = xN[max_index] - x[max_index]
		h_max_y = yN[max_index] - y[max_index]
		h_max_z = zN[max_index] - z[max_index]

		h2_max = h_max_x*h_max_x + h_max_y*h_max_y + h_max_z*h_max_z

		if h2 > h2_max
			max_index = t
		end


		h_p_x = xN_p[t] - x_p[t]
		h_p_y = yN_p[t] - y_p[t]
		h_p_z = zN_p[t] - z_p[t]

		if p.Robust
			S, S_p = ErrorMetric(p,
				(hint_x, hint_y, hint_z),
				(h_x, h_y, h_z),
				(h_p_x, h_p_y, h_p_z),
			)

			S_x[t]   = S[1];
			S_y[t]   = S[2];
			S_z[t]   = S[3]
			S_p_x[t] = S_p[1];
			S_p_y[t] = S_p[2];
			S_p_z[t] = S_p[3]

			# Kinematikus blokk (VSSM)
			xDes_p[t] = xN_p[t] + p.Λ * h_x + p.K_VSSM * tanh(S_x[t] / p.w)
			yDes_p[t] = yN_p[t] + p.Λ * h_y + p.K_VSSM * tanh(S_y[t] / p.w)
			zDes_p[t] = zN_p[t] + p.Λ * h_z + p.K_VSSM * tanh(S_z[t] / p.w)

			desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
		else
			xDes_p[t] = p.Λ^2 * hint_x + 2 * p.Λ * h_x + xN_p[t]
			yDes_p[t] = p.Λ^2 * hint_y + 2 * p.Λ * h_y + yN_p[t]
			zDes_p[t] = p.Λ^2 * hint_z + 2 * p.Λ * h_z + zN_p[t]

			desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
		end

		# Deformáció (ACC blokk)
		if p.Adaptive && t > 4
			past_input, past_input_p = G_MIMO(p, past_input, past_input_p, past_response, desired, error_limit)
		else
			past_input_p .= past_input
			past_input .= desired
		end

		past_inputs[t, :] .= past_input
		past_inputs_p[t, :] .= past_input_p

		xDef_p[t] = past_input[1];
		yDef_p[t] = past_input[2];
		zDef_p[t] = past_input[3]
		xDef_pp[t] = past_input_p[1];
		yDef_pp[t] = past_input_p[2];
		zDef_pp[t] = past_input_p[3]

		# Irányítójelek (u)
		u_x[t] = xDef_p[t] - p.σa * (y[t] - x[t])
		u_y[t] = yDef_p[t] - x[t] * (p.ρa - z[t]) + y[t]
		u_z[t] = zDef_p[t] - x[t] * y[t] + p.βa * z[t]

		# Rendszer
		x_p[t] = p.σe * (y[t] - x[t]) + u_x[t]
		y_p[t] = x[t] * (p.ρe - z[t]) - y[t] + u_y[t]
		z_p[t] = x[t] * y[t] - p.βe * z[t] + u_z[t]

		past_response .= [x_p[t], y_p[t], z_p[t]]
		past_responses[t, :] .= past_response

		# Gyorsulás
		x_pp[t] = (p.σe - p.σa) * (y_p[t] - x_p[t]) + xDef_pp[t]
		y_pp[t] = (p.ρe - p.ρa) * x_p[t] + yDef_pp[t]
		z_pp[t] = (p.βa - p.βe) * z_p[t] + zDef_pp[t]
		accelerations[t, :] .= [x_pp[t], y_pp[t], z_pp[t]]

		# Integrálás
		x[t+1] = x[t] + δt * x_p[t]
		y[t+1] = y[t] + δt * y_p[t]
		z[t+1] = z[t] + δt * z_p[t]

		hint_x += δt * h_x
		hint_y += δt * h_y
		hint_z += δt * h_z
	end

	# Plot
	if do_plot
		# Nominális és megvalósult pályák
		fig_caption = "nominal_realized_trajectories"
		fig = figure(fig_caption)
		title("Nominális és Megvalósult Pályák")
		subplot(311)
		ax1 = gca()
		grid(true)
		ylabel("X [m]")
		plot(time_mem[1:l], xN[1:l], color = "#D55E00", linewidth = 1.5, label = "Nominális", alpha = 0.8)
		plot(time_mem[1:l], x[1:l], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)

		subplot(312, sharex = ax1)
		ax2 = gca()
		grid(true)
		ylabel("Y [m]")
		plot(time_mem[1:l], yN[1:l], color = "#D55E00", linewidth = 1.5, label = "Nominális")
		plot(time_mem[1:l], y[1:l], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)

		subplot(313, sharex = ax2)
		ax3 = gca()
		grid(true)
		xlabel("Idő [s]")
		ylabel("Z [m]")
		plot(time_mem[1:l], zN[1:l], color = "#D55E00", linewidth = 1.5, label = "Nominális")
		plot(time_mem[1:l], z[1:l], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)

		handles, labels = ax1.get_legend_handles_labels()
		fig.legend(handles, labels, loc = "lower center", bbox_to_anchor = (0.5, -0.01), ncol = 2)
		subplots_adjust(hspace = 0.2, bottom = 0.2)
		tight_layout()
		fig[:canvas][:draw]()
		if (p.save_pdf)
			savefig("allplots_lorenz.pdf")
		end

		# 3D pályák
		fig = figure("Trajectory_tracking_3D")
		title("Nominális és Megvalósult Pályák 3D")
		grid(true)
		plot3D(xN[1:l], yN[1:l], zN[1:l], color = "red", label = "Nominális")
		plot3D(x[1:l], y[1:l], z[1:l], color = "green", label = "Megvalósult", linestyle = "--")
		legend(loc = 1, borderaxespad = 0)

		# Sebességek
		fig_caption = "velocities"
		fig = figure(fig_caption)
		grid(true)
		title("Nominális és Megvalósult Sebességek")

		subplot(311)
		ax1 = gca()
		grid(true)
		ylabel("X [m/s]")
		plot(time_mem[1:l], past_responses[1:l, 1], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
		plot(time_mem[1:l], xN_p[1:l], color = "red", linewidth = 2, label = "Nominális")

		subplot(312, sharex = ax1)
		ax2 = gca()
		grid(true)
		ylabel("Y [m/s]")
		plot(time_mem[1:l], yN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
		plot(time_mem[1:l], past_responses[1:l, 2], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")

		subplot(313, sharex = ax2)
		ax3 = gca()
		grid(true)
		xlabel("Idő [s]")
		ylabel("Z [m/s]")
		plot(time_mem[1:l], zN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
		plot(time_mem[1:l], past_responses[1:l, 3], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")

		tight_layout()
		subplots_adjust(hspace = 0.0)
		fig[:canvas][:draw]()
		if (p.save_pdf)
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
		end

		# Gyorsulások
		fig_caption = "accelerations"
		fig = figure(fig_caption)
		grid(true)
		title("Nominális és Megvalósult Gyorsulások")

		subplot(311)
		ax1 = gca()
		grid(true)
		ylabel("X [m/s²]")
		plot(time_mem[1:l], xN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
		plot(time_mem[1:l], accelerations[1:l, 1], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")

		subplot(312, sharex = ax1)
		ax2 = gca()
		grid(true)
		ylabel("Y [m/s²]")
		plot(time_mem[1:l], yN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
		plot(time_mem[1:l], accelerations[1:l, 2], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")

		subplot(313, sharex = ax2)
		ax3 = gca()
		grid(true)
		xlabel("Idő [s]")
		ylabel("Z [m/s²]")
		plot(time_mem[1:l], zN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
		plot(time_mem[1:l], accelerations[1:l, 3], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")

		tight_layout()
		subplots_adjust(hspace = 0.0)
		fig[:canvas][:draw]()
		if (p.save_pdf)
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
		end

		# Hibák
		fig_caption = "tracking_error"
		fig = figure(fig_caption)
		subplot(311)
		ax1 = gca()
		grid(true)
		title("Követési Hibák az Idő Függvényében")
		ylabel("X [m]")
		plot(time_mem[1:l], xN[1:l] .- x[1:l], color = "red", linewidth = 2)

		subplot(312, sharex = ax1)
		ax2 = gca()
		grid(true)
		ylabel("Y [m]")
		plot(time_mem[1:l], yN[1:l] .- y[1:l], color = "red", linewidth = 2)

		subplot(313, sharex = ax2)
		ax3 = gca()
		grid(true)
		xlabel("Idő [s]")
		ylabel("Z [m]")
		plot(time_mem[1:l], zN[1:l] .- z[1:l], color = "red", linewidth = 2)

		tight_layout()
		subplots_adjust(hspace = 0.190)
		fig[:canvas][:draw]()
		if (p.save_pdf)
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
		end

		# Irányítójelek
		fig_caption = "control_signal"
		figure(fig_caption)
		grid(true)
		title("Irányítójelek az Idő Függvényében")
		xlabel("Idő [s]")
		ylabel("Irányítójel [N]")
		plot(time_mem[1:l], u_x[1:l], color = "red", label = "\$u_x\$")
		plot(time_mem[1:l], u_y[1:l], color = "green", label = "\$u_y\$", linestyle = "--")
		plot(time_mem[1:l], u_z[1:l], color = "blue", label = "\$u_z\$")
		legend(loc = "lower left", fancybox = "True")
		tight_layout()
		gcf()[:canvas][:draw]()
		if (p.save_pdf)
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
		end
		# Fázisterek
		fig_caption = "phase_trajectories_x"
		figure(fig_caption)
		grid(true)
		title("Fázis tér X irányban")
		xlabel("Pozíció [m]")
		ylabel("Sebesség [m/s]")
		plot(x_p[1:l], x[1:l], color = "green", linewidth = 2.5, label = "Megvalósult", linestyle = "--")
		plot(xN_p[1:l], xN[1:l], color = "red", linewidth = 2, label = "Nominális")
		legend(loc = "lower left", fancybox = "True")
		tight_layout()
		if (p.save_pdf)
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
		end

		fig_caption = "phase_trajectories_y"
		figure(fig_caption)
		grid(true)
		title("Fázis tér Y irányban")
		xlabel("Pozíció [m]")
		ylabel("Sebesség [m/s]")
		plot(y_p[1:l], y[1:l], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
		plot(yN_p[1:l], yN[1:l], color = "red", linewidth = 2, label = "Nominális")
		legend(loc = "lower left", fancybox = "True")
		tight_layout()

		if (p.save_pdf)
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
		end

		fig_caption = "phase_trajectories_z"
		figure(fig_caption)
		grid(true)
		title("Fázis tér Z irányban")
		xlabel("Pozíció [m]")
		ylabel("Sebesség [m/s]")
		plot(z_p[1:l], z[1:l], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
		plot(zN_p[1:l], zN[1:l], color = "red", linewidth = 2, label = "Nominális")
		legend(loc = "lower left", fancybox = "True")
		tight_layout()

		if p.save_pdf
			savefig("temp_lorenz.pdf")
			append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

			# combine and save a timestamped PDF
			try
				ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
				fname = joinpath(p.pdf_dir, "lorenz_$(ts).pdf")
				# ensure directory
				isdir(p.pdf_dir) || mkpath(p.pdf_dir)
				savefig(fname)
			catch e
				@warn "Could not save PDF: $e"
			end
		end

		show()
	end

	h_x = xN[max_index] - x[max_index]
	h_y = yN[max_index] - y[max_index]
	h_z = zN[max_index] - z[max_index]

	return sqrt(h_x*h_x + h_y*h_y + h_z*h_z)
end

end # module
