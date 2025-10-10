################################
# The Lorenz System            #
# Controller type: RFPT o VSSM #
# MIMO System                  #
################################

##########################
# The equation of motion #
##########################

##############
# ẋ=σ(y-x)   #
# ẏ=x(ρ-z)-y #
# ż=xy-βz    #
##############

using LinearAlgebra
using PyPlot

Adaptive = 1  #RFPT
Robust = 1    #VSSM
Ploting = 1
SingleRun = 1
#################
# Time variable #
#################
δt = 1e-2
LONG = Int(2e4)
# LONG_ZOOM = Int(1.5e3 * 0.75)
l = LONG - 1

######################
# Control Parameters #
######################
K = 1e2
B = -1
A = 1.97e-3

########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ = 1
K_VSSM = 50
w = 1

#############################################
# Parameters to design a Nominal Trajectory #
#############################################
A₁ = 2
ω₁ = 0.5
A₂ = 3
ω₂ = 0.7
A₃ = 1
ω₃ = 1

##########################
# Exact model parameters #
##########################
βₑ = 8 / 3
σₑ = 5
ρₑ = 40

############################
# Approx. model parameters #
############################
βₐ = 7 / 3
σₐ = 4
ρₐ = 36

time_mem = zeros(LONG)

# Nominal Trajectory
xN = zeros(LONG)
xN_p = zeros(LONG)
xN_pp = zeros(LONG)

yN = zeros(LONG)
yN_p = zeros(LONG)
yN_pp = zeros(LONG)

zN = zeros(LONG)
zN_p = zeros(LONG)
zN_pp = zeros(LONG)

# Desired
xDes_p = zeros(LONG)
yDes_p = zeros(LONG)
zDes_p = zeros(LONG)

# xDes_pp = zeros(LONG)
# yDes_pp = zeros(LONG)
# zDes_pp = zeros(LONG)

# Deformed
xDef_p = zeros(LONG)
yDef_p = zeros(LONG)
zDef_p = zeros(LONG)

xDef_pp = zeros(LONG)
yDef_pp = zeros(LONG)
zDef_pp = zeros(LONG)

# Control Signal
u_x = zeros(LONG)
u_y = zeros(LONG)
u_z = zeros(LONG)


x_p = zeros(LONG)
y_p = zeros(LONG)
z_p = zeros(LONG)

x_pp = zeros(LONG)
y_pp = zeros(LONG)
z_pp = zeros(LONG)

x = zeros(LONG)
y = zeros(LONG)
z = zeros(LONG)

S_x = zeros(LONG)
S_p_x = zeros(LONG)
S_y = zeros(LONG)
S_p_y = zeros(LONG)
S_z = zeros(LONG)
S_p_z = zeros(LONG)

past_input = zeros(3)
past_input_p = zeros(3)
past_inputs = zeros(Float64, LONG, 3)
past_inputs_p = zeros(Float64, LONG, 3)
past_response = zeros(3)
past_responses = zeros(Float64, LONG, 3)
accelerations = zeros(Float64, LONG, 3)

function ErrorMetric(hint, h, h_p) # = Λ * hint + h
	S = Λ * hint + h
	S_p = Λ * h + h_p
	return S, S_p
end
###########
# ACC itt #
###########
function G_MIMO(r_prev, r_prev_last, f_prev, des_now, err_limit, K, B, A)
	Amatr_h = (f_prev - des_now)
	error_norm = norm(Amatr_h)
	if error_norm > err_limit
		e_direction = Amatr_h / error_norm
		B_factor = B * tanh(A * error_norm)
		G = (1 + B_factor) * r_prev + B_factor * K * e_direction
	else
		G = r_prev
	end
	G_p = (G - r_prev_last)
	return G, G_p
end

function log()
	date_string = Dates.format(now(), "mm-dd_HH-MM")
	fileName = "./Plots/Lorenz/" * date_string * ".pdf"
	mv("allplots_lorenz.pdf", fileName)
	file = open("./Plots/Lorenz/log.txt", "a")
	line_breaker = "\n####################################################\n"
	file_text = string(line_breaker, date_string, " Következö paraméterekkel volt használva:\nControl Params:\n",
		"K= ", K, "\tB= ", B, "\tA= ", A,
		"\nTime Parameters:\nLONG=", LONG, "\tδt=", δt,
		"\nKinematik Block parameters:\nΛ=", Λ, "\tK_VSSM=", K_VSSM, "\tw=", w,
		"\nApproximate Model Parameters:\nβₐ=", βₐ, "\tσₐ=", σₐ, "\tρₐ=", ρₐ,
		"\nExact Model Parameters:\nβₑ=", βₑ, "\tσₑ=", σₑ, "\tρₑ=", ρₑ,
		"\nNominal Trajectory Parameters:\nω=", [ω₁, ω₂, ω₃], "\nAmp=", [A₁, A₂, A₃], line_breaker)
	write(file, file_text)
	close(file)
end

function singlerun(hint_x, hint_y, hint_z, past_input, past_input_p, past_response, error_limit, idx, lastPlot = true)
	print('.')
	if (Ploting == 1 && SingleRun == 0 && lastPlot)
		close("all")
	end
	for t ∈ 1:l
		time_mem[t] = δt * t

		#Nominal trajectory for the actual time frame
		xN[t] = A₁ * sin(ω₁ * time_mem[t])
		xN_p[t] = A₁ * ω₁ * cos(ω₁ * time_mem[t])
		xN_pp[t] = -A₁ * ω₁^2 * sin(ω₁ * time_mem[t])

		yN[t] = A₂ * sin(ω₂ * time_mem[t])
		yN_p[t] = A₂ * ω₂ * cos(ω₂ * time_mem[t])
		yN_pp[t] = -A₂ * ω₂^2 * sin(ω₂ * time_mem[t])

		zN[t] = A₃ * sin(ω₃ * time_mem[t])
		zN_p[t] = A₃ * ω₃ * cos(ω₃ * time_mem[t])
		zN_pp[t] = -A₃ * ω₃^2 * sin(ω₃ * time_mem[t])

		#Errors
		h_x = xN[t] - x[t]
		h_y = yN[t] - y[t]
		h_z = zN[t] - z[t]

		h_p_x = xN_p[t] - x_p[t]
		h_p_y = yN_p[t] - y_p[t]
		h_p_z = zN_p[t] - z_p[t]

		if Robust == 1
			#the error metric
			S, S_p = ErrorMetric([hint_x, hint_y, hint_z], [h_x, h_y, h_z], [h_p_x, h_p_y, h_p_z])
			S_x[t] = S[1]
			S_y[t] = S[2]
			S_z[t] = S[3]

			S_p_x[t] = S_p[1]
			S_p_y[t] = S_p[2]
			S_p_z[t] = S_p[3]

			#kinblock
			xDes_p[t] = xN_p[t] + Λ * h_x + K_VSSM * tanh(S_x[t] / w)
			yDes_p[t] = yN_p[t] + Λ * h_y + K_VSSM * tanh(S_y[t] / w)
			zDes_p[t] = zN_p[t] + Λ * h_z + K_VSSM * tanh(S_z[t] / w)

			# xDes_pp[t] = xN_pp[t] + Λ * h_p_x + K_VSSM * (1 / w * sech(S_x[t] / w)^2 * S_p_x[t])
			# yDes_pp[t] = yN_pp[t] + Λ * h_p_y + K_VSSM * (1 / w * sech(S_y[t] / w)^2 * S_p_y[t])
			# zDes_pp[t] = zN_pp[t] + Λ * h_p_z + K_VSSM * (1 / w * sech(S_z[t] / w)^2 * S_p_z[t])

			desired = [xDes_p[t], yDes_p[t], zDes_p[t]]

		else

			#kinblock
			xDes_p[t] = Λ^2 * hint_x + 2 * Λ * h_x + xN_p[t]
			yDes_p[t] = Λ^2 * hint_y + 2 * Λ * h_y + yN_p[t]
			zDes_p[t] = Λ^2 * hint_z + 2 * Λ * h_z + zN_p[t]

			desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
		end

		# Deformation
		if Adaptive == 1 && t > 4
			###########
			# ACC itt #
			###########
			past_input, past_input_p = G_MIMO(past_input, past_input_p, past_response, desired, error_limit, K, B, A)
		else
			past_input_p = past_input
			past_input = desired
		end

		past_inputs[t, :] = past_input
		past_inputs_p[t, :] = past_input_p

		xDef_p[t] = past_input[1]
		yDef_p[t] = past_input[2]
		zDef_p[t] = past_input[3]

		###########
		# ACC itt #
		###########
		xDef_pp[t] = past_input_p[1]
		yDef_pp[t] = past_input_p[2]
		zDef_pp[t] = past_input_p[3]

		#Control Signal
		u_x[t] = xDef_p[t] - σₐ * (y[t] - x[t])
		u_y[t] = yDef_p[t] - x[t] * (ρₐ - z[t]) + y[t]
		u_z[t] = zDef_p[t] - x[t] * y[t] + βₐ * z[t]

		#System
		x_p[t] = σₑ * (y[t] - x[t]) + u_x[t]
		y_p[t] = x[t] * (ρₑ - z[t]) - y[t] + u_y[t]
		z_p[t] = x[t] * y[t] - βₑ * z[t] + u_z[t]
		past_response = [x_p[t], y_p[t], z_p[t]]
		past_responses[t, :] = past_response

		#Acceleration
		x_pp[t] = (σₑ - σₐ) * (y_p[t] - x_p[t]) + xDef_pp[t]
		y_pp[t] = (ρₑ - ρₐ) * x_p[t] + yDef_pp[t]
		z_pp[t] = (βₐ - βₑ) * z_p[t] + zDef_pp[t]
		accelerations[t, :] = [x_pp[t], y_pp[t], z_pp[t]]


		#Integrals
		x[t+1] = x[t] + δt * x_p[t]
		y[t+1] = y[t] + δt * y_p[t]
		z[t+1] = z[t] + δt * z_p[t]

		hint_x = hint_x + δt * h_x
		hint_y = hint_y + δt * h_y
		hint_z = hint_z + δt * h_z
	end



	##############################################
	# Nominal and Realized Trajectories Plotting #
	##############################################

	fig_caption = "nominal_realized_trajectories"
	fig = figure(fig_caption)
	title("Nominális és Megvalósult Pályák")
	subplot(311)
	ax1 = gca()
	grid1 = grid(true)

	ylabel("X [m]")

	plot(time_mem[1:l], xN[1:l], color = "#D55E00", linewidth = 1.5, label = "Nominális", alpha = 0.8)
	plot(time_mem[1:l], x[1:l], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)
	# fill_between(time_mem[1:l], xN[1:l], x[1:l], color = "gray", alpha = 0.3, label = "Error Band X")

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m]")
	plot(time_mem[1:l], yN[1:l], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:l], y[1:l], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)
	# fill_between(time_mem[1:l], yN[1:l], x[1:l], color = "gray", alpha = 0.3, label = "Error Band Y ")

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m]")
	plot(time_mem[1:l], zN[1:l], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:l], z[1:l], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)
	# fill_between(time_mem[1:l], zN[1:l], x[1:l], color = "gray", alpha = 0.3, label = "Error Band Z")

	handles, labels = ax1.get_legend_handles_labels()
	fig.legend(handles, labels, loc = "lower center", bbox_to_anchor = (0.5, -0.01), ncol = 2)
	subplots_adjust(hspace = 0.2, bottom = 0.2)

	tight_layout()
	fig[:canvas][:draw]()
	savefig("allplots_lorenz.pdf")

	fig = figure("Trajectory_tracking_3D")
	title("Nominális és Megvalósult Pályák 3D")
	grid(true)
	plot3D(xN[1:l], yN[1:l], zN[1:l], color = "red", label = "Nominális")
	plot3D(x[1:l], y[1:l], z[1:l], color = "green", label = "Megvalósult", linestyle = "--")
	legend(loc = 1, borderaxespad = 0)

	#######################
	# Velocities Plotting #
	#######################
	fig_caption = "velocities"
	fig = figure(fig_caption)
	grid(true)
	title("Nominális és Megvalósult Sebességek")

	subplot(311)
	ax1 = gca()
	grid1 = grid(true)
	ylabel("X [m/s]")

	plot(time_mem[1:l], past_responses[1:l, 1], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	plot(time_mem[1:l], xN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
	# plot(time_mem[1:l], past_inputs[1:l, 1], color = "blue", label = L"\dot{x}^{Des}", linestyle = "-.")

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m/s]")

	plot(time_mem[1:l], yN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:l], past_responses[1:l, 2], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs[1:l, 2], color = "blue", label = L"\dot{y}^{Des}", linestyle = "-.")

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m/s]")
	plot(time_mem[1:l], zN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:l], past_responses[1:l, 3], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs[1:l, 3], color = "blue", label = "Desired", linestyle = "-.")
	legend(loc = "lower left", fancybox = "True")

	tight_layout()
	subplots_adjust(hspace = 0.0)
	fig[:canvas][:draw]()

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

	#######################
	# Acceleration Plotting #
	#######################
	fig_caption = "accelerations"
	fig = figure(fig_caption)
	grid(true)
	title("Nominális és Megvalósult Gyorsulások")

	subplot(311)
	ax1 = gca()
	grid1 = grid(true)
	ylabel("X [m/s²]")


	plot(time_mem[1:l], xN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:l], accelerations[1:l, 1], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs_p[1:l, 1], color = "blue", label = L"\ddot{x}^{Des}", linestyle = "-.")

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m/s²]")

	plot(time_mem[1:l], yN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:l], accelerations[1:l, 2], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs_p[1:l, 2], color = "blue", label = L"\ddot{y}^{Des}", linestyle = "-.")

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m/s²]")
	plot(time_mem[1:l], zN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:l], accelerations[1:l, 3], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs_p[1:l, 3], color = "blue", label = "Desired", linestyle = "-.")
	handles, labels = ax1.get_legend_handles_labels()
	fig.legend(handles, labels, loc = "lower center", bbox_to_anchor = (0.5, -0.01), ncol = 2)
	subplots_adjust(hspace = 0.2, bottom = 0.2)

	tight_layout()
	subplots_adjust(hspace = 0.0)
	fig[:canvas][:draw]()

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

	############################
	# Tracking Errors Plotting #
	############################

	fig_caption = "tracking_error"
	fig = figure(fig_caption)

	subplot(311)
	ax1 = gca()
	grid1 = grid(true)
	title("Követési Hibák az Idő Függvényében")
	ylabel("X [m]")
	plot(time_mem[1:l], xN[1:l] - x[1:l], color = "red", linewidth = 2)
	# fill_between(time_mem[1:l], xN[1:l], x[1:l], color = "gray", alpha = 0.3)

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m]")
	plot(time_mem[1:l], yN[1:l] - y[1:l], color = "red", linewidth = 2)

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m]")
	plot(time_mem[1:l], zN[1:l] - z[1:l], color = "red", linewidth = 2)

	# legend(loc = "lower left", fancybox = "True")
	tight_layout()
	subplots_adjust(hspace = 0.190)
	fig[:canvas][:draw]()

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

	###########################
	# Control signal Plotting #
	###########################

	fig_caption = "control_signal"
	figure(fig_caption)
	grid(true)
	title("Irányítójelek az Idő Függvényében")
	xlabel("Idő [s]")
	ylabel("Irányítójel [N]")

	plot(time_mem[1:l], u_x[1:l], color = "red", label = L"$u_x$")
	plot(time_mem[1:l], u_y[1:l], color = "green", label = L"$u_y$", linestyle = "--")
	plot(time_mem[1:l], u_z[1:l], color = "blue", label = L"$u_z$")

	legend(loc = "lower left", fancybox = "True")
	legend()
	tight_layout()
	fig[:canvas][:draw]()

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

	###############################
	# Phase Trajectories Plotting #
	###############################

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

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

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


	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)

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

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)


end

function run_simulation(q0, q_p0, q_pp0, init_Amp, init_ω, t_range)

	if (SingleRun == 1)
		x[1] = q0[1]
		y[1] = q0[2]
		z[1] = q0[3]
		x_p[1] = q_p0[1]
		y_p[1] = q_p0[2]
		z_p[1] = q_p0[3]
		x_pp[1] = q_pp0[1]
		y_pp[1] = q_pp0[2]
		z_pp[1] = q_pp0[3]
		singlerun(hint_x, hint_y, hint_z, past_input, past_input_p, past_response, error_limit, 1)
	else
		# Define initial conditions: (q0, q_p0)
		condition_q = init_Amp .* sin.(init_ω .* t_range)
		condition_qp = init_Amp * init_ω .* cos.(init_ω .* t_range)
		for (idx, (q0, q_p0)) in enumerate(zip(condition_q, condition_qp))
			# Initialize memory arrays again for each run
			x[1] = q0
			y[1] = q0
			z[1] = q0
			x_p[1] = q_p0
			y_p[1] = q_p0
			z_p[1] = q_p0

			singlerun(hint_x, hint_y, hint_z, past_input, past_input_p, past_response, error_limit, idx)
		end
	end
	log()
	if (Ploting == 1 && SingleRun == 1)
		show()
	end
end

#initial conditions
hint_x = 0
hint_y = 0
hint_z = 0

error_limit = 1e-3

x[1] = A₁ * sin(ω₁ * δt)
y[1] = A₂ * sin(ω₂ * δt)
z[1] = A₃ * sin(ω₃ * δt)
x_p[1] = A₁ * ω₁ * cos(ω₁ * δt)
y_p[1] = A₂ * ω₂ * cos(ω₂ * δt)
z_p[1] = A₃ * ω₃ * cos(ω₃ * δt)
x_pp[1] = -A₁ * ω₁^2 * sin(ω₁ * δt)
y_pp[1] = -A₂ * ω₂^2 * sin(ω₂ * δt)
z_pp[1] = -A₃ * ω₃^2 * sin(ω₃ * δt)
# Setting initial condition
t_max = 20.0
init_ω = 0.5
init_Amp = 2

# időlépések
δranget = 1
t_range = 0:δranget:t_max

using PDFmerger
using Dates
run_simulation((x[1], y[1], z[1]), (x_p[1], y_p[1], z_p[1]), (x_pp[1], y_pp[1], z_pp[1]), init_Amp, init_ω, t_range)
show()
