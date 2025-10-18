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

adaptive = 1  # RFPT
robust = 1    # VSSM
plotting = 1  # was Ploting
# always single run; multi-run logic removed
#################
# Time variable #
#################
δt = 1e-2
N = Int(2e4)
last_idx = N - 1

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

time_mem = zeros(N)

"""Nominal trajectory arrays (position, velocity, acceleration)"""
x_nom_pos = zeros(N)
x_nom_vel = zeros(N)
x_nom_acc = zeros(N)

y_nom_pos = zeros(N)
y_nom_vel = zeros(N)
y_nom_acc = zeros(N)

z_nom_pos = zeros(N)
z_nom_vel = zeros(N)
z_nom_acc = zeros(N)

"""Desired velocities (first derivative)"""
x_des_vel = zeros(N)
y_des_vel = zeros(N)
z_des_vel = zeros(N)

"""Deformed (adaptive) velocities and accelerations"""
x_def_vel = zeros(N)
y_def_vel = zeros(N)
z_def_vel = zeros(N)

x_def_acc = zeros(N)
y_def_acc = zeros(N)
z_def_acc = zeros(N)

"""Control signals"""
u_x = zeros(N)
u_y = zeros(N)
u_z = zeros(N)


x_pos = zeros(N)
y_pos = zeros(N)
z_pos = zeros(N)

x_vel = zeros(N)
y_vel = zeros(N)
z_vel = zeros(N)

x_acc = zeros(N)
y_acc = zeros(N)
z_acc = zeros(N)

S_x = zeros(N)
S_p_x = zeros(N)
S_y = zeros(N)
S_p_y = zeros(N)
S_z = zeros(N)
S_p_z = zeros(N)

past_inputs = zeros(Float64, N, 3)
past_inputs_p = zeros(Float64, N, 3)
past_responses = zeros(Float64, N, 3)
past_input = zeros(3)
past_input_p = zeros(3)
past_inputs = zeros(Float64, N, 3)
past_inputs_p = zeros(Float64, N, 3)
past_response = zeros(3)
past_responses = zeros(Float64, N, 3)
realized_acc = zeros(Float64, N, 3)

function error_metrics(h_int, err_pos, err_vel) # robust metric
	S = Λ * h_int + err_pos
	S_p = Λ * err_pos + err_vel
	return S, S_p
end
###########
# ACC itt #
###########
function g_mimo(r_prev, r_prev_last, f_prev, des_now, err_limit, K, B, A)
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
	fileName = "./data/Lorenz/" * date_string * ".pdf"
	mv("allplots_lorenz.pdf", fileName)
	file = open("./data/Lorenz/log.txt", "a")
	line_breaker = "\n####################################################\n"
	file_text = string(line_breaker, date_string, " Parameters:\nControl Params:\n",
		"K= ", K, "\tB= ", B, "\tA= ", A,
		"\nTime:\nn=", N, "\tδt=", δt,
		"\nKinematic Block:\nΛ=", Λ, "\tK_VSSM=", K_VSSM, "\tw=", w,
		"\nApproximate Model Parameters:\nβₐ=", βₐ, "\tσₐ=", σₐ, "\tρₐ=", ρₐ,
		"\nExact Model Parameters:\nβₑ=", βₑ, "\tσₑ=", σₑ, "\tρₑ=", ρₑ,
		"\nNominal Trajectory Parameters:\nω=", [ω₁, ω₂, ω₃], "\nAmp=", [A₁, A₂, A₃], line_breaker)
	write(file, file_text)
	close(file)
end

function lorenz_single_run(h_int_x, h_int_y, h_int_z, past_input, past_input_p, past_response, error_limit)
	print('.')
	for t ∈ 1:last_idx
		time_mem[t] = δt * t

		#Nominal trajectory for the actual time frame
		x_nom_pos[t] = A₁ * sin(ω₁ * time_mem[t])
		x_nom_vel[t] = A₁ * ω₁ * cos(ω₁ * time_mem[t])
		x_nom_acc[t] = -A₁ * ω₁^2 * sin(ω₁ * time_mem[t])

		y_nom_pos[t] = A₂ * sin(ω₂ * time_mem[t])
		y_nom_vel[t] = A₂ * ω₂ * cos(ω₂ * time_mem[t])
		y_nom_acc[t] = -A₂ * ω₂^2 * sin(ω₂ * time_mem[t])

		z_nom_pos[t] = A₃ * sin(ω₃ * time_mem[t])
		z_nom_vel[t] = A₃ * ω₃ * cos(ω₃ * time_mem[t])
		z_nom_acc[t] = -A₃ * ω₃^2 * sin(ω₃ * time_mem[t])

		#Errors
		err_x_pos = x_nom_pos[t] - x_pos[t]
		err_y_pos = y_nom_pos[t] - y_pos[t]
		err_z_pos = z_nom_pos[t] - z_pos[t]

		err_x_vel = x_nom_vel[t] - x_vel[t]
		err_y_vel = y_nom_vel[t] - y_vel[t]
		err_z_vel = z_nom_vel[t] - z_vel[t]

		if robust == 1
			#the error metric
			S, S_p = error_metrics([h_int_x, h_int_y, h_int_z], [err_x_pos, err_y_pos, err_z_pos], [err_x_vel, err_y_vel, err_z_vel])
			S_x[t] = S[1]
			S_y[t] = S[2]
			S_z[t] = S[3]

			S_p_x[t] = S_p[1]
			S_p_y[t] = S_p[2]
			S_p_z[t] = S_p[3]

			#kinblock
			x_des_vel[t] = x_nom_vel[t] + Λ * err_x_pos + K_VSSM * tanh(S_x[t] / w)
			y_des_vel[t] = y_nom_vel[t] + Λ * err_y_pos + K_VSSM * tanh(S_y[t] / w)
			z_des_vel[t] = z_nom_vel[t] + Λ * err_z_pos + K_VSSM * tanh(S_z[t] / w)

			# xDes_pp[t] = xN_pp[t] + Λ * h_p_x + K_VSSM * (1 / w * sech(S_x[t] / w)^2 * S_p_x[t])
			# yDes_pp[t] = yN_pp[t] + Λ * h_p_y + K_VSSM * (1 / w * sech(S_y[t] / w)^2 * S_p_y[t])
			# zDes_pp[t] = zN_pp[t] + Λ * h_p_z + K_VSSM * (1 / w * sech(S_z[t] / w)^2 * S_p_z[t])

			desired = [x_des_vel[t], y_des_vel[t], z_des_vel[t]]

		else

			#kinblock
			x_des_vel[t] = Λ^2 * h_int_x + 2 * Λ * err_x_pos + x_nom_vel[t]
			y_des_vel[t] = Λ^2 * h_int_y + 2 * Λ * err_y_pos + y_nom_vel[t]
			z_des_vel[t] = Λ^2 * h_int_z + 2 * Λ * err_z_pos + z_nom_vel[t]

			desired = [x_des_vel[t], y_des_vel[t], z_des_vel[t]]
		end

		# Deformation
		if adaptive == 1 && t > 4
			###########
			# ACC itt #
			###########
			past_input, past_input_p = g_mimo(past_input, past_input_p, past_response, desired, error_limit, K, B, A)
		else
			past_input_p = past_input
			past_input = desired
		end

		past_inputs[t, :] = past_input
		past_inputs_p[t, :] = past_input_p

		x_def_vel[t] = past_input[1]
		y_def_vel[t] = past_input[2]
		z_def_vel[t] = past_input[3]

		###########
		# ACC itt #
		###########
		x_def_acc[t] = past_input_p[1]
		y_def_acc[t] = past_input_p[2]
		z_def_acc[t] = past_input_p[3]

		#Control Signal
		u_x[t] = x_def_vel[t] - σₐ * (y_pos[t] - x_pos[t])
		u_y[t] = y_def_vel[t] - x_pos[t] * (ρₐ - z_pos[t]) + y_pos[t]
		u_z[t] = z_def_vel[t] - x_pos[t] * y_pos[t] + βₐ * z_pos[t]

		#System
		x_vel[t] = σₑ * (y_pos[t] - x_pos[t]) + u_x[t]
		y_vel[t] = x_pos[t] * (ρₑ - z_pos[t]) - y_pos[t] + u_y[t]
		z_vel[t] = x_pos[t] * y_pos[t] - βₑ * z_pos[t] + u_z[t]
		past_response = [x_vel[t], y_vel[t], z_vel[t]]
		past_responses[t, :] = past_response

		#Acceleration
		x_acc[t] = (σₑ - σₐ) * (y_vel[t] - x_vel[t]) + x_def_acc[t]
		y_acc[t] = (ρₑ - ρₐ) * x_vel[t] + y_def_acc[t]
		z_acc[t] = (βₐ - βₑ) * z_vel[t] + z_def_acc[t]
		realized_acc[t, :] = [x_acc[t], y_acc[t], z_acc[t]]


		#Integrals
		x_pos[t+1] = x_pos[t] + δt * x_vel[t]
		y_pos[t+1] = y_pos[t] + δt * y_vel[t]
		z_pos[t+1] = z_pos[t] + δt * z_vel[t]

		h_int_x = h_int_x + δt * err_x_pos
		h_int_y = h_int_y + δt * err_y_pos
		h_int_z = h_int_z + δt * err_z_pos
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

	plot(time_mem[1:last_idx], x_nom_pos[1:last_idx], color = "#D55E00", linewidth = 1.5, label = "Nominális", alpha = 0.8)
	plot(time_mem[1:last_idx], x_pos[1:last_idx], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)
	# fill_between(time_mem[1:l], xN[1:l], x[1:l], color = "gray", alpha = 0.3, label = "Error Band X")

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m]")
	plot(time_mem[1:last_idx], y_nom_pos[1:last_idx], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:last_idx], y_pos[1:last_idx], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)
	# fill_between(time_mem[1:l], yN[1:l], x[1:l], color = "gray", alpha = 0.3, label = "Error Band Y ")

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m]")
	plot(time_mem[1:last_idx], z_nom_pos[1:last_idx], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:last_idx], z_pos[1:last_idx], linestyle = "--", color = "#009E73", linewidth = 2.5, label = "Megvalósult", alpha = 0.8)
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
	plot3D(x_nom_pos[1:last_idx], y_nom_pos[1:last_idx], z_nom_pos[1:last_idx], color = "red", label = "Nominális")
	plot3D(x_pos[1:last_idx], y_pos[1:last_idx], z_pos[1:last_idx], color = "green", label = "Megvalósult", linestyle = "--")
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

	plot(time_mem[1:last_idx], past_responses[1:last_idx, 1], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	plot(time_mem[1:last_idx], x_nom_vel[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	# plot(time_mem[1:l], past_inputs[1:l, 1], color = "blue", label = L"\dot{x}^{Des}", linestyle = "-.")

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m/s]")

	plot(time_mem[1:last_idx], y_nom_vel[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], past_responses[1:last_idx, 2], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs[1:l, 2], color = "blue", label = L"\dot{y}^{Des}", linestyle = "-.")

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m/s]")
	plot(time_mem[1:last_idx], z_nom_vel[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], past_responses[1:last_idx, 3], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
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


	plot(time_mem[1:last_idx], x_nom_acc[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], realized_acc[1:last_idx, 1], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs_p[1:l, 1], color = "blue", label = L"\ddot{x}^{Des}", linestyle = "-.")

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m/s²]")

	plot(time_mem[1:last_idx], y_nom_acc[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], realized_acc[1:last_idx, 2], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
	# plot(time_mem[1:l], past_inputs_p[1:l, 2], color = "blue", label = L"\ddot{y}^{Des}", linestyle = "-.")

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m/s²]")
	plot(time_mem[1:last_idx], z_nom_acc[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], realized_acc[1:last_idx, 3], color = "green", linewidth = 3, label = "Megvalósult", linestyle = "--")
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
	plot(time_mem[1:last_idx], x_nom_pos[1:last_idx] - x_pos[1:last_idx], color = "red", linewidth = 2)
	# fill_between(time_mem[1:l], xN[1:l], x[1:l], color = "gray", alpha = 0.3)

	subplot(312, sharex = ax1)
	ax2 = gca()
	grid(true)
	ylabel("Y [m]")
	plot(time_mem[1:last_idx], y_nom_pos[1:last_idx] - y_pos[1:last_idx], color = "red", linewidth = 2)

	subplot(313, sharex = ax2)
	ax3 = gca()
	grid(true)
	xlabel("Idő [s]")
	ylabel("Z [m]")
	plot(time_mem[1:last_idx], z_nom_pos[1:last_idx] - z_pos[1:last_idx], color = "red", linewidth = 2)

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

	plot(time_mem[1:last_idx], u_x[1:last_idx], color = "red", label = L"$u_x$")
	plot(time_mem[1:last_idx], u_y[1:last_idx], color = "green", label = L"$u_y$", linestyle = "--")
	plot(time_mem[1:last_idx], u_z[1:last_idx], color = "blue", label = L"$u_z$")

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
	plot(x_vel[1:last_idx], x_pos[1:last_idx], color = "green", linewidth = 2.5, label = "Megvalósult", linestyle = "--")
	plot(x_nom_vel[1:last_idx], x_nom_pos[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
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
	plot(y_vel[1:last_idx], y_pos[1:last_idx], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
	plot(y_nom_vel[1:last_idx], y_nom_pos[1:last_idx], color = "red", linewidth = 2, label = "Nominális")

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
	plot(z_vel[1:last_idx], z_pos[1:last_idx], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
	plot(z_nom_vel[1:last_idx], z_nom_pos[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	legend(loc = "lower left", fancybox = "True")
	tight_layout()

	savefig("temp_lorenz.pdf")
	append_pdf!("allplots_lorenz.pdf", "temp_lorenz.pdf", cleanup = true)
	date_string = Dates.format(now(), "mm-dd_HH-MM");
	fileName = "./data/Lorenz/" * date_string * ".pdf";
	mv("allplots_lorenz.pdf", fileName)

end

function simulate_lorenz()
	# initial already set before call
	lorenz_single_run(h_int_x, h_int_y, h_int_z, past_input, past_input_p, past_response, error_limit)
	log()
	if (plotting == 1)
		show()
	end
end

#initial conditions
h_int_x = 0
h_int_y = 0
h_int_z = 0

error_limit = 1e-3

x_pos[1] = A₁ * sin(ω₁ * δt)
y_pos[1] = A₂ * sin(ω₂ * δt)
z_pos[1] = A₃ * sin(ω₃ * δt)
x_vel[1] = A₁ * ω₁ * cos(ω₁ * δt)
y_vel[1] = A₂ * ω₂ * cos(ω₂ * δt)
z_vel[1] = A₃ * ω₃ * cos(ω₃ * δt)
x_acc[1] = -A₁ * ω₁^2 * sin(ω₁ * δt)
y_acc[1] = -A₂ * ω₂^2 * sin(ω₂ * δt)
z_acc[1] = -A₃ * ω₃^2 * sin(ω₃ * δt)
# Setting initial condition
t_max = 20.0
init_ω = 0.5
init_Amp = 2

# időlépések
δranget = 1
t_range = 0:δranget:t_max

using PDFmerger
using Dates
# simulate_lorenz()
show()
