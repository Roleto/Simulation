################################
# The Rössler System            #
# Controller type: RFPT o VSSM #
# MIMO System                  #
################################

####################
# The equation of motion #
####################

##########
# ẋ=-(y+z)    #
# ẏ=x+ay #
# ż=b+z(x-c)    #

##########

using LinearAlgebra
using PyPlot
using PDFmerger
using Dates

# Flags (lowercase)
adaptive = 1  # RFPT
robust = 1    # VSSM
plotting = 1

# Time variables
Δt = 1e-3
N = Int(2e4)
last_idx = N - 1

# Control Parameters
K = 1e2
B = -1
A = 1.97e-3

########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ = 1
K_VSSM = 50
w = 1

# Nominal trajectory parameters
A₁ = 2;
ω₁ = 0.5
A₂ = 3;
ω₂ = 0.7
A₃ = 1;
ω₃ = 1.0

# Exact model parameters
aₑ = 0.01
bₑ = 0.2
cₑ = 5.7

# Approximate model parameters
aₐ = 0.1
bₐ = 0.3
cₐ = 5.5

error_limit = 1e-3

# Adaptive deformation function
function g_mimo(prev_input, prev_response, desired, err_limit, K, B, A)
	Δ = (prev_response - desired)
	nrm = norm(Δ)
	if nrm > err_limit
		dir = Δ / nrm
		Bf = B * tanh(A * nrm)
		return (1 + Bf) * prev_input + Bf * K * dir
	else
		return prev_input
	end
end

# Robust error metric (integral + position error)
error_metric(h_int, h) = Λ * h_int + h

# Kinematic block mapping desired velocities (first derivative version used in original)
function kin_block(S, h, q_nom_vel)
	s_vec = K * tanh.(S ./ w)
	qλ = q_nom_vel * Λ
	return [qλ[1] * h[1] * s_vec[1], qλ[2] * h[2] * s_vec[2], qλ[3] * h[3] * s_vec[3]]
end

function rossler_single_run()
	# Allocate arrays
	time_mem = zeros(N)
	x_pos = zeros(N);
	y_pos = zeros(N);
	z_pos = zeros(N)
	x_vel = zeros(N);
	y_vel = zeros(N);
	z_vel = zeros(N)
	# nominal
	x_nom_pos = zeros(N);
	y_nom_pos = zeros(N);
	z_nom_pos = zeros(N)
	x_nom_vel = zeros(N);
	y_nom_vel = zeros(N);
	z_nom_vel = zeros(N)
	# desired & deformed velocities
	x_des_vel = zeros(N);
	y_des_vel = zeros(N);
	z_des_vel = zeros(N)
	x_def_vel = zeros(N);
	y_def_vel = zeros(N);
	z_def_vel = zeros(N)
	# control signals
	u_ctrl_x = zeros(N);
	u_ctrl_y = zeros(N);
	u_ctrl_z = zeros(N)
	# sliding variables
	S_x = zeros(N);
	S_y = zeros(N);
	S_z = zeros(N)
	# integral errors
	h_int_x = 0.0;
	h_int_y = 0.0;
	h_int_z = 0.0

	# initial conditions from nominal at Δt
	x_pos[1] = A₁ * sin(ω₁ * Δt);
	y_pos[1] = A₂ * sin(ω₂ * Δt);
	z_pos[1] = A₃ * sin(ω₃ * Δt)
	x_vel[1] = A₁ * ω₁ * cos(ω₁ * Δt);
	y_vel[1] = A₂ * ω₂ * cos(ω₂ * Δt);
	z_vel[1] = A₃ * ω₃ * cos(ω₃ * Δt)

	past_input = zeros(3);
	past_response = zeros(3)

	for t in 1:last_idx
		time_mem[t] = Δt * t
		# nominal trajectory
		x_nom_pos[t] = A₁ * sin(ω₁ * time_mem[t]);
		x_nom_vel[t] = A₁ * ω₁ * cos(ω₁ * time_mem[t])
		y_nom_pos[t] = A₂ * sin(ω₂ * time_mem[t]);
		y_nom_vel[t] = A₂ * ω₂ * cos(ω₂ * time_mem[t])
		z_nom_pos[t] = A₃ * sin(ω₃ * time_mem[t]);
		z_nom_vel[t] = A₃ * ω₃ * cos(ω₃ * time_mem[t])
		# errors
		h_x = x_nom_pos[t] - x_pos[t]
		h_y = y_nom_pos[t] - y_pos[t]
		h_z = z_nom_pos[t] - z_pos[t]
		# desired velocity
		if robust == 1
			S_vec = error_metric([h_int_x, h_int_y, h_int_z], [h_x, h_y, h_z])
			S_x[t] = S_vec[1];
			S_y[t] = S_vec[2];
			S_z[t] = S_vec[3]
			x_des_vel[t] = x_nom_vel[t] + Λ * h_x + K_VSSM * tanh(S_x[t] / w)
			y_des_vel[t] = y_nom_vel[t] + Λ * h_y + K_VSSM * tanh(S_y[t] / w)
			z_des_vel[t] = z_nom_vel[t] + Λ * h_z + K_VSSM * tanh(S_z[t] / w)
		else
			x_des_vel[t] = x_nom_vel[t] + Λ^2 * h_int_x + 2 * Λ * h_x
			y_des_vel[t] = y_nom_vel[t] + Λ^2 * h_int_y + 2 * Λ * h_y
			z_des_vel[t] = z_nom_vel[t] + Λ^2 * h_int_z + 2 * Λ * h_z
		end
		desired = [x_des_vel[t], y_des_vel[t], z_des_vel[t]]
		# deformation
		if adaptive == 1 && t > 3
			past_input = g_mimo(past_input, past_response, desired, error_limit, K, B, A)
		else
			past_input = desired
		end
		x_def_vel[t] = past_input[1];
		y_def_vel[t] = past_input[2];
		z_def_vel[t] = past_input[3]
		# control signals (approximate inverse of model)
		u_ctrl_x[t] = x_def_vel[t] + y_pos[t] + z_pos[t]
		u_ctrl_y[t] = y_def_vel[t] - x_pos[t] - aₐ * y_pos[t]
		u_ctrl_z[t] = z_def_vel[t] - bₐ - z_pos[t] * (x_pos[t] - cₐ)
		# system dynamics (exact parameters)
		x_vel[t] = -y_pos[t] - z_pos[t] + u_ctrl_x[t]
		y_vel[t] = x_pos[t] + aₑ * y_pos[t] + u_ctrl_y[t]
		z_vel[t] = bₑ + z_pos[t] * (x_pos[t] - cₑ) + u_ctrl_z[t]
		past_response = [x_vel[t], y_vel[t], z_vel[t]]
		# integrate
		x_pos[t+1] = x_pos[t] + Δt * x_vel[t]
		y_pos[t+1] = y_pos[t] + Δt * y_vel[t]
		z_pos[t+1] = z_pos[t] + Δt * z_vel[t]
		# update integrals
		h_int_x += Δt * h_x;
		h_int_y += Δt * h_y;
		h_int_z += Δt * h_z
	end
	# Plotting
	fig = figure("nominal_realized_trajectories")
	title("Nominális és Megvalósult Pályák")
	subplot(311);
	ylabel("X [m]");
	grid(true)
	plot(time_mem[1:last_idx], x_nom_pos[1:last_idx], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:last_idx], x_pos[1:last_idx], color = "#009E73", linestyle = "--", linewidth = 2.5, label = "Megvalósult")
	subplot(312);
	ylabel("Y [m]");
	grid(true)
	plot(time_mem[1:last_idx], y_nom_pos[1:last_idx], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:last_idx], y_pos[1:last_idx], color = "#009E73", linestyle = "--", linewidth = 2.5, label = "Megvalósult")
	subplot(313);
	ylabel("Z [m]");
	xlabel("Idő [s]");
	grid(true)
	plot(time_mem[1:last_idx], z_nom_pos[1:last_idx], color = "#D55E00", linewidth = 1.5, label = "Nominális")
	plot(time_mem[1:last_idx], z_pos[1:last_idx], color = "#009E73", linestyle = "--", linewidth = 2.5, label = "Megvalósult")
	legend(loc = "lower right")
	tight_layout();
	savefig("allplots_rossler.pdf")

	fig = figure("Trajectory_tracking_3D");
	title("Nominális és Megvalósult Pályák 3D");
	grid(true)
	plot3D(x_nom_pos[1:last_idx], y_nom_pos[1:last_idx], z_nom_pos[1:last_idx], color = "red", label = "Nominális")
	plot3D(x_pos[1:last_idx], y_pos[1:last_idx], z_pos[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1)

	# velocities
	fig = figure("velocities");
	title("Sebességek");
	grid(true)
	subplot(311);
	ylabel("X [m/s]");
	grid(true)
	plot(time_mem[1:last_idx], x_nom_vel[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], x_vel[1:last_idx], color = "green", linewidth = 3, linestyle = "--", label = "Megvalósult")
	subplot(312);
	ylabel("Y [m/s]");
	grid(true)
	plot(time_mem[1:last_idx], y_nom_vel[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], y_vel[1:last_idx], color = "green", linewidth = 3, linestyle = "--", label = "Megvalósult")
	subplot(313);
	ylabel("Z [m/s]");
	xlabel("Idő [s]");
	grid(true)
	plot(time_mem[1:last_idx], z_nom_vel[1:last_idx], color = "red", linewidth = 2, label = "Nominális")
	plot(time_mem[1:last_idx], z_vel[1:last_idx], color = "green", linewidth = 3, linestyle = "--", label = "Megvalósult")
	legend(loc = "lower left");
	tight_layout();
	savefig("temp_rossler.pdf");
	append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

	# tracking error
	fig = figure("tracking_error");
	title("Követési Hibák");
	grid(true)
	subplot(311);
	ylabel("X [m]");
	grid(true);
	plot(time_mem[1:last_idx], x_nom_pos[1:last_idx] - x_pos[1:last_idx], color = "red")
	subplot(312);
	ylabel("Y [m]");
	grid(true);
	plot(time_mem[1:last_idx], y_nom_pos[1:last_idx] - y_pos[1:last_idx], color = "red")
	subplot(313);
	ylabel("Z [m]");
	xlabel("Idő [s]");
	grid(true);
	plot(time_mem[1:last_idx], z_nom_pos[1:last_idx] - z_pos[1:last_idx], color = "red")
	tight_layout();
	savefig("temp_rossler.pdf");
	append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

	# control signals
	fig = figure("control_signal");
	title("Irányítójelek");
	grid(true)
	xlabel("Idő [s]");
	ylabel("Irányítójel [N]")
	plot(time_mem[1:last_idx], u_ctrl_x[1:last_idx], color = "red", label = "u_x")
	plot(time_mem[1:last_idx], u_ctrl_y[1:last_idx], color = "green", linestyle = "--", label = "u_y")
	plot(time_mem[1:last_idx], u_ctrl_z[1:last_idx], color = "blue", label = "u_z")
	legend(loc = "lower left");
	tight_layout();
	savefig("temp_rossler.pdf");
	append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

	# phase plots
	fig = figure("phase_x");
	title("Fázistér X");
	grid(true);
	xlabel("Pozíció [m]");
	ylabel("Sebesség [m/s]")
	plot(x_nom_pos[1:last_idx], x_nom_vel[1:last_idx], color = "red", label = "Nominális")
	plot(x_pos[1:last_idx], x_vel[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend();
	tight_layout();
	savefig("temp_rossler.pdf");
	append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

	fig = figure("phase_y");
	title("Fázistér Y");
	grid(true);
	xlabel("Pozíció [m]");
	ylabel("Sebesség [m/s]")
	plot(y_nom_pos[1:last_idx], y_nom_vel[1:last_idx], color = "red", label = "Nominális")
	plot(y_pos[1:last_idx], y_vel[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend();
	tight_layout();
	savefig("temp_rossler.pdf");
	append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

	fig = figure("phase_z");
	title("Fázistér Z");
	grid(true);
	xlabel("Pozíció [m]");
	ylabel("Sebesség [m/s]")
	plot(z_nom_pos[1:last_idx], z_nom_vel[1:last_idx], color = "red", label = "Nominális")
	plot(z_pos[1:last_idx], z_vel[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend();
	tight_layout();
	savefig("temp_rossler.pdf");
	append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

	date_string = Dates.format(now(), "mm-dd_HH-MM")
	fileName = "./data/Rossler/" * date_string * ".pdf"
	mv("allplots_rossler.pdf", fileName)
	return (x_pos = x_pos, y_pos = y_pos, z_pos = z_pos, x_vel = x_vel, y_vel = y_vel, z_vel = z_vel,
		x_nom_pos = x_nom_pos, y_nom_pos = y_nom_pos, z_nom_pos = z_nom_pos)
end

function log()
	file = open("./data/Rossler/log.txt", "a")
	date_string = Dates.format(now(), "mm-dd_HH-MM")
	line_breaker = "\n####################################################\n"
	file_text = string(line_breaker, date_string, " Parameters:\nControl:\nK=", K, "\tB=", B, "\tA=", A, "\tw=", w, "\tΛ=", Λ,
		"\nApproximate Model Parameters:\naₐ=", aₐ, "\tbₐ=", bₐ, "\tcₐ=", cₐ,
		"\nExact Model Parameters:\naₑ=", aₑ, "\tbₑ=", bₑ, "\tcₑ=", cₑ,
		"\nNominal Trajectory Parameters:\nω=[", ω₁, ",", ω₂, ",", ω₃, "] Amp=[", A₁, ",", A₂, ",", A₃, "]",
		line_breaker)
	write(file, file_text);
	close(file)
end

function simulate_rossler()
	res = rossler_single_run()
	log()
	if plotting == 1
		show()
	end
	return res
end




