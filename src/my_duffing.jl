#############################################
# Duffing Oscillator Controlled by RFPT o VSSM #
#############################################

######################
# Control Parameters #
######################
adaptive = 1
robust = 1
plotting = 1
#########
# Time  #
#########
δt = 1e-3
N = Int(2e4)
last_idx = N - 1

######################
# Control Parameters #
######################
K = 1e5
B = -1
A = 1e-5


########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ = 1
w = 1
K_VSSM = 500

####################################
# Parameters for Nominal Trajectory#
####################################
ω = 0.5
Amp = 2

##########################
# Exact Model Parameters #
##########################
αₑ = 1
δₑ = 0.2
βₑ = 1

################################
# Approximate Model Parameters #
################################
αₐ = 0.8
δₐ = 0.1
βₐ = 0.9
βₐ = 0.9

Exact(q, q_p, u) = αₑ * q + δₑ * q_p + βₑ * q^3 + u

Approx(q, q_p, q_Np) = q_Np - αₐ * q - δₐ * q_p - βₐ * q^3 #mₐ*q_Np[t]#-μₐ*(1-q^2)*q_p+ωₐ^2*q+αₐ*q^3+λₐ*q^5
ErrorMetrics(h_int, h, h_p) = Λ^2 * h_int + 2 * Λ * h + h_p

KinBlock(S, h, h_p, qN_pp) = K_VSSM * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + qN_pp

G(past_input, past_response, xDnow) = (K + past_input) * (1 + B * tanh(A * (past_response - xDnow))) - K

function nominal_traj(t)
	q_nom_pos_val = Amp * sin(ω * t)
	q_nom_vel_val = ω * Amp * cos(ω * t)
	q_nom_acc_val = -Amp * ω^2 * sin(ω * t)
	return q_nom_pos_val, q_nom_vel_val, q_nom_acc_val
end

function log()
	date_string = Dates.format(now(), "mm-dd_HH-MM")
	line_breaker = "\n####################################################\n"
	is_adaptive = Bool(1 == adaptive)
	is_robust = Bool(1 == robust)
	file = open("./data/Duffing/log.txt", "a")
	file_text = string(line_breaker, date_string, " Parameters:\nadaptive=", is_adaptive, " robust=", is_robust,
		"\nControl Params:\n",
		"K=", K, "\tB=", B, "\tA=", A, "\tw=", w, "\tΛ=", Λ,
		"\nTime:\nδt=", δt, "\tN=", N, "\n",
		"\nApproximate Model Parameters:\nαₐ=", αₐ, "\tδₐ=", δₐ, "\tβₐ=", βₐ,
		"\nExact Model Parameters:\nαₑ=", αₑ, "\tδₑ=", δₑ, "\tβₑ=", βₑ,
		"\nNominal Trajectory Parameters:\nω=", ω, "\tAmp=", Amp, line_breaker)
	write(file, file_text)
	close(file)
end

function duffing_single_run()
	# allocate fresh arrays
	time_mem = zeros(N)
	q_pos = zeros(N);
	q_vel = zeros(N);
	q_acc = zeros(N)
	q_nom_pos = zeros(N);
	q_nom_vel = zeros(N);
	q_nom_acc = zeros(N)
	q_des_acc = zeros(N);
	q_def_acc = zeros(N);
	u_ctrl = zeros(N)
	err_int = 0.0
	# initial conditions
	q_pos[1] = Amp * sin(ω * δt)
	q_vel[1] = Amp * ω * cos(ω * δt)
	q_acc[1] = -Amp * ω^2 * sin(ω * δt)
	for t in 1:last_idx
		time_mem[t] = t * δt

		# Nominal trajectory
		q_nom_pos[t], q_nom_vel[t], q_nom_acc[t] = nominal_traj(time_mem[t])

		# Errors
		err_pos = q_nom_pos[t] - q_pos[t]
		err_vel = q_nom_vel[t] - q_vel[t]

		# Desired acc
		if robust == 1
			S = ErrorMetrics(err_int, err_pos, err_vel)
			q_des_acc[t] = KinBlock(S, err_pos, err_vel, q_nom_acc[t])
		else
			q_des_acc[t] = q_nom_acc[t] + Λ^3 * err_int + 3 * Λ^2 * err_pos + 3 * Λ * err_vel
		end

		# Deformation
		if adaptive == 1 && t > 3
			q_def_acc[t] = G(q_def_acc[t-1], q_acc[t-1], q_des_acc[t])
		else
			q_def_acc[t] = q_des_acc[t]
		end

		u_ctrl[t] = Approx(q_pos[t], q_vel[t], q_def_acc[t])

		# Exact system response
		q_acc[t] = Exact(q_pos[t], q_vel[t], u_ctrl[t])

		# Euler integration
		q_vel[t+1] = q_vel[t] + δt * q_acc[t]
		q_pos[t+1] = q_pos[t] + δt * q_vel[t]
		err_int += δt * err_pos
	end

	# Plots
	figure("Trajectory_tracking");
	grid(true);
	title("Pályakövetés az idő függvényében");
	xlabel("Idő [s]");
	ylabel("Pozíció [m]")
	plot(time_mem[1:last_idx], q_nom_pos[1:last_idx], color = "red", label = "Nominális")
	plot(time_mem[1:last_idx], q_pos[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult");
	legend(loc = 1);
	savefig("allplots_duffing.pdf")
	figure("Acceleration");
	grid(true);
	title("Gyorsulások");
	xlabel("Idő [s]");
	ylabel("Gyorsulás [m/s²]")
	plot(time_mem[1:last_idx], q_nom_acc[1:last_idx], color = "#7684FF", label = "Névleges")
	plot(time_mem[1:last_idx], q_acc[1:last_idx], color = "#FFAA41", linestyle = "--", label = "Realizált")
	plot(time_mem[1:last_idx], q_des_acc[1:last_idx], color = "#2A40FF", linestyle = "-.", label = "Desired");
	legend(loc = 1);
	savefig("temp_duffing.pdf");
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)
	figure("Tracking_Error");
	grid(true);
	title("Követési hiba");
	xlabel("Idő [s]");
	ylabel("Hiba [m]")
	plot(time_mem[1:last_idx], q_nom_pos[1:last_idx] - q_pos[1:last_idx], color = "red");
	savefig("temp_duffing.pdf");
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)
	figure("Control_Signal");
	grid(true);
	title("Irányítójel");
	xlabel("Idő [s]");
	ylabel("Irányítójel [N]")
	plot(time_mem[1:last_idx], u_ctrl[1:last_idx], color = "red");
	savefig("temp_duffing.pdf");
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)
	figure("Phase_Space");
	grid(true);
	title("Fázistér");
	xlabel("Pozíció [m]");
	ylabel("Sebesség [m/s]")
	plot(q_nom_pos[1:last_idx], q_nom_vel[1:last_idx], color = "red", label = "Nominális")
	plot(q_pos[1:last_idx], q_vel[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult");
	legend(loc = 1);
	savefig("temp_duffing.pdf");
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)
	date_string = Dates.format(now(), "mm-dd_HH-MM");
	fileName = "./data/Duffing/" * date_string * ".pdf";
	mv("allplots_duffing.pdf", fileName)
	log();
	if plotting == 1
		show()
	end
	return (q_pos = q_pos, q_vel = q_vel, q_acc = q_acc, q_nom_pos = q_nom_pos)
end

function simulate_duffing()
	return duffing_single_run()
end

##############################
# Define arrays for Plotting #
##############################

time_mem = zeros(N) # t

q_pos = zeros(N) # megvalósult pálya (position)
q_vel = zeros(N) # velocity
q_acc = zeros(N) # acceleration

q_nom_pos = zeros(N) # nominális pálya
q_nom_vel = zeros(N)
q_nom_acc = zeros(N)

u_ctrl = zeros(N)  # szabályozó jel

q_des_acc = zeros(N) # desired acceleration
q_def_acc = zeros(N) # deformed acceleration (adaptive)

err_int = 0
q_pos[1] = Amp * sin(ω * δt)
q_vel[1] = Amp * ω * cos(ω * δt)
q_acc[1] = -Amp * ω^2 * sin(ω * δt)

# Setting initial condition
t_max = 20.0
init_ω = 0.5
init_Amp = 2

# időlépések
δranget = 1
t_range = 0:δranget:t_max

using PyPlot
using PDFmerger
using Dates
# simulate_duffing()
# run_simulation(1, 5, 2, init_Amp, init_ω, t_range)



show()
