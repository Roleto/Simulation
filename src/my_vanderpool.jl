####################################################
# Van der Pol Oscillator Controlled by RFPT o VSSM #
####################################################
adaptive = 1
robust = 1
plotting = 1  # formerly Ploting

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
ω_nom = 0.5
Amp_nom = 2

##########################
# Exact Model Parameters #
##########################
mₑ = 1
μₑ = 0.4
ωₑ = 0.46
αₑ = 1
λₑ = 0.1

################################
# Approximate Model Parameters #
################################
mₐ = 0.8
μₐ = 0.5
ωₐ = 0.42
αₐ = 0.9
λₐ = 0.09
ll = 0.015

error_metric(h_int, h, h_p) = Λ^2 * h_int + 2 * Λ * h + h_p
kin_block(S, h, h_p, q_nom_acc) = K * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + q_nom_acc

function exact_acc(u_ctrl, q_pos, q_vel)
	q_acc = (u_ctrl + μₑ * (1 - q_pos^2) * q_vel - ωₑ^2 * q_pos - αₑ * q_pos^3 - λₑ * q_pos^5) / mₑ
	return q_acc
end

# Approximate inverse (kept minimal per original intent)
function approx_ctrl(q_acc)
	u_ctrl = mₐ * q_acc  # simplified (nonlinear terms omitted as in original)
	return u_ctrl * ll
end

function adapt_g(past_input, past_response, desired_now)
	out = (K + past_input) * (1 + B * tanh(A * (past_response - desired_now))) - K
	return out
end

function vanderpol_nominal_traj(t)
	q_nom_pos = Amp_nom * sin(ω_nom * t * Δt)
	q_nom_vel = ω_nom * Amp_nom * cos(ω_nom * t * Δt)
	q_nom_acc = -Amp_nom * ω_nom^2 * sin(ω_nom * t * Δt)
	return q_nom_pos, q_nom_vel, q_nom_acc
end
function log()
	file = open("./data/VanDerPol/log.txt", "a")
	line_breaker = "\n####################################################\n"
	date_string = Dates.format(now(), "mm-dd_HH-MM")
	file_text = string(line_breaker, date_string,
		" Initial: q0=", q_pos[1], " q_p0=", q_vel[1], " q_pp0=", q_acc[1],
		"\n Parameters:\nadaptive=", Bool(adaptive==1), " robust=", Bool(robust==1), "\nControl Params:\n",
		"K=", K, "\tB=", B, "\tA=", A, "\tw=", w, "\tΛ=", Λ,
		"\nTime:\nΔt=", Δt, "\tN=", N, "\n",
		"\nApproximate Model Parameters:\nμₐ=", μₐ, "\tωₐ=", ωₐ, "\tαₐ=", αₐ, "\tλₐ=", λₐ, "\tmₐ=", mₐ, "\tll=", ll,
		"\nExact Model Parameters:\nμₑ=", μₑ, "\tωₑ=", ωₑ, "\tαₑ=", αₑ, "\tλₑ=", λₑ, "\tmₑ=", mₑ,
		"\nNominal Trajectory Parameters:\nω_nom=", ω_nom, "\tAmp_nom=", Amp_nom, line_breaker)
	write(file, file_text)
	close(file)
end

function vanderpol_single_run()
	# allocate
	time_mem = zeros(N)
	q_pos = zeros(N)
	q_vel = zeros(N)
	q_acc = zeros(N)
	q_nom_pos = zeros(N)
	q_nom_vel = zeros(N)
	q_nom_acc = zeros(N)
	q_des_acc = zeros(N)
	q_def_acc = zeros(N)
	u_ctrl = zeros(N)
	q_err_int = 0.0
	# initial conditions from nominal at first step
	q_pos[1] = Amp_nom * sin(ω_nom * Δt)
	q_vel[1] = Amp_nom * ω_nom * cos(ω_nom * Δt)
	q_acc[1] = -Amp_nom * ω_nom^2 * sin(ω_nom * Δt)
	for t in 1:last_idx
		time_mem[t] = t * Δt
		(q_nom_pos[t], q_nom_vel[t], q_nom_acc[t]) = vanderpol_nominal_traj(t)
		# errors
		q_err_pos = q_nom_pos[t] - q_pos[t]
		q_err_vel = q_nom_vel[t] - q_vel[t]
		# desired acceleration via robust or non-robust law
		if robust == 1
			S = Λ^2 * q_err_int + 2 * Λ * q_err_pos + q_err_vel
			q_des_acc[t] = kin_block(S, q_err_pos, q_err_vel, q_nom_acc[t])
		else
			q_des_acc[t] = q_nom_acc[t] + Λ^3 * q_err_int + 3 * Λ^2 * q_err_pos + 3 * Λ * q_err_vel
		end
		# adaptive deformation
		if adaptive == 1 && t > 10
			q_def_acc[t] = adapt_g(q_def_acc[t-1], q_acc[t-1], q_des_acc[t])
		else
			q_def_acc[t] = q_des_acc[t]
		end
		# control
		u_ctrl[t] = approx_ctrl(q_def_acc[t])
		# exact system response
		q_acc[t] = exact_acc(u_ctrl[t], q_pos[t], q_vel[t])
		# integrate
		q_vel[t+1] = q_vel[t] + Δt * q_acc[t]
		q_pos[t+1] = q_pos[t] + Δt * q_vel[t]
		q_err_int += Δt * q_err_pos
	end
	# Plotting 

	figure("Trajectory_tracking")
	grid(true)
	title("Pályakövetés az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Pozíció [m]")
	plot(time_mem[1:last_idx], q_nom_pos[1:last_idx], color = "red", label = "Nominális")
	plot(time_mem[1:last_idx], q_pos[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("allplots_vander.pdf")

	figure("Velocities")
	grid(true)
	title("Sebesség az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Sebesség [m/s]")
	plot(time_mem[1:last_idx], q_nom_vel[1:last_idx], color = "red", label = "Nominális")
	plot(time_mem[1:last_idx], q_vel[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Accelerations")
	grid(true)
	title("Gyorsulás az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Gyorsulás [m/s²]")
	plot(time_mem[1:last_idx], q_nom_acc[1:last_idx], color = "red", label = "Nominális")
	plot(time_mem[1:last_idx], q_acc[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Tracking_Error")
	title("Követési hiba az idő függvényében")
	grid(true)
	xlabel("Idő [s]")
	ylabel("Követési hiba [m]")
	plot(time_mem[1:last_idx], q_nom_pos[1:last_idx] - q_pos[1:last_idx], color = "red")
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Control_Signal")
	title("Irányítójel az idő függvényében")
	grid(true)
	xlabel("Idő [s]")
	ylabel("Irányítójel [N]")
	plot(time_mem[1:last_idx], u_ctrl[1:last_idx], color = "red")
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Phase_Space")
	title("Fázistér")
	xlabel("Pozíció [m]")
	ylabel("Sebesség [m/s]")
	grid(true)
	plot(q_nom_vel[1:last_idx], q_nom_pos[1:last_idx], color = "red", label = "Nominális")
	plot(q_vel[1:last_idx], q_pos[1:last_idx], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	date_string = Dates.format(now(), "mm-dd_HH-MM")
	fileName = "./data/VanDerPol/" * date_string * ".pdf"
	mv("allplots_vander.pdf", fileName)
	return (q_pos = q_pos, q_vel = q_vel, q_acc = q_acc, q_nom_pos = q_nom_pos, q_nom_vel = q_nom_vel, q_nom_acc = q_nom_acc)
end

# Required packages
using PyPlot
using PDFmerger
using Dates

show()
