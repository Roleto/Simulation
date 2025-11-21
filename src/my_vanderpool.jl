using PyPlot
using Pkg
PyPlot.matplotlib.use("TkAgg")
Pkg.activate(".")
####################################################
# Van der Pol Oscillator Controlled by RFPT o VSSM #
####################################################
Adaptive = 1
Robust = 1
Ploting = 1
SingleRun = 1

#########
# Time  #
#########
δt = 1e-3
LONG = Int(2e4)
l = LONG - 1

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

ErrorMetrics(h_int, h, h_p) = Λ^2 * h_int + 2 * Λ * h + h_p
KinBlock(S, h, h_p, qN_pp) = K * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + qN_pp

function Exact(u, q, q_p)
	q_pp = (u + μₑ * (1 - q^2) * q_p - ωₑ^2 * q - αₑ * q^3 - λₑ * q^5) / mₑ
	return q_pp
end

# function Approx(q, q_p, q_pp)
function Approx(q_pp)
	u = mₐ * q_pp# - μₐ * (1 - q^2) * q_p + ωₐ^2 * q + αₐ * q^3 + λₐ * q^5
	return u * ll
end

function sigmoid(x)
	s = x / (1 + abs(x))
	return s
end

function G(past_input, past_response, xDnow)
	out = (K + past_input) * (1 + B * tanh(A * (past_response - xDnow))) - K
	return out
end

function nominalTraj(t)
	qN = Amp * sin(ω * t * δt)
	q_pN = ω * Amp * cos(ω * t * δt)
	q_ppN = -Amp * ω^2 * sin(ω * t * δt)
	return [qN, q_pN, q_ppN]
end
function log()
	file = open("./Plots/Vanderpool/log.txt", "a")
	line_breaker = "\n####################################################\n"
	isAdaptiv = Bool(1 == Adaptive)
	isRobust = Bool(1 == Robust)
	date_string = Dates.format(now(), "mm-dd_HH-MM")
	file_text = string(line_breaker, date_string, "Kezdőértékek: q0=", q_mem[1], "q_p0=", q_p_mem[1], "q_pp0=", q_pp_mem[1],
		"\n Következö paraméterekkel volt használva:
		\nAdaptiv=", isAdaptiv, "Robust=", isRobust, "\nControl Params:\n",
		"K= ", K, "\tB= ", B, "\tA= ", A, "\tw= ", w, "\tΛ= ", Λ,
		"\nTime variable :\nδt=", δt, "\t=LONG", LONG, "\n",
		"\nApproximate Model Parameters:\nμₐ=", μₐ, "\tωₐ=", ωₐ, "\tαₐ=", αₐ, "\tλₐ=", λₐ, "\tmₐ=", mₐ, "\tll=", ll,
		"\nExact Model Parameters:\nμₑ=", μₑ, "\tωₑ=", ωₑ, "\tαₑ=", αₑ, "\tλₑ=", λₑ, "\tmₑ=", mₑ,
		"\nNominal Trajectory Parameters:\nω=", ω, "\tAmp=", Amp, line_breaker)
	write(file, file_text)
	close(file)
end
function singlerun(h_int, q, q_p, q_pp, qN_p, qN_pp, u, idx, lastPlot = true)
	print('.')
	if (Ploting == 1 && SingleRun == 0 && lastPlot)
		close("all")
	end
	for t ∈ 1:l
		time_mem[t] = t * δt

		#define the nominal trajectory
		(qN_mem[t], qN_p[t], qN_pp[t]) = nominalTraj(t)

		#Compute the Error
		h = qN_mem[t] - q[t]
		h_p = qN_p[t] - q_p[t]

		#the kinematic block
		if Robust == 1
			S = Λ^2 * h_int + 2 * Λ * h + h_p
			qDes_pp_mem[t] = KinBlock(S, h, h_p, qN_pp[t])
			# K_VSSM * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + qN_pp[i]
		else
			qDes_pp_mem[t] = qN_pp[t] + Λ^3 * h_int + 3 * Λ^2 * h + 3 * Λ * h_p
		end

		#Deformation
		if Adaptive == 1 && t > 10
			qDef_pp[t] = G(qDef_pp[t-1], q_pp[t-1], qDes_pp_mem[t])
		else
			qDef_pp[t] = qDes_pp_mem[t]
		end

		#Compute control signal
		u[t] = Approx(qDef_pp[t])

		# Compute the exact systems's respons
		q_pp[t] = Exact(u[t], q[t], q_p[t])

		#Integrate back with Euler's method
		q_p[t+1] = q_p[t] + δt * q_pp[t]
		q[t+1] = q[t] + δt * q_p[t]
		q_mem[t] = q[t]
		h_int = h_int + δt * h
	end
	# Plotting 

	figure("Trajectory_tracking")
	grid(true)
	title("Pályakövetés az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Pozíció [m]")
	plot(time_mem[1:l], qN_mem[1:l], color = "red", label = "Nominális")
	plot(time_mem[1:l], q_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("allplots_vander.pdf")

	figure("Velocities")
	grid(true)
	title("Sebesség az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Sebesség [m/s]")
	plot(time_mem[1:l], qN_p[1:l], color = "red", label = "Nominális")
	plot(time_mem[1:l], q_p_mem[1:l], color = "green", label = "Megvalósult", linestyle = "--")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Accelerations")
	grid(true)
	title("Gyorsulás az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Gyorsulás [m/s²]")
	plot(time_mem[1:l], qN_pp[1:l], color = "red", label = "Nominális")
	plot(time_mem[1:l], q_pp_mem[1:l], color = "green", label = "Megvalósult", linestyle = "--")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Tracking_Error")
	title("Követési hiba az idő függvényében")
	grid(true)
	xlabel("Idő [s]")
	ylabel("Követési hiba [m]")
	plot(time_mem[1:l], qN_mem[1:l] - q_mem[1:l], color = "red")
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Control_Signal")
	title("Irányítójel az idő függvényében")
	grid(true)
	xlabel("Idő [s]")
	ylabel("Irányítójel [N]")
	plot(time_mem[1:l], u[1:l], color = "red")
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	figure("Phase_Space")
	title("Fázistér")
	xlabel("Pozíció [m]")
	ylabel("Sebesség [m/s]")
	grid(true)
	plot(qN_p[1:l], qN_mem[1:l], color = "red", label = "Nominális")
	plot(q_p_mem[1:l], q_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_vander.pdf")
	append_pdf!("allplots_vander.pdf", "temp_vander.pdf", cleanup = true)

	# date_string = Dates.format(now(), "mm-dd_HH-MM")
	# fileName = "./Plots/Vanderpool/" * date_string * "_$idx.pdf"
	# mv("allplots_vander.pdf", fileName)
	# if (Ploting == 0 && lastPlot)
	# 	close("all")
	# 	figure("trajectoria")
	# 	grid(true)
	# 	title("Trajectory  Run 1 Robust = $Robust \n (q₀=$(q[1]), q̇₀=$(q_p[1]))")
	# 	xlabel("Idő [s]")
	# 	ylabel("position")
	# 	plot(time_mem[1:l], qN_mem[1:l], color = "red", label = "Nominális")
	# 	plot(time_mem[1:l], q_mem[1:l], color = "green", linestyle = "--", label = "Megvalósult")
	# 	legend(loc = 1, borderaxespad = 0)
	# 	savefig("./Plots/Vanderpool/TEMP/traj_run_$idx.png")
	# 	show()
	# end
end

function run_simulation(q0, q_p0, q_pp0, init_Amp, init_ω, t_range)
	global time_mem = zeros(LONG)
	global q_mem = zeros(LONG)
	global q_p_mem = zeros(LONG)
	global q_pp_mem = zeros(LONG)
	global qN_mem = zeros(LONG)
	global qN_p_mem = zeros(LONG)
	global qN_pp_mem = zeros(LONG)
	global h_int = 0
	global u_mem = zeros(LONG)
	if (SingleRun == 1)
		q_mem[1] = q0
		q_p_mem[1] = q_p0
		q_pp_mem[1] = q_pp0
		singlerun(h_int, q_mem, q_p_mem, q_pp_mem, qN_p_mem, qN_pp_mem, u_mem, 1)
	else
		# Define initial conditions: (q0, q_p0)
		# condition_q = (+-)6.3:-0.1:(-+)6.3 ezeken belül stabil a basic beálításokkal
		condition_q = init_Amp .* sin.(init_ω .* t_range)
		condition_qp = init_Amp * init_ω .* cos.(init_ω .* t_range)
		for (idx, (q0, q_p0)) in enumerate(zip(condition_q, condition_qp))
			# Initialize memory arrays again for each run
			q_mem[1] = q0
			q_p_mem[1] = q_p0
			singlerun(h_int, q_mem, q_p_mem, q_pp_mem, qN_p_mem, qN_pp_mem, u_mem, idx, (idx == length(condition_q)))
		end
	end
	# log()
	if (Ploting == 1 && SingleRun == 1)
		show()
	end
end

##############################
# Define arrays for Plotting #
##############################

time_mem = zeros(LONG) # t

q_mem = zeros(LONG) #[q] A megvalósult pálya
q_p_mem = zeros(LONG)
q_pp_mem = zeros(LONG)

qN_mem = zeros(LONG) #[qN] A nominális pálya
qN_p_mem = zeros(LONG)
qN_pp_mem = zeros(LONG)

u_mem = zeros(LONG)  #[u] A szabályozó jel elmentése

qDes_pp_mem = zeros(LONG) #[qDes_pp] A PID korrekciós adatok
qDef_pp = zeros(LONG) #Az Adaptívan torzított jelre

h_int = 0 #A követési hiba integrálja
q_mem[1] = Amp * sin(ω * δt) # Kezdeti pozíció
q_p_mem[1] = Amp * ω * cos(ω * δt) # Kezdeti sebesség
q_pp_mem[1] = -Amp * ω^2 * sin(ω * δt)


# Setting initial condition
t_max = 20.0
init_ω = 0.5
init_Amp = 2

# időlépések
δranget = 1
t_range = 0:δranget:t_max

##############
# Simulation #
##############
using PyPlot
using PDFmerger
using Dates
run_simulation(q_mem[1], q_p_mem[1], q_pp_mem[1], init_Amp, init_ω, t_range)
# run_simulation(1, 5, 2, init_Amp, init_ω, t_range)

show()
