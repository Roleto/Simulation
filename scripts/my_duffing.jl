#############################################
# Duffing Oscillator Controlled by RFPT o VSSM #
#############################################
using Pkg
# Pkg.add("Tk")
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("Qt5Agg")
######################
# Control Parameters #
######################
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
αₑ = 1
δₑ = 0.2
βₑ = 1

################################
# Approximate Model Parameters #
################################
αₐ = 0.8
δₐ = 0.1
βₐ = 0.9

Exact(q, q_p, u) = αₑ * q + δₑ * q_p + βₑ * q^3 + u

Approx(q, q_p, q_Np) = q_Np - αₐ * q - δₐ * q_p - βₐ * q^3 #mₐ*q_Np[t]#-μₐ*(1-q^2)*q_p+ωₐ^2*q+αₐ*q^3+λₐ*q^5
ErrorMetrics(h_int, h, h_p) = Λ^2 * h_int + 2 * Λ * h + h_p

KinBlock(S, h, h_p, qN_pp) = K_VSSM * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + qN_pp

G(past_input, past_response, xDnow) = (K + past_input) * (1 + B * tanh(A * (past_response - xDnow))) - K

function nominalTraj(t)
	qN = Amp * sin(ω * t)
	q_pN = ω * Amp * cos(ω * t)
	q_ppN = -Amp * ω^2 * sin(ω * t)
	return qN, q_pN, q_ppN
end

function log()
	date_string = Dates.format(now(), "mm-dd_HH-MM")
	line_breaker = "\n####################################################\n"
	isAdaptiv = Bool(1 == Adaptive)
	isRobust = Bool(1 == Robust)

	file = open("./Plots/Duffing/log.txt", "a")
	file_text = string(line_breaker, date_string, " Következö paraméterekkel volt használva:\nAdaptiv=", isAdaptiv, "Robust=", isRobust,
		"\nControl Params:\n",
		"K=", K, "\tB=", B, "\tA=", A, "\tw=", w, "\tΛ=", Λ,
		"\nTime variable :\nδt=", δt, "\t=LONG", LONG, "\n",
		"\nApproximate Model Parameters:\nαₐ=", αₐ, "\tδₐ=", δₐ, "\tβₐ=", βₐ,
		"\nExact Model Parameters:\nαₑ=", αₑ, "\tδₑ=", δₑ, "\tβₑ=", βₑ,
		"\nNominal Trajectory Parameters:\nω=", ω, "\tAmp=", Amp, line_breaker,
	)

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
		qN_mem[t], qN_p_mem[t], qN_pp_mem[t] = nominalTraj(time_mem[t])

		# Compute the Errors
		h = qN_mem[t] - q_mem[t]
		h_p = qN_p_mem[t] - q_p_mem[t]

		if Robust == 1
			S = ErrorMetrics(h_int, h, h_p)
			qDes_pp_mem[t] = KinBlock(S, h, h_p, qN_pp_mem[t])
		else
			qDes_pp_mem[t] = qN_pp_mem[t] + Λ^3 * h_int + 3 * Λ^2 * h + 3 * Λ * h_p
		end

		if Adaptive == 1 && t > 3
			qDef_pp_mem[t] = G(qDef_pp_mem[t-1], q_pp_mem[t-1], qDes_pp_mem[t])
		else
			qDef_pp_mem[t] = qDes_pp_mem[t]
		end

		u_mem[t] = Approx(q_mem[t], q_p_mem[t], qDef_pp_mem[t])

		# Compute the exact systems's respons
		q_pp_mem[t] = Exact(q_mem[t], q_p_mem[t], u_mem[t])
		q_pp = αₑ * q_mem[t] + δₑ * q_p_mem[t] + βₑ * q_mem[t]^3 + u_mem[t]

		#Integrate back with Euler's method
		q_p_mem[t+1] = q_p_mem[t] + δt * q_pp_mem[t]
		q_mem[t+1] = q_mem[t] + δt * q_p_mem[t]
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
	savefig("allplots_duffing.pdf")

	#velocity gyorsan kell

	figure("Acceleration")
	grid(true)
	title("Accelerations vs Time")
	xlabel("Time [s]")
	ylabel("Acceleration [m/s²]")
	# ylabel(L"Acceleration $\frac{m}{s^2}$")
	plot(time_mem[1:l], qN_pp_mem[1:l], color = "#7684FF", label = "Névleges")
	plot(time_mem[1:l], q_pp_mem[1:l], color = "#FFAA41", label = "Realizált", linestyle = "--")
	plot(time_mem[1:l], qDes_pp_mem[1:l], color = "#2A40FF", label = "Desired", linestyle = "-.")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_duffing.pdf")
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)

	figure("Tracking_Error")
	title("Követési hiba az idő függvényében")
	grid(true)
	xlabel("Idő [s]")
	ylabel("Követési hiba [m]")
	plot(time_mem[1:l], qN_mem[1:l] - q_mem[1:l], color = "red")
	savefig("temp_duffing.pdf")
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)

	figure("Control_Signal")
	grid(true)
	title("Irányítójel az idő függvényében")
	xlabel("Idő [s]")
	ylabel("Irányítójel [N]")
	plot(time_mem[1:l], u_mem[1:l], color = "red")
	savefig("temp_duffing.pdf")
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)

	figure("Phase_Space")
	grid(true)
	title("Fázistér")
	xlabel("Pozíció [m]")
	ylabel("Sebesség [m/s]")
	plot(qN_mem[1:l, :], qN_p_mem[1:l, :], color = "red", label = "Nominális")
	plot(q_mem[1:l, :], q_p_mem[1:l, :], color = "green", linestyle = "--", label = "Megvalósult")
	legend(loc = 1, borderaxespad = 0)
	savefig("temp_duffing.pdf")
	append_pdf!("allplots_duffing.pdf", "temp_duffing.pdf", cleanup = true)

	date_string = Dates.format(now(), "mm-dd_HH-MM")
	fileName = "./Plots/Duffing/" * date_string * "_$idx.pdf"
	mv("allplots_duffing.pdf", fileName)
	if (lastPlot && Ploting == 0)
		close("all")
		figure("trajectoria")
		grid(true)
		title("Trajectory  Run 1 Robust = $Robust \n (q₀=$(q[1]), q̇₀=$(q_p[1]))")
		xlabel("Time [s]")
		ylabel("position")
		plot(time_mem[1:l], qN_mem[1:l], color = "red", label = "Névleges")
		plot(time_mem[1:l], q_mem[1:l], color = "green", linestyle = "--", label = "Realizált")
		legend(loc = 1, borderaxespad = 0)
		savefig("./Plots/Duffing/TEMP/traj_run_$idx.png")
	end
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
		condition_q = init_Amp .* sin.(init_ω .* t_range)
		condition_qp = init_Amp * init_ω .* cos.(init_ω .* t_range)
		for (idx, (q0, q_p0)) in enumerate(zip(condition_q, condition_qp))
			# Initialize memory arrays again for each run
			q_mem[1] = q0
			q_p_mem[1] = q_p0
			singlerun(h_int, q_mem, q_p_mem, q_pp_mem, qN_p_mem, qN_pp_mem, u_mem, idx, (idx == length(condition_q)))
		end
	end
	log()
	if (Ploting == 1 && SingleRun == 1)
		show()
	end
end


##############################
# Define arrays for Plotting #
##############################

time_mem = zeros(LONG) #t

q_mem = zeros(LONG) #A megvalósult pálya
q_p_mem = zeros(LONG)
q_pp_mem = zeros(LONG)

qN_mem = zeros(LONG) #A nominális pálya
qN_p_mem = zeros(LONG)
qN_pp_mem = zeros(LONG)

u_mem = zeros(LONG)  #A szabályozó jel

qDes_pp_mem = zeros(LONG)
qDef_pp_mem = zeros(LONG)

h_int = 0
q_mem[1] = 10.0
q_p_mem[1] = 1.0
q_pp_mem[1] = -Amp * ω^2 * sin(ω * δt)

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
run_simulation(q_mem[1], q_p_mem[1], q_pp_mem[1], init_Amp, init_ω, t_range)
# run_simulation(1, 5, 2, init_Amp, init_ω, t_range)



show()
