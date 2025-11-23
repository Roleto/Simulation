#############################################
# Duffing Oscillator Controlled by RFPT o VSSM #
#############################################
using PDFmerger
using Dates

######################
# Control Parameters #
######################
Adaptive = 1
Robust = 1
SingleRun = 1
#########
# Time  #
#########
δt = 1e-3
N = Int(2e4)
LONG = N - 1

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

function singlerun(q, q_p, ploting = true)
	q_mem[1] = q
	q_p_mem[1] = q_p
	max_index = 1
	h_int = 0
	for t ∈ 1:(LONG-1)
		time_mem[t] = t * δt

		#define the nominal trajectory
		qN_mem[t], qN_p_mem[t], qN_pp_mem[t] = nominalTraj(time_mem[t])

		# Compute the Errors
		h = qN_mem[t] - q_mem[t]
		h_p = qN_p_mem[t] - q_p_mem[t]

		# Finding max value
		h_max = qN_mem[max_index] - q_mem[max_index]
		if (h*h) < (h_max*h_max)
			max_index = t
		end

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
	if ploting
		# Plotting 
		l=LONG - l
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
	end

	date_string = Dates.format(now(), "mm-dd_HH-MM")
	random_hash = rand(1:100)
	fileName = "../Data/Duffing/" * date_string * "_$random_hash" * ".pdf"
	mv("allplots_duffing.pdf", fileName)
	return qN_mem[max_index] - q_mem[max_index]
end

"""
	init_duffing(;plot=true)

Allocate/initialize globals for Duffing simulation. Returns initial (q0, q_p0).
Set `Ploting` externally if desired. Use returned values to call `duffing_single_run`.
"""
function init_duffing(q0, q_p0; plot = true)
	global time_mem = zeros(LONG) #t
	global q_mem = zeros(LONG) # realized trajectory
	global q_p_mem = zeros(LONG)
	global q_pp_mem = zeros(LONG)
	global qN_mem = zeros(LONG) # nominal trajectory
	global qN_p_mem = zeros(LONG)
	global qN_pp_mem = zeros(LONG)
	global u_mem = zeros(LONG)  # control signal
	global qDes_pp_mem = zeros(LONG)
	global qDef_pp_mem = zeros(LONG)
	global h_int = 0
	# initial state from nominal trajectory at first time step
	# q_mem[1] = q0
	# q_p_mem[1] = q_p0
	q_mem[1] = Amp * sin(ω * δt)
	q_p_mem[1] = Amp * ω * cos(ω * δt)
	q_pp_mem[1] = -Amp * ω^2 * sin(ω * δt)
end

# Convenience wrapper to call original singlerun after init
function duffing_single_run(q0, q_p0; plot = true)
	init_duffing(q0, q_p0; plot = plot)
	max_error = singlerun(q0, q_p0, plot)
	if plot
		show()
	end
	return max_error
end

function duffing_single_run(; plot = true)
	q0=0
	q_p0=0
	return duffing_single_run(q0, q_p0; plot = plot)
end

# ---------------------------------------------------------------
# Archived original multi-run grid search (run_simulation)
# Keeping for future reference; can be re-enabled or refactored.
#
# function run_simulation(q_param, q_p_param, q_pp_param)
#   global time_mem = zeros(LONG)
#   global q_mem = zeros(LONG)
#   global q_p_mem = zeros(LONG)
#   global q_pp_mem = zeros(LONG)
#   global qN_mem = zeros(LONG)
#   global qN_p_mem = zeros(LONG)
#   global qN_pp_mem = zeros(LONG)
#   global h_int = 0
#   global u_mem = zeros(LONG)
#   position = zeros(441)
#   velocity = zeros(441)
#   error_max = zeros(441)
#   index = 1
#   for q0 in -1:0.1:1
#       for q_p0 in -1:0.1:1
#           position[index] = q0
#           velocity[index] = q_p0
#           error_max[index] = singlerun(q0, q_p0, false)
#           index += 1
#       end
#   end
#   figure("Error_Max")
#   grid(true)
#   title("Hiba Tér")
#   xlabel("Pozíció [m]")
#   ylabel("Sebesség [m/s]")
#   plot(position[1:441], error_max[1:441], color = "red", label = "Hiba")
#   legend(loc = 1, borderaxespad = 0)
# end
# ---------------------------------------------------------------

##############################
# Define arrays for Plotting #
##############################

# Initialization moved to init_duffing()

# Setting initial condition
t_max = 20.0
init_ω = 0.5
init_Amp = 2

# időlépések
δranget = 1
t_range = 0:δranget:t_max
