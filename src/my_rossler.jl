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
using PyPlot
using Pkg
PyPlot.matplotlib.use("TkAgg")
Pkg.activate(".")
using LinearAlgebra
using PyPlot

Adaptive = 1  #RFPT
Robust = 1    #VSSM

#################
# Time variable #
#################
δt = 1e-3
LONG = Int(2e4)
l = LONG - 1

#################
# Control Parameters #
#################
K = 1e2
B = -1
A = 1.97e-3

########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ = 1
K_VSSM = 50
w = 1

####################################
# Parameters for Nominal Trajectory#
####################################
A₁ = 2
ω₁ = 0.5
A₂ = 3
ω₂ = 0.7
A₃ = 1
ω₃ = 1

##########################
# Exact model parameters #
##########################
aₑ = 0.01
bₑ = 0.2
cₑ = 5.7

############################
# Approx. model parameters #
############################
aₐ = 0.1
bₐ = 0.3
cₐ = 5.5

function G_MIMO(past_input, past_response, desired, err_limit, K, B, A)
	Amatr_h = (past_response - desired)
	error_norm = norm(Amatr_h, 2)
	if error_norm > err_limit
		e_direction = Amatr_h / error_norm
		B_factor = B * tanh(A * error_norm)
		G = (1 + B_factor) * past_input + B_factor * K * e_direction
	else
		G = past_input
	end
	return G
end

ErrorMetric(hint, h) = Λ * hint + h

function KinBlock(S, h, qN_p)
	s = K * [tanh(S[1] / w), tanh(S[2] / w), tanh(S[3] / w)]
	qn = qN_p * Λ
	x = qn[1] * h[1] * s[1]
	y = qn[2] * h[2] * s[2]
	z = qn[3] * h[3] * s[3]
	return [x, y, z]
end



time_mem = zeros(LONG)

# Nominal Trajectory
# qN = zeros(Float64, LONG, 3)
# qN_p = zeros(Float64, LONG, 3)

xN = zeros(LONG)
xN_p = zeros(LONG)

yN = zeros(LONG)
yN_p = zeros(LONG)

zN = zeros(LONG)
zN_p = zeros(LONG)

# Desired
# qDes_p = zeros(Float64, LONG, 3)
xDes_p = zeros(LONG)
yDes_p = zeros(LONG)
zDes_p = zeros(LONG)

# Deformed
# qDef_p = zeros(Float64, LONG, 3)
xDef_p = zeros(LONG)
yDef_p = zeros(LONG)
zDef_p = zeros(LONG)


# Control Signal
# u = zeros(Float64, LONG, 3)
u_x = zeros(LONG)
u_y = zeros(LONG)
u_z = zeros(LONG)

# q_p = zeros(Float64, LONG, 3)
x_p = zeros(LONG)
y_p = zeros(LONG)
z_p = zeros(LONG)

# qA = zeros(Float64, LONG, 3)
x = zeros(LONG)
y = zeros(LONG)
z = zeros(LONG)

S_x = zeros(LONG)
S_p_x = zeros(LONG)
S_y = zeros(LONG)
S_p_y = zeros(LONG)
S_z = zeros(LONG)
S_p_z = zeros(LONG)

#initial conditions
hint_x = 0
hint_y = 0
hint_z = 0

past_input = zeros(3)
past_response = zeros(3)
past_responses = zeros(Float64, LONG, 3)

# errors = zeros(Float64, LONG, 3)
error_limit = 1e-3

x[1] = A₁ * sin(ω₁ * δt)
y[1] = A₂ * sin(ω₂ * δt)
z[1] = A₃ * sin(ω₃ * δt)
x_p[1] = A₁ * ω₁ * cos(ω₁ * δt)
y_p[1] = A₂ * ω₂ * cos(ω₂ * δt)
z_p[1] = A₃ * ω₃ * cos(ω₃ * δt)
# x_pp[1] = -A₁ * ω₁^2 * sin(ω₁ * δt)
# y_pp[1] = -A₂ * ω₂^2 * sin(ω₂ * δt)
# z_pp[1] = -A₃ * ω₃^2 * sin(ω₃ * δt)

println("Simulation started")
for t ∈ 1:l
	global hint_x
	global hint_y
	global hint_z
	global past_input
	global past_response
	global error_limit

	time_mem[t] = δt * t
	#Nominal trajectory for the actual time frame
	# qN[t, :] = [Amp[j] * sin(ω[j] * time_mem[t]) for j ∈ 1:3]
	xN[t] = A₁ * sin(ω₁ * time_mem[t])
	xN_p[t] = A₁ * ω₁ * cos(ω₁ * time_mem[t])

	yN[t] = A₂ * sin(ω₂ * time_mem[t])
	yN_p[t] = A₂ * ω₂ * cos(ω₂ * time_mem[t])

	zN[t] = A₃ * sin(ω₃ * time_mem[t])
	zN_p[t] = A₃ * ω₃ * cos(ω₃ * time_mem[t])

	h_x = xN[t] - x[t]
	h_y = yN[t] - y[t]
	h_z = zN[t] - z[t]

	if Robust == 1
		S = ErrorMetric([hint_x, hint_y, hint_z], [h_x, h_y, h_z])

		S_x[t] = S[1]
		S_y[t] = S[2]
		S_z[t] = S[3]

		xDes_p[t] = xN_p[t] + Λ * h_x + K_VSSM * tanh(S_x[t] / w)
		yDes_p[t] = yN_p[t] + Λ * h_y + K_VSSM * tanh(S_y[t] / w)
		zDes_p[t] = zN_p[t] + Λ * h_z + K_VSSM * tanh(S_z[t] / w)


		desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
	else
		xDes_p[t] = Λ^2 * hint_x + 2 * Λ * h_x + xN_p[t]
		yDes_p[t] = Λ^2 * hint_y + 2 * Λ * h_y + yN_p[t]
		zDes_p[t] = Λ^2 * hint_z + 2 * Λ * h_z + zN_p[t]

		desired = [xDes_p[t], yDes_p[t], zDes_p[t]]
	end

	# Deformation
	if Adaptive == 1 && t > 3
		past_input = G_MIMO(past_input, past_response, desired, error_limit, K, B, A)
	else
		past_input = desired
	end

	xDef_p[t] = past_input[1]
	yDef_p[t] = past_input[2]
	zDef_p[t] = past_input[3]

	#Control Signal                                                                                                                                      

	u_x[t] = xDef_p[t] + y[t] + z[t]
	u_y[t] = yDef_p[t] - x[t] - aₐ * y[t]
	u_z[t] = zDef_p[t] - bₐ - z[t] * (x[t] - cₐ)

	#modell
	x_p[t] = -y[t] - z[t] + u_x[t]
	y_p[t] = x[t] + aₑ * y[t] + u_y[t]
	z_p[t] = bₑ + z[t] * (x[t] - cₑ) + u_z[t]

	past_response = [x_p[t], y_p[t], z_p[t]]
	past_responses[t, :] = past_response

	#Integrals
	x[t+1] = x[t] + δt * x_p[t]
	y[t+1] = y[t] + δt * y_p[t]
	z[t+1] = z[t] + δt * z_p[t]

	hint_x = hint_x + δt * h_x
	hint_y = hint_y + δt * h_y
	hint_z = hint_z + δt * h_z
end
println("Simulation Ended")

############
# Plotting #
############

using PDFmerger


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

legend(loc = "lower right", fancybox = "True")
tight_layout()
subplots_adjust(hspace = 0.2)
fig[:canvas][:draw]()
savefig("allplots_rossler.pdf")

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
title(L"$1^{st} Time Derivatives$")

subplot(311)
ax1 = gca()
grid1 = grid(true)
ylabel(L"$\dot{x}$")

plot(time_mem[1:l], xN_p[1:l], color = "red", linewidth = 2, label = L"\dot{x}^{N}")
plot(time_mem[1:l], past_responses[1:l, 1], color = "green", linewidth = 3, label = L"\dot{x}", linestyle = "--")

subplot(312, sharex = ax1)
ax2 = gca()
grid(true)
ylabel(L"$\dot{y}$")

plot(time_mem[1:l], yN_p[1:l], color = "red", linewidth = 2, label = L"\dot{y}^{N}")
plot(time_mem[1:l], past_responses[1:l, 2], color = "green", linewidth = 3, label = L"\dot{y}", linestyle = "--")

subplot(313, sharex = ax2)
ax3 = gca()
grid(true)
xlabel("t, [s]")
ylabel(L"$\dot{z}$")
plot(time_mem[1:l], zN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
plot(time_mem[1:l], past_responses[1:l, 3], color = "green", linewidth = 3, label = "Realized", linestyle = "--")
legend(loc = "lower left", fancybox = "True")

tight_layout()
subplots_adjust(hspace = 0.0)
fig[:canvas][:draw]()

savefig("temp_rossler.pdf")
append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)


#######################
# Acceleration Plotting #
#######################
# fig_caption = "accelerations"
# fig = figure(fig_caption)
# grid(true)
# title(L"$2^{st} Time Derivatives$")

# subplot(311)
# ax1 = gca()
# grid1 = grid(true)
# ylabel(L"$\ddot{x}$")

# plot(time_mem[1:l], xN_pp[1:l], color = "red", linewidth = 2, label = L"\ddot{x}^{N}")
# plot(time_mem[1:l], accelerations[1:l, 1], color = "green", linewidth = 3, label = L"\ddot{x}", linestyle = "--")

# subplot(312, sharex = ax1)
# ax2 = gca()
# grid(true)
# ylabel(L"$\dot{y}$")

# plot(time_mem[1:l], yN_pp[1:l], color = "red", linewidth = 2, label = L"\ddot{y}^{N}")
# plot(time_mem[1:l], accelerations[1:l, 2], color = "green", linewidth = 3, label = L"\ddot{y}", linestyle = "--")

# subplot(313, sharex = ax2)
# ax3 = gca()
# grid(true)
# xlabel("t, [s]")
# ylabel(L"$\ddot{z}$")
# plot(time_mem[1:l], zN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
# plot(time_mem[1:l], accelerations[1:l, 3], color = "green", linewidth = 3, label = "Realized", linestyle = "--")
# plot(time_mem[1:l], past_inputs_p[1:l, 3], color = "blue", label = "Desired", linestyle = "-.")
# legend(loc = "lower left", fancybox = "True")

# tight_layout()
# subplots_adjust(hspace = 0.0)
# fig[:canvas][:draw]()
# savefig("temp_rossler.pdf")
# append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)


# #######################
# # Acceleration Plotting #
# #######################
# fig_caption = "accelerations"
# fig = figure(fig_caption)
# grid(true)
# title(L"$2^{st} Time Derivatives$")

# subplot(311)
# ax1 = gca()
# grid1 = grid(true)
# ylabel(L"$\ddot{x}$")

# plot(time_mem[1:l], xN_pp[1:l], color = "red", linewidth = 2, label = L"\ddot{x}^{N}")
# plot(time_mem[1:l], accelerations[1:l, 1], color = "green", linewidth = 3, label = L"\ddot{x}", linestyle = "--")

# subplot(312, sharex = ax1)
# ax2 = gca()
# grid(true)
# ylabel(L"$\dot{y}$")

# plot(time_mem[1:l], yN_pp[1:l], color = "red", linewidth = 2, label = L"\ddot{y}^{N}")
# plot(time_mem[1:l], accelerations[1:l, 2], color = "green", linewidth = 3, label = L"\ddot{y}", linestyle = "--")

# subplot(313, sharex = ax2)
# ax3 = gca()
# grid(true)
# xlabel("t, [s]")
# ylabel(L"$\ddot{z}$")
# plot(time_mem[1:l], zN_pp[1:l], color = "red", linewidth = 2, label = "Nominális")
# plot(time_mem[1:l], accelerations[1:l, 3], color = "green", linewidth = 3, label = "Realized", linestyle = "--")

# legend(loc = "lower left", fancybox = "True")

# tight_layout()
# subplots_adjust(hspace = 0.0)
# fig[:canvas][:draw]()
# savefig("temp_rossler.pdf")
# append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)


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

savefig("temp_rossler.pdf")
append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

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
tight_layout()
fig[:canvas][:draw]()
savefig("temp_rossler.pdf")
append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

###############################
# Phase Trajectories Plotting #
###############################

fig_caption = "phase_trajectories_x"
figure(fig_caption)
grid(true)
title("Fázis tér X irányban")
xlabel("Pozíció [m]")
ylabel("Sebesség [m/s]")
plot(xN[1:l], xN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
plot(x[1:l], x_p[1:l], color = "green", linestyle = "--", linewidth = 2.5, label = "Megvalósult")
legend(loc = "lower left", fancybox = "True")
tight_layout()

savefig("temp_rossler.pdf")
append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

fig_caption = "phase_trajectories_y"
figure(fig_caption)
grid(true)
title("Fázis tér Y irányban")
xlabel("Pozíció [m]")
ylabel("Sebesség [m/s]")
plot(yN[1:l], yN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
plot(y[1:l], y_p[1:l], color = "green", linewidth = 2.75, linestyle = "--", label = "Megvalósult")
legend(loc = "upper left", fancybox = "True")
tight_layout()
savefig("temp_rossler.pdf")
append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

fig_caption = "phase_trajectories_z"
figure(fig_caption)
grid(true)
title("Fázis tér Z irányban")
xlabel("Pozíció [m]")
ylabel("Sebesség [m/s]")
plot(zN[1:l], zN_p[1:l], color = "red", linewidth = 2, label = "Nominális")
plot(z[1:l], z_p[1:l], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
legend(loc = "upper left", fancybox = "True")
tight_layout()

savefig("temp_rossler.pdf")
append_pdf!("allplots_rossler.pdf", "temp_rossler.pdf", cleanup = true)

############
# Log file #
############
using Dates

# date_string = Dates.format(now(), "mm-dd_HH-MM")
# fileName = "./Plots/Rossler/" * date_string * ".pdf"
# mv("allplots_rossler.pdf", fileName)
# file = open("./Plots/Rossler/log.txt", "a")
# line_breaker = "\n####################################################\n"

# file_text = string(line_breaker, date_string, " Következö paraméterekkel volt használva:\nControl Params:\n",
# 	"K= ", K, "\tB= ", B, "\tA= ", A, "\tw= ", w, "\tΛ= ", Λ,
# 	"\nApproximate Model Parameters:\naₐ=", aₐ, "\tbₐ=", bₐ, "\tcₐ=", cₐ,
# 	"\nTime variable :\nδt=", δt, "\t=LONG", LONG, "\n",
# 	"\nExact Model Parameters:\naₑ=", aₑ, "\tbₑ=", bₑ, "\tcₑ=", cₑ,
# 	"\nNominal Trajectory Parameters:\nω=", [ω₁, ω₂, ω₃], "\nAmp=", [A₁, A₂, A₃], line_breaker)
# write(file, file_text)
# close(file)


# fig_caption = "nominal_realized_trajectories"
# fig = figure(fig_caption)
# grid("True")
# title("Nominal and Realized Trajectories")

# subplot(311)
# ax1 = gca()
# grid1 = grid("True")
# ylabel(L"x")
# plot(time_mem[3:LONG-1], qN[3:LONG-1, 1], color = "red", linewidth = 2, label = "Nominális")
# plot(time_mem[3:LONG-1], qA[3:LONG-1, 1], color = "green", linewidth = 2.5, label = "Megvalósult", linestyle = "--")
# plot(time_mem[3:LONG-1], qA_pid[3:LONG-1, 1], color = "blue", label = "Realized Non-adaptive", linestyle = "-.")

# subplot(312, sharex = ax1)
# ax2 = gca()
# grid("True")
# ylabel(L"y")
# plot(time_mem[3:LONG-1], qN[3:LONG-1, 2], color = "red", linewidth = 2, label = "Nominális")
# plot(time_mem[3:LONG-1], qA[3:LONG-1, 2], color = "green", linewidth = 2.5, label = "Megvalósult", linestyle = "--")
# plot(time_mem[3:LONG-1], qA_pid[3:LONG-1, 2], color = "blue", label = "Realized Non-adaptive", linestyle = "-.")

# subplot(313, sharex = ax2)
# ax3 = gca()
# grid("True")
# xlabel(L"t, [$s$]")
# ylabel(L"z")
# plot(time_mem[3:LONG-1], qN[3:LONG-1, 3], color = "red", linewidth = 2, label = "Nominális")
# plot(time_mem[3:LONG-1], qA[3:LONG-1, 3], color = "green", linewidth = 2.5, label = "Megvalósult", linestyle = "--")
# plot(time_mem[3:LONG-1], qA_pid[3:LONG-1, 3], color = "blue", label = "Realized Non-adaptive", linestyle = "-.")

# legend()
# tight_layout()
# subplots_adjust(hspace = 0.0)
# fig[:canvas][:draw]()
# savefig("allplots.pdf")

# fig_caption = "tracking_error"
# fig = figure(fig_caption)
# grid("True")
# title("Tracking Errors vs Time")

# subplot(311)
# ax1 = gca()
# grid1 = grid("True")
# ylabel(L"$x^N-x$")
# plot(time_mem[3:LONG-1], errors[3:LONG-1, 1], color = "red", linewidth = 2, label = "Adaptive")
# plot(time_mem[3:LONG-1], errors_pid[3:LONG-1, 1], color = "green", linewidth = 2.5, linestyle = "--", label = "Non-adaptive")

# subplot(312, sharex = ax1)
# ax2 = gca()
# grid("True")
# ylabel(L"$y^N-y$")
# plot(time_mem[3:LONG-1], errors[3:LONG-1, 2], color = "red", linewidth = 2, label = "Adaptive")
# plot(time_mem[3:LONG-1], errors_pid[3:LONG-1, 2], color = "green", linewidth = 2.5, linestyle = "--", label = "Non-adaptive")

# subplot(313, sharex = ax2)
# ax3 = gca()
# grid("True")
# xlabel(L"t, [$s$]")
# ylabel(L"$z^N-z$")
# plot(time_mem[3:LONG-1], errors[3:LONG-1, 3], color = "red", linewidth = 2, label = "Adaptive")
# plot(time_mem[3:LONG-1], errors_pid[3:LONG-1, 3], color = "green", linewidth = 2.5, linestyle = "--", label = "Non-adaptive")

# legend()
# tight_layout()
# subplots_adjust(hspace = 0.190)
# fig[:canvas][:draw]()

# fig_caption = "phase_trajectories_x"
# figure(fig_caption)
# grid("True")
# title("Phase Trajectories")
# xlabel(L"x")
# ylabel(L"$\dot{x}$")
# plot(qN[3:LONG-1, 1], qN_p[3:LONG-1, 1], color = "red", linewidth = 2, label = "Nominális")
# plot(qA[3:LONG-1, 1], past_responses[3:LONG-1, 1], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
# legend(loc = "lower left", fancybox = "True")
# tight_layout()

# savefig("temp.pdf")
# append_pdf!("allplots.pdf", "temp.pdf", cleanup = true)

# fig_caption = "phase_trajectories_y"
# figure(fig_caption)
# grid("True")
# title("Phase Trajectories")
# xlabel(L"y")
# ylabel(L"$\dot{y}$")
# plot(qN[3:LONG-1, 2], qN_p[3:LONG-1, 2], color = "red", linewidth = 2, label = "Nominális")
# plot(qA[3:LONG-1, 2], past_responses[3:LONG-1, 2], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
# legend(loc = "upper left", fancybox = "True")
# tight_layout()

# fig_caption = "phase_trajectories_z"
# figure(fig_caption)
# grid("True")
# title("Phase Trajectories")
# xlabel(L"z")
# ylabel(L"$\dot{z}$")
# plot(qN[3:LONG-1, 3], qN_p[3:LONG-1, 3], color = "red", linewidth = 2, label = "Nominális")
# plot(qA[3:LONG-1, 3], past_responses[3:LONG-1, 3], color = "green", linewidth = 2.5, linestyle = "--", label = "Megvalósult")
# legend(loc = "upper left", fancybox = "True")
# tight_layout()

# fig_caption = "control_signal"
# figure(fig_caption)
# grid("True")
# title("Control Signals vs Time")
# xlabel(L"t, $[s]$")
# ylabel(L"u")
# plot(time_mem[3:LONG-1], u[3:LONG-1, 1], color = "red", label = L"$u_1$")
# plot(time_mem[3:LONG-1], u[3:LONG-1, 2], color = "green", label = L"$u_1$")
# plot(time_mem[3:LONG-1], u[3:LONG-1, 3], color = "blue", label = L"$u_3$")
# legend(loc = "lower left", fancybox = "True")
# tight_layout()


# fig_caption = "velocities"
# fig = figure(fig_caption)
# grid("True")
# title(L"$1^{st} Time Derivatives$")

# subplot(311)
# ax1 = gca()
# grid1 = grid("True")
# ylabel(L"$\dot{x}$")
# plot(time_mem[3:LONG-1], qN_p[3:LONG-1, 1], color = "red", linewidth = 2, label = L"\dot{x}^{N}")
# plot(time_mem[3:LONG-1], past_responses[3:LONG-1, 1], color = "green", linewidth = 3, label = L"\dot{x}", linestyle = "--")
# plot(time_mem[3:LONG-1], past_inputs[3:LONG-1, 1], color = "blue", label = L"\dot{x}^{Des}", linestyle = "-.")

# subplot(312, sharex = ax1)
# ax2 = gca()
# grid("True")
# ylabel(L"$\dot{y}$")
# plot(time_mem[3:LONG-1], qN_p[3:LONG-1, 2], color = "red", linewidth = 2, label = L"\dot{y}^{N}")
# plot(time_mem[3:LONG-1], past_responses[3:LONG-1, 2], color = "green", linewidth = 3, label = L"\dot{y}", linestyle = "--")
# plot(time_mem[3:LONG-1], past_inputs[3:LONG-1, 2], color = "blue", label = L"\dot{y}^{Des}", linestyle = "-.")

# subplot(313, sharex = ax2)
# ax3 = gca()
# grid("True")
# xlabel("t, [s]")
# ylabel(L"$\dot{z}$")
# plot(time_mem[3:LONG-1], qN_p[3:LONG-1, 3], color = "red", linewidth = 2, label = "Nominális")
# plot(time_mem[3:LONG-1], past_responses[3:LONG-1, 3], color = "green", linewidth = 3, label = "Realized", linestyle = "--")
# plot(time_mem[3:LONG-1], past_inputs[3:LONG-1, 3], color = "blue", label = "Desired", linestyle = "-.")
# legend(loc = "lower left", fancybox = "True")

# tight_layout()
# subplots_adjust(hspace = 0.0)
# fig[:canvas][:draw]()

show()


