using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("Qt5Agg")
using LinearAlgebra
Adaptive=1  #RFPT
Robust=1    #VSSM
######################
# Control Parameters #
######################
K=1e2
B=-1
A=1.97e-3


#########
# Time  #
#########
δt=1e-3
LONG=Int(2e4)
l=LONG-1

########################################
# Kinematic Block Parameter (2nd order)#
########################################
Λ=1
K_VSSM=50
w=1
####################################
# Parameters for Nominal Trajectory#
####################################


A₁=2
ω₁=0.5
A₂=3
ω₂=0.7
A₃=1
ω₃=1

##########################
# Exact Model Parameters #
##########################

βₑ=8/3
σₑ=5
ρₑ=40

################################
# Approximate Model Parameters #
################################

βₐ=7/3
σₐ=4
ρₐ=36

# The Original Deformation function for MIMO case
function G_MIMO(past_input, past_response, desired, error_limit, Kc, Bc, Ac) # Inputs and Outputs are vectors
	#  need control parameters  K B A
	error_norm=norm(past_response - desired, 2)
	# If the norm of the error is greater then the limnit compute the deformation
	# (it is not near the Fixed Point)
	if error_norm>error_limit
		e_direction=(past_response-desired)/error_norm
		B_factor = Bc * tanh(Ac * error_norm)
		out=(1+B_factor)*past_input + B_factor*Kc*e_direction
	else
		out=past_input # Almost in the Fixed Point
	end
	return out
end


# Time
t=zeros(LONG)

# Nominal Trajectory
xN=zeros(LONG)
xN_p=zeros(LONG)

yN=zeros(LONG)
yN_p=zeros(LONG)

zN=zeros(LONG)
zN_p=zeros(LONG)

# Desired
xDes_p=zeros(LONG)
yDes_p=zeros(LONG)
zDes_p=zeros(LONG)

# Deformed
xDef_p=zeros(LONG)
yDef_p=zeros(LONG)
zDef_p=zeros(LONG)

# Control Signal
u_x=zeros(LONG)
u_y=zeros(LONG)
u_z=zeros(LONG)


x_p=zeros(LONG)
y_p=zeros(LONG)
z_p=zeros(LONG)

x=zeros(LONG)
y=zeros(LONG)
z=zeros(LONG)

S_x=zeros(LONG)
S_p_x=zeros(LONG)
S_y=zeros(LONG)
S_p_y=zeros(LONG)
S_z=zeros(LONG)
S_p_z=zeros(LONG)

#initial conditions
hint_x=0
hint_y=0
hint_z=0

past_input=[0, 0, 0]
past_response=[0, 0, 0]
error_limit=1e-3

x[1]=A₁*sin(ω₁*δt)
y[1]=A₂*sin(ω₂*δt)
z[1]=A₃*sin(ω₃*δt)
x[1]=0
y[1]=0
z[1]=0
for i ∈ 1:l
	global hint_x
	global hint_y
	global hint_z
	global past_input
	global past_response
	global error_limit
	#Compute the time in seconds
	t[i]=δt*i
	# Compute the Nominal trajectory.
	xN[i]=A₁*sin(ω₁*t[i])
	xN_p[i]=A₁*ω₁*cos(ω₁*t[i])

	yN[i]=A₂*sin(ω₂*t[i])
	yN_p[i]=A₂*ω₂*cos(ω₂*t[i])

	zN[i]=A₃*sin(ω₃*t[i])
	zN_p[i]=A₃*ω₃*cos(ω₃*t[i])

	# Compute the Error.
	h_x=xN[i]-x[i]
	h_y=yN[i]-y[i]
	h_z=zN[i]-z[i]

	if Robust==1
		S_x[i]=Λ*hint_x+h_x
		S_y[i]=Λ*hint_y+h_y
		S_z[i]=Λ*hint_z+h_z
		xDes_p[i]=xN_p[i]+Λ*h_x+K_VSSM*tanh(S_x[i]/w)
		yDes_p[i]=yN_p[i]+Λ*h_y+K_VSSM*tanh(S_y[i]/w)
		zDes_p[i]=zN_p[i]+Λ*h_z+K_VSSM*tanh(S_z[i]/w)

		desired=[xDes_p[i], yDes_p[i], zDes_p[i]]
	else
		xDes_p[i]=Λ^2*hint_x+2*Λ*h_x+xN_p[i]
		yDes_p[i]=Λ^2*hint_y+2*Λ*h_y+yN_p[i]
		zDes_p[i]=Λ^2*hint_z+2*Λ*h_z+zN_p[i]

		desired=[xDes_p[i], yDes_p[i], zDes_p[i]]
	end

	#Deformation
	if Adaptive==1 && i>3
		past_input=G_MIMO(past_input, past_response, desired, error_limit, K, B, A)
	else
		past_input=desired
	end #if

	xDef_p[i]=past_input[1]
	yDef_p[i]=past_input[2]
	zDef_p[i]=past_input[3]
	#Control Signal
	u_x[i]=xDef_p[i] - σₐ * (y[i] - x[i])
	u_y[i]=yDef_p[i] - x[i] * (ρₐ - z[i]) + y[i]
	u_z[i]=zDef_p[i] - x[i] * y[i] + βₐ * z[i]

	#System
	x_p[i] = σₑ * (y[i] - x[i]) + u_x[i]
	y_p[i] = x[i] * (ρₑ - z[i]) - y[i] + u_y[i]
	z_p[i] = x[i] * y[i] - βₑ * z[i] + u_z[i]
	past_response = [x_p[i], y_p[i], z_p[i]]

	#Integrals
	x[i+1]=x[i]+δt*x_p[i]
	y[i+1]=y[i]+δt*y_p[i]
	z[i+1]=z[i]+δt*z_p[i]

	hint_x=hint_x+δt*h_x
	hint_y=hint_y+δt*h_y
	hint_z=hint_z+δt*h_z
end# for

figure(1)
grid(true)
title("traj. tracking")
plot(t[1:l], xN[1:l])
plot(t[1:l], yN[1:l])
plot(t[1:l], zN[1:l])
plot(t[1:l], x[1:l], "r--")
plot(t[1:l], y[1:l], "r--")
plot(t[1:l], z[1:l], "r--")

figure(2)
title("Tracking error")
grid(true)
plot(t[1:l], xN[1:l]-x[1:l])
plot(t[1:l], yN[1:l]-y[1:l])
plot(t[1:l], zN[1:l]-z[1:l])

figure(3)
title("Phase Space - X")
grid(true)
plot(xN[1:l], xN_p[1:l])
plot(x[1:l], x_p[1:l], "r--")


figure(4)
title("Phase Space - Y")
grid(true)
plot(yN[1:l], yN_p[1:l])
plot(y[1:l], y_p[1:l], "r--")

figure(5)
title("Phase Space - Z")
grid(true)
plot(zN[1:l], zN_p[1:l])
plot(z[1:l], z_p[1:l], "r--")

function my_max(arr)
	maxval = arr[1]
	for x in arr
		if x > maxval
			maxval = x
		end
	end
	return maxval
end

maxhiba_x = my_max(abs.(xN[1:l] .- x[1:l]))
maxhiba_y = my_max(abs.(yN[1:l] .- y[1:l]))
maxhiba_z = my_max(abs.(zN[1:l] .- z[1:l]))
maxhiba = my_max([maxhiba_x, maxhiba_y, maxhiba_z])
println("maxhiba_x :", maxhiba_x)
println("maxhiba_y :", maxhiba_y)
println("maxhiba_z :", maxhiba_z)
println("maxhiba   :", maxhiba)

show()