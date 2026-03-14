module VanDerPolModule
using PyPlot

export run_vanderpol

function run_vanderpol(q0=nothing, q_p0=nothing; do_plot=false)
#############################################
# Van der Pol Oscillator Controlled by RFPT #
#############################################
Adaptive=1
Robust=1
######################
# Control Parameters #
######################
K=1e5
B=-1
A=1e-5

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
K_VSSM=500
w=1
####################################
# Parameters for Nominal Trajectory#
####################################

ω=0.5
Amp=2
##########################
# Exact Model Parameters #
##########################

μₑ=0.4
ωₑ=0.46
αₑ=1
λₑ=0.1
mₑ=1

################################
# Approximate Model Parameters #
################################

μₐ=0.5
ωₐ=0.42
αₐ=0.9
λₐ=0.09
mₐ=0.8


##############################
# Define arrays for Plotting #
##############################
#time
t=zeros(LONG)

q=zeros(LONG) #A megvalósult pálya
q_p=zeros(LONG)
q_pp=zeros(LONG)

qN=zeros(LONG) #A nominális pálya
qN_p=zeros(LONG)
qN_pp=zeros(LONG)

u=zeros(LONG) #A szanályozó jel elmentése

qDes_pp=zeros(LONG)#A PID korrekciós adatok
qDef_pp=zeros(LONG) #Az Adaptívan torzított jelre

hint=0 #A követési hiba integrálja
q[1]   = (q0   === nothing) ? Amp*sin(ω*δt)   : q0
q_p[1] = (q_p0 === nothing) ? Amp*ω*cos(ω*δt) : q_p0

for i=1:l
      #Compute the time in seconds
      t[i]=δt*i
      # Compute the Nominal trajectory.
      qN[i]=Amp*sin(ω*t[i])
      qN_p[i]=Amp*ω*cos(ω*t[i]) #_p is for d/dt 
      qN_pp[i]=-Amp*ω^2*sin(ω*t[i])
      # Compute the Error.
      local h=qN[i]-q[i]
      local h_p=qN_p[i]-q_p[i]
    
      if Robust==1
        S=Λ^2*hint+2*Λ*h+h_p
        qDes_pp[i]=K_VSSM*tanh(S/w)+Λ^2*h+2*Λ*h_p+qN_pp[i]
      else
        qDes_pp[i]=qN_pp[i]+Λ^3*hint+3*Λ^2*h+3*Λ*h_p
      end
      
      #Deformation
      if Adaptive==1 && i>3
          qDef_pp[i]=(qDef_pp[i-1]+K)*(1+B*tanh(A*(q_pp[i-1]-qDes_pp[i])))-K
      else
          qDef_pp[i]=qDes_pp[i]
      end #if
      #Compute the control signal
      u[i]=mₐ*qDef_pp[i]#-μₐ*(1-q[i]^2)*q_p[i]+ωₐ^2*q[i]+αₐ*q[i]^3+λₐ*q[i]^5
      # Compute the exact systems's respons
      q_pp[i]=(u[i]+μₑ*(1-q[i]^2)*q_p[i]-ωₑ^2*q[i]-αₑ*q[i]^3-λₑ*q[i]^5)/mₑ
      #Integrate back with Euler's method
      q_p[i+1]=q_p[i]+δt*q_pp[i]
      q[i+1]=q[i]+δt*q_p[i]
      hint=hint+δt*h
end# for

    if do_plot
        figure()
        grid(true)
        title("VanDerPol – Trajectory Tracking")
        xlabel("time")
        ylabel("position")
        plot(t[1:l], qN[1:l], label="nominal")
        plot(t[1:l], q[1:l], "r--", label="actual")
        legend()

        figure()
        grid(true)
        title("VanDerPol – Tracking Error")
        xlabel("time")
        ylabel("error")
        plot(t[1:l], qN[1:l] .- q[1:l])

        figure()
        grid(true)
        title("VanDerPol – Phase Space")
        xlabel("q")
        ylabel("q_p")
        plot(qN[1:l], qN_p[1:l])
        plot(q[1:l], q_p[1:l], "r--")

        figure()
        grid(true)
        title("VanDerPol – Control Signal")
        xlabel("time")
        ylabel("signal")
        plot(t[1:l], u[1:l])
    end

    maxhiba = maximum(abs.(qN[1:l] .- q[1:l]))
    println("maxhiba :", maxhiba)
    return maxhiba
end # run_vanderpol

end # module VanDerPolModule