module DuffingModule
using PyPlot

export run_duffing

function run_duffing(q0=10.0, q_p0=10.0; do_plot=false)
    #############################################
    # Duffing Oscillator Controlled by RFPT #
    #############################################
    Adaptive = 1
    Robust = 1
    ######################
    # Control Parameters #
    ######################
    K = 1e5
    B = -1
    A = 1e-5

    #########
    # Time  #
    #########
    δt = 1e-3
    LONG = Int(2e4)
    l = LONG - 1

    ########################################
    # Kinematic Block Parameter (2nd order)#
    ########################################
    Λ = 1
    K_VSSM = 500
    w = 1
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


    ##############################
    # Define arrays for Plotting #
    ##############################
    #time
    t = zeros(LONG)

    q = zeros(LONG) #A megvalósult pálya
    q_p = zeros(LONG)
    q_pp = zeros(LONG)

    qN = zeros(LONG) #A nominális pálya
    qN_p = zeros(LONG)
    qN_pp = zeros(LONG)
    h = zeros(LONG)

    u = zeros(LONG) #A szanályozó jel elmentése

    qDes_pp = zeros(LONG)#A PID korrekciós adatok
    qDef_pp = zeros(LONG) #Az Adaptívan torzított jelre

    hint = 0 #A követési hiba integrálja
    #q[1]=Amp*sin(ω*δt)
    #q_p[1]=Amp*ω*cos(ω*δt)

    q[1] = q0
    q_p[1] = q_p0

    for i = 1:l
        #Compute the time in seconds
        t[i] = δt * i
        # Compute the Nominal trajectory.
        qN[i] = Amp * sin(ω * t[i])
        qN_p[i] = Amp * ω * cos(ω * t[i]) #_p is for d/dt 
        qN_pp[i] = -Amp * ω^2 * sin(ω * t[i])
        # Compute the Error.
        h[i] = qN[i] - q[i]
        h_p = qN_p[i] - q_p[i]

        if Robust == 1
            S = Λ^2 * hint + 2 * Λ * h[i] + h_p
            qDes_pp[i] = K_VSSM * tanh(S / w) + Λ^2 * h[i] + 2 * Λ * h_p + qN_pp[i]
        else
            qDes_pp[i] = qN_pp[i] + Λ^3 * hint + 3 * Λ^2 * h[i] + 3 * Λ * h_p
        end

        #Deformation
        if Adaptive == 1 && i > 3
            qDef_pp[i] = (qDef_pp[i-1] + K) * (1 + B * tanh(A * (q_pp[i-1] - qDes_pp[i]))) - K
        else
            qDef_pp[i] = qDes_pp[i]
        end #if
        #Compute the control signal
        u[i] = qDef_pp[i] - αₐ * q[i]
        # Compute the exact systems's respons
        q_pp[i] = αₑ * q[i] + δₑ * q_p[i] - βₑ * q[i]^3 + u[i]
        #Integrate back with Euler's method
        q_p[i+1] = q_p[i] + δt * q_pp[i]
        q[i+1] = q[i] + δt * q_p[i]
        hint = hint + δt * h[i]
    end# for

    if do_plot
        figure()
        grid(true)
        title("Duffing – Trajectory Tracking")
        xlabel("time")
        ylabel("position")
        plot(t[1:l], qN[1:l], label="nominal")
        plot(t[1:l], q[1:l], "r--", label="actual")
        legend()

        figure()
        grid(true)
        title("Duffing – Tracking Error")
        xlabel("time")
        ylabel("error")
        plot(t[1:l], qN[1:l] .- q[1:l])

        figure()
        grid(true)
        title("Duffing – Phase Space")
        xlabel("q")
        ylabel("q_p")
        plot(qN[1:l], qN_p[1:l])
        plot(q[1:l], q_p[1:l], "r--")

        figure()
        grid(true)
        title("Duffing – Control Signal")
        xlabel("time")
        ylabel("signal")
        plot(t[1:l], u[1:l])
    end

    maxhiba = maximum(abs.(h[1:l]))
    println("maxhiba :", maxhiba)
    return maxhiba
end # run_duffing

end # module DuffingModule
