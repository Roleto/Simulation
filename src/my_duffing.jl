module DuffingModule
using PyPlot

export run_duffing

function run_duffing(q0=2.0 * sin(0.5 * 1e-3), qp0=2.0 * 0.5 * cos(0.5 * 1e-3);
    Adaptive=1,
    Robust=1,
    K=1e5,
    B=-1.0,
    A=1e-5,
    δt=1e-3,
    LONG=Int(2e4),
    Λ=1,
    K_VSSM=500,
    w=1,
    ω=0.5,
    Amp=2,
    αₑ=1.0,
    δₑ=0.2,
    βₑ=1.0,
    αₐ=0.8,
    δₐ=0.1,
    βₐ=0.9,
    do_plot=false)

    l = LONG - 1
    t = zeros(LONG)
    q = zeros(LONG)
    q_p = zeros(LONG)
    q_pp = zeros(LONG)
    qN = zeros(LONG)
    qN_p = zeros(LONG)
    qN_pp = zeros(LONG)
    h = zeros(LONG)
    qDes_pp = zeros(LONG)
    qDef_pp = zeros(LONG)
    u = zeros(LONG)

    q[1] = q0
    q_p[1] = qp0
    hint = 0.0

    for i = 1:l
        t[i] = δt * i
        qN[i] = Amp * sin(ω * t[i])
        qN_p[i] = Amp * ω * cos(ω * t[i])
        qN_pp[i] = -Amp * ω^2 * sin(ω * t[i])

        h[i] = qN[i] - q[i]
        h_p = qN_p[i] - q_p[i]

        if Robust == 1
            S = Λ^2 * hint + 2 * Λ * h[i] + h_p
            qDes_pp[i] = K_VSSM * tanh(S / w) + Λ^2 * h[i] + 2 * Λ * h_p + qN_pp[i]
        else
            qDes_pp[i] = qN_pp[i] + Λ^3 * hint + 3 * Λ^2 * h[i] + 3 * Λ * h_p
        end

        if Adaptive == 1 && i > 3
            qDef_pp[i] = (qDef_pp[i-1] + K) * (1 + B * tanh(A * (q_pp[i-1] - qDes_pp[i]))) - K
        else
            qDef_pp[i] = qDes_pp[i]
        end

        u[i] = qDef_pp[i] - αₐ * q[i]
        q_pp[i] = αₑ * q[i] + δₑ * q_p[i] - βₑ * q[i]^3 + u[i]

        q_p[i+1] = q_p[i] + δt * q_pp[i]
        q[i+1] = q[i] + δt * q_p[i]
        hint = hint + δt * h[i]
    end

    if do_plot
        figure("Duffing - Trajectory")
        grid(true)
        title("Duffing – Trajectory Tracking")
        xlabel("time")
        ylabel("position")
        plot(t[1:l], qN[1:l], color="red", label="nominal")
        plot(t[1:l], q[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("Duffing - Velocity")
        grid(true)
        title("Duffing – Velocity")
        xlabel("time")
        ylabel("velocity")
        plot(t[1:l], qN_p[1:l], color="red", label="nominal")
        plot(t[1:l], q_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("Duffing - Acceleration")
        grid(true)
        title("Duffing – Accelerations")
        xlabel("time")
        ylabel("acceleration")
        plot(t[1:l], qN_pp[1:l], color="red", label="nominal")
        plot(t[1:l], q_pp[1:l], color="green", linestyle="--", label="actual")
        plot(t[1:l], qDes_pp[1:l], color="blue", linestyle="-.", label="desired")
        legend()

        figure("Duffing - Error")
        grid(true)
        title("Duffing – Tracking Error")
        xlabel("time")
        ylabel("error")
        plot(t[1:l], h[1:l], color="red")

        figure("Duffing - Control Signal")
        grid(true)
        title("Duffing – Control Signal")
        xlabel("time")
        ylabel("signal")
        plot(t[1:l], u[1:l], color="red")

        figure("Duffing - Phase Space")
        grid(true)
        title("Duffing – Phase Space")
        xlabel("position")
        ylabel("velocity")
        plot(qN[1:l], qN_p[1:l], color="red", label="nominal")
        plot(q[1:l], q_p[1:l], color="green", linestyle="--", label="actual")
        legend()
    end

    return maximum(abs.(h[1:l]))
end

end # module DuffingModule
