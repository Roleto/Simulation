module VanDerPolModule
using PyPlot

export run_vanderpol

function run_vanderpol(q0=2.0 * sin(0.5 * 1e-3), qp0=2.0 * 0.5 * cos(0.5 * 1e-3);
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
    μₑ=0.4,
    ωₑ=0.46,
    αₑ=1.0,
    λₑ=0.1,
    mₑ=1.0,
    μₐ=0.5,
    ωₐ=0.42,
    αₐ=0.9,
    λₐ=0.09,
    mₐ=0.8,
    do_plot=false)

    l = LONG - 1
    t = zeros(LONG)
    q = zeros(LONG)
    q_p = zeros(LONG)
    q_pp = zeros(LONG)
    qN = zeros(LONG)
    qN_p = zeros(LONG)
    qN_pp = zeros(LONG)
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

        h = qN[i] - q[i]
        h_p = qN_p[i] - q_p[i]

        if Robust == 1
            S = Λ^2 * hint + 2 * Λ * h + h_p
            qDes_pp[i] = K_VSSM * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + qN_pp[i]
        else
            qDes_pp[i] = qN_pp[i] + Λ^3 * hint + 3 * Λ^2 * h + 3 * Λ * h_p
        end

        if Adaptive == 1 && i > 3
            qDef_pp[i] = (qDef_pp[i-1] + K) * (1 + B * tanh(A * (q_pp[i-1] - qDes_pp[i]))) - K
        else
            qDef_pp[i] = qDes_pp[i]
        end

        u[i] = mₐ * qDef_pp[i]
        q_pp[i] = (u[i] + μₑ * (1 - q[i]^2) * q_p[i] - ωₑ^2 * q[i] - αₑ * q[i]^3 - λₑ * q[i]^5) / mₑ

        q_p[i+1] = q_p[i] + δt * q_pp[i]
        q[i+1] = q[i] + δt * q_p[i]
        hint = hint + δt * h
    end

    if do_plot
        figure("VanDerPol - Trajectory")
        grid(true)
        title("Van der Pol – Trajectory Tracking")
        xlabel("time")
        ylabel("position")
        plot(t[1:l], qN[1:l], color="red", label="nominal")
        plot(t[1:l], q[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("VanDerPol - Velocity")
        grid(true)
        title("Van der Pol – Velocity")
        xlabel("time")
        ylabel("velocity")
        plot(t[1:l], qN_p[1:l], color="red", label="nominal")
        plot(t[1:l], q_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("VanDerPol - Acceleration")
        grid(true)
        title("Van der Pol – Acceleration")
        xlabel("time")
        ylabel("acceleration")
        plot(t[1:l], qN_pp[1:l], color="red", label="nominal")
        plot(t[1:l], q_pp[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("VanDerPol - Error")
        grid(true)
        title("Van der Pol – Tracking Error")
        xlabel("time")
        ylabel("error")
        plot(t[1:l], qN[1:l] .- q[1:l], color="red")

        figure("VanDerPol - Control Signal")
        grid(true)
        title("Van der Pol – Control Signal")
        xlabel("time")
        ylabel("signal")
        plot(t[1:l], u[1:l], color="red")

        figure("VanDerPol - Phase Space")
        grid(true)
        title("Van der Pol – Phase Space")
        xlabel("q")
        ylabel("q̇")
        plot(qN[1:l], qN_p[1:l], color="red", label="nominal")
        plot(q[1:l], q_p[1:l], color="green", linestyle="--", label="actual")
        legend()
        show()
    end

    return maximum(abs.(qN[1:l] .- q[1:l]))
end

end # module VanDerPolModule