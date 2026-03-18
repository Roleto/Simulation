# ============================================================
# main.jl  –  RFPT + VSSM Simulation Runner  (Fresh Start)
# ============================================================
#
# HOW TO USE:
#   julia main.jl                  – runs the single-simulation demo
#   include("main.jl") in REPL    – loads all functions interactively
#
# STRUCTURE:
#   Section 1 – Simulation functions (one per system)
#               Each accepts initial conditions + optional keyword
#               parameter overrides and returns max tracking error.
#   Section 2 – Run a single simulation (default OR custom params)
#   Section 3 – Parameter sweep over (position × velocity) grid
#   Section 4 – 3D plotting  (x=position, y=velocity, z=max_error)
#
# Original backup saved in: main_original_backup.jl
# ============================================================

using PyPlot
using LinearAlgebra

# ============================================================
# Section 1:  Simulation Functions
# ============================================================

# ------------------------------------------------------------
# 1a. Duffing Oscillator  (2nd order, SISO)
#     Positional arguments:  q0   – initial position
#                            qp0  – initial velocity
#     Returns: maximum absolute tracking error
# ------------------------------------------------------------
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
        show()
    end

    return maximum(abs.(h[1:l]))
end

# ------------------------------------------------------------
# 1b. Van der Pol Oscillator  (2nd order, SISO)
#     Positional arguments:  q0   – initial position
#                            qp0  – initial velocity
#     Returns: maximum absolute tracking error
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# MIMO deformation helper  (shared by Lorenz and Rossler)
# ------------------------------------------------------------
function G_MIMO(past_input, past_response, desired, error_limit, Kc, Bc, Ac)
    error_norm = norm(past_response - desired, 2)
    if error_norm > error_limit
        e_direction = (past_response - desired) / error_norm
        B_factor = Bc * tanh(Ac * error_norm)
        return (1 + B_factor) * past_input + B_factor * Kc * e_direction
    else
        return past_input
    end
end

# ------------------------------------------------------------
# 1c. Lorenz System  (3D, MIMO)
#     Positional arguments:  x0, y0, z0  – initial conditions
#     Returns: combined maximum tracking error  sqrt(ex²+ey²+ez²)
# ------------------------------------------------------------
function run_lorenz(x0=0.0, y0=0.0, z0=0.0;
    Adaptive=1,
    Robust=1,
    K=1e2,
    B=-1.0,
    A=1.97e-3,
    δt=1e-3,
    LONG=Int(2e4),
    Λ=1,
    K_VSSM=50,
    w=1,
    A₁=2, ω₁=0.5,
    A₂=3, ω₂=0.7,
    A₃=1, ω₃=1,
    βₑ=8 / 3, σₑ=5, ρₑ=40,
    βₐ=7 / 3, σₐ=4, ρₐ=36,
    error_limit=1e-3,
    do_plot=false)

    l = LONG - 1
    t = zeros(LONG)
    xN = zeros(LONG)
    xN_p = zeros(LONG)
    yN = zeros(LONG)
    yN_p = zeros(LONG)
    zN = zeros(LONG)
    zN_p = zeros(LONG)
    xDes = zeros(LONG)
    xDef = zeros(LONG)
    yDes = zeros(LONG)
    yDef = zeros(LONG)
    zDes = zeros(LONG)
    zDef = zeros(LONG)
    x_p = zeros(LONG)
    y_p = zeros(LONG)
    z_p = zeros(LONG)
    x = zeros(LONG)
    y = zeros(LONG)
    z = zeros(LONG)
    u_x = zeros(LONG)
    u_y = zeros(LONG)
    u_z = zeros(LONG)

    x[1] = x0
    y[1] = y0
    z[1] = z0
    hint_x = 0.0
    hint_y = 0.0
    hint_z = 0.0
    past_input = [0.0, 0.0, 0.0]
    past_response = [0.0, 0.0, 0.0]

    for i = 1:l
        t[i] = δt * i
        xN[i] = A₁ * sin(ω₁ * t[i])
        xN_p[i] = A₁ * ω₁ * cos(ω₁ * t[i])
        yN[i] = A₂ * sin(ω₂ * t[i])
        yN_p[i] = A₂ * ω₂ * cos(ω₂ * t[i])
        zN[i] = A₃ * sin(ω₃ * t[i])
        zN_p[i] = A₃ * ω₃ * cos(ω₃ * t[i])

        h_x = xN[i] - x[i]
        h_y = yN[i] - y[i]
        h_z = zN[i] - z[i]

        if Robust == 1
            S_x = Λ * hint_x + h_x
            S_y = Λ * hint_y + h_y
            S_z = Λ * hint_z + h_z
            xDes[i] = xN_p[i] + Λ * h_x + K_VSSM * tanh(S_x / w)
            yDes[i] = yN_p[i] + Λ * h_y + K_VSSM * tanh(S_y / w)
            zDes[i] = zN_p[i] + Λ * h_z + K_VSSM * tanh(S_z / w)
        else
            xDes[i] = xN_p[i] + Λ^2 * hint_x + 2 * Λ * h_x
            yDes[i] = yN_p[i] + Λ^2 * hint_y + 2 * Λ * h_y
            zDes[i] = zN_p[i] + Λ^2 * hint_z + 2 * Λ * h_z
        end

        desired = [xDes[i], yDes[i], zDes[i]]

        if Adaptive == 1 && i > 3
            past_input = G_MIMO(past_input, past_response, desired, error_limit, K, B, A)
        else
            past_input = desired
        end

        xDef[i] = past_input[1]
        yDef[i] = past_input[2]
        zDef[i] = past_input[3]

        u_x[i] = xDef[i] - σₐ * (y[i] - x[i])
        u_y[i] = yDef[i] - x[i] * (ρₐ - z[i]) + y[i]
        u_z[i] = zDef[i] - x[i] * y[i] + βₐ * z[i]

        x_p[i] = σₑ * (y[i] - x[i]) + u_x[i]
        y_p[i] = x[i] * (ρₑ - z[i]) - y[i] + u_y[i]
        z_p[i] = x[i] * y[i] - βₑ * z[i] + u_z[i]
        past_response = [x_p[i], y_p[i], z_p[i]]

        x[i+1] = x[i] + δt * x_p[i]
        y[i+1] = y[i] + δt * y_p[i]
        z[i+1] = z[i] + δt * z_p[i]

        hint_x += δt * h_x
        hint_y += δt * h_y
        hint_z += δt * h_z
    end

    ex = maximum(abs.(xN[1:l] .- x[1:l]))
    ey = maximum(abs.(yN[1:l] .- y[1:l]))
    ez = maximum(abs.(zN[1:l] .- z[1:l]))

    # Numerical accelerations (2nd derivative via finite difference)
    x_pp  = vcat(0.0, diff(x_p[1:l])  ./ δt)
    y_pp  = vcat(0.0, diff(y_p[1:l])  ./ δt)
    z_pp  = vcat(0.0, diff(z_p[1:l])  ./ δt)
    xN_pp = vcat(0.0, diff(xN_p[1:l]) ./ δt)
    yN_pp = vcat(0.0, diff(yN_p[1:l]) ./ δt)
    zN_pp = vcat(0.0, diff(zN_p[1:l]) ./ δt)

    if do_plot
        # -- Trajectory (3 subplots) --
        fig1 = figure("Lorenz - Trajectory")
        fig1.suptitle("Lorenz – Trajectory Tracking")
        subplot(311)
        ax1 = gca()
        grid(true)
        ylabel("X")
        plot(t[1:l], xN[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1)
        grid(true)
        ylabel("Y")
        plot(t[1:l], yN[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1)
        grid(true)
        ylabel("Z")
        xlabel("time")
        plot(t[1:l], zN[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

        # -- 3D Trajectory --
        fig2 = figure("Lorenz - Trajectory 3D")
        ax3d = fig2.add_subplot(111, projection="3d")
        ax3d.plot(xN[1:l], yN[1:l], zN[1:l], color="red", label="nominal")
        ax3d.plot(x[1:l], y[1:l], z[1:l], color="green", linestyle="--", label="actual")
        ax3d.set_xlabel("X")
        ax3d.set_ylabel("Y")
        ax3d.set_zlabel("Z")
        ax3d.set_title("Lorenz – 3D Trajectory")
        ax3d.legend()

        # -- Velocities (3 subplots) --
        fig3 = figure("Lorenz - Velocities")
        fig3.suptitle("Lorenz – Velocities")
        subplot(311)
        ax1 = gca()
        grid(true)
        ylabel("X [/s]")
        plot(t[1:l], xN_p[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x_p[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1)
        grid(true)
        ylabel("Y [/s]")
        plot(t[1:l], yN_p[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y_p[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1)
        grid(true)
        ylabel("Z [/s]")
        xlabel("time")
        plot(t[1:l], zN_p[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z_p[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

        # -- Accelerations (3 subplots) --
        fig_acc = figure("Lorenz - Accelerations")
        fig_acc.suptitle("Lorenz – Accelerations")
        subplot(311); ax1 = gca(); grid(true); ylabel("X [/s²]")
        plot(t[1:l], xN_pp[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x_pp[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1); grid(true); ylabel("Y [/s²]")
        plot(t[1:l], yN_pp[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y_pp[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1); grid(true); ylabel("Z [/s²]"); xlabel("time")
        plot(t[1:l], zN_pp[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z_pp[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

        # -- Tracking Errors (3 subplots) --
        fig4 = figure("Lorenz - Tracking Error")
        fig4.suptitle("Lorenz – Tracking Errors")
        subplot(311)
        ax1 = gca()
        grid(true)
        ylabel("X error")
        plot(t[1:l], xN[1:l] .- x[1:l], color="red", linewidth=1.5)
        subplot(312, sharex=ax1)
        grid(true)
        ylabel("Y error")
        plot(t[1:l], yN[1:l] .- y[1:l], color="red", linewidth=1.5)
        subplot(313, sharex=ax1)
        grid(true)
        ylabel("Z error")
        xlabel("time")
        plot(t[1:l], zN[1:l] .- z[1:l], color="red", linewidth=1.5)
        tight_layout()

        # -- Control Signals --
        figure("Lorenz - Control Signals")
        grid(true)
        title("Lorenz – Control Signals")
        xlabel("time")
        ylabel("signal")
        plot(t[1:l], u_x[1:l], color="red", label="u_x")
        plot(t[1:l], u_y[1:l], color="green", linestyle="--", label="u_y")
        plot(t[1:l], u_z[1:l], color="blue", label="u_z")
        legend()
        tight_layout()

        # -- Phase Spaces --
        figure("Lorenz - Phase X")
        grid(true)
        title("Lorenz – Phase Space X")
        xlabel("X")
        ylabel("Ẋ")
        plot(xN[1:l], xN_p[1:l], color="red", label="nominal")
        plot(x[1:l], x_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("Lorenz - Phase Y")
        grid(true)
        title("Lorenz – Phase Space Y")
        xlabel("Y")
        ylabel("Ẏ")
        plot(yN[1:l], yN_p[1:l], color="red", label="nominal")
        plot(y[1:l], y_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("Lorenz - Phase Z")
        grid(true)
        title("Lorenz – Phase Space Z")
        xlabel("Z")
        ylabel("Ż")
        plot(zN[1:l], zN_p[1:l], color="red", label="nominal")
        plot(z[1:l], z_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        show()
    end

    return sqrt(ex^2 + ey^2 + ez^2)
end

# ------------------------------------------------------------
# 1d. Rossler System  (3D, MIMO)
#     Positional arguments:  x0, y0, z0  – initial conditions
#     Returns: combined maximum tracking error  sqrt(ex²+ey²+ez²)
# ------------------------------------------------------------
function run_rossler(x0=2.0 * sin(0.5 * 1e-3), y0=3.0 * sin(0.7 * 1e-3), z0=1.0 * sin(1.0 * 1e-3);
    Adaptive=1,
    Robust=1,
    K=1e2,
    B=-1.0,
    A=1.97e-3,
    δt=1e-3,
    LONG=Int(2e4),
    Λ=1,
    K_VSSM=50,
    w=1,
    A₁=2, ω₁=0.5,
    A₂=3, ω₂=0.7,
    A₃=1, ω₃=1,
    aₑ=0.01, bₑ=0.2, cₑ=5.7,
    aₐ=0.1, bₐ=0.3, cₐ=5.5,
    error_limit=1e-3,
    do_plot=false)

    l = LONG - 1
    t = zeros(LONG)
    xN = zeros(LONG)
    xN_p = zeros(LONG)
    yN = zeros(LONG)
    yN_p = zeros(LONG)
    zN = zeros(LONG)
    zN_p = zeros(LONG)
    xDes = zeros(LONG)
    xDef = zeros(LONG)
    yDes = zeros(LONG)
    yDef = zeros(LONG)
    zDes = zeros(LONG)
    zDef = zeros(LONG)
    x_p = zeros(LONG)
    y_p = zeros(LONG)
    z_p = zeros(LONG)
    x = zeros(LONG)
    y = zeros(LONG)
    z = zeros(LONG)
    u_x = zeros(LONG)
    u_y = zeros(LONG)
    u_z = zeros(LONG)

    x[1] = x0
    y[1] = y0
    z[1] = z0
    hint_x = 0.0
    hint_y = 0.0
    hint_z = 0.0
    past_input = [0.0, 0.0, 0.0]
    past_response = [0.0, 0.0, 0.0]

    for i = 1:l
        t[i] = δt * i
        xN[i] = A₁ * sin(ω₁ * t[i])
        xN_p[i] = A₁ * ω₁ * cos(ω₁ * t[i])
        yN[i] = A₂ * sin(ω₂ * t[i])
        yN_p[i] = A₂ * ω₂ * cos(ω₂ * t[i])
        zN[i] = A₃ * sin(ω₃ * t[i])
        zN_p[i] = A₃ * ω₃ * cos(ω₃ * t[i])

        h_x = xN[i] - x[i]
        h_y = yN[i] - y[i]
        h_z = zN[i] - z[i]

        if Robust == 1
            S_x = Λ * hint_x + h_x
            S_y = Λ * hint_y + h_y
            S_z = Λ * hint_z + h_z
            xDes[i] = xN_p[i] + Λ * h_x + K_VSSM * tanh(S_x / w)
            yDes[i] = yN_p[i] + Λ * h_y + K_VSSM * tanh(S_y / w)
            zDes[i] = zN_p[i] + Λ * h_z + K_VSSM * tanh(S_z / w)
        else
            xDes[i] = xN_p[i] + Λ^2 * hint_x + 2 * Λ * h_x
            yDes[i] = yN_p[i] + Λ^2 * hint_y + 2 * Λ * h_y
            zDes[i] = zN_p[i] + Λ^2 * hint_z + 2 * Λ * h_z
        end

        desired = [xDes[i], yDes[i], zDes[i]]

        if Adaptive == 1 && i > 3
            past_input = G_MIMO(past_input, past_response, desired, error_limit, K, B, A)
        else
            past_input = desired
        end

        xDef[i] = past_input[1]
        yDef[i] = past_input[2]
        zDef[i] = past_input[3]

        u_x[i] = xDef[i] + y[i] + z[i]
        u_y[i] = yDef[i] - x[i] - aₐ * y[i]
        u_z[i] = zDef[i] - bₐ - z[i] * (x[i] - cₐ)

        x_p[i] = -y[i] - z[i] + u_x[i]
        y_p[i] = x[i] + aₑ * y[i] + u_y[i]
        z_p[i] = bₑ + z[i] * (x[i] - cₑ) + u_z[i]
        past_response = [x_p[i], y_p[i], z_p[i]]

        x[i+1] = x[i] + δt * x_p[i]
        y[i+1] = y[i] + δt * y_p[i]
        z[i+1] = z[i] + δt * z_p[i]

        hint_x += δt * h_x
        hint_y += δt * h_y
        hint_z += δt * h_z
    end

    ex = maximum(abs.(xN[1:l] .- x[1:l]))
    ey = maximum(abs.(yN[1:l] .- y[1:l]))
    ez = maximum(abs.(zN[1:l] .- z[1:l]))

    # Numerical accelerations (2nd derivative via finite difference)
    x_pp  = vcat(0.0, diff(x_p[1:l])  ./ δt)
    y_pp  = vcat(0.0, diff(y_p[1:l])  ./ δt)
    z_pp  = vcat(0.0, diff(z_p[1:l])  ./ δt)
    xN_pp = vcat(0.0, diff(xN_p[1:l]) ./ δt)
    yN_pp = vcat(0.0, diff(yN_p[1:l]) ./ δt)
    zN_pp = vcat(0.0, diff(zN_p[1:l]) ./ δt)

    if do_plot
        # -- Trajectory (3 subplots) --
        fig1 = figure("Rossler - Trajectory")
        fig1.suptitle("Rössler – Trajectory Tracking")
        subplot(311)
        ax1 = gca()
        grid(true)
        ylabel("X")
        plot(t[1:l], xN[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1)
        grid(true)
        ylabel("Y")
        plot(t[1:l], yN[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1)
        grid(true)
        ylabel("Z")
        xlabel("time")
        plot(t[1:l], zN[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

        # -- Velocities (3 subplots) --
        fig_vel = figure("Rossler - Velocities")
        fig_vel.suptitle("Rössler – Velocities")
        subplot(311); ax1 = gca(); grid(true); ylabel("X [/s]")
        plot(t[1:l], xN_p[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x_p[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1); grid(true); ylabel("Y [/s]")
        plot(t[1:l], yN_p[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y_p[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1); grid(true); ylabel("Z [/s]"); xlabel("time")
        plot(t[1:l], zN_p[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z_p[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

        # -- Accelerations (3 subplots) --
        fig_acc = figure("Rossler - Accelerations")
        fig_acc.suptitle("Rössler – Accelerations")
        subplot(311); ax1 = gca(); grid(true); ylabel("X [/s²]")
        plot(t[1:l], xN_pp[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x_pp[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1); grid(true); ylabel("Y [/s²]")
        plot(t[1:l], yN_pp[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y_pp[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1); grid(true); ylabel("Z [/s²]"); xlabel("time")
        plot(t[1:l], zN_pp[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z_pp[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

        # -- Tracking Errors (3 subplots) --
        fig2 = figure("Rossler - Tracking Error")
        fig2.suptitle("Rössler – Tracking Errors")
        subplot(311)
        ax1 = gca()
        grid(true)
        ylabel("X error")
        plot(t[1:l], xN[1:l] .- x[1:l], color="red", linewidth=1.5)
        subplot(312, sharex=ax1)
        grid(true)
        ylabel("Y error")
        plot(t[1:l], yN[1:l] .- y[1:l], color="red", linewidth=1.5)
        subplot(313, sharex=ax1)
        grid(true)
        ylabel("Z error")
        xlabel("time")
        plot(t[1:l], zN[1:l] .- z[1:l], color="red", linewidth=1.5)
        tight_layout()

        # -- Control Signals --
        figure("Rossler - Control Signals")
        grid(true)
        title("Rössler – Control Signals")
        xlabel("time")
        ylabel("signal")
        plot(t[1:l], u_x[1:l], color="red", label="u_x")
        plot(t[1:l], u_y[1:l], color="green", linestyle="--", label="u_y")
        plot(t[1:l], u_z[1:l], color="blue", label="u_z")
        legend()
        tight_layout()

        # -- Phase Spaces --
        figure("Rossler - Phase X")
        grid(true)
        title("Rössler – Phase Space X")
        xlabel("X")
        ylabel("Ẋ")
        plot(xN[1:l], xN_p[1:l], color="red", label="nominal")
        plot(x[1:l], x_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("Rossler - Phase Y")
        grid(true)
        title("Rössler – Phase Space Y")
        xlabel("Y")
        ylabel("Ẏ")
        plot(yN[1:l], yN_p[1:l], color="red", label="nominal")
        plot(y[1:l], y_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        figure("Rossler - Phase Z")
        grid(true)
        title("Rössler – Phase Space Z")
        xlabel("Z")
        ylabel("Ż")
        plot(zN[1:l], zN_p[1:l], color="red", label="nominal")
        plot(z[1:l], z_p[1:l], color="green", linestyle="--", label="actual")
        legend()

        show()
    end

    return sqrt(ex^2 + ey^2 + ez^2)
end


# ============================================================
# Section 2:  Run a Single Simulation
# ============================================================
# Two ways to call a simulation:
#
#   Method 1 – default parameters (just call with no arguments):
#       e = run_duffing()
#
#   Method 2 – custom initial conditions (positional args):
#       e = run_duffing(5.0, 3.0)
#
#   Method 2b – also override any named parameter:
#       e = run_duffing(5.0, 3.0; K=2e5, Amp=3)
# ============================================================

println("=== Duffing  (default params) ===")
e_duffing = run_duffing()
println("  max error: ", e_duffing)

#println("=== Duffing  (custom: q0=5, qp0=3) ===")
#e_duffing_custom = run_duffing(5.0, 3.0)
#println("  max error: ", e_duffing_custom)

#println("=== VanDerPol  (default) ===")
#e_vdp = run_vanderpol()
#println("  max error: ", e_vdp)

#println("=== Lorenz  (default) ===")
#e_lorenz = run_lorenz()
#println("  max error: ", e_lorenz)

#println("=== Rossler  (default) ===")
#e_rossler = run_rossler()
#println("  max error: ", e_rossler)


# ============================================================
# Section 3:  Parameter Sweep  (position × velocity grid)
# ============================================================
# Choose which simulation to sweep and set ranges.
# For 3D systems (Lorenz / Rossler) we sweep x0 × y0 and fix z0.
#
# Uncomment this block to run the sweep.
# ============================================================

#= ---- Sweep for Duffing / VanDerPol ----
positions  = range(-15.0, 15.0, length=20)
velocities = range(-15.0, 15.0, length=20)
errors     = zeros(length(positions), length(velocities))

for (i, q0) in enumerate(positions)
    for (j, qp0) in enumerate(velocities)
        errors[i, j] = run_duffing(q0, qp0)
        # swap with run_vanderpol(q0, qp0) to sweep VanDerPol
    end
end

pos_grid = [q0  for q0  in positions,  qp0 in velocities]
vel_grid = [qp0 for q0  in positions,  qp0 in velocities]
=#


# ============================================================
# Section 4:  3D Plotting
#             x = initial position
#             y = initial velocity
#             z = maximum tracking error
# ============================================================
# Uncomment after the sweep in Section 3 is done.
# Both a scatter plot and a surface plot are provided –
# decide which one you prefer and remove the other.
# ============================================================

#= ---- Scatter plot ----
fig1 = figure("Scatter")
ax1  = fig1.add_subplot(111, projection="3d")
ax1.scatter(vec(pos_grid), vec(vel_grid), vec(errors),
            c=vec(errors), cmap="viridis")
ax1.set_xlabel("Initial position  q₀")
ax1.set_ylabel("Initial velocity  q̇₀")
ax1.set_zlabel("Max tracking error")
ax1.set_title("Duffing – Max error vs initial conditions (scatter)")
tight_layout()
=#

#= ---- Surface plot ----
fig2 = figure("Surface")
ax2  = fig2.add_subplot(111, projection="3d")
ax2.plot_surface(pos_grid, vel_grid, errors, cmap="viridis", alpha=0.85)
ax2.set_xlabel("Initial position  q₀")
ax2.set_ylabel("Initial velocity  q̇₀")
ax2.set_zlabel("Max tracking error")
ax2.set_title("Duffing – Max error vs initial conditions (surface)")
tight_layout()
=#

#show()
