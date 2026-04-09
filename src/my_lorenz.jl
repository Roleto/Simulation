module LorenzModule
using PyPlot
using LinearAlgebra

export run_lorenz

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

    x_pp = vcat(0.0, diff(x_p[1:l]) ./ δt)
    y_pp = vcat(0.0, diff(y_p[1:l]) ./ δt)
    z_pp = vcat(0.0, diff(z_p[1:l]) ./ δt)
    xN_pp = vcat(0.0, diff(xN_p[1:l]) ./ δt)
    yN_pp = vcat(0.0, diff(yN_p[1:l]) ./ δt)
    zN_pp = vcat(0.0, diff(zN_p[1:l]) ./ δt)

    if do_plot
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

        fig2 = figure("Lorenz - Trajectory 3D")
        ax3d = fig2.add_subplot(111, projection="3d")
        ax3d.plot(xN[1:l], yN[1:l], zN[1:l], color="red", label="nominal")
        ax3d.plot(x[1:l], y[1:l], z[1:l], color="green", linestyle="--", label="actual")
        ax3d.set_xlabel("X")
        ax3d.set_ylabel("Y")
        ax3d.set_zlabel("Z")
        ax3d.set_title("Lorenz – 3D Trajectory")
        ax3d.legend()

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

        fig_acc = figure("Lorenz - Accelerations")
        fig_acc.suptitle("Lorenz – Accelerations")
        subplot(311)
        ax1 = gca()
        grid(true)
        ylabel("X [/s²]")
        plot(t[1:l], xN_pp[1:l], color="red", linewidth=1.5, label="nominal")
        plot(t[1:l], x_pp[1:l], color="green", linestyle="--", linewidth=1.5, label="actual")
        legend(loc="upper right")
        subplot(312, sharex=ax1)
        grid(true)
        ylabel("Y [/s²]")
        plot(t[1:l], yN_pp[1:l], color="red", linewidth=1.5)
        plot(t[1:l], y_pp[1:l], color="green", linestyle="--", linewidth=1.5)
        subplot(313, sharex=ax1)
        grid(true)
        ylabel("Z [/s²]")
        xlabel("time")
        plot(t[1:l], zN_pp[1:l], color="red", linewidth=1.5)
        plot(t[1:l], z_pp[1:l], color="green", linestyle="--", linewidth=1.5)
        tight_layout()

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

        # show()
    end

    return sqrt(ex^2 + ey^2 + ez^2)
end

end # module LorenzModule