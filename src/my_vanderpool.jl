module VanDerPolModule
using PyPlot
using Random

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
    do_plot=false,
    metrics=false,
    tail=0.2,
    tag="actual",        # jelmagyarázat-felirat → több futás egy ábrán összehasonlítható
    acc_clip=true,       # gyorsulás ábra: az eleji RFPT-tüske ne nyomja össze a skálát
    mismatch=nothing,    # robusztusság: ε szám → a becsült paraméterek pₐ = (1+ε)·pₑ (nothing = a fenti ₐ értékek)
    noise=0.0,           # robusztusság: mérési zaj szórása (a szabályozó által mért állapot és válasz zajos)
    seed=1)              # a zaj véletlenszám-magja (ismételhető futás)

    if mismatch !== nothing
        μₐ, ωₐ, αₐ, λₐ, mₐ = (1 + mismatch) .* (μₑ, ωₑ, αₑ, λₑ, mₑ)
    end
    rng = MersenneTwister(seed)
    nz() = noise > 0 ? noise * randn(rng) : 0.0    # noise = 0 → pontosan a régi viselkedés

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
    h_arr = zeros(LONG)      # követési hiba (metrikákhoz)
    S_arr = zeros(LONG)      # csúszófelület (csak ábrához)

    q[1] = q0
    q_p[1] = qp0
    hint = 0.0

    for i = 1:l
        t[i] = δt * i
        qN[i] = Amp * sin(ω * t[i])
        qN_p[i] = Amp * ω * cos(ω * t[i])
        qN_pp[i] = -Amp * ω^2 * sin(ω * t[i])

        # mért állapot (zajjal) – ezt látja a szabályozó; a valódi hiba h_arr[i] az ábrákhoz
        qm = q[i] + nz()
        qpm = q_p[i] + nz()
        h_arr[i] = qN[i] - q[i]
        h = qN[i] - qm
        h_p = qN_p[i] - qpm

        if Robust == 1
            S = Λ^2 * hint + 2 * Λ * h + h_p
            S_arr[i] = S
            qDes_pp[i] = K_VSSM * tanh(S / w) + Λ^2 * h + 2 * Λ * h_p + qN_pp[i]
        else
            qDes_pp[i] = qN_pp[i] + Λ^3 * hint + 3 * Λ^2 * h + 3 * Λ * h_p
        end

        if Adaptive == 1 && i > 3
            qDef_pp[i] = (qDef_pp[i-1] + K) * (1 + B * tanh(A * (q_pp[i-1] + nz() - qDes_pp[i]))) - K
        else
            qDef_pp[i] = qDes_pp[i]
        end

        u[i] = mₐ * qDef_pp[i] - μₐ * (1 - qm^2) * qpm + ωₐ^2 * qm + αₐ * qm^3 + λₐ * qm^5
        q_pp[i] = (u[i] + μₑ * (1 - q[i]^2) * q_p[i] - ωₑ^2 * q[i] - αₑ * q[i]^3 - λₑ * q[i]^5) / mₑ

        q_p[i+1] = q_p[i] + δt * q_pp[i]
        q[i+1] = q[i] + δt * q_p[i]
        hint = hint + δt * h
    end

    if do_plot
        # több futás egy ábrán: a nominális csak egyszer kerül a jelmagyarázatba,
        # egyedi tag esetén a szín automatikus (különböző futások különböző színnel)
        nl() = isempty(gca().lines) ? "nominal" : "_nolegend_"
        c = tag == "actual" ? "green" : nothing
        c2 = tag == "actual" ? "red" : nothing

        figure("VanDerPol - Trajectory")
        grid(true)
        title("Van der Pol – Trajectory Tracking")
        xlabel("time")
        ylabel("position")
        plot(t[1:l], qN[1:l], color="red", label=nl())
        plot(t[1:l], q[1:l], color=c, linestyle="--", label=tag)
        legend()

        figure("VanDerPol - Velocity")
        grid(true)
        title("Van der Pol – Velocity")
        xlabel("time")
        ylabel("velocity")
        plot(t[1:l], qN_p[1:l], color="red", label=nl())
        plot(t[1:l], q_p[1:l], color=c, linestyle="--", label=tag)
        legend()

        figure("VanDerPol - Acceleration")
        grid(true)
        title("Van der Pol – Accelerations")
        xlabel("time")
        ylabel("acceleration")
        plot(t[1:l], qN_pp[1:l], color="red", label=nl())
        plot(t[1:l], q_pp[1:l], color=c, linestyle="--", label=tag)
        plot(t[1:l], qDes_pp[1:l], linestyle="-.", label="desired ($tag)")
        if acc_clip      # skála az eleji tranziens (első 5%) nélküli értékekből
            r0 = max(4, round(Int, 0.05 * l)):l
            lim = 1.2 * maximum(abs, vcat(qN_pp[r0], q_pp[r0], qDes_pp[r0]))
            isfinite(lim) && lim > 0 && ylim(-lim, lim)
        end
        legend()

        figure("VanDerPol - Error")
        grid(true)
        title("Van der Pol – Tracking Error")
        xlabel("time")
        ylabel("error")
        axhline(0, color="black", linewidth=0.8)      # cél: h = 0
        plot(t[1:l], h_arr[1:l], color=c2, label=tag)
        legend()

        figure("VanDerPol - Error (log)")
        grid(true, which="both")
        title("Van der Pol – |Tracking Error| (log)")
        xlabel("time")
        ylabel("|error|")
        semilogy(t[1:l], abs.(h_arr[1:l]) .+ 1e-16, color=c2, label=tag)
        legend()

        figure("VanDerPol - Control Signal")
        grid(true)
        title("Van der Pol – Control Signal")
        xlabel("time")
        ylabel("signal")
        # ideális vezérlés: ami a PONTOS modellel a nominális pályán tartaná a rendszert
        uN = mₑ .* qN_pp[1:l] .- μₑ .* (1 .- qN[1:l] .^ 2) .* qN_p[1:l] .+ ωₑ^2 .* qN[1:l] .+
             αₑ .* qN[1:l] .^ 3 .+ λₑ .* qN[1:l] .^ 5
        plot(t[1:l], uN, color="black", linewidth=0.8, label=(isempty(gca().lines) ? "nominal (ideális)" : "_nolegend_"))
        plot(t[1:l], u[1:l], color=c2, linestyle="--", label=tag)
        legend()

        figure("VanDerPol - Sliding Surface")
        grid(true)
        title("Van der Pol – Sliding Surface S(t)")
        xlabel("time")
        ylabel("S")
        axhline(0, color="black", linewidth=0.8)
        plot(t[1:l], S_arr[1:l], color=c2, label=tag)
        legend()

        figure("VanDerPol - Deformed Acceleration")
        grid(true)
        title("Van der Pol – RFPT: desired vs deformed")
        xlabel("time")
        ylabel("acceleration")
        plot(t[1:l], qDes_pp[1:l], linestyle="-.", label="desired ($tag)")
        plot(t[1:l], qDef_pp[1:l], label="deformed ($tag)")
        legend()

        figure("VanDerPol - Phase Space")
        grid(true)
        title("Van der Pol – Phase Space")
        xlabel("q")
        ylabel("q̇")
        plot(qN[1:l], qN_p[1:l], color="red", label=nl())
        plot(q[1:l], q_p[1:l], color=c, linestyle="--", label=tag)
        legend()
    end

    if metrics
        r = (l-round(Int, tail * l)+1):l     # állandósult ablak: utolsó `tail` rész
        ok = all(isfinite, q[1:l]) && maximum(abs, q[1:l]) < 1e3 && maximum(abs, u[1:l]) < 100
        return (rms_ss=sqrt(sum(abs2, h_arr[r]) / length(r)), maxe_ss=maximum(abs, h_arr[r]),
            rms_all=sqrt(sum(abs2, h_arr[1:l]) / l), maxe_all=maximum(abs, h_arr[1:l]),
            maxu=maximum(abs, u[1:l]), ok=ok)
    end

    return maximum(abs.(qN[1:l] .- q[1:l]))
end

end # module VanDerPolModule