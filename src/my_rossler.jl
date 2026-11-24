module RosslerModule
using PyPlot
using LinearAlgebra
using Random

export run_rossler

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
  do_plot=false,
  metrics=false,
  tail=0.2,
  tag="actual",        # jelmagyarázat-felirat → több futás egy ábrán összehasonlítható
  vel_clip=true,       # sebesség ábrák: az eleji tüske ne nyomja össze a skálát
  mismatch=nothing,    # robusztusság: ε szám → a becsült paraméterek pₐ = (1+ε)·pₑ (nothing = a fenti ₐ értékek)
  noise=0.0,           # robusztusság: mérési zaj szórása (a szabályozó által mért állapot és válasz zajos)
  seed=1)              # a zaj véletlenszám-magja (ismételhető futás)

  if mismatch !== nothing
    aₐ, bₐ, cₐ = (1 + mismatch) .* (aₑ, bₑ, cₑ)
  end
  rng = MersenneTwister(seed)
  nz() = noise > 0 ? noise * randn(rng) : 0.0    # noise = 0 → pontosan a régi viselkedés

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
  Hx = zeros(LONG); Hy = zeros(LONG); Hz = zeros(LONG)   # követési hibák (metrikákhoz, ábrákhoz)
  Sx = zeros(LONG); Sy = zeros(LONG); Sz = zeros(LONG)   # csúszófelületek (csak ábrához)

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

    # mért állapot (zajjal) – ezt látja a szabályozó; a valódi hiba Hx/Hy/Hz az ábrákhoz
    xm = x[i] + nz(); ym = y[i] + nz(); zm = z[i] + nz()
    Hx[i] = xN[i] - x[i]; Hy[i] = yN[i] - y[i]; Hz[i] = zN[i] - z[i]
    h_x = xN[i] - xm
    h_y = yN[i] - ym
    h_z = zN[i] - zm

    if Robust == 1
      S_x = Λ * hint_x + h_x
      S_y = Λ * hint_y + h_y
      S_z = Λ * hint_z + h_z
      Sx[i] = S_x; Sy[i] = S_y; Sz[i] = S_z
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

    u_x[i] = xDef[i] + ym + zm
    u_y[i] = yDef[i] - xm - aₐ * ym
    u_z[i] = zDef[i] - bₐ - zm * (xm - cₐ)

    x_p[i] = -y[i] - z[i] + u_x[i]
    y_p[i] = x[i] + aₑ * y[i] + u_y[i]
    z_p[i] = bₑ + z[i] * (x[i] - cₑ) + u_z[i]
    past_response = [x_p[i] + nz(), y_p[i] + nz(), z_p[i] + nz()]   # mért válasz az RFPT-nek

    x[i+1] = x[i] + δt * x_p[i]
    y[i+1] = y[i] + δt * y_p[i]
    z[i+1] = z[i] + δt * z_p[i]

    # korai kilépés: numerikus robbanás vagy állapot-explózió
    if !isfinite(x[i+1]) || !isfinite(y[i+1]) || !isfinite(z[i+1]) ||
       abs(x[i+1]) > 1e4  || abs(y[i+1]) > 1e4  || abs(z[i+1]) > 1e4
        return metrics ? (rms_ss=Inf, maxe_ss=Inf, rms_all=Inf, maxe_all=Inf, maxu=Inf, ok=false) : Inf
    end

    hint_x += δt * h_x
    hint_y += δt * h_y
    hint_z += δt * h_z
  end

  ex = maximum(abs.(xN[1:l] .- x[1:l]))
  ey = maximum(abs.(yN[1:l] .- y[1:l]))
  ez = maximum(abs.(zN[1:l] .- z[1:l]))

  if do_plot
    T = t[1:l]
    c = tag == "actual" ? "green" : nothing     # egyedi tag → automatikus szín
    c2 = tag == "actual" ? "red" : nothing
    names = ("X", "Y", "Z")

    # ideális vezérlés: ami a PONTOS modellel a nominális pályán tartaná a rendszert
    uNx = xN_p[1:l] .+ yN[1:l] .+ zN[1:l]
    uNy = yN_p[1:l] .- xN[1:l] .- aₑ .* yN[1:l]
    uNz = zN_p[1:l] .- bₑ .- zN[1:l] .* (xN[1:l] .- cₑ)

    # Pozíció: nominal, actual
    axs = axes3("Rossler - Trajectory", "Rössler – Trajectory Tracking")
    for (k, (n, a)) in enumerate(((xN, x), (yN, y), (zN, z)))
      axs[k].plot(T, n[1:l], color="red", linewidth=1.5, label=nlab(axs[k]))
      axs[k].plot(T, a[1:l], color=c, linestyle="--", linewidth=1.5, label=tag)
      axs[k].set_ylabel(names[k])
    end
    finish3(axs)

    # 3D pálya
    fig3d = figure("Rossler - Trajectory 3D")
    ax3d = isempty(fig3d.axes) ? fig3d.add_subplot(111, projection="3d") : fig3d.axes[1]
    ax3d.plot(xN[1:l], yN[1:l], zN[1:l], color="red", label=nlab(ax3d))
    ax3d.plot(x[1:l], y[1:l], z[1:l], color=c, linestyle="--", label=tag)
    ax3d.set_xlabel("X"); ax3d.set_ylabel("Y"); ax3d.set_zlabel("Z")
    ax3d.set_title("Rössler – 3D Trajectory")
    ax3d.legend()

    # Sebesség: nominal, actual, desired (elsőrendű rendszer → a desired itt van)
    axs = axes3("Rossler - Velocities", "Rössler – Velocities")
    for (k, (n, a, d)) in enumerate(((xN_p, x_p, xDes), (yN_p, y_p, yDes), (zN_p, z_p, zDes)))
      axs[k].plot(T, n[1:l], color="red", linewidth=1.5, label=nlab(axs[k]))
      axs[k].plot(T, a[1:l], color=c, linestyle="--", linewidth=1.5, label=tag)
      axs[k].plot(T, d[1:l], linestyle="-.", linewidth=1.2, label="desired ($tag)")
      axs[k].set_ylabel("$(names[k]) [/s]")
      if vel_clip      # skála az eleji tranziens (első 5%) nélkül
        r0 = max(4, round(Int, 0.05 * l)):l
        lim = 1.2 * maximum(abs, vcat(n[r0], a[r0], d[r0]))
        isfinite(lim) && lim > 0 && axs[k].set_ylim(-lim, lim)
      end
    end
    finish3(axs)

    # RFPT: desired vs deformed
    axs = axes3("Rossler - Deformed Velocity", "Rössler – RFPT: desired vs deformed")
    for (k, (d, f)) in enumerate(((xDes, xDef), (yDes, yDef), (zDes, zDef)))
      axs[k].plot(T, d[1:l], linestyle="-.", label="desired ($tag)")
      axs[k].plot(T, f[1:l], label="deformed ($tag)")
      axs[k].set_ylabel("$(names[k]) [/s]")
    end
    finish3(axs)

    # Hiba komponensenként, cél: 0
    axs = axes3("Rossler - Tracking Error", "Rössler – Tracking Errors")
    for (k, H) in enumerate((Hx, Hy, Hz))
      isempty(axs[k].lines) && axs[k].axhline(0, color="black", linewidth=0.8)
      axs[k].plot(T, H[1:l], color=c2, linewidth=1.5, label=tag)
      axs[k].set_ylabel("$(names[k]) error")
    end
    finish3(axs)

    # Hiba normája log skálán
    figure("Rossler - Error (log)")
    grid(true, which="both")
    title("Rössler – ‖Tracking Error‖ (log)")
    xlabel("time")
    ylabel("‖error‖")
    semilogy(T, sqrt.(Hx[1:l] .^ 2 .+ Hy[1:l] .^ 2 .+ Hz[1:l] .^ 2) .+ 1e-16, color=c2, label=tag)
    legend()

    # Vezérlőjel: ideális (nominal) vs actual
    axs = axes3("Rossler - Control Signals", "Rössler – Control Signals")
    for (k, (un, ua)) in enumerate(((uNx, u_x), (uNy, u_y), (uNz, u_z)))
      axs[k].plot(T, un, color="black", linewidth=0.8, label=(isempty(axs[k].lines) ? "nominal (ideális)" : "_nolegend_"))
      axs[k].plot(T, ua[1:l], color=c2, linestyle="--", label=tag)
      axs[k].set_ylabel("u_$(lowercase(names[k]))")
    end
    finish3(axs)

    # Csúszófelületek
    axs = axes3("Rossler - Sliding Surface", "Rössler – Sliding Surfaces S(t)")
    for (k, S) in enumerate((Sx, Sy, Sz))
      isempty(axs[k].lines) && axs[k].axhline(0, color="black", linewidth=0.8)
      axs[k].plot(T, S[1:l], color=c2, label=tag)
      axs[k].set_ylabel("S_$(lowercase(names[k]))")
    end
    finish3(axs)

    # Fázisterek
    for (k, (n, np_, a, ap)) in enumerate(((xN, xN_p, x, x_p), (yN, yN_p, y, y_p), (zN, zN_p, z, z_p)))
      figure("Rossler - Phase $(names[k])")
      grid(true)
      title("Rössler – Phase Space $(names[k])")
      xlabel(names[k])
      ylabel("$(names[k])̇")
      plot(n[1:l], np_[1:l], color="red", label=nlab(gca()))
      plot(a[1:l], ap[1:l], color=c, linestyle="--", label=tag)
      legend()
    end
  end

  if metrics
    # hibanorma időben: ‖h(t)‖ = sqrt(hx² + hy² + hz²)
    en = sqrt.(Hx[1:l] .^ 2 .+ Hy[1:l] .^ 2 .+ Hz[1:l] .^ 2)
    r = (l-round(Int, tail * l)+1):l     # állandósult ablak: utolsó `tail` rész
    um = max(maximum(abs, u_x[1:l]), maximum(abs, u_y[1:l]), maximum(abs, u_z[1:l]))
    return (rms_ss=sqrt(sum(abs2, en[r]) / length(r)), maxe_ss=maximum(en[r]),
      rms_all=sqrt(sum(abs2, en) / l), maxe_all=maximum(en),
      maxu=um, ok=all(isfinite, en))
  end

  return sqrt(ex^2 + ey^2 + ez^2)
end

# ── Ábra-segédfüggvények ──────────────────────────────────────
# 3 egymás alatti tengely (X, Y, Z); ha az ábra már létezik, ugyanazokat adja vissza
# → több futás egymásra rajzolható
function axes3(name, ttl)
  fig = figure(name)
  if isempty(fig.axes)
    fig.subplots(3, 1, sharex=true)
    fig.suptitle(ttl)
  end
  return fig.axes
end

# a nominális csak egyszer kerül a jelmagyarázatba
nlab(ax) = isempty(ax.lines) ? "nominal" : "_nolegend_"

function finish3(axs)
  for ax in axs
    ax.grid(true)
    ax.legend(loc="upper right", fontsize=8)
  end
  axs[end].set_xlabel("time")
  tight_layout()
end

end # module RosslerModule