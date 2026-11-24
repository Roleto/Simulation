# Szakdolgozat – szimulációs jegyzetek (3. fejezet)

RFPT + VSSM hibrid szabályozó, Julia. Folyamatosan bővülő jegyzet: eredmények, magyarázatok, döntések.

---

## 1. Futtatás

```julia
# a Simulation mappából, julia --project=.
using Revise
includet("scripts/initial_condition_sweeps.jl")

ic_run(:duffing)                                   # default kezdőérték + default paraméterek
ic_run(:duffing, 7.3, -4.1)                        # saját kezdőérték
ic_run(:duffing, 7.3, -4.1; Λ=100, K_VSSM=1000)    # + paraméter felülírás
ic_run(:lorenz, 1.0, 0.5, 0.2)                     # 3D modell: 3 kezdőérték

ic_sweep(:duffing, -20:2:20, -20:2:20)             # kezdeti feltétel sweep, maxe_ss, log színskála
ic_sweep(:vanderpol, -5:0.5:5, -5:0.5:5; metric=:rms_ss, Λ=56, K_VSSM=1000)

# robusztusság (abláció: hibrid / csak RFPT / csak VSSM)
robust_sweep(:duffing)                                       # modellhiba ε = -0.9 … 2.0
robust_sweep(:duffing, [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]; param=:noise)

# összehasonlítás egy ábrán (tag = felirat)
close("all")
ic_run(:vanderpol, 3.2, -1.7; Λ=1,  tag="Λ=1")
ic_run(:vanderpol, 3.2, -1.7; Λ=56, tag="Λ=56")

# abláció
ic_run(:duffing, 5.0, 0.0; Robust=0,   tag="csak RFPT")
ic_run(:duffing, 5.0, 0.0; Adaptive=0, tag="csak VSSM")
ic_run(:duffing, 5.0, 0.0;             tag="hibrid")
```

Mutatók (`metrics=true`): `maxe_ss`, `rms_ss` (utolsó 20%, állandósult), `maxe_all`, `rms_all`, `maxu`.
Lorenz / Rössler: a hiba a ‖h(t)‖ = √(hx²+hy²+hz²) norma.

---

## 2. Ábrák modellenként

| Ábra | Tartalom | Mit mutat |
|---|---|---|
| Pozíció | nominal, actual | követés |
| Sebesség | nominal, actual (+ desired, deformed – Lorenz/Rössler) | követés; elsőrendű rendszernél itt dolgozik az RFPT |
| Gyorsulás (Duffing, VdP) | nominal, actual, desired | az RFPT a kívánt gyorsulást valósítja meg |
| RFPT: desired vs deformed | qDes, qDef | mennyit „torzít” az RFPT a modellhiba miatt |
| Hiba, hiba log skálán | h(t), 0-s cél | konvergencia, elért pontosság |
| Vezérlőjel | ideális (nominal) u, actual u | beavatkozás; actual → ideális = az RFPT megtanulta a modellhibát |
| Csúszófelület S(t) | S, 0-s cél | elérési fázis, csúszó mód, chattering |
| Fázistér | nominal, actual | beállás a referenciára |
| 3D pálya (Lorenz, Rössler) | nominal, actual | — |

- **Desired** csak ott létezik, ahol a szabályozó előírja: másodrendű rendszernél (Duffing, VdP) gyorsulás, elsőrendűnél (Lorenz, Rössler) sebesség szinten. Pozícióra nincs desired.
- **Hiba**: a h = nominal − actual maga az eltérés, ezért nincs külön nominal/desired görbéje; a cél 0.
- **Ideális (nominális) vezérlés**: az a u, amely a *pontos* modellel a rendszert hibátlanul a nominális pályán tartaná.
  - Duffing: `uN = q̈N − αₑ·qN + δₑ·q̇N + βₑ·qN³`
  - VdP: `uN = mₑ·q̈N − μₑ(1−qN²)q̇N + ωₑ²qN + αₑqN³ + λₑqN⁵`
  - Lorenz: `uNx = ẋN − σₑ(yN−xN)`, `uNy = ẏN − xN(ρₑ−zN) + yN`, `uNz = żN − xN·yN + βₑzN`
  - Rössler: `uNx = ẋN + yN + zN`, `uNy = ẏN − xN − aₑ·yN`, `uNz = żN − bₑ − zN(xN − cₑ)`
- Lorenz/Rössler: a numerikus deriváltból számolt „gyorsulás” ábra el lett hagyva – elsőrendű rendszernél nem szabályozott mennyiség, nincs tartalma.

---

## 3. Eredmények, megfigyelések

### 3.1 A teljes időre vett max hiba = a kezdeti hiba
A default sweep hőtérképeken (`maximum(abs.(h))` az egész futásra) a hiba ≈ |q₀| (Duffing, VdP), ill. ≈ √(x₀²+y₀²) (Lorenz, Rössler). A maximum mindig t=0-ban van → a szabályozóról nem mond semmit.
**Döntés:** kezdeti feltétel vizsgálatnál az állandósult hibát (`maxe_ss`, utolsó 20%) nézzük.
Default paraméterekkel minden kezdőérték stabil volt (nincs NaN / divergált cella).

### 3.2 Duffing – kezdeti feltétel sweep (default: Λ=1, K_VSSM=500)
- Rács: q₀, q̇₀ ∈ −20…20, lépés 2.
- `maxe_ss` ≈ 1e-4 mindenhol (log10: −3.96 … −4.03) → **gyakorlatilag független a kezdőértéktől.**
- q₀ ≈ −20 körül ~15%-kal nagyobb (valószínűleg lassabb lecsengés).
- Az 1e-4 szint numerikus alsó korlát, lásd 3.3.

### 3.3 Numerikus hiba-padló: a hiba ∝ δt
VdP, default (Λ=1, K_VSSM=500), default kezdőérték – az állandósult hiba a lépésközzel arányos:

| δt | max hiba (2. félidő) |
|---|---|
| 1e-3 | 1.0e-4 |
| 5e-4 | 5.0e-5 |
| 2.5e-4 | 2.5e-5 |

- Ok: a rendszer explicit Euler-lépéssel halad, a nominális sebesség analitikus → δt/2 fáziskésés (a h' ≠ ḣ eltérés ≈ δt/2·q̈N ≈ 2.5e-4).
- S ≈ 0 közben is ~1e-4 a hiba → nem a szabályozó hibája, hanem a diszkretizációé.
- Nagyobb Λ elnyomja (Λ² erősítés) → 3.4.
- Ez magyarázza a Duffing sweep ≈1e-4 szintjét is.

### 3.4 Λ hatása – állandósult hiba (VdP, default kezdőérték, K_VSSM=500)

| Λ | maxe_ss |
|---|---|
| 1 | 9.5e-5 |
| 5 | 5.1e-6 |
| 10 | 1.3e-6 |
| 20 | 3.4e-7 |
| 56 | 4.9e-8 |
| 100 | 1.8e-8 |

- ≈ 1/Λ² skálázás; Λ tízszerese ≈ két nagyságrend.
- Julia futás megerősítette (log hiba ábra: ~1e-4 / ~1.3e-6 / ~5e-8 Λ = 1 / 10 / 56).
- Default kezdőértéknél a többi ábrán (pozíció, sebesség, u, fázistér) nincs látható különbség – csak a hibán.

### 3.5 S(t) – nincs chattering
- VdP default: S sima, ±2.5e-6, a referencia periódusával hullámzik → a tanh határrétegén belüli csúszó mód, **nem chattering** (az gyors, lépésről lépésre ugráló jel lenne).
- t=0 tüske (~2e-4): az első 3 lépés, amikor az RFPT még nem aktív – rendben van.
- **Figyelem:** S definíciója tartalmazza Λ-t (`S = Λ²∫h + 2Λh + ḣ`), ezért különböző Λ-k S-értékei közvetlenül nem összehasonlíthatók (nagyobb Λ → nagyobb S ugyanakkora hibához). Λ=56-nál S ~5.6e-5, de sima.

### 3.6 RFPT működése (VdP)
- **Vezérlőjel:** az actual u pontosan ráfekszik az ideális u-ra → az RFPT teljesen kompenzálja a modellhibát (pl. mₐ=0.8 vs mₑ=1). Erős ábra a dolgozatba.
- **Desired vs deformed:** a desired ≈ nominális gyorsulás (kis hiba), a deformed jelentősen eltér tőle – ennyit kell az RFPT-nek torzítania, hogy a valós gyorsulás a desired legyen.

### 3.7 Tranziens: pontosság ↔ gyorsaság ↔ beavatkozás (VdP, kezdőérték 3.2, −1.7)

| Λ | K_VSSM | beállás (\|h\|<1e-2) | maxe_ss | max \|u\| |
|---|---|---|---|---|
| 1 | 500 | 7.7 s | 9.5e-5 | 422 |
| 2 | 500 | 3.9 s | 3.0e-5 | 443 |
| 3 | 500 | 2.6 s | 1.4e-5 | 448 |
| 5 | 500 | 1.6 s | 5.1e-6 | 479 |
| 10 | 500 | 0.9 s | 1.3e-6 | 662 |
| 56 | 500 | 0.8 s | 4.9e-8 | 8131 |
| 1 | 50 | 7.7 s | 1.0e-4 | 32 |
| 3 | 50 | 2.9 s | 1.6e-5 | 43 |
| 5 | 50 | 2.0 s | 6.2e-6 | 44 |
| 10 | 50 | 1.8 s | 1.8e-6 | 237 |

(Python-replika számolta, a Julia ábrák egyeznek vele.)

Megfigyelések a Julia ábrákon:
- **Λ=1:** sima, „szép” beállás, de lassú (~7 s) és a végső hiba a legnagyobb.
- **Λ=10:** gyors, de nagy túllövés (hiba 0.9-ig, pálya −0.7-ig).
- **Λ=56 (opt.):** leggyorsabb, legpontosabb, de hirtelen; a hibán lapos plató 0.2–0.75 s között = a tanh telítődése (elérési fázis, a VSSM teljes erővel tol). Óriási u-csúcs.
- **Λ=3, K_VSSM=50:** nincs nagy kezdeti ugrás az u-ban (~43), ~1.5–2 s alatt ráfekszik az ideális u-ra, sima → **legkiegyensúlyozottabb.**
- **A u-csúcsot főleg a K_VSSM okozza, nem a Λ.**

**Következtetés:** nincs egyetlen „legjobb” paraméter → kompromisszum. Javaslat: modellenként két paraméterkészlet:
- **„pontos”** – pl. VdP: Λ=56, K_VSSM=1000 (min. hiba, nagy beavatkozási csúcs)
- **„kíméletes”** – pl. VdP: Λ=3, K_VSSM=50 (kis, sima vezérlőjel, gyors beállás, hiba ~1e-5)

### 3.8 Robusztusság – Duffing (default: Λ=1, K_VSSM=500, default kezdőérték) ✅ Julia
Ellenőrzés: `ic_run(:duffing)` → maxe_ss = 9.437e-5, ugyanaz, mint a módosítás előtt (a `mismatch`/`noise` alapértéken nem változtat semmit).

**Modellhiba** – a becsült paraméterek pₐ = (1+ε)·pₑ (ε=0 pontos modell, ε=−0.9 a valós 10%-a, ε=2 a háromszorosa), `maxe_ss`:

| ε | hibrid | csak RFPT | csak VSSM |
|---|---|---|---|
| −0.9 | 9.45e-5 | 1.32e-3 | 3.49e-3 |
| −0.5 | 9.44e-5 | 7.47e-4 | 1.92e-3 |
| −0.2 | 9.43e-5 | 3.30e-4 | 7.26e-4 |
| 0 | 9.43e-5 | 2.49e-4 | 9.43e-5 |
| 0.5 | 9.42e-5 | 7.18e-4 | 2.07e-3 |
| 1.0 | 9.41e-5 | 1.44e-3 | 4.09e-3 |
| 2.0 | 9.39e-5 | 2.88e-3 | 8.18e-3 |

(teljes táblázat 0.1-es lépéssel: `robust_sweep(:duffing)`; a Python-előrejelzéssel pontosan egyezik)

- **Hibrid:** a hiba a teljes −90% … +200% tartományban **gyakorlatilag állandó** (9.39–9.45e-5, a numerikus padló) → a dolgozat fő állítása.
- **Csak VSSM:** pontos modellnél ugyanolyan jó, mint a hibrid, de a hiba **lineárisan nő |ε|-nal** (ε=2-nél ~87×).
- **Csak RFPT** (PID kinematikai blokkal): pontos modellnél is ~2.6× rosszabb, és szintén lineárisan romlik (ε=2-nél ~30×).
- Egyik eset sem divergált a tartományban.

**Mérési zaj** – szórás a mért q, q̇, q̈ értékekre (a régi, rögzített αₐ, δₐ, βₐ modellhibával), `maxe_ss`:

| zaj | hibrid | csak RFPT | csak VSSM |
|---|---|---|---|
| 0 | 9.44e-5 | 3.08e-4 | 2.77e-4 |
| 1e-6 | 9.44e-5 | 3.08e-4 | 2.77e-4 |
| 1e-5 | 9.43e-5 | 3.07e-4 | 2.77e-4 |
| 1e-4 | 9.52e-5 | 3.01e-4 | 2.74e-4 |
| 1e-3 | 1.27e-4 | 2.60e-4 | 2.82e-4 |
| 1e-2 | 1.10e-3 | 1.05e-3 | 1.34e-3 |

- **Hibrid:** 1e-4 zajszintig érzéketlen; 1e-3-nál +35%; 1e-2-nél ~10× (≈1e-3) – ekkor a zaj dominál, mindhárom eset ~1e-3.
- A hiba ~0.1× a zaj szórása → a szabályozó szűri a zajt (1e-2 zaj = a ±2 amplitúdó 0.5%-a).
- A csak RFPT enyhe javulása 1e-3-ig valószínűleg véletlen (egyetlen zajminta, `seed=1`) – nem érdemes értelmezni; ha kell, több seed átlaga.
- `ic_run(:duffing; noise=1e-3)`: maxe_ss = 1.27e-4, maxu 5.5 → 9.9.

**Zajos futás ábrái (noise=1e-3) – zajerősítés, nem chattering:**
- Gyorsulás, desired, deformed, u: ±2–4 zaj. Ok: S ≈ ±0.005 (a mért hibából), S ≪ w=1 → a tanh lineáris → a VSSM tag ≈ 500·S ≈ ±2.5. A K_VSSM=500 nagy erősítés, a mérési zajt felnagyítja.
- Nem chattering: az a kapcsolásból eredő határciklus lenne; ez véletlen zaj egy lineáris erősítőn át.
- Sebesség: kis zaj (±0.002); pozíció, hiba: ±1.5e-4 – a rendszer kétszer integrál → aluláteresztő szűrő, a követés jó marad.

K_VSSM hatása zajjal (Python, noise=1e-3, utolsó 20%):

| K_VSSM | maxe_ss zaj nélkül | maxe_ss zajjal | u-zaj szórása |
|---|---|---|---|
| 500 | 9.4e-5 | 1.2e-4 | 1.28 |
| 100 | 9.6e-5 | 1.2e-4 | 0.23 |
| 50 | 9.9e-5 | 1.2e-4 | 0.11 |
| 10 | 1.2e-4 | 1.2e-4 | 0.02 |

→ K_VSSM=50: ~12× kisebb vezérlőjel-zaj, gyakorlatilag azonos hiba. Megerősíti a „kíméletes” paraméterkészletet (3.7).
Ellenőrzés Juliában: `ic_run(:duffing; noise=1e-3, tag="K=500")`, `ic_run(:duffing; noise=1e-3, K_VSSM=50, tag="K=50")`.

**Fogalmak:**
- `mismatch=ε`: csak a szabályozó *becsült* paramétereit tolja el (αₐ, δₐ, βₐ = (1+ε)·(αₑ, δₑ, βₑ)); a valós rendszer, az erősítések (Λ, K_VSSM, K, A, B), a kezdőérték és a referencia változatlan. Default (nothing): a kézi 0.8 / 0.1 / 0.9 (−20% / −50% / −10%).
- `noise=σ`: minden lépésben új, 0 várható értékű, σ szórású véletlen szám a *mért* q, q̇ és a mért q̈ (RFPT) értékére. A valódi állapotot nem rontja közvetlenül – a szabályozó számol rosszabb u-t, és ez hat a rendszerre. Az ábrák/metrikák a valódi hibát mutatják. `seed`: ugyanaz a seed → ugyanaz a zajsorozat.
- `robust_sweep(model, levels; param, metric, ic, tail, …)`: minden szintre 3 futás (hibrid / csak RFPT / csak VSSM), a mutató log skálán a szint függvényében.
- `tail`: az állandósult ablak hossza (default 0.2 = utolsó 20%); `tail=0.5` is megadható bármelyik függvénynél. Figyelni: a tranziens ne lógjon bele.

### 3.9 Robusztusság – VdP, Lorenz, Rössler (default paraméterek, default kezdőérték) ✅ Julia
Ellenőrzés (alapértéken változatlan): VdP maxe_ss = 9.447e-5, Lorenz 2.599e-4, Rössler 9.807e-5.

**Modellhiba** (`robust_sweep(:modell)`, pₐ = (1+ε)·pₑ minden paraméterre), `maxe_ss`, NaN = divergált:

| ε | VdP hibrid | VdP RFPT | VdP VSSM | Lorenz hibrid | Lorenz RFPT | Lorenz VSSM | Rössler hibrid | Rössler RFPT | Rössler VSSM |
|---|---|---|---|---|---|---|---|---|---|
| −0.9 | 1.08e-4 | 2.21e-2 | 6.04e-2 | 3.02e-3 | 9.82e-2 | NaN | 3.74e-4 | 1.34e-2 | 6.90e-2 |
| −0.5 | 9.38e-5 | 2.54e-3 | 7.32e-3 | 1.11e-3 | 4.42e-2 | 0.399 | 2.21e-4 | 7.92e-3 | 3.92e-2 |
| 0 | 9.43e-5 | 2.49e-4 | 9.43e-5 | 9.67e-5 | 3.82e-3 | 1.06e-5 | 9.67e-5 | 3.82e-3 | 1.06e-5 |
| 0.5 | 9.45e-5 | 8.41e-4 | 2.65e-3 | 1.17e-3 | 4.36e-2 | 0.413 | 2.25e-4 | 8.14e-3 | 4.15e-2 |
| 0.6 | 2.49e-4 | 9.49e-4 | 2.99e-3 | 1.51e-3 | 5.35e-2 | 0.881 | 2.61e-4 | 9.42e-3 | 5.01e-2 |
| 0.8 | 4.65e-3 | 1.13e-3 | 3.55e-3 | 2.48e-3 | 7.88e-2 | NaN | 3.36e-4 | 1.21e-2 | 6.77e-2 |
| 1.0 | 3.57e-2 | NaN | 4.02e-3 | 4.46e-3 | 0.123 | NaN | 4.13e-4 | 1.49e-2 | 8.57e-2 |
| 1.1 | NaN | NaN | 4.22e-3 | 6.94e-3 | 0.217 | NaN | 4.51e-4 | 1.62e-2 | 9.49e-2 |
| 1.2 | NaN | NaN | 4.40e-3 | NaN | NaN | NaN | 4.90e-4 | 1.76e-2 | 0.104 |
| 2.0 | NaN | NaN | 5.50e-3 | NaN | NaN | NaN | 7.97e-4 | 2.90e-2 | 0.184 |

Megfigyelések:
- **Rössler:** a hibrid a teljes tartományban stabil és messze a legjobb (VSSM-nél ~200×, RFPT-nél ~35× kisebb hiba); a hiba enyhén, lineárisan nő |ε|-nal.
- **Lorenz:** a hibrid −0.9 … 1.1 között stabil, 30–100× jobb az RFPT-nél; a csak VSSM már ε=−0.8 és ε=0.7-nél divergál. **ε ≥ 1.2: a hibrid is szétesik** – oka még nincs ellenőrizve (hipotézis: a szükséges korrekció mértéke meghaladja azt, amit az RFPT lépésenként követni tud, K=1e2; próba: `robust_sweep(:lorenz, 1.0:0.1:2.0; K=1e3)`).
  - **K=1e3 eredménye (A változatlan 1.97e-3):** a csak RFPT ε=2-ig stabil lett és pontosabb (ε=1: 8.2e-3 vs 0.123), de a **hibrid már ε=1.0-tól divergál** (K=1e2-vel ε=1.0 még 4.5e-3 volt). → A hipotézis így nem igazolható: a K növelése az A·K-t is 0.197-ről 1.97-re növelte, ami a kontrakciós határ (2) közelébe viszi az RFPT-t; a VSSM nagy erősítésű desired jelével együtt ez instabil.
  - A Lorenznél a bemeneti erősítés 1 (nincs tömegparaméter), ezért az ε nem a VdP-s gain-feltételen keresztül hat; valószínűbb ok: az állapotfüggő kompenzációs hiba (pl. ε·ρ·x, ε·x·y) nagy, a kaotikus csatolás miatt gyorsan nő.
  - **K=1e3, A=1.97e-4 (A·K = 0.197 változatlan) ✅:** a hibrid **ε=2-ig stabil**, és pontosabb is:

    | ε | hibrid (K=1e2) | hibrid (K=1e3, A·K azonos) | csak RFPT (K=1e3) |
    |---|---|---|---|
    | 1.0 | 4.46e-3 | 1.96e-3 | 8.63e-2 |
    | 1.1 | 6.94e-3 | 2.15e-3 | 9.47e-2 |
    | 1.2 | NaN | 2.35e-3 | 0.103 |
    | 1.5 | NaN | 2.93e-3 | 0.128 |
    | 2.0 | NaN | 3.94e-3 | 0.169 |

  - **Magyarázat:** a K az RFPT „munkatartománya” – a deformált jel nagyságrendjének a K alatt kell maradnia (a leképezés (x+K)-val skáláz). Nagy ε-nál a Lorenz állapotfüggő kompenzációs hibája (ε·ρ·x, ε·x·y, x ~ 20) a szükséges korrekciót K=100 fölé viszi → szétesik. Nagyobb K-val (az A·K szorzat megtartásával) ez a korlát kitolódik.
  - **Két külön szerepű paraméter:** A·K → konvergencia (kontrakció, VdP-n igazolva); K → a kompenzálható modellhiba nagysága (Lorenzen igazolva). A korábbi K=1e3, A=1.97e-3 próba azért esett szét, mert az A·K is 1.97-re nőtt (a kontrakciós határ közelébe).
- **VdP:** a hibrid ε ≤ 0.5-ig a numerikus padlón, **ε=0.6-tól romlik, ε ≥ 1.1 divergál** → RFPT konvergenciahatár (lásd lent). A csak VSSM (RFPT nélkül) végig stabil marad → ott a VSSM a „biztonsági háló”.
- ε=0-nál a Lorenz és a Rössler azonos értékeket ad (VSSM 1.06e-5, RFPT 3.82e-3): pontos modellnél a kompenzáció után a hibadinamika mindkettőnél ugyanaz, és a referencia is azonos → konzisztencia-ellenőrzésnek jó.
- ε=0-nál a csak VSSM a legjobb (pontos modell: nincs mit tanulni, az RFPT csak „zajt” ad) – de ez a valóságban nem fordul elő.

**VdP szétesés oka – RFPT konvergenciafeltétel (Python-ellenőrzés):**

| ε | minden paraméter | csak mₐ | minden más (mₐ pontos) |
|---|---|---|---|
| 0.5 | 9.45e-5 | 9.43e-5 | 9.46e-5 |
| 0.6 | 2.49e-4 | 1.43e-4 | 9.47e-5 |
| 0.8 | 4.65e-3 | 1.04e-3 | 9.48e-5 |
| 1.0 | 3.57e-2 | 2.62e-3 | 9.49e-5 |
| 1.5 | NaN | NaN | 9.53e-5 |
| 2.0 | NaN | NaN | 9.57e-5 |

Csak mₐ/mₑ arány: 1.5 → 9.4e-5, 1.6 → 1.4e-4, 1.8 → 1.0e-3, 2.0 → 2.6e-3, **2.1 → divergál**.

- **Kizárólag a tömeg becslése (mₐ) okozza**; a többi paraméter akár 3× tévedése sem számít.
- Magyarázat: az RFPT iteráció a fixpont közelében ≈ lineáris, a kontrakciós tényező ≈ 1 − A·K·(mₐ/mₑ). Itt A·K = 1e-5·1e5 = 1 → konvergencia, ha |1 − mₐ/mₑ| < 1, azaz **0 < mₐ/mₑ < 2**. A szimuláció pontosan ezt adja (2.1-nél divergál); 1.6 fölött a lassú, váltakozó előjelű konvergencia miatt már romlik a pontosság (a desired minden lépésben változik).
- A Duffingban nincs tömegparaméter (a bemenet erősítése a modellben és a valóságban is 1) → a tényező ≈ 0 → ezért nem függ a hibrid a modellhibától.
- Következmény: az RFPT robusztus a „lassú” paraméterekre, de a **bemeneti erősítés (tömeg) becslésére korlátos** – ez az elméleti konvergenciafeltétel szimulációs igazolása.
- **Julia-val igazolva:** `robust_sweep(:vanderpol; A=5e-6)` (A·K = 0.5 → elméleti határ mₐ/mₑ < 4): a hibrid a **teljes −0.9 … 2.0 tartományban stabil**, ~9.4e-5 (a numerikus padló). Egyedül ε=−0.9-nél nagyobb (1.7e-4): ott mₐ/mₑ = 0.1 → kontrakciós tényező ≈ 1 − 0.05 = 0.95, lassú konvergencia. A csak RFPT is végig stabil lett.
- **Tervezési szabály:** az A·K szorzatot a bemeneti erősítés becsült tévedéséhez kell választani: A·K·(mₐ/mₑ) ∈ (0, 2). Kisebb A·K → szélesebb robusztussági tartomány, de lassabb tanulás (kis mₐ/mₑ-nél).

**Mérési zaj** (`param=:noise`), `maxe_ss`:

| σ | VdP hibrid | VdP RFPT | VdP VSSM | Lorenz hibrid | Lorenz RFPT | Lorenz VSSM | Rössler hibrid | Rössler RFPT | Rössler VSSM |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 9.45e-5 | 3.32e-4 | 8.39e-4 | 2.60e-4 | 1.12e-2 | 8.53e-2 | 9.81e-5 | 3.78e-3 | 3.95e-3 |
| 1e-5 | 9.44e-5 | 3.33e-4 | 8.39e-4 | 2.64e-4 | 1.12e-2 | 8.53e-2 | 1.01e-4 | 3.78e-3 | 3.95e-3 |
| 1e-4 | 9.53e-5 | 3.41e-4 | 8.35e-4 | 3.12e-4 | 1.12e-2 | 8.53e-2 | 1.36e-4 | 3.79e-3 | 3.98e-3 |
| 1e-3 | 1.27e-4 | 4.21e-4 | 8.04e-4 | 7.98e-4 | 1.15e-2 | 8.54e-2 | 6.28e-4 | 3.82e-3 | 4.40e-3 |
| 1e-2 | 1.11e-3 | 1.22e-3 | 1.57e-3 | 7.04e-3 | 1.84e-2 | 8.82e-2 | 6.45e-3 | 4.40e-3 | 9.15e-3 |

- A hibrid 1e-4 zajig ~érzéketlen (VdP), Lorenz/Rössler-nél 1e-4-től már nő; nagyobb zajnál a hiba ≈ 0.6–0.7 × σ (Lorenz, Rössler) ill. ≈ 0.1 × σ (Duffing, VdP – a másodrendű rendszer kétszer integrál, jobban szűr).
- Rössler σ=1e-2: a csak RFPT (4.4e-3) jobb, mint a hibrid (6.5e-3) → a VSSM tag zajerősítése (K_VSSM·tanh, lásd 3.8); kisebb K_VSSM-mel valószínűleg javul.
- Divergálás zaj miatt nem volt.

---

## 4. Modellek és paraméterek

| Modell | Default | „Optimális” (sweep) | max error (opt.) |
|---|---|---|---|
| Duffing | Λ=1, K_VSSM=500, w=1 | Λ=100, K_VSSM=1000 | 4.887e-7 (GAS) |
| Van der Pol | Λ=1, K_VSSM=500, w=1 | Λ=56, K_VSSM=1000 | 1.257e-6 (GAS) |
| Lorenz | Λ=1, K_VSSM=50, w=1 | Λ=10, K_VSSM=100 | 2.532e-3 (UUB, attraktor korlát) |
| Rössler | Λ=1, K_VSSM=50, w=1 | Λ=10, K_VSSM=1000 | 1.002e-4 (GAS) |

RFPT: Duffing/VdP K=1e5, A=1e-5; Lorenz/Rössler K=1e2, A=1.97e-3; B=−1.

Default (kézi) modellhibák (ₐ vs ₑ):
- Duffing: α 0.8/1.0, δ 0.1/0.2, β 0.9/1.0
- VdP: μ 0.5/0.4, ω 0.42/0.46, α 0.9/1.0, λ 0.09/0.1, m 0.8/1.0
- Lorenz: σ 4/5, ρ 36/40, β (7/3)/(8/3)
- Rössler: a 0.1/0.01 (**10×!**), b 0.3/0.2, c 5.5/5.7

⚠️ **Ellenőrizendő:** a régi `src/main.jl` saját, eltérő modellmásolatokat tartalmaz (nem a `my_*.jl`-t használja):
- Duffing: `u = qDef_pp − αₐq` (a `my_duffing.jl`-ben teljes kompenzáció), és a rendszerben `+δₑq̇` vs `−δₑq̇`.
- VdP: `u = mₐ·qDef_pp` (a `my_vanderpool.jl`-ben teljes kompenzáció).
- Lorenz, Rössler: azonos.

Ha az „optimális” max error értékek a `main.jl`-ből jöttek, a Duffing és VdP értékeket a `my_*.jl` modellel újra kell mérni. A `main.jl` nincs használatban.

---

## 5. A szimulációs fejezet terve

1. **Alapműködés** – modellenként optimális futás, minden ábra (2. pont).
2. **Abláció** – csak RFPT / csak VSSM / hibrid → a hibrid indoklása.
3. **Paraméter-érzékenység** – Λ × K_VSSM, A·K (megvan, `run_all_sweeps`).
4. **Kezdeti feltétel vizsgálat** – `maxe_ss` + beállási idő a kezdőérték függvényében.
5. **Pontosság ↔ beavatkozás kompromisszum** – „pontos” vs „kíméletes” paraméterkészlet (3.7).
6. **Robusztusság** – modellhiba skálázása (pl. 0/20/50/100%), mérési zaj, külső zavarás.
7. **Összesítő táblázat** – modellenként: maxe_ss, RMS, beállási idő, max|u|, GAS/UUB.

---

## 6. Állapot

- [x] Duffing – ábrák, metrics, ic_run / ic_sweep
- [x] Van der Pol – ábrák, metrics, Λ-összehasonlítás, tranziens kompromisszum
- [x] Lorenz – ábrák (desired a sebességen, 3D, S, ideális u), metrics – tesztelve, jó
- [ ] Rössler – ugyanaz, mint Lorenz → tesztelésre vár
- [ ] Duffing: „pontos” vs „kíméletes” paraméterkészlet
- [ ] Optimális értékek újramérése a `my_*.jl` modellekkel (4. pont)
- [ ] Abláció, összesítő táblázat
- [x] Robusztusság: Duffing (`mismatch`, `noise`) – Julia-val lefuttatva
- [x] Robusztusság: VdP, Lorenz, Rössler – Julia-val lefuttatva (3.9)
- [x] Lorenz ε ≥ 1.2 szétesés oka: az RFPT K munkatartománya; K=1e3, A=1.97e-4 → ε=2-ig stabil
- [x] VdP: RFPT-határ igazolva kisebb A·K-val (`A=5e-6` → ε=2-ig stabil)

---

## 7. Kódváltozások (napló)

- `my_duffing.jl`, `my_vanderpool.jl`: `tag`, `acc_clip`, hiba log ábra, ideális u, S(t), desired vs deformed; VdP gyorsulás ábra visszakapcsolva; VdP `metrics`.
- `my_lorenz.jl`: `metrics`, `tag`, `vel_clip`; desired + deformed a sebesség szinten; ideális u; S(t); hibanorma log; numerikus gyorsulás ábra elhagyva; ábrák egymásra rajzolhatók (`axes3`).
- `my_rossler.jl`: ugyanazok a változások, mint a Lorenznél.
- `my_duffing.jl`: `mismatch`, `noise`, `seed` (robusztusság); a szabályozó a mért (zajos) állapotot használja, az ábrák/metrikák a valódi hibát.
- `my_vanderpool.jl`, `my_lorenz.jl`, `my_rossler.jl`: `mismatch`, `noise`, `seed` – mint a Duffingnál; Lorenz/Rössler: a mért x, y, z és a mért sebesség (RFPT válasz) zajos.
- `runner.jl`: `plot_heatmap_2d(...; xlab, ylab)` – Lorenz/Rössler sweepnél x₀, y₀ tengelyfelirat.
- `runner.jl`: `plot_heatmap_2d(...; logscale, label)`.
- `scripts/initial_condition_sweeps.jl`: `ic_run`, `ic_sweep`, `pstr`, `MODELS`; a `rundefault()` az új függvényeket használja.
