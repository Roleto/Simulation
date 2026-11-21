# Paraméter-sweepek – leírás

Cél: a hibrid szabályozó (RFPT + csúszó módú) paramétereinek megválasztása és a stabilitási
határ numerikus ellenőrzése mind a négy rendszerre (Duffing, van der Pol, Lorenz, Rössler).
Minden futásból egy számot jegyzünk fel: a maximális követési hibát. A divergált futás NaN.

## Futtatás

```julia
# julia --project=.   (a Simulation mappából)
using Revise
includet("src/runner.jl")
includet("scripts/run_all_sweeps.jl")
run_all_sweeps()                    # mind a 4 modell
run_all_sweeps(models=["Lorenz"])   # egy modell
run_all_sweeps(cap=1.0)             # szigorúbb divergálás-küszöb (alap: 1e3)
```

Egyedi sweep:

```julia
P, E = sweep_AK_1d(run_duffing, 10 .^ range(-3, 0.5, length=25), 1e5)
save_AK_1d("Duffing", P, E)

a, b, E = sweep_params(run_lorenz, :Λ, 10 .^ range(0, 2, length=9), :K_VSSM, 10 .^ range(0, 4, length=9))
save_sweep_2d("Lorenz", a, b, E, :Λ, :K_VSSM)
```

## Mit csinálnak a sweepek

| Sweep | Paraméterek | Kimenet | Mit mutat |
|---|---|---|---|
| `A·K` (1D) | `A = P/K`, K fix | `AK_1d.csv/.png` | Az adaptáció lépésköze: túl kicsi → lassú; plató → jó; túl nagy → instabil. A stabilitási feltétel (`0 < A·K·f′ < 2`) numerikus ellenőrzése. |
| `Λ × K_VSSM` (2D) | Λ, K_VSSM | `sweep_Λ_K_VSSM.csv/.png` | A csúszó módú rész két paramétere együtt; hol van a jó tartomány és a törés (chattering). |
| `K_VSSM × w` (2D) | K_VSSM, w | `sweep_K_VSSM_w.csv/.png` | A tanh simítás (w) és a K_VSSM kölcsönhatása. |

Az `A·K` egydimenziós, mert az RFPT lépésköze kis hibánál ≈ `A·K` (a tanh linearizálható), így csak a
szorzat számít.

## Rácsok

- `A·K`: 10^(−3) … 10^(0,5), 25 pont (logaritmikus)
- Λ: 10^0 … 10^2, 9 pont
- K_VSSM: 10^0 … 10^4, 9 pont
- w: 10^(−2) … 10^1, 7 pont

## Alapértékek (a modellfájlokból)

| | K | A | Λ | K_VSSM | w |
|---|---|---|---|---|---|
| Duffing | 1e5 | 1e-5 | 1 | 500 | 1 |
| van der Pol | 1e5 | 1e-5 | 1 | 500 | 1 |
| Lorenz | 1e2 | 1,97e-3 | 1 | 50 | 1 |
| Rössler | 1e2 | 1,97e-3 | 1 | 50 | 1 |

B = −1 fix, K (=K_c) rögzített, A hangolandó.

## Kimenet

`data/Sweeps/<modell>/` (modell: `Duffing`, `VanDerPol`, `Lorenz`, `Rossler`). A CSV-kben a 2D mátrix
sorai a második, oszlopai az első paraméter értékei; üres/NaN cella = divergált vagy `cap` fölötti hiba.
A hőtérképen a szín a hiba tízes alapú logaritmusa (sötét = jó), a szürke cella a használhatatlan futás.
Az `A·K` ábra tengelyei logaritmikusak, a feliratok a tényleges értékek.

## Hibamutató és korlátok

- Mindenhol a maximális követési hiba (Lorenz/Rössler: a három koordináta maximumainak euklideszi normája).
- A divergálás-szűrő a Lorenznél és a Rösslernél gyenge: a `cap=1e3` küszöb a 0,1–1 nagyságrendű, de
  százszoros hibát nem fogja meg. Szigorítás: `cap=1.0`.
- A Duffingnál a részletes mutatók (`rms_ss`, `maxe_ss`, `rms_all`, `maxe_all`, `maxu`, `tail` ablak)
  a `scripts/sweeps.jl`-ben vannak (`sw_duffing_L_KVSSM()`, `sw_duffing_KVSSM_w()`), és az `ok` feltétel
  a vezérlőjelet is ellenőrzi (`maxu < 100`).

## Eddigi eredmények (Duffing, `scripts/sweeps.jl`)

- `A·K` plató ≈ 0,03–1,5, hiba ≈ 1e-4; a stabilitási határ ≈ 1,6–2.
- `Λ × K_VSSM`: K_VSSM ≤ 1000-ig a hiba Λ-val és K_VSSM-mel monoton csökken (az optimum a rács szélén);
  K_VSSM ≥ 3162 felett chattering (`maxu` ≈ 2800 fölött). Javasolt: Λ ≈ 10–30, K_VSSM ≈ 100–300.
- `A·K` javasolt: ≈ 0,1–0,5 (a határtól távol).
