"""
═══════════════════════════════════════════════════════════════════
  OFFICIAL VALIDATION BENCHMARK  ─  AD-HTC Thermodynamic Solver
═══════════════════════════════════════════════════════════════════
  Uses the EXACT same solver equations found in app.py / script.js,
  but implemented inline so we avoid tkinter/GUI imports.

  Benchmark inputs:
    T1 = 300 K,  rp = 12,  ηc = 0.85,  ηt = 0.88
    TIT (T3) = 1450 K,  ṁ = 40 kg/s
  ──────────────────────────────────────────────────────────────
  Imports ONLY: thermodynamics.py  (no tkinter, no matplotlib)
═══════════════════════════════════════════════════════════════════
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd                             # already required by app.py deps
from thermodynamics import (GasTable, SteamTable,
                            calculate_hrsg, validate_results)

# ── Terminal colours ───────────────────────────────────────────────
GREEN = "\033[92m"; RED = "\033[91m"; YELLOW = "\033[93m"
CYAN  = "\033[96m"; BOLD = "\033[1m"; RESET  = "\033[0m"

passed = [0]; failed = [0]

def ok(msg):     print(f"  {GREEN}✔ PASS{RESET}  {msg}")
def fail(msg):   print(f"  {RED}✗ FAIL{RESET}  {msg}")
def header(msg): print(f"\n{BOLD}{CYAN}{'═'*62}{RESET}\n{BOLD}  {msg}{RESET}\n{'─'*62}")

def check(label, value, lo, hi, unit=""):
    tag = f"{label} = {value:.3f}{unit}  [expected {lo}–{hi}{unit}]"
    if lo <= value <= hi:
        ok(tag); passed[0] += 1
    else:
        fail(tag); failed[0] += 1

def check_bool(label, ok_cond, detail=""):
    if ok_cond:
        ok(label + (f"  ({detail})" if detail else "")); passed[0] += 1
    else:
        fail(label + (f"  ({detail})" if detail else "")); failed[0] += 1


# ═══════════════════════════════════════════════════════════════════
#  SOLVER — mirrors calculate_gas_cycle() in app.py exactly
# ═══════════════════════════════════════════════════════════════════

def run_gas_cycle(T1, P1, rp, TIT, eta_c, eta_t, eta_cc, m_air, LHV):
    states = []
    h1 = GasTable.h(T1);  s1 = GasTable.s(T1, P1)
    states.append({'State':'1','Location':'Air Inlet',
                   'T(K)':T1,'P(kPa)':P1,'h(kJ/kg)':h1,'s(kJ/kgK)':s1})

    P2 = P1 * rp
    T2s = GasTable.T_isentropic(T1, P1, P2)
    h2s = GasTable.h(T2s)
    h2  = h1 + (h2s - h1) / eta_c
    T2  = T2s + (T2s - T1) * (1.0/eta_c - 1.0)
    s2  = GasTable.s(T2, P2)
    states.append({'State':'2','Location':'Compressor Exit',
                   'T(K)':T2,'P(kPa)':P2,'h(kJ/kg)':h2,'s(kJ/kgK)':s2})

    P3 = P2
    states.append({'State':'3','Location':'Combustor Inlet',
                   'T(K)':T2,'P(kPa)':P3,'h(kJ/kg)':h2,'s(kJ/kgK)':s2})

    T4 = TIT;  P4 = P2 * 0.96
    h4 = GasTable.h(T4);  s4 = GasTable.s(T4, P4)
    states.append({'State':'4','Location':'Combustor Exit (TIT)',
                   'T(K)':T4,'P(kPa)':P4,'h(kJ/kg)':h4,'s(kJ/kgK)':s4})

    P5  = P1 * 1.02
    T5s = GasTable.T_isentropic(T4, P4, P5)
    h5s = GasTable.h(T5s)
    h5  = h4 - eta_t * (h4 - h5s)
    T5  = T5s + (T4 - T5s) * (1.0 - eta_t)
    s5  = GasTable.s(T5, P5)
    states.append({'State':'5','Location':'Turbine Exhaust',
                   'T(K)':T5,'P(kPa)':P5,'h(kJ/kg)':h5,'s(kJ/kgK)':s5})

    w_c = h2 - h1
    w_t = h4 - h5
    q_in_ideal = h4 - h2
    q_in = q_in_ideal / eta_cc
    w_net = w_t - w_c
    eta   = w_net / q_in * 100 if q_in > 0 else 0
    bwr   = w_c / w_t  * 100 if w_t > 0 else 0
    m_fuel = (q_in * m_air) / LHV
    return dict(states=states, df=pd.DataFrame(states),
                w_c=w_c, w_t=w_t, q_in=q_in, w_net=w_net,
                eta=eta, bwr=bwr, m_fuel=m_fuel, T_exhaust=T5,
                T2=T2, m_air=m_air, LHV=LHV,
                cp_avg=GasTable.cp_avg(T1,T4),
                gamma_avg=GasTable.gamma_avg(T1,T4))


def run_steam_cycle(P_boiler, T_steam, P_cond, eta_st, eta_fp):
    states = []
    Pa, Ta, Pb = P_boiler, T_steam, P_cond

    ha = SteamTable.h_super(Pa, Ta);  sa = SteamTable.s_super(Pa, Ta)
    states.append({'State':'a','Location':'Boiler Exit (SH)',
                   'T(K)':Ta,'P(kPa)':Pa,'h(kJ/kg)':ha,'s(kJ/kgK)':sa,'x':'-'})

    sb_s = sa
    sg_c = SteamTable.sg(Pb);  sf_c = SteamTable.sf(Pb)
    if sb_s < sg_c:
        xb_s  = SteamTable.x_from_s(Pb, sb_s)
        hb_s  = SteamTable.h_from_x(Pb, xb_s)
        Tb    = SteamTable.T_sat(Pb)
    else:
        Tb    = SteamTable.T_from_s_super(Pb, sb_s)
        hb_s  = SteamTable.h_super(Pb, Tb)

    hb    = ha - eta_st * (ha - hb_s)
    hf_c  = SteamTable.hf(Pb);  hfg_c = SteamTable.hfg(Pb)
    xb    = (hb - hf_c) / hfg_c if hfg_c > 0 else 1.0
    if xb <= 1.0:
        sb    = sf_c + xb * (sg_c - sf_c)
        Tb    = SteamTable.T_sat(Pb);  xb_str = f'{xb:.3f}'
    else:
        sb    = SteamTable.s_super(Pb, Tb);   xb_str = '-'
    states.append({'State':'b','Location':'ST Exit',
                   'T(K)':Tb,'P(kPa)':Pb,'h(kJ/kg)':hb,'s(kJ/kgK)':sb,'x':xb_str})

    Tc = SteamTable.T_sat(Pb);  hc = SteamTable.hf(Pb);  sc = SteamTable.sf(Pb)
    states.append({'State':'c','Location':'Condenser Exit',
                   'T(K)':Tc,'P(kPa)':Pc,'h(kJ/kg)':hc,'s(kJ/kgK)':sc,'x':'0.000'})

    vf_val = SteamTable.vf(Pb);  wp_s = vf_val * (Pa - Pb)
    wp = wp_s / eta_fp;  hd = hc + wp;  sd = sc + 0.001;  Td = Tc + wp / 4.18
    states.append({'State':'d','Location':'Feed Pump Exit',
                   'T(K)':Td,'P(kPa)':Pa,'h(kJ/kg)':hd,'s(kJ/kgK)':sd,'x':'-'})

    w_st = ha - hb;  w_fp = wp;  q_boiler = ha - hd;  q_cond = hb - hc
    w_net = w_st - w_fp
    eta   = w_net / q_boiler * 100 if q_boiler > 0 else 0
    return dict(states=states, df=pd.DataFrame(states),
                w_st=w_st, w_fp=w_fp, q_boiler=q_boiler,
                q_cond=q_cond, w_net=w_net, eta=eta)


def run_ad_htc(m_biomass, moisture_split, ad_yield, LHV, m_fuel):
    m_rich    = m_biomass * moisture_split
    m_lean    = m_biomass * (1 - moisture_split)
    biogas_vol = m_rich * ad_yield        # m³/s
    rho_biogas = 1.15
    m_biogas   = biogas_vol * rho_biogas  # kg/s
    E_biogas   = m_biogas * LHV           # kJ/s = kW
    E_demand   = m_fuel * LHV
    surplus    = m_biogas - m_fuel
    return dict(biogas_vol=biogas_vol, m_biogas=m_biogas,
                E_biogas=E_biogas, E_demand=E_demand, surplus=surplus)


# ═══════════════════════════════════════════════════════════════════
#  FIX TYPO in steam function (Pc not defined — use Pb)
# ═══════════════════════════════════════════════════════════════════
# Patch the typo inline:
import types as _types
_orig = run_steam_cycle
def run_steam_cycle(P_boiler, T_steam, P_cond, eta_st, eta_fp):
    states = []
    Pa, Ta, Pb = P_boiler, T_steam, P_cond
    ha = SteamTable.h_super(Pa, Ta);  sa = SteamTable.s_super(Pa, Ta)
    states.append({'State':'a','Location':'Boiler Exit (SH)',
                   'T(K)':Ta,'P(kPa)':Pa,'h(kJ/kg)':ha,'s(kJ/kgK)':sa,'x':'-'})
    sb_s = sa
    sg_c = SteamTable.sg(Pb);  sf_c = SteamTable.sf(Pb)
    if sb_s < sg_c:
        xb_s = SteamTable.x_from_s(Pb, sb_s);  hb_s = SteamTable.h_from_x(Pb, xb_s)
        Tb   = SteamTable.T_sat(Pb)
    else:
        Tb   = SteamTable.T_from_s_super(Pb, sb_s);  hb_s = SteamTable.h_super(Pb, Tb)
    hb = ha - eta_st * (ha - hb_s)
    hf_c = SteamTable.hf(Pb);  hfg_c = SteamTable.hfg(Pb)
    xb = (hb - hf_c) / hfg_c if hfg_c > 0 else 1.0
    if xb <= 1.0:
        sb = sf_c + xb*(sg_c-sf_c);  Tb = SteamTable.T_sat(Pb);  xb_str = f'{xb:.3f}'
    else:
        sb = SteamTable.s_super(Pb, Tb);  xb_str = '-'
    states.append({'State':'b','Location':'ST Exit',
                   'T(K)':Tb,'P(kPa)':Pb,'h(kJ/kg)':hb,'s(kJ/kgK)':sb,'x':xb_str})
    Tc = SteamTable.T_sat(Pb);  hc = SteamTable.hf(Pb);  sc = SteamTable.sf(Pb)
    states.append({'State':'c','Location':'Condenser Exit',
                   'T(K)':Tc,'P(kPa)':Pb,'h(kJ/kg)':hc,'s(kJ/kgK)':sc,'x':'0.000'})
    vf_val = SteamTable.vf(Pb);  wp_s = vf_val*(Pa-Pb)
    wp = wp_s/eta_fp;  hd = hc+wp;  sd = sc+0.001;  Td = Tc+wp/4.18
    states.append({'State':'d','Location':'Feed Pump Exit',
                   'T(K)':Td,'P(kPa)':Pa,'h(kJ/kg)':hd,'s(kJ/kgK)':sd,'x':'-'})
    w_st = ha-hb;  w_fp = wp;  q_boiler = ha-hd;  q_cond = hb-hc
    w_net = w_st-w_fp;  eta = w_net/q_boiler*100 if q_boiler > 0 else 0
    return dict(states=states, df=pd.DataFrame(states),
                w_st=w_st, w_fp=w_fp, q_boiler=q_boiler,
                q_cond=q_cond, w_net=w_net, eta=eta)


# ═══════════════════════════════════════════════════════════════════
#  RUN BENCHMARK
# ═══════════════════════════════════════════════════════════════════

# ─── Fixed benchmark inputs ───────────────────────────────────────
T1, P1, rp, TIT   = 300, 101.325, 12, 1450
eta_c, eta_t       = 0.85, 0.88
eta_cc             = 1.0
m_air              = 40        # kg/s
LHV_fuel           = 50000     # kJ/kg (nat gas approx)

P_boiler = 4000; T_steam = 673; P_cond = 10
eta_st = 0.85;  eta_fp = 0.90

m_biomass = 6; moisture_split = 0.5; ad_yield = 0.45

# ─── Run solvers ─────────────────────────────────────────────────
gas   = run_gas_cycle(T1, P1, rp, TIT, eta_c, eta_t, eta_cc, m_air, LHV_fuel)
steam = run_steam_cycle(P_boiler, T_steam, P_cond, eta_st, eta_fp)

cp_exh  = GasTable.cp_avg(gas['T_exhaust'], T1)
hrsg    = calculate_hrsg(gas['T_exhaust'], 420, m_air, cp_exh, 0.85, 15, steam['q_boiler'])
m_steam = hrsg['m_steam']

ad = run_ad_htc(m_biomass, moisture_split, ad_yield, LHV_fuel, gas['m_fuel'])


# ────────────────────────────────────────────────────────────────
header("PART 1 — BRAYTON (Gas) CYCLE")

T2    = gas['T2']
T5    = gas['T_exhaust']
Wc    = gas['w_c']  * m_air / 1e3   # MW
Wt    = gas['w_t']  * m_air / 1e3   # MW
Wn    = gas['w_net']* m_air / 1e3   # MW
Qin   = gas['q_in'] * m_air / 1e3   # MW
eta_B = gas['eta']

print(f"""
  T1  = {T1} K          (input)
  T2  = {T2:.1f} K       (compressor exit)
  TIT = {TIT} K          (turbine inlet - input)
  T5  = {T5:.1f} K       (turbine exhaust)
  W_compressor = {Wc:.2f} MW
  W_turbine    = {Wt:.2f} MW
  W_net (gas)  = {Wn:.2f} MW
  Q_in (heat)  = {Qin:.2f} MW
  η_Brayton    = {eta_B:.1f} %
  Cp_avg       = {gas['cp_avg']:.4f} kJ/kg·K
  γ_avg        = {gas['gamma_avg']:.4f}
""")

check("T2 — compressor exit",   T2,    655,  675,  " K")
check("Wc — compressor work",   Wc,     14,   15,  " MW")
check("T5 — turbine exhaust",   T5,    820,  860,  " K")
check("Wt — turbine work",      Wt,     24,   26,  " MW")
check("Wnet — gas net power",   Wn,      9,   11,  " MW")
check("Qin — heat input",       Qin,    30,   32,  " MW")
check("η_Brayton",              eta_B,  30,   36,  " %")


# ────────────────────────────────────────────────────────────────
header("PART 2 — HRSG HEAT RECOVERY")

Q_avail = hrsg['Q_available'] / 1e3   # MW
Q_recov = hrsg['Q_recovered'] / 1e3   # MW

print(f"""
  T_exhaust    = {T5:.1f} K
  T_stack      = {hrsg['T_stack_actual']:.1f} K
  Q_available  = {Q_avail:.2f} MW
  Q_recovered  = {Q_recov:.2f} MW  (ε = 0.85)
  m_steam      = {m_steam:.2f} kg/s
  Pinch OK     = {hrsg['pinch_ok']}
""")

check("Q_available (exhaust heat)", Q_avail, 15, 17,  " MW")
check("Q_recovered (HRSG output)",  Q_recov, 13, 14,  " MW")
check_bool("HRSG output < Gas Q_in (physical sanity)",
           Q_recov < Qin,
           f"{Q_recov:.2f} < {Qin:.2f} MW")


# ────────────────────────────────────────────────────────────────
header("PART 3 — STEAM (RANKINE) CYCLE")

Wst   = steam['w_st']  * m_steam / 1e3   # MW
Wst_kg = steam['w_net']                   # kJ/kg
eta_R = steam['eta']

print(f"""
  w_st (per kg steam) = {steam['w_st']:.1f} kJ/kg
  w_net (per kg)      = {Wst_kg:.1f} kJ/kg
  m_steam             = {m_steam:.2f} kg/s
  W_ST (MW)           = {Wst:.2f} MW
  η_Rankine           = {eta_R:.1f} %
""")

check("W_net,steam",  Wst,    4,   6,  " MW")
check("η_Rankine",    eta_R, 25,  35,  " %")
check_bool("W_steam < Q_recovered (physical sanity)",
           Wst < Q_recov,
           f"{Wst:.2f} < {Q_recov:.2f} MW")


# ────────────────────────────────────────────────────────────────
header("PART 4 — AD BIOGAS ENERGY SUPPLY")

E_biogas_MW = ad['E_biogas'] / 1e3
E_demand_MW = ad['E_demand'] / 1e3

print(f"""
  Biomass feed    = {m_biomass} kg/s
  moisture_split  = {moisture_split}  →  m_rich = {m_biomass*moisture_split} kg/s
  ad_yield        = {ad_yield} m³/kg
  Biogas volume   = {ad['biogas_vol']:.3f} m³/s  (expected ≈ 2.7)
  m_biogas        = {ad['m_biogas']:.3f} kg/s   (expected ≈ 3.24)
  LHV used        = {LHV_fuel/1000:.0f} MJ/kg (nat gas; benchmark uses 21 MJ/kg for 68 MW)
  E_biogas        = {E_biogas_MW:.1f} MW
  Fuel demand     = {gas['m_fuel']:.3f} kg/s
  E_demand        = {E_demand_MW:.1f} MW
  Surplus (mass)  = {ad['surplus']:.3f} kg/s
""")

check("Biogas volume",   ad['biogas_vol'],  2.5,  2.9,  " m³/s")
check("m_biogas",        ad['m_biogas'],    2.9,  3.6,  " kg/s")
check_bool("AD supplies > turbine fuel demand",
           ad['E_biogas'] > gas['q_in'] * m_air,
           f"E_biogas={E_biogas_MW:.1f} MW vs Q_in={Qin:.1f} MW")


# ────────────────────────────────────────────────────────────────
header("PART 5 — COMBINED CYCLE")

W_comb     = Wn + Wst
eta_comb   = W_comb / Qin * 100 if Qin > 0 else 0

print(f"""
  W_net,gas    = {Wn:.2f} MW
  W_net,steam  = {Wst:.2f} MW
  W_combined   = {W_comb:.2f} MW
  η_combined   = {eta_comb:.1f} %
""")

check("W_net,gas",    Wn,       9,  11,  " MW")
check("W_net,steam",  Wst,      4,   6,  " MW")
check("W_combined",   W_comb,  14,  16,  " MW")
check("η_combined",   eta_comb, 45,  58,  " %")


# ════════════════════════════════════════════════════════════════
header("VALIDATION SUMMARY")
total = passed[0] + failed[0]
print(f"  Total checks : {total}")
print(f"  {GREEN}Passed{RESET}       : {passed[0]}")
print(f"  {RED}Failed{RESET}       : {failed[0]}")
if failed[0] == 0:
    print(f"\n  {GREEN}{BOLD}✔ ALL {total} CHECKS PASSED — Solver is structurally CORRECT!{RESET}")
else:
    pct = passed[0]/total*100
    print(f"\n  {YELLOW}{BOLD}  {passed[0]}/{total} ({pct:.0f}%) checks passed{RESET}")
    print(f"  {RED}{BOLD}✗ {failed[0]} check(s) failed — see details above{RESET}")
print(f"{'═'*62}\n")
