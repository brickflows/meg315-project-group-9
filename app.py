"""
AD-HTC Fuel-Enhanced Power Gas Cycle — Main Application
========================================================
Tkinter GUI with embedded matplotlib charts and schematic.
Generates H-S chart for HTC steam cycle and T-Ḣ chart for
the gas power cycle on clicking 'Analyse'.

© 2025 Energhx Research Group — Faculty of Engineering, University of Lagos
"""

import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, ArrowStyle
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

from thermodynamics import GasTable, SteamTable


# ═══════════════════════════════════════════════════════
#  CYCLE CALCULATIONS
# ═══════════════════════════════════════════════════════

def calculate_gas_cycle(params):
    """
    Brayton gas cycle: Air → Compressor → Biogas Combustor → Turbine → Exhaust.
    Returns dict with state points and performance.
    """
    T1, P1 = params['T1'], params['P1']
    rp = params['rp']
    TIT = params['TIT']
    eta_c, eta_t = params['eta_c'], params['eta_t']

    states = []

    # State 1: Air inlet
    h1 = GasTable.h(T1)
    s1 = GasTable.s(T1, P1)
    states.append({'State': '①', 'Location': 'Air Inlet',
                   'T(K)': T1, 'P(kPa)': P1, 'h(kJ/kg)': h1, 's(kJ/kg·K)': s1})

    # State 2: Compressor exit
    P2 = P1 * rp
    T2s = GasTable.T_isentropic(T1, P1, P2)
    h2s = GasTable.h(T2s)
    h2 = h1 + (h2s - h1) / eta_c
    T2 = T2s + (T2s - T1) * (1.0 / eta_c - 1.0)
    s2 = GasTable.s(T2, P2)
    states.append({'State': '②', 'Location': 'Compressor Exit',
                   'T(K)': T2, 'P(kPa)': P2, 'h(kJ/kg)': h2, 's(kJ/kg·K)': s2})

    # State 3: Combustor inlet (same as state 2)
    P3 = P2
    states.append({'State': '③', 'Location': 'Combustor Inlet',
                   'T(K)': T2, 'P(kPa)': P3, 'h(kJ/kg)': h2, 's(kJ/kg·K)': s2})

    # State 4: Combustor exit (TIT)
    T4 = TIT
    P4 = P2 * 0.96  # ~4% pressure drop
    h4 = GasTable.h(T4)
    s4 = GasTable.s(T4, P4)
    states.append({'State': '④', 'Location': 'Combustor Exit (TIT)',
                   'T(K)': T4, 'P(kPa)': P4, 'h(kJ/kg)': h4, 's(kJ/kg·K)': s4})

    # State 5: Turbine exhaust
    P5 = P1 * 1.02
    T5s = GasTable.T_isentropic(T4, P4, P5)
    h5s = GasTable.h(T5s)
    h5 = h4 - eta_t * (h4 - h5s)
    T5 = T5s + (T4 - T5s) * (1.0 - eta_t)
    s5 = GasTable.s(T5, P5)
    states.append({'State': '⑤', 'Location': 'Turbine Exhaust',
                   'T(K)': T5, 'P(kPa)': P5, 'h(kJ/kg)': h5, 's(kJ/kg·K)': s5})

    w_c = h2 - h1
    w_t = h4 - h5
    q_in = h4 - h2
    w_net = w_t - w_c
    eta = w_net / q_in * 100 if q_in > 0 else 0
    bwr = w_c / w_t * 100 if w_t > 0 else 0
    m_fuel = (q_in * params['m_air']) / params['LHV'] if params['LHV'] > 0 else 0

    return {
        'states': states,
        'df': pd.DataFrame(states),
        'w_c': w_c, 'w_t': w_t, 'q_in': q_in, 'w_net': w_net,
        'eta': eta, 'bwr': bwr, 'm_fuel': m_fuel, 'T_exhaust': T5,
    }


def calculate_steam_cycle(params):
    """
    HTC Steam (Rankine) cycle: Boiler → Steam Turbine → Condenser → Feed Pump.
    """
    Pa = params['P_boiler']  # kPa
    Ta = params['T_steam']   # K
    Pb = params['P_cond']    # kPa
    eta_st = params['eta_st']
    eta_fp = params['eta_fp']

    states = []

    # State a: Boiler exit (superheated)
    ha = SteamTable.h_super(Pa, Ta)
    sa = SteamTable.s_super(Pa, Ta)
    states.append({'State': 'ⓐ', 'Location': 'Boiler Exit (Superheated)',
                   'T(K)': Ta, 'P(kPa)': Pa, 'h(kJ/kg)': ha, 's(kJ/kg·K)': sa, 'x': '—'})

    # State b: ST exhaust (isentropic + efficiency)
    sb_s = sa
    sg_c = SteamTable.sg(Pb)
    sf_c = SteamTable.sf(Pb)
    if sb_s < sg_c:
        xb_s = SteamTable.x_from_s(Pb, sb_s)
        hb_s = SteamTable.h_from_x(Pb, xb_s)
        Tb = SteamTable.T_sat(Pb)
    else:
        Tb = SteamTable.T_from_s_super(Pb, sb_s)
        hb_s = SteamTable.h_super(Pb, Tb)

    hb = ha - eta_st * (ha - hb_s)
    hf_c = SteamTable.hf(Pb)
    hfg_c = SteamTable.hfg(Pb)
    xb = (hb - hf_c) / hfg_c if hfg_c > 0 else 1.0
    if xb <= 1.0:
        sb = sf_c + xb * (sg_c - sf_c)
        Tb = SteamTable.T_sat(Pb)
        xb_str = f'{xb:.3f}'
    else:
        sb = SteamTable.s_super(Pb, Tb)
        xb_str = '—'

    states.append({'State': 'ⓑ', 'Location': 'ST Exit',
                   'T(K)': Tb, 'P(kPa)': Pb, 'h(kJ/kg)': hb, 's(kJ/kg·K)': sb, 'x': xb_str})

    # State c: Condenser exit (saturated liquid)
    Tc = SteamTable.T_sat(Pb)
    hc = SteamTable.hf(Pb)
    sc = SteamTable.sf(Pb)
    states.append({'State': 'ⓒ', 'Location': 'Condenser Exit (Sat. Liq.)',
                   'T(K)': Tc, 'P(kPa)': Pb, 'h(kJ/kg)': hc, 's(kJ/kg·K)': sc, 'x': '0.000'})

    # State d: Feed pump exit
    vf_val = SteamTable.vf(Pb)
    wp_s = vf_val * (Pa - Pb)
    wp = wp_s / eta_fp
    hd = hc + wp
    sd = sc + 0.001
    Td = Tc + wp / 4.18
    states.append({'State': 'ⓓ', 'Location': 'Feed Pump Exit',
                   'T(K)': Td, 'P(kPa)': Pa, 'h(kJ/kg)': hd, 's(kJ/kg·K)': sd, 'x': '—'})

    w_st = ha - hb
    w_fp = wp
    q_boiler = ha - hd
    q_cond = hb - hc
    w_net = w_st - w_fp
    eta = w_net / q_boiler * 100 if q_boiler > 0 else 0

    return {
        'states': states,
        'df': pd.DataFrame(states),
        'w_st': w_st, 'w_fp': w_fp, 'q_boiler': q_boiler,
        'q_cond': q_cond, 'w_net': w_net, 'eta': eta,
    }


def calculate_ad_htc(params):
    """AD-HTC mass and energy balance."""
    m_total = params['m_biomass']
    m_rich = m_total * params['moisture_split']
    m_lean = m_total * (1 - params['moisture_split'])
    biogas_vol = m_rich * params['ad_yield']
    biogas_density = 1.15  # kg/m³
    m_biogas = biogas_vol * biogas_density
    htc_energy = m_lean * 1.5 * (params['htc_temp'] - 298)
    return {
        'm_total': m_total, 'm_rich': m_rich, 'm_lean': m_lean,
        'biogas_vol': biogas_vol, 'm_biogas': m_biogas,
        'htc_energy': htc_energy,
    }


# ═══════════════════════════════════════════════════════
#  PLOTTING FUNCTIONS
# ═══════════════════════════════════════════════════════

DARK_BG = '#0d0d1a'
CARD_BG = '#141428'
TEXT_COL = '#e2e8f0'
GRID_COL = '#1e293b'
ACCENT_CYAN = '#38bdf8'
ACCENT_PURPLE = '#a78bfa'
ACCENT_ORANGE = '#fb923c'
ACCENT_GREEN = '#4ade80'
ACCENT_RED = '#f87171'
ACCENT_BLUE = '#60a5fa'
ACCENT_YELLOW = '#fbbf24'
ACCENT_TEAL = '#14b8a6'


def draw_schematic(fig):
    """Draw the AD-HTC cycle schematic on a matplotlib figure."""
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 9)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_facecolor(DARK_BG)
    fig.patch.set_facecolor(DARK_BG)

    box_kw = dict(boxstyle='round,pad=0.3', linewidth=1.5)
    text_kw = dict(ha='center', va='center', fontsize=8, fontweight='bold',
                   color='white', family='sans-serif')
    label_kw = dict(ha='center', va='center', fontsize=6.5, color='#94a3b8',
                    family='sans-serif')
    arrow_kw = dict(arrowstyle='->', mutation_scale=14, linewidth=1.8)

    # ── Gas cycle section box ──
    gas_box = FancyBboxPatch((0.5, 0.3), 11, 2.8, boxstyle='round,pad=0.2',
                              facecolor='#0f172a', edgecolor='#1e3a5f',
                              linewidth=1, linestyle='--', alpha=0.5)
    ax.add_patch(gas_box)
    ax.text(6, 0.5, 'AD-HTC FUEL-ENHANCED GAS POWER CYCLE (BRAYTON)',
            ha='center', fontsize=7, color=ACCENT_BLUE, alpha=0.6, fontweight='bold')

    # ── Compressor ──
    comp = FancyBboxPatch((1.5, 1.2), 2, 1.2, facecolor='#2563eb', edgecolor='#3b82f6', **box_kw)
    ax.add_patch(comp)
    ax.text(2.5, 1.85, 'Compressor', **text_kw)
    ax.text(2.5, 1.55, 'C', fontsize=7, ha='center', va='center', color='#bfdbfe')

    # Air inlet
    ax.annotate('', xy=(1.5, 1.8), xytext=(0.5, 1.8),
                arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=1.5))
    ax.text(0.6, 1.45, 'Air', **label_kw)

    # State ① and ②
    ax.plot(1.5, 1.8, 'o', color=ACCENT_CYAN, markersize=6, zorder=5)
    ax.text(1.3, 2.1, '①', fontsize=8, color=ACCENT_CYAN, fontweight='bold')
    ax.plot(3.5, 1.8, 'o', color=ACCENT_CYAN, markersize=6, zorder=5)
    ax.text(3.65, 2.1, '②', fontsize=8, color=ACCENT_CYAN, fontweight='bold')

    # Shaft dashed line
    ax.plot([3.5, 8.5], [1.8, 1.8], '--', color='#64748b', linewidth=2)
    ax.text(6, 1.55, '← Shaft →', fontsize=6, ha='center', color='#64748b')

    # ── Biogas Combustion Chamber ──
    comb = FancyBboxPatch((4.5, 3.8), 3, 1.2, facecolor='#d97706', edgecolor='#f59e0b', **box_kw)
    ax.add_patch(comb)
    ax.text(6, 4.45, 'Biogas Combustion', **text_kw)
    ax.text(6, 4.15, 'Chamber', fontsize=7, ha='center', va='center', color='#fef3c7')

    # Compressor → Combustor (up)
    ax.annotate('', xy=(4.5, 4.4), xytext=(3.5, 1.8),
                arrowprops=dict(arrowstyle='->', color=ACCENT_BLUE, lw=2,
                                connectionstyle='arc3,rad=-0.2'))
    ax.plot(4.5, 4.4, 'o', color=ACCENT_YELLOW, markersize=6, zorder=5)
    ax.text(4.2, 4.6, '③', fontsize=8, color=ACCENT_YELLOW, fontweight='bold')

    # Combustor → Turbine (down-right)
    ax.annotate('', xy=(8.5, 1.8), xytext=(7.5, 4.4),
                arrowprops=dict(arrowstyle='->', color=ACCENT_ORANGE, lw=2,
                                connectionstyle='arc3,rad=0.2'))
    ax.plot(7.5, 4.4, 'o', color=ACCENT_YELLOW, markersize=6, zorder=5)
    ax.text(7.65, 4.65, '④', fontsize=8, color=ACCENT_YELLOW, fontweight='bold')

    # ── Turbine ──
    turb = FancyBboxPatch((8.5, 1.2), 2, 1.2, facecolor='#dc2626', edgecolor='#ef4444', **box_kw)
    ax.add_patch(turb)
    ax.text(9.5, 1.85, 'Turbine', **text_kw)
    ax.text(9.5, 1.55, 'GT', fontsize=7, ha='center', va='center', color='#fecaca')

    # Exhaust
    ax.annotate('', xy=(11.5, 1.8), xytext=(10.5, 1.8),
                arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=1.5))
    ax.text(11.2, 1.45, 'Exhaust\nGases', **label_kw)
    ax.plot(10.5, 1.8, 'o', color=ACCENT_RED, markersize=6, zorder=5)
    ax.text(10.65, 2.1, '⑤', fontsize=8, color=ACCENT_RED, fontweight='bold')

    # ── HTC Steam Cycle section ──
    htc_box = FancyBboxPatch((1.2, 5.0), 5, 2.5, boxstyle='round,pad=0.2',
                              facecolor='#0f0f2a', edgecolor='#312e81',
                              linewidth=1, linestyle='--', alpha=0.4)
    ax.add_patch(htc_box)
    ax.text(3.7, 7.2, 'HTC STEAM CYCLE', fontsize=7, ha='center',
            color=ACCENT_PURPLE, alpha=0.6, fontweight='bold')

    # ── Boiler ──
    boiler = FancyBboxPatch((3.5, 6.3), 2, 0.8, facecolor='#be123c', edgecolor='#e11d48', **box_kw)
    ax.add_patch(boiler)
    ax.text(4.5, 6.7, 'Boiler', **text_kw)

    # ── HTC Reactor ──
    reactor = FancyBboxPatch((3.5, 5.2), 2, 0.8, facecolor='#ea580c', edgecolor='#f97316', **box_kw)
    ax.add_patch(reactor)
    ax.text(4.5, 5.6, 'Reactor', **text_kw)
    ax.text(4.5, 5.35, 'HTC', fontsize=6.5, ha='center', va='center', color='#fed7aa')

    # Steam turbine symbol (small triangle)
    st_pts = np.array([[2.3, 5.8], [2.8, 6.3], [2.3, 6.8]])
    st_tri = plt.Polygon(st_pts, facecolor='#7c3aed', edgecolor='#a78bfa', linewidth=1.5)
    ax.add_patch(st_tri)
    ax.text(2.4, 6.3, 'ST', fontsize=6.5, ha='center', va='center', color='white', fontweight='bold')

    # State points for steam cycle
    ax.plot(3.5, 6.7, 'o', color=ACCENT_PURPLE, markersize=5, zorder=5)
    ax.text(3.3, 6.9, 'ⓑ', fontsize=7, color=ACCENT_PURPLE, fontweight='bold')
    ax.plot(5.5, 6.7, 'o', color=ACCENT_PURPLE, markersize=5, zorder=5)
    ax.text(5.6, 6.9, 'ⓒ', fontsize=7, color=ACCENT_PURPLE, fontweight='bold')
    ax.plot(3.5, 5.6, 'o', color=ACCENT_ORANGE, markersize=5, zorder=5)
    ax.text(3.25, 5.4, 'ⓐ', fontsize=7, color=ACCENT_ORANGE, fontweight='bold')

    # Boiler → ST
    ax.annotate('', xy=(2.8, 6.5), xytext=(3.5, 6.7),
                arrowprops=dict(arrowstyle='->', color=ACCENT_PURPLE, lw=1.5))
    # Reactor ↔ Boiler
    ax.annotate('', xy=(4.5, 6.3), xytext=(4.5, 6.0),
                arrowprops=dict(arrowstyle='<->', color='#94a3b8', lw=1.2))
    # Condensate return
    ax.plot([2.3, 1.8, 1.8, 2.3], [5.8, 5.8, 6.8, 6.8], '--',
            color=ACCENT_PURPLE, linewidth=1, alpha=0.4)
    ax.text(1.6, 6.3, 'Cond.\nReturn', fontsize=5.5, ha='center', color=ACCENT_PURPLE, alpha=0.5)

    # Volatile output from reactor
    ax.annotate('', xy=(4.5, 4.6), xytext=(4.5, 5.2),
                arrowprops=dict(arrowstyle='->', color='#94a3b8', lw=1.2))
    ax.text(4.5, 4.8, 'Volatile Matters\n& Feedstock Waste', fontsize=5.5,
            ha='center', color='#94a3b8')

    # ── Biomass Feedstock Homogenizer (top-left) ──
    homo = FancyBboxPatch((1.5, 7.8), 2.5, 0.8, facecolor='#7c3aed', edgecolor='#a78bfa', **box_kw)
    ax.add_patch(homo)
    ax.text(2.75, 8.25, 'Biomass Feedstock', fontsize=7, ha='center',
            va='center', color='white', fontweight='bold')
    ax.text(2.75, 8.0, 'Homogenizer', fontsize=6.5, ha='center', va='center', color='#ddd6fe')

    # Biomass input
    ax.annotate('', xy=(1.5, 8.2), xytext=(0.5, 8.2),
                arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=1.8))
    ax.text(0.4, 8.5, 'Biomass\nFeedstock', fontsize=6, color='#94a3b8')

    # Homogenizer → Moisture-lean → HTC Reactor
    ax.annotate('', xy=(3.5, 5.6), xytext=(2.75, 7.8),
                arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=1.5,
                                connectionstyle='arc3,rad=-0.1'))
    ax.text(2.3, 7.0, 'Moisture\n-lean', fontsize=5.5, color=ACCENT_GREEN)

    # ── AD (Anaerobic Digestion) (top-right) ──
    ad = FancyBboxPatch((7.5, 7.3), 1.2, 0.8, facecolor='#0891b2', edgecolor='#06b6d4', **box_kw)
    ax.add_patch(ad)
    ax.text(8.1, 7.7, 'AD', **text_kw)

    # Homogenizer → AD (moisture-rich)
    ax.annotate('', xy=(7.5, 8.0), xytext=(4.0, 8.2),
                arrowprops=dict(arrowstyle='->', color='#22d3ee', lw=1.5))
    ax.text(5.8, 8.5, 'Moisture-rich Biomass', fontsize=6, ha='center', color='#22d3ee')

    # ── Enhanced Biogas Collector ──
    coll = FancyBboxPatch((7.0, 5.8), 2.2, 0.8, facecolor='#0d9488', edgecolor='#14b8a6', **box_kw)
    ax.add_patch(coll)
    ax.text(8.1, 6.25, 'Enhanced Biogas', fontsize=7, ha='center',
            va='center', color='white', fontweight='bold')
    ax.text(8.1, 6.0, 'Collector', fontsize=6.5, ha='center', va='center', color='#99f6e4')

    # AD → Collector
    ax.annotate('', xy=(8.1, 6.6), xytext=(8.1, 7.3),
                arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=1.5))

    # Collector → Combustion Chamber (biogas fuel)
    ax.annotate('', xy=(6, 5.0), xytext=(7.5, 6.0),
                arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=2,
                                connectionstyle='arc3,rad=0.2'))
    ax.text(7.0, 5.3, 'Biogas\nFuel', fontsize=6, color=ACCENT_GREEN, fontweight='bold')

    # Collector → Biogas Distribution
    ax.annotate('', xy=(10.5, 6.2), xytext=(9.2, 6.2),
                arrowprops=dict(arrowstyle='->', color=ACCENT_GREEN, lw=1.5))
    ax.text(10.6, 6.5, 'Biogas Distribution', fontsize=6, color=ACCENT_GREEN)
    ax.text(10.6, 6.2, 'to Building Envelopes', fontsize=5.5, color=ACCENT_GREEN, alpha=0.7)

    # ── Footer ──
    ax.text(6, 0.1, '©2025 Energhx Research Group — Faculty of Engineering, University of Lagos',
            ha='center', fontsize=5.5, color='#475569')

    fig.tight_layout(pad=0.5)


def plot_hs_chart(fig, steam_result):
    """H-S diagram for HTC steam (Rankine) cycle."""
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_facecolor(DARK_BG)
    fig.patch.set_facecolor(DARK_BG)

    states = steam_result['states']
    sm = {s['State']: s for s in states}

    # Saturation dome
    sf_arr, hf_arr, sg_arr, hg_arr = SteamTable.saturation_dome(60)
    dome_s = np.concatenate([sf_arr, sg_arr[::-1], [sf_arr[0]]])
    dome_h = np.concatenate([hf_arr, hg_arr[::-1], [hf_arr[0]]])
    ax.fill(dome_s, dome_h, alpha=0.05, color='#94a3b8')
    ax.plot(dome_s, dome_h, '--', color='#94a3b8', alpha=0.35, linewidth=1.2, label='Saturation Dome')

    # Process lines
    processes = [
        ('Boiler (ⓓ→ⓐ)', 'ⓓ', 'ⓐ', ACCENT_ORANGE),
        ('Steam Turbine (ⓐ→ⓑ)', 'ⓐ', 'ⓑ', ACCENT_PURPLE),
        ('Condenser (ⓑ→ⓒ)', 'ⓑ', 'ⓒ', ACCENT_CYAN),
        ('Feed Pump (ⓒ→ⓓ)', 'ⓒ', 'ⓓ', ACCENT_BLUE),
    ]
    for name, id_from, id_to, color in processes:
        s1, s2 = sm[id_from], sm[id_to]
        s_vals = np.linspace(s1['s(kJ/kg·K)'], s2['s(kJ/kg·K)'], 25)
        h_vals = np.linspace(s1['h(kJ/kg)'], s2['h(kJ/kg)'], 25)
        ax.plot(s_vals, h_vals, '-', color=color, linewidth=2.5, label=name)

    # State markers
    for sid in ['ⓓ', 'ⓐ', 'ⓑ', 'ⓒ']:
        st = sm[sid]
        ax.plot(st['s(kJ/kg·K)'], st['h(kJ/kg)'], 'o', color=ACCENT_CYAN,
                markersize=10, markeredgecolor='white', markeredgewidth=1.5, zorder=5)
        ax.annotate(f" {sid} {st['Location']}", (st['s(kJ/kg·K)'], st['h(kJ/kg)']),
                    fontsize=7, color=TEXT_COL, fontweight='bold')

    ax.set_xlabel('Entropy s [kJ/(kg·K)]', color=TEXT_COL, fontsize=11, fontweight='600')
    ax.set_ylabel('Enthalpy h [kJ/kg]', color=TEXT_COL, fontsize=11, fontweight='600')
    ax.set_title('H-S Diagram — HTC Steam (Rankine) Cycle', color=TEXT_COL,
                 fontsize=13, fontweight='bold', pad=12)
    ax.legend(fontsize=7.5, facecolor=CARD_BG, edgecolor='#334155',
              labelcolor=TEXT_COL, loc='best')
    ax.tick_params(colors='#64748b')
    ax.grid(True, alpha=0.1, color=ACCENT_CYAN)
    for spine in ax.spines.values():
        spine.set_color('#1e293b')
    fig.tight_layout(pad=1.0)


def plot_th_chart(fig, gas_result):
    """T-Ḣ diagram for gas (Brayton) cycle processes."""
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_facecolor(DARK_BG)
    fig.patch.set_facecolor(DARK_BG)

    states = gas_result['states']
    procs = [
        ('Compression (①→②)', 0, 1, ACCENT_BLUE),
        ('Combustion (③→④)', 2, 3, ACCENT_ORANGE),
        ('Expansion (④→⑤)', 3, 4, ACCENT_RED),
    ]

    cum_h = 0
    markers_x, markers_y, markers_labels = [0], [states[0]['T(K)']], [states[0]['State']]

    for name, i_from, i_to, color in procs:
        s1, s2 = states[i_from], states[i_to]
        dh = abs(s2['h(kJ/kg)'] - s1['h(kJ/kg)'])
        h_vals = np.linspace(cum_h, cum_h + dh, 30)
        T_vals = np.linspace(s1['T(K)'], s2['T(K)'], 30)
        ax.plot(h_vals, T_vals, '-', color=color, linewidth=2.5, label=name)
        cum_h += dh
        markers_x.append(cum_h)
        markers_y.append(s2['T(K)'])
        markers_labels.append(s2['State'])

    # State markers
    ax.plot(markers_x, markers_y, 'D', color=ACCENT_ORANGE, markersize=8,
            markeredgecolor='white', markeredgewidth=1.5, zorder=5, label='State Points')
    for x, y, lbl in zip(markers_x, markers_y, markers_labels):
        ax.annotate(f' {lbl}', (x, y), fontsize=8, color=TEXT_COL, fontweight='bold')

    ax.set_xlabel('Cumulative Heat / Work Rate Ḣ [kJ/kg]', color=TEXT_COL,
                  fontsize=11, fontweight='600')
    ax.set_ylabel('Temperature T [K]', color=TEXT_COL, fontsize=11, fontweight='600')
    ax.set_title('T-Ḣ Diagram — Gas Power Cycle (Brayton)', color=TEXT_COL,
                 fontsize=13, fontweight='bold', pad=12)
    ax.legend(fontsize=7.5, facecolor=CARD_BG, edgecolor='#334155',
              labelcolor=TEXT_COL, loc='best')
    ax.tick_params(colors='#64748b')
    ax.grid(True, alpha=0.1, color=ACCENT_ORANGE)
    for spine in ax.spines.values():
        spine.set_color('#1e293b')
    fig.tight_layout(pad=1.0)


# ═══════════════════════════════════════════════════════
#  TKINTER GUI APPLICATION
# ═══════════════════════════════════════════════════════

class ADHTCApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('AD-HTC Fuel-Enhanced Power Gas Cycle Analyser')
        self.geometry('1400x850')
        self.configure(bg=DARK_BG)
        self.minsize(1100, 700)

        # Style
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure('.', background=DARK_BG, foreground=TEXT_COL)
        self.style.configure('TFrame', background=DARK_BG)
        self.style.configure('TLabel', background=DARK_BG, foreground=TEXT_COL,
                             font=('Segoe UI', 9))
        self.style.configure('TLabelframe', background=CARD_BG, foreground=ACCENT_CYAN,
                             font=('Segoe UI', 10, 'bold'))
        self.style.configure('TLabelframe.Label', background=CARD_BG, foreground=ACCENT_CYAN)
        self.style.configure('TEntry', fieldbackground='#1a1a2e', foreground=TEXT_COL,
                             font=('Consolas', 9))
        self.style.configure('Accent.TButton', background='#2563eb', foreground='white',
                             font=('Segoe UI', 11, 'bold'), padding=(20, 10))
        self.style.map('Accent.TButton',
                       background=[('active', '#3b82f6'), ('pressed', '#1d4ed8')])
        self.style.configure('Treeview', background='#1a1a2e', foreground=TEXT_COL,
                             fieldbackground='#1a1a2e', font=('Consolas', 8),
                             rowheight=22)
        self.style.configure('Treeview.Heading', background=CARD_BG,
                             foreground=ACCENT_CYAN, font=('Segoe UI', 8, 'bold'))

        self._build_ui()
        self._draw_initial_schematic()

    def _build_ui(self):
        # ── Header ──
        header = tk.Frame(self, bg='#080816', height=70)
        header.pack(fill='x', side='top')
        header.pack_propagate(False)
        tk.Label(header, text='AD-HTC Fuel-Enhanced Power Gas Cycle',
                 font=('Segoe UI', 16, 'bold'), fg=ACCENT_CYAN, bg='#080816').pack(pady=(12, 0))
        tk.Label(header, text='Anaerobic Digestion + Hydrothermal Carbonization — Energhx Research Group',
                 font=('Segoe UI', 9), fg='#64748b', bg='#080816').pack()

        # ── Main paned window ──
        self.pw = tk.PanedWindow(self, orient='horizontal', bg=DARK_BG,
                                  sashwidth=4, sashrelief='flat')
        self.pw.pack(fill='both', expand=True, padx=8, pady=8)

        # Left: Schematic + Charts (notebook)
        self.left_frame = tk.Frame(self.pw, bg=DARK_BG)
        self.pw.add(self.left_frame, stretch='always', width=900)

        self.notebook = ttk.Notebook(self.left_frame)
        self.notebook.pack(fill='both', expand=True)

        # Tab 1: Schematic
        self.schematic_frame = tk.Frame(self.notebook, bg=DARK_BG)
        self.notebook.add(self.schematic_frame, text='  Schematic  ')
        self.fig_schematic = plt.Figure(figsize=(9, 6), facecolor=DARK_BG)
        self.canvas_schematic = FigureCanvasTkAgg(self.fig_schematic, self.schematic_frame)
        self.canvas_schematic.get_tk_widget().pack(fill='both', expand=True)

        # Tab 2: H-S Chart
        self.hs_frame = tk.Frame(self.notebook, bg=DARK_BG)
        self.notebook.add(self.hs_frame, text='  H-S Chart (Steam)  ')
        self.fig_hs = plt.Figure(figsize=(8, 5), facecolor=DARK_BG)
        self.canvas_hs = FigureCanvasTkAgg(self.fig_hs, self.hs_frame)
        self.canvas_hs.get_tk_widget().pack(fill='both', expand=True)
        self.toolbar_hs = NavigationToolbar2Tk(self.canvas_hs, self.hs_frame)
        self.toolbar_hs.update()

        # Tab 3: T-Ḣ Chart
        self.th_frame = tk.Frame(self.notebook, bg=DARK_BG)
        self.notebook.add(self.th_frame, text='  T-Ḣ Chart (Gas)  ')
        self.fig_th = plt.Figure(figsize=(8, 5), facecolor=DARK_BG)
        self.canvas_th = FigureCanvasTkAgg(self.fig_th, self.th_frame)
        self.canvas_th.get_tk_widget().pack(fill='both', expand=True)
        self.toolbar_th = NavigationToolbar2Tk(self.canvas_th, self.th_frame)
        self.toolbar_th.update()

        # Tab 4: Data Tables
        self.tables_frame = tk.Frame(self.notebook, bg=DARK_BG)
        self.notebook.add(self.tables_frame, text='  State Tables  ')

        # Right: Parameters panel
        self.right_frame = tk.Frame(self.pw, bg=CARD_BG)
        self.pw.add(self.right_frame, stretch='never', width=320)

        self._build_params_panel()

    def _build_params_panel(self):
        """Build the parameter inputs on the right side."""
        canvas = tk.Canvas(self.right_frame, bg=CARD_BG, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.right_frame, orient='vertical', command=canvas.yview)
        scroll_frame = tk.Frame(canvas, bg=CARD_BG)

        scroll_frame.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=scroll_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)

        # Enable mouse wheel scrolling
        def _on_mouse_wheel(event):
            canvas.yview_scroll(-1 * (event.delta // 120), 'units')
        canvas.bind_all('<MouseWheel>', _on_mouse_wheel)

        self.entries = {}

        def add_section(parent, title):
            lf = tk.LabelFrame(parent, text=title, bg=CARD_BG, fg=ACCENT_CYAN,
                               font=('Segoe UI', 9, 'bold'), padx=10, pady=5)
            lf.pack(fill='x', padx=10, pady=(8, 4))
            return lf

        def add_entry(parent, key, label, default):
            f = tk.Frame(parent, bg=CARD_BG)
            f.pack(fill='x', pady=2)
            tk.Label(f, text=label, bg=CARD_BG, fg='#94a3b8',
                     font=('Consolas', 7.5), anchor='w').pack(fill='x')
            e = tk.Entry(f, bg='#1a1a2e', fg=TEXT_COL, insertbackground=TEXT_COL,
                         font=('Consolas', 9), relief='flat', bd=0,
                         highlightthickness=1, highlightcolor=ACCENT_CYAN,
                         highlightbackground='#334155')
            e.insert(0, str(default))
            e.pack(fill='x', ipady=3)
            self.entries[key] = e

        # ── Gas Cycle ──
        sec = add_section(scroll_frame, '🔥 Gas Power Cycle (Brayton)')
        add_entry(sec, 'T1', 'Ambient Temperature T₁ (K)', 298)
        add_entry(sec, 'P1', 'Ambient Pressure P₁ (kPa)', 101.325)
        add_entry(sec, 'rp', 'Pressure Ratio rp', 10)
        add_entry(sec, 'TIT', 'Turbine Inlet Temp T₃ (K)', 1400)
        add_entry(sec, 'eta_c', 'Compressor Eff η_c', 0.86)
        add_entry(sec, 'eta_t', 'Turbine Eff η_t', 0.89)
        add_entry(sec, 'LHV', 'Biogas LHV (kJ/kg)', 20000)
        add_entry(sec, 'm_air', 'Air Mass Flow (kg/s)', 50)

        # ── Steam Cycle ──
        sec = add_section(scroll_frame, '💨 HTC Steam Cycle (Rankine)')
        add_entry(sec, 'P_boiler', 'Boiler Pressure (kPa)', 4000)
        add_entry(sec, 'T_steam', 'Steam Temperature (K)', 673)
        add_entry(sec, 'P_cond', 'Condenser Pressure (kPa)', 10)
        add_entry(sec, 'eta_st', 'Steam Turbine Eff η_st', 0.85)
        add_entry(sec, 'eta_fp', 'Feed Pump Eff η_fp', 0.80)

        # ── AD / Biomass ──
        sec = add_section(scroll_frame, '🌿 Biomass / AD Parameters')
        add_entry(sec, 'm_biomass', 'Biomass Feed Rate (kg/s)', 5)
        add_entry(sec, 'moisture_split', 'Moisture-Rich Fraction', 0.6)
        add_entry(sec, 'ad_yield', 'AD Biogas Yield (m³/kg)', 0.4)
        add_entry(sec, 'htc_temp', 'HTC Reactor Temp (K)', 523)

        # ── Analyse Button ──
        btn_frame = tk.Frame(scroll_frame, bg=CARD_BG)
        btn_frame.pack(fill='x', padx=10, pady=15)
        self.analyse_btn = tk.Button(
            btn_frame, text='📊  Analyse', font=('Segoe UI', 12, 'bold'),
            bg='#2563eb', fg='white', activebackground='#3b82f6',
            activeforeground='white', relief='flat', cursor='hand2',
            bd=0, padx=20, pady=10, command=self._run_analysis
        )
        self.analyse_btn.pack(fill='x')

    def _get_params(self):
        p = {}
        for key, entry in self.entries.items():
            try:
                p[key] = float(entry.get())
            except ValueError:
                messagebox.showerror('Input Error', f'Invalid value for {key}')
                return None
        return p

    def _draw_initial_schematic(self):
        draw_schematic(self.fig_schematic)
        self.canvas_schematic.draw()

    def _run_analysis(self):
        params = self._get_params()
        if params is None:
            return

        try:
            gas = calculate_gas_cycle(params)
            steam = calculate_steam_cycle(params)
            ad = calculate_ad_htc(params)
        except Exception as e:
            messagebox.showerror('Calculation Error', str(e))
            return

        # Plot H-S chart
        plot_hs_chart(self.fig_hs, steam)
        self.canvas_hs.draw()

        # Plot T-Ḣ chart
        plot_th_chart(self.fig_th, gas)
        self.canvas_th.draw()

        # Populate tables
        self._populate_tables(gas, steam, ad, params)

        # Switch to H-S tab
        self.notebook.select(1)

    def _populate_tables(self, gas, steam, ad, params):
        """Fill the State Tables tab with data."""
        for widget in self.tables_frame.winfo_children():
            widget.destroy()

        main_canvas = tk.Canvas(self.tables_frame, bg=DARK_BG, highlightthickness=0)
        main_sb = ttk.Scrollbar(self.tables_frame, orient='vertical', command=main_canvas.yview)
        inner = tk.Frame(main_canvas, bg=DARK_BG)
        inner.bind('<Configure>', lambda e: main_canvas.configure(scrollregion=main_canvas.bbox('all')))
        main_canvas.create_window((0, 0), window=inner, anchor='nw')
        main_canvas.configure(yscrollcommand=main_sb.set)
        main_sb.pack(side='right', fill='y')
        main_canvas.pack(side='left', fill='both', expand=True)

        # ── Performance Summary ──
        perf_lf = tk.LabelFrame(inner, text='⚡ Cycle Performance Summary',
                                 bg=CARD_BG, fg=ACCENT_CYAN,
                                 font=('Segoe UI', 10, 'bold'), padx=10, pady=8)
        perf_lf.pack(fill='x', padx=10, pady=8)

        metrics = [
            ('Gas Net Work', f"{gas['w_net']:.1f} kJ/kg"),
            ('Steam Net Work', f"{steam['w_net']:.1f} kJ/kg"),
            ('Gas Cycle η', f"{gas['eta']:.1f} %"),
            ('Steam Cycle η', f"{steam['eta']:.1f} %"),
            ('Back Work Ratio', f"{gas['bwr']:.1f} %"),
            ('Heat Input (CC)', f"{gas['q_in']:.1f} kJ/kg"),
            ('Gas Cycle Power', f"{gas['w_net'] * params['m_air'] / 1000:.2f} MW"),
            ('GT Exhaust Temp', f"{gas['T_exhaust']:.0f} K"),
        ]
        perf_frame = tk.Frame(perf_lf, bg=CARD_BG)
        perf_frame.pack(fill='x')
        for i, (label, val) in enumerate(metrics):
            r, c = divmod(i, 4)
            f = tk.Frame(perf_frame, bg='#1a1a2e', padx=8, pady=6)
            f.grid(row=r, column=c, padx=4, pady=3, sticky='nsew')
            perf_frame.columnconfigure(c, weight=1)
            tk.Label(f, text=val, font=('Consolas', 12, 'bold'),
                     fg=ACCENT_CYAN, bg='#1a1a2e').pack()
            tk.Label(f, text=label, font=('Segoe UI', 7),
                     fg='#94a3b8', bg='#1a1a2e').pack()

        # ── AD-HTC Balance ──
        ad_lf = tk.LabelFrame(inner, text='🌿 AD-HTC Mass & Energy Balance',
                               bg=CARD_BG, fg=ACCENT_GREEN,
                               font=('Segoe UI', 10, 'bold'), padx=10, pady=8)
        ad_lf.pack(fill='x', padx=10, pady=8)

        ad_metrics = [
            ('Total Biomass', f"{ad['m_total']:.2f} kg/s"),
            ('Moisture-Rich (→AD)', f"{ad['m_rich']:.2f} kg/s"),
            ('Moisture-Lean (→HTC)', f"{ad['m_lean']:.2f} kg/s"),
            ('Biogas Production', f"{ad['biogas_vol']:.2f} m³/s"),
            ('Biogas Mass Flow', f"{ad['m_biogas']:.2f} kg/s"),
            ('Fuel Demand', f"{gas['m_fuel']:.2f} kg/s"),
            ('HTC Thermal Energy', f"{ad['htc_energy'] / 1000:.2f} MW"),
            ('HTC Reactor Temp', f"{params['htc_temp']:.0f} K"),
        ]
        ad_frame = tk.Frame(ad_lf, bg=CARD_BG)
        ad_frame.pack(fill='x')
        for i, (label, val) in enumerate(ad_metrics):
            r, c = divmod(i, 4)
            f = tk.Frame(ad_frame, bg='#1a1a2e', padx=8, pady=6)
            f.grid(row=r, column=c, padx=4, pady=3, sticky='nsew')
            ad_frame.columnconfigure(c, weight=1)
            tk.Label(f, text=val, font=('Consolas', 11, 'bold'),
                     fg=ACCENT_GREEN, bg='#1a1a2e').pack()
            tk.Label(f, text=label, font=('Segoe UI', 7),
                     fg='#94a3b8', bg='#1a1a2e').pack()

        # ── Gas Cycle Table ──
        self._add_treeview_table(inner, '🔥 Gas Cycle State Points', gas['df'],
                                 ['State', 'Location', 'T(K)', 'P(kPa)', 'h(kJ/kg)', 's(kJ/kg·K)'])

        # ── Steam Cycle Table ──
        self._add_treeview_table(inner, '💧 HTC Steam Cycle State Points', steam['df'],
                                 ['State', 'Location', 'T(K)', 'P(kPa)', 'h(kJ/kg)', 's(kJ/kg·K)', 'x'])

    def _add_treeview_table(self, parent, title, df, columns):
        lf = tk.LabelFrame(parent, text=title, bg=CARD_BG, fg=ACCENT_CYAN,
                            font=('Segoe UI', 10, 'bold'), padx=10, pady=8)
        lf.pack(fill='x', padx=10, pady=8)

        tree = ttk.Treeview(lf, columns=columns, show='headings', height=len(df))
        for col in columns:
            tree.heading(col, text=col)
            tree.column(col, width=100, anchor='center')
        tree.column('Location', width=200, anchor='w')

        for _, row in df.iterrows():
            vals = []
            for col in columns:
                v = row.get(col, '')
                if isinstance(v, float):
                    vals.append(f'{v:.4f}' if 's(' in col else f'{v:.1f}')
                else:
                    vals.append(str(v))
            tree.insert('', 'end', values=vals)

        tree.pack(fill='x')


# ═══════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════

if __name__ == '__main__':
    app = ADHTCApp()
    app.mainloop()
