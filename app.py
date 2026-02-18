"""
AD-HTC Fuel-Enhanced Power Gas Cycle Analyser — Python GUI
===========================================================
Tkinter + matplotlib desktop application with:
  • HRSG-coupled gas+steam combined cycle
  • Second-law (exergy) analysis
  • AD energy balance with fuel supply/demand
  • Engineering validation warnings
  • H-s and T-h diagrams with correct labelling
"""
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

from thermodynamics import (GasTable, SteamTable, STEAM_SOURCE,
                            calculate_hrsg, exergy_flow_gas,
                            exergy_destruction_component, second_law_efficiency,
                            fuel_exergy, validate_inputs, validate_results)

# ── Colour palette ──
DARK_BG = '#0d0d1a'; CARD_BG = '#141428'; TEXT_COL = '#e2e8f0'
ACCENT_CYAN = '#38bdf8'; ACCENT_PURPLE = '#a78bfa'; ACCENT_ORANGE = '#fb923c'
ACCENT_GREEN = '#4ade80'; ACCENT_RED = '#f87171'; ACCENT_BLUE = '#60a5fa'
ACCENT_YELLOW = '#fbbf24'; ACCENT_TEAL = '#14b8a6'

# ═══════════════════════════════════════════════
#  CYCLE CALCULATIONS
# ═══════════════════════════════════════════════

def calculate_gas_cycle(p):
    T1, P1, rp, TIT = p['T1'], p['P1'], p['rp'], p['TIT']
    eta_c, eta_t, eta_cc = p['eta_c'], p['eta_t'], p['eta_cc']
    states = []
    h1, s1 = GasTable.h(T1), GasTable.s(T1, P1)
    states.append({'State':'1','Location':'Air Inlet','T(K)':T1,'P(kPa)':P1,'h(kJ/kg)':h1,'s(kJ/kgK)':s1})
    P2 = P1 * rp
    T2s = GasTable.T_isentropic(T1, P1, P2)
    h2s = GasTable.h(T2s)
    h2 = h1 + (h2s - h1) / eta_c
    T2 = T2s + (T2s - T1) * (1.0/eta_c - 1.0)
    s2 = GasTable.s(T2, P2)
    states.append({'State':'2','Location':'Compressor Exit','T(K)':T2,'P(kPa)':P2,'h(kJ/kg)':h2,'s(kJ/kgK)':s2})
    P3 = P2
    states.append({'State':'3','Location':'Combustor Inlet','T(K)':T2,'P(kPa)':P3,'h(kJ/kg)':h2,'s(kJ/kgK)':s2})
    T4 = TIT; P4 = P2 * 0.96
    h4 = GasTable.h(T4); s4 = GasTable.s(T4, P4)
    states.append({'State':'4','Location':'Combustor Exit (TIT)','T(K)':T4,'P(kPa)':P4,'h(kJ/kg)':h4,'s(kJ/kgK)':s4})
    P5 = P1 * 1.02
    T5s = GasTable.T_isentropic(T4, P4, P5)
    h5s = GasTable.h(T5s)
    h5 = h4 - eta_t * (h4 - h5s)
    T5 = T5s + (T4 - T5s) * (1.0 - eta_t)
    s5 = GasTable.s(T5, P5)
    states.append({'State':'5','Location':'Turbine Exhaust','T(K)':T5,'P(kPa)':P5,'h(kJ/kg)':h5,'s(kJ/kgK)':s5})
    w_c = h2 - h1; w_t = h4 - h5
    q_in_ideal = h4 - h2
    q_in = q_in_ideal / eta_cc
    w_net = w_t - w_c
    eta = w_net / q_in * 100 if q_in > 0 else 0
    bwr = w_c / w_t * 100 if w_t > 0 else 0
    m_fuel = (q_in * p['m_air']) / p['LHV'] if p['LHV'] > 0 else 0
    f_ratio = m_fuel / p['m_air'] if p['m_air'] > 0 else 0
    cp_avg = GasTable.cp_avg(T1, T4)
    gamma_avg = GasTable.gamma_avg(T1, T4)
    return {
        'states': states, 'df': pd.DataFrame(states),
        'w_c':w_c,'w_t':w_t,'q_in':q_in,'w_net':w_net,'eta':eta,'bwr':bwr,
        'm_fuel':m_fuel,'f_ratio':f_ratio,'T_exhaust':T5,
        'cp_avg':cp_avg,'gamma_avg':gamma_avg,
    }

def calculate_steam_cycle(p):
    Pa, Ta, Pb = p['P_boiler'], p['T_steam'], p['P_cond']
    eta_st, eta_fp = p['eta_st'], p['eta_fp']
    states = []
    ha = SteamTable.h_super(Pa, Ta); sa = SteamTable.s_super(Pa, Ta)
    states.append({'State':'a','Location':'Boiler Exit (SH)','T(K)':Ta,'P(kPa)':Pa,'h(kJ/kg)':ha,'s(kJ/kgK)':sa,'x':'-'})
    sb_s = sa; sg_c = SteamTable.sg(Pb); sf_c = SteamTable.sf(Pb)
    if sb_s < sg_c:
        xb_s = SteamTable.x_from_s(Pb, sb_s)
        hb_s = SteamTable.h_from_x(Pb, xb_s)
        Tb = SteamTable.T_sat(Pb)
    else:
        Tb = SteamTable.T_from_s_super(Pb, sb_s)
        hb_s = SteamTable.h_super(Pb, Tb)
    hb = ha - eta_st * (ha - hb_s)
    hf_c = SteamTable.hf(Pb); hfg_c = SteamTable.hfg(Pb)
    xb = (hb - hf_c) / hfg_c if hfg_c > 0 else 1.0
    if xb <= 1.0:
        sb = sf_c + xb * (sg_c - sf_c); Tb = SteamTable.T_sat(Pb); xb_s = f'{xb:.3f}'
    else:
        sb = SteamTable.s_super(Pb, Tb); xb_s = '-'
    states.append({'State':'b','Location':'ST Exit','T(K)':Tb,'P(kPa)':Pb,'h(kJ/kg)':hb,'s(kJ/kgK)':sb,'x':xb_s})
    Tc = SteamTable.T_sat(Pb); hc = SteamTable.hf(Pb); sc = SteamTable.sf(Pb)
    states.append({'State':'c','Location':'Condenser Exit','T(K)':Tc,'P(kPa)':Pb,'h(kJ/kg)':hc,'s(kJ/kgK)':sc,'x':'0.000'})
    vf_val = SteamTable.vf(Pb); wp_s = vf_val * (Pa - Pb)
    wp = wp_s / eta_fp; hd = hc + wp; sd = sc + 0.001; Td = Tc + wp / 4.18
    states.append({'State':'d','Location':'Feed Pump Exit','T(K)':Td,'P(kPa)':Pa,'h(kJ/kg)':hd,'s(kJ/kgK)':sd,'x':'-'})
    w_st = ha - hb; w_fp = wp; q_boiler = ha - hd; q_cond = hb - hc
    w_net = w_st - w_fp; eta = w_net / q_boiler * 100 if q_boiler > 0 else 0
    return {
        'states':states,'df':pd.DataFrame(states),
        'w_st':w_st,'w_fp':w_fp,'q_boiler':q_boiler,'q_cond':q_cond,
        'w_net':w_net,'eta':eta,
    }

def calculate_ad_htc(p, gas):
    m_total = p['m_biomass']
    m_rich = m_total * p['moisture_split']
    m_lean = m_total * (1 - p['moisture_split'])
    biogas_vol = m_rich * p['ad_yield']
    rho_biogas = 1.15
    m_biogas = biogas_vol * rho_biogas
    E_biogas = m_biogas * p['LHV']
    E_demand = gas['m_fuel'] * p['LHV']
    renewable_frac = min(E_biogas / E_demand, 1.0) * 100 if E_demand > 0 else 0
    htc_energy = m_lean * 1.5 * (p['htc_temp'] - 298)
    surplus = m_biogas - gas['m_fuel']
    return {
        'm_total':m_total,'m_rich':m_rich,'m_lean':m_lean,
        'biogas_vol':biogas_vol,'m_biogas':m_biogas,
        'E_biogas':E_biogas,'E_demand':E_demand,
        'renewable_frac':renewable_frac,'htc_energy':htc_energy,
        'surplus':surplus,
    }

def calculate_exergy(gas, steam, p):
    T0 = p['T1']; h0 = GasTable.h(T0); s0 = GasTable.s(T0, p['P1'])
    gs = gas['states']
    e = [exergy_flow_gas(s['h(kJ/kg)'], s['s(kJ/kgK)'], h0, s0, T0) for s in gs]
    m = p['m_air']
    I_comp = exergy_destruction_component(m, gs[1]['s(kJ/kgK)'], gs[0]['s(kJ/kgK)'], T0=T0)
    I_comb = exergy_destruction_component(m, gs[3]['s(kJ/kgK)'], gs[2]['s(kJ/kgK)'],
                                          Q=m*gas['q_in'], T_boundary=p['TIT'], T0=T0)
    I_turb = exergy_destruction_component(m, gs[4]['s(kJ/kgK)'], gs[3]['s(kJ/kgK)'], T0=T0)
    I_total = I_comp + I_comb + I_turb
    E_fuel = fuel_exergy(gas['m_fuel'], p['LHV'])
    W_net_total = gas['w_net'] * m
    eta_II = second_law_efficiency(W_net_total, E_fuel)
    S_gen_total = I_total / T0 if T0 > 0 else 0
    return {
        'e_states':e, 'I_comp':I_comp, 'I_comb':I_comb, 'I_turb':I_turb,
        'I_total':I_total, 'E_fuel':E_fuel, 'eta_II':eta_II,
        'S_gen_total':S_gen_total,
    }

# ═══════════════════════════════════════════════
#  PLOTTING
# ═══════════════════════════════════════════════

def _style_ax(ax, fig):
    ax.set_facecolor(DARK_BG); fig.patch.set_facecolor(DARK_BG)
    ax.tick_params(colors='#64748b')
    for sp in ax.spines.values(): sp.set_color('#1e293b')

def draw_schematic(fig):
    fig.clear(); ax = fig.add_subplot(111)
    ax.set_xlim(0,12); ax.set_ylim(0,9); ax.set_aspect('equal'); ax.axis('off')
    ax.set_facecolor(DARK_BG); fig.patch.set_facecolor(DARK_BG)
    bkw = dict(boxstyle='round,pad=0.3', linewidth=1.5)
    tkw = dict(ha='center',va='center',fontsize=8,fontweight='bold',color='white',family='sans-serif')
    lkw = dict(ha='center',va='center',fontsize=7,color='#94a3b8',family='sans-serif')
    # Gas cycle box
    ax.add_patch(FancyBboxPatch((0.5,0.3),11,2.8,boxstyle='round,pad=0.2',facecolor='#0f172a',edgecolor='#1e3a5f',lw=1,ls='--',alpha=0.5))
    ax.text(6,0.5,'AD-HTC FUEL-ENHANCED GAS POWER CYCLE (BRAYTON)',ha='center',fontsize=7,color=ACCENT_BLUE,alpha=0.6,fontweight='bold')
    # Compressor
    ax.add_patch(FancyBboxPatch((1.5,1.2),2,1.2,facecolor='#2563eb',edgecolor='#3b82f6',**bkw))
    ax.text(2.5,1.85,'Compressor',**tkw); ax.text(2.5,1.55,'C',fontsize=7,ha='center',va='center',color='#bfdbfe')
    ax.annotate('',xy=(1.5,1.8),xytext=(0.5,1.8),arrowprops=dict(arrowstyle='->',color='#94a3b8',lw=1.5))
    ax.text(0.6,1.45,'Air',**lkw)
    ax.plot(1.5,1.8,'o',color=ACCENT_CYAN,ms=6,zorder=5); ax.text(1.3,2.1,'1',fontsize=8,color=ACCENT_CYAN,fontweight='bold')
    ax.plot(3.5,1.8,'o',color=ACCENT_CYAN,ms=6,zorder=5); ax.text(3.65,2.1,'2',fontsize=8,color=ACCENT_CYAN,fontweight='bold')
    ax.plot([3.5,8.5],[1.8,1.8],'--',color='#64748b',lw=2)
    ax.text(6,1.55,'<- Shaft ->',fontsize=6,ha='center',color='#64748b')
    # Combustion chamber
    ax.add_patch(FancyBboxPatch((4.5,3.8),3,1.2,facecolor='#d97706',edgecolor='#f59e0b',**bkw))
    ax.text(6,4.45,'Biogas Combustion',**tkw); ax.text(6,4.15,'Chamber',fontsize=7,ha='center',va='center',color='#fef3c7')
    ax.annotate('',xy=(4.5,4.4),xytext=(3.5,1.8),arrowprops=dict(arrowstyle='->',color=ACCENT_BLUE,lw=2,connectionstyle='arc3,rad=-0.2'))
    ax.plot(4.5,4.4,'o',color=ACCENT_YELLOW,ms=6,zorder=5); ax.text(4.2,4.6,'3',fontsize=8,color=ACCENT_YELLOW,fontweight='bold')
    ax.annotate('',xy=(8.5,1.8),xytext=(7.5,4.4),arrowprops=dict(arrowstyle='->',color=ACCENT_ORANGE,lw=2,connectionstyle='arc3,rad=0.2'))
    ax.plot(7.5,4.4,'o',color=ACCENT_YELLOW,ms=6,zorder=5); ax.text(7.65,4.65,'4',fontsize=8,color=ACCENT_YELLOW,fontweight='bold')
    # Turbine
    ax.add_patch(FancyBboxPatch((8.5,1.2),2,1.2,facecolor='#dc2626',edgecolor='#ef4444',**bkw))
    ax.text(9.5,1.85,'Turbine',**tkw); ax.text(9.5,1.55,'GT',fontsize=7,ha='center',va='center',color='#fecaca')
    ax.annotate('',xy=(11.5,1.8),xytext=(10.5,1.8),arrowprops=dict(arrowstyle='->',color='#94a3b8',lw=1.5))
    ax.text(11.2,1.45,'Exhaust',**lkw)
    ax.plot(10.5,1.8,'o',color=ACCENT_RED,ms=6,zorder=5); ax.text(10.65,2.1,'5',fontsize=8,color=ACCENT_RED,fontweight='bold')
    # HRSG box
    ax.add_patch(FancyBboxPatch((9.2,3.2),2.2,0.8,facecolor='#7c2d12',edgecolor='#ea580c',**bkw))
    ax.text(10.3,3.65,'HRSG',**tkw)
    ax.annotate('',xy=(10.3,3.2),xytext=(10.3,2.4),arrowprops=dict(arrowstyle='->',color=ACCENT_ORANGE,lw=1.5))
    ax.annotate('',xy=(6.2,5.6),xytext=(9.2,3.6),arrowprops=dict(arrowstyle='->',color=ACCENT_RED,lw=1.5,connectionstyle='arc3,rad=0.3'))
    ax.text(8.0,4.5,'Steam\nHeat',fontsize=6,color=ACCENT_RED)
    # HTC Steam Cycle
    ax.add_patch(FancyBboxPatch((1.2,5.0),5,2.5,boxstyle='round,pad=0.2',facecolor='#0f0f2a',edgecolor='#312e81',lw=1,ls='--',alpha=0.4))
    ax.text(3.7,7.2,'HTC STEAM CYCLE (RANKINE)',fontsize=7,ha='center',color=ACCENT_PURPLE,alpha=0.6,fontweight='bold')
    # Boiler
    ax.add_patch(FancyBboxPatch((3.5,6.3),2,0.8,facecolor='#be123c',edgecolor='#e11d48',**bkw))
    ax.text(4.5,6.7,'Boiler',**tkw)
    # HTC Reactor
    ax.add_patch(FancyBboxPatch((3.5,5.2),2,0.8,facecolor='#ea580c',edgecolor='#f97316',**bkw))
    ax.text(4.5,5.6,'HTC Reactor',**tkw)
    # ST
    st_pts = np.array([[2.3,5.8],[2.8,6.3],[2.3,6.8]])
    ax.add_patch(plt.Polygon(st_pts,facecolor='#7c3aed',edgecolor='#a78bfa',lw=1.5))
    ax.text(2.4,6.3,'ST',fontsize=7,ha='center',va='center',color='white',fontweight='bold')
    ax.plot(3.5,6.7,'o',color=ACCENT_PURPLE,ms=5,zorder=5); ax.text(3.3,6.9,'b',fontsize=7,color=ACCENT_PURPLE,fontweight='bold')
    ax.plot(5.5,6.7,'o',color=ACCENT_PURPLE,ms=5,zorder=5); ax.text(5.6,6.9,'a',fontsize=7,color=ACCENT_PURPLE,fontweight='bold')
    ax.annotate('',xy=(2.8,6.5),xytext=(3.5,6.7),arrowprops=dict(arrowstyle='->',color=ACCENT_PURPLE,lw=1.5))
    ax.annotate('',xy=(4.5,6.3),xytext=(4.5,6.0),arrowprops=dict(arrowstyle='<->',color='#94a3b8',lw=1.2))
    ax.plot([2.3,1.8,1.8,2.3],[5.8,5.8,6.8,6.8],'--',color=ACCENT_PURPLE,lw=1,alpha=0.4)
    ax.text(1.6,6.3,'Cond.',fontsize=6,ha='center',color=ACCENT_PURPLE,alpha=0.5)
    # Biomass Homogenizer
    ax.add_patch(FancyBboxPatch((1.5,7.8),2.5,0.8,facecolor='#7c3aed',edgecolor='#a78bfa',**bkw))
    ax.text(2.75,8.25,'Biomass Feedstock',fontsize=7,ha='center',va='center',color='white',fontweight='bold')
    ax.text(2.75,8.0,'Homogenizer',fontsize=7,ha='center',va='center',color='#ddd6fe')
    ax.annotate('',xy=(1.5,8.2),xytext=(0.5,8.2),arrowprops=dict(arrowstyle='->',color=ACCENT_GREEN,lw=1.8))
    ax.text(0.4,8.5,'Biomass',fontsize=6,color='#94a3b8')
    ax.annotate('',xy=(3.5,5.6),xytext=(2.75,7.8),arrowprops=dict(arrowstyle='->',color=ACCENT_GREEN,lw=1.5,connectionstyle='arc3,rad=-0.1'))
    ax.text(2.3,7.0,'Moisture\n-lean',fontsize=6,color=ACCENT_GREEN)
    # AD
    ax.add_patch(FancyBboxPatch((7.5,7.3),1.2,0.8,facecolor='#0891b2',edgecolor='#06b6d4',**bkw))
    ax.text(8.1,7.7,'AD',**tkw)
    ax.annotate('',xy=(7.5,8.0),xytext=(4.0,8.2),arrowprops=dict(arrowstyle='->',color='#22d3ee',lw=1.5))
    ax.text(5.8,8.5,'Moisture-rich',fontsize=6,ha='center',color='#22d3ee')
    # Biogas Collector
    ax.add_patch(FancyBboxPatch((7.0,5.8),2.2,0.8,facecolor='#0d9488',edgecolor='#14b8a6',**bkw))
    ax.text(8.1,6.25,'Enhanced Biogas',fontsize=7,ha='center',va='center',color='white',fontweight='bold')
    ax.text(8.1,6.0,'Collector',fontsize=7,ha='center',va='center',color='#99f6e4')
    ax.annotate('',xy=(8.1,6.6),xytext=(8.1,7.3),arrowprops=dict(arrowstyle='->',color=ACCENT_GREEN,lw=1.5))
    ax.annotate('',xy=(6,5.0),xytext=(7.5,6.0),arrowprops=dict(arrowstyle='->',color=ACCENT_GREEN,lw=2,connectionstyle='arc3,rad=0.2'))
    ax.text(7.0,5.3,'Biogas\nFuel',fontsize=6,color=ACCENT_GREEN,fontweight='bold')
    ax.annotate('',xy=(10.5,6.2),xytext=(9.2,6.2),arrowprops=dict(arrowstyle='->',color=ACCENT_GREEN,lw=1.5))
    ax.text(10.6,6.5,'Biogas Dist.',fontsize=6,color=ACCENT_GREEN)
    ax.text(6,0.1,'Energhx Research Group - University of Lagos',ha='center',fontsize=6,color='#475569')
    fig.tight_layout(pad=0.5)

def plot_hs_chart(fig, steam):
    fig.clear(); ax = fig.add_subplot(111); _style_ax(ax, fig)
    sm = {s['State']:s for s in steam['states']}
    sf_a, hf_a, sg_a, hg_a = SteamTable.saturation_dome(60)
    dome_s = np.concatenate([sf_a, sg_a[::-1], [sf_a[0]]])
    dome_h = np.concatenate([hf_a, hg_a[::-1], [hf_a[0]]])
    ax.fill(dome_s, dome_h, alpha=0.05, color='#94a3b8')
    ax.plot(dome_s, dome_h, '--', color='#94a3b8', alpha=0.35, lw=1.2, label='Saturation Dome')
    procs = [('Boiler (d->a)','d','a',ACCENT_ORANGE),('ST (a->b)','a','b',ACCENT_PURPLE),
             ('Condenser (b->c)','b','c',ACCENT_CYAN),('Pump (c->d)','c','d',ACCENT_BLUE)]
    for name,f,t,col in procs:
        s1,s2 = sm[f],sm[t]
        ax.plot(np.linspace(s1['s(kJ/kgK)'],s2['s(kJ/kgK)'],25),
                np.linspace(s1['h(kJ/kg)'],s2['h(kJ/kg)'],25),'-',color=col,lw=2.5,label=name)
    for sid in ['d','a','b','c']:
        st = sm[sid]
        ax.plot(st['s(kJ/kgK)'],st['h(kJ/kg)'],'o',color=ACCENT_CYAN,ms=10,mec='white',mew=1.5,zorder=5)
        ax.annotate(f" {sid} {st['Location']}",xy=(st['s(kJ/kgK)'],st['h(kJ/kg)']),fontsize=7,color=TEXT_COL,fontweight='bold')
    ax.set_xlabel('Entropy s [kJ/(kg K)]',color=TEXT_COL,fontsize=11)
    ax.set_ylabel('Enthalpy h [kJ/kg]',color=TEXT_COL,fontsize=11)
    ax.set_title('h-s Diagram - HTC Steam (Rankine) Cycle',color=TEXT_COL,fontsize=13,fontweight='bold',pad=12)
    ax.legend(fontsize=8,facecolor=CARD_BG,edgecolor='#334155',labelcolor=TEXT_COL,loc='best')
    ax.grid(True,alpha=0.1,color=ACCENT_CYAN); fig.tight_layout(pad=1.0)

def plot_th_chart(fig, gas):
    fig.clear(); ax = fig.add_subplot(111); _style_ax(ax, fig)
    states = gas['states']
    procs = [('Compression (1->2)',0,1,ACCENT_BLUE),('Combustion (3->4)',2,3,ACCENT_ORANGE),('Expansion (4->5)',3,4,ACCENT_RED)]
    cum_h = 0
    mx, my, ml = [0],[states[0]['T(K)']],[states[0]['State']]
    for name,i,j,col in procs:
        s1,s2 = states[i],states[j]
        dh = abs(s2['h(kJ/kg)']-s1['h(kJ/kg)'])
        h_v = np.linspace(cum_h,cum_h+dh,30); T_v = np.linspace(s1['T(K)'],s2['T(K)'],30)
        ax.plot(h_v,T_v,'-',color=col,lw=2.5,label=name)
        cum_h += dh; mx.append(cum_h); my.append(s2['T(K)']); ml.append(s2['State'])
    ax.plot(mx,my,'D',color=ACCENT_ORANGE,ms=8,mec='white',mew=1.5,zorder=5,label='State Points')
    for x,y,l in zip(mx,my,ml):
        ax.annotate(f' {l}',(x,y),fontsize=8,color=TEXT_COL,fontweight='bold')
    ax.set_xlabel('Specific Enthalpy Change [kJ/kg]',color=TEXT_COL,fontsize=11)
    ax.set_ylabel('Temperature T [K]',color=TEXT_COL,fontsize=11)
    ax.set_title('T-h Diagram - Gas Power Cycle (Brayton)',color=TEXT_COL,fontsize=13,fontweight='bold',pad=12)
    ax.legend(fontsize=8,facecolor=CARD_BG,edgecolor='#334155',labelcolor=TEXT_COL,loc='best')
    ax.grid(True,alpha=0.1,color=ACCENT_ORANGE); fig.tight_layout(pad=1.0)

# ═══════════════════════════════════════════════
#  TKINTER GUI
# ═══════════════════════════════════════════════

class ADHTCApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('AD-HTC Fuel-Enhanced Power Gas Cycle Analyser')
        self.geometry('1400x900'); self.configure(bg=DARK_BG); self.minsize(1100,700)
        st = ttk.Style(); st.theme_use('clam')
        st.configure('.',background=DARK_BG,foreground=TEXT_COL)
        st.configure('TFrame',background=DARK_BG); st.configure('TLabel',background=DARK_BG,foreground=TEXT_COL,font=('Segoe UI',9))
        st.configure('TLabelframe',background=CARD_BG,foreground=ACCENT_CYAN,font=('Segoe UI',10,'bold'))
        st.configure('TLabelframe.Label',background=CARD_BG,foreground=ACCENT_CYAN)
        st.configure('TEntry',fieldbackground='#1a1a2e',foreground=TEXT_COL,font=('Consolas',9))
        st.configure('Treeview',background='#1a1a2e',foreground=TEXT_COL,fieldbackground='#1a1a2e',font=('Consolas',8),rowheight=22)
        st.configure('Treeview.Heading',background=CARD_BG,foreground=ACCENT_CYAN,font=('Segoe UI',8,'bold'))
        self._build_ui(); self._draw_initial_schematic()

    def _build_ui(self):
        hdr = tk.Frame(self,bg='#080816',height=60); hdr.pack(fill='x',side='top'); hdr.pack_propagate(False)
        tk.Label(hdr,text='AD-HTC Fuel-Enhanced Power Gas Cycle',font=('Segoe UI',15,'bold'),fg=ACCENT_CYAN,bg='#080816').pack(pady=(8,0))
        tk.Label(hdr,text=f'Energhx Research Group | Gas: Polynomial Cp(T) | Steam: {STEAM_SOURCE}',
                 font=('Segoe UI',8),fg='#64748b',bg='#080816').pack()
        self.pw = tk.PanedWindow(self,orient='horizontal',bg=DARK_BG,sashwidth=4,sashrelief='flat')
        self.pw.pack(fill='both',expand=True,padx=8,pady=8)
        self.left = tk.Frame(self.pw,bg=DARK_BG); self.pw.add(self.left,stretch='always',width=900)
        self.nb = ttk.Notebook(self.left); self.nb.pack(fill='both',expand=True)
        # Tabs
        self.sch_f = tk.Frame(self.nb,bg=DARK_BG); self.nb.add(self.sch_f,text='  Schematic  ')
        self.fig_sch = plt.Figure(figsize=(9,6),facecolor=DARK_BG)
        self.cv_sch = FigureCanvasTkAgg(self.fig_sch,self.sch_f); self.cv_sch.get_tk_widget().pack(fill='both',expand=True)
        self.hs_f = tk.Frame(self.nb,bg=DARK_BG); self.nb.add(self.hs_f,text='  h-s Chart (Steam)  ')
        self.fig_hs = plt.Figure(figsize=(8,5),facecolor=DARK_BG)
        self.cv_hs = FigureCanvasTkAgg(self.fig_hs,self.hs_f); self.cv_hs.get_tk_widget().pack(fill='both',expand=True)
        NavigationToolbar2Tk(self.cv_hs,self.hs_f).update()
        self.th_f = tk.Frame(self.nb,bg=DARK_BG); self.nb.add(self.th_f,text='  T-h Chart (Gas)  ')
        self.fig_th = plt.Figure(figsize=(8,5),facecolor=DARK_BG)
        self.cv_th = FigureCanvasTkAgg(self.fig_th,self.th_f); self.cv_th.get_tk_widget().pack(fill='both',expand=True)
        NavigationToolbar2Tk(self.cv_th,self.th_f).update()
        self.tbl_f = tk.Frame(self.nb,bg=DARK_BG); self.nb.add(self.tbl_f,text='  Results & Tables  ')
        self.right = tk.Frame(self.pw,bg=CARD_BG); self.pw.add(self.right,stretch='never',width=330)
        self._build_params()

    def _build_params(self):
        cv = tk.Canvas(self.right,bg=CARD_BG,highlightthickness=0)
        sb = ttk.Scrollbar(self.right,orient='vertical',command=cv.yview)
        sf = tk.Frame(cv,bg=CARD_BG)
        sf.bind('<Configure>',lambda e:cv.configure(scrollregion=cv.bbox('all')))
        cv.create_window((0,0),window=sf,anchor='nw'); cv.configure(yscrollcommand=sb.set)
        sb.pack(side='right',fill='y'); cv.pack(side='left',fill='both',expand=True)
        cv.bind_all('<MouseWheel>',lambda e:cv.yview_scroll(-1*(e.delta//120),'units'))
        self.entries = {}
        def sec(parent,title):
            lf = tk.LabelFrame(parent,text=title,bg=CARD_BG,fg=ACCENT_CYAN,font=('Segoe UI',9,'bold'),padx=10,pady=5)
            lf.pack(fill='x',padx=10,pady=(6,3)); return lf
        def ent(parent,key,label,default):
            f = tk.Frame(parent,bg=CARD_BG); f.pack(fill='x',pady=1)
            tk.Label(f,text=label,bg=CARD_BG,fg='#94a3b8',font=('Segoe UI',8),anchor='w').pack(fill='x')
            e = tk.Entry(f,bg='#1a1a2e',fg=TEXT_COL,insertbackground=TEXT_COL,font=('Consolas',9),
                         relief='flat',bd=0,highlightthickness=1,highlightcolor=ACCENT_CYAN,highlightbackground='#334155')
            e.insert(0,str(default)); e.pack(fill='x',ipady=3); self.entries[key] = e
        s = sec(sf,'Gas Power Cycle (Brayton)')
        ent(s,'T1','Ambient Temp T1 (K)',298); ent(s,'P1','Ambient Press P1 (kPa)',101.325)
        ent(s,'rp','Pressure Ratio rp',10); ent(s,'TIT','Turbine Inlet Temp TIT (K)',1400)
        ent(s,'eta_c','Compressor Eff eta_c',0.86); ent(s,'eta_t','Turbine Eff eta_t',0.89)
        ent(s,'eta_cc','Combustion Eff eta_cc',0.98)
        ent(s,'LHV','Biogas LHV (kJ/kg)',20000); ent(s,'m_air','Air Mass Flow (kg/s)',50)
        s = sec(sf,'HTC Steam Cycle (Rankine)')
        ent(s,'P_boiler','Boiler Pressure (kPa)',4000); ent(s,'T_steam','Steam Temp (K)',673)
        ent(s,'P_cond','Condenser Press (kPa)',10)
        ent(s,'eta_st','ST Eff eta_st',0.85); ent(s,'eta_fp','Pump Eff eta_fp',0.80)
        s = sec(sf,'HRSG Coupling')
        ent(s,'T_stack','Stack Temp (K)',420); ent(s,'eta_hrsg','HRSG Effectiveness',0.85)
        ent(s,'pinch_dT','Pinch Delta-T (K)',15)
        s = sec(sf,'Biomass / AD Parameters')
        ent(s,'m_biomass','Biomass Feed (kg/s)',5); ent(s,'moisture_split','Moisture-Rich Frac',0.6)
        ent(s,'ad_yield','AD Yield (m3/kg)',0.4); ent(s,'htc_temp','HTC Reactor Temp (K)',523)
        bf = tk.Frame(sf,bg=CARD_BG); bf.pack(fill='x',padx=10,pady=10)
        tk.Button(bf,text='Analyse',font=('Segoe UI',12,'bold'),bg='#2563eb',fg='white',
                  activebackground='#3b82f6',activeforeground='white',relief='flat',cursor='hand2',
                  bd=0,padx=20,pady=10,command=self._run).pack(fill='x')

    def _get_p(self):
        p = {}
        for k,e in self.entries.items():
            try: p[k] = float(e.get())
            except ValueError: messagebox.showerror('Input Error',f'Invalid: {k}'); return None
        return p

    def _draw_initial_schematic(self):
        draw_schematic(self.fig_sch); self.cv_sch.draw()

    def _run(self):
        p = self._get_p()
        if not p: return
        # Input validation
        warns_in = validate_inputs(p)
        try:
            gas = calculate_gas_cycle(p)
            steam = calculate_steam_cycle(p)
            ad = calculate_ad_htc(p, gas)
            # HRSG coupling
            cp_exh = GasTable.cp_avg(gas['T_exhaust'], p['T1'])
            hrsg = calculate_hrsg(gas['T_exhaust'], p['T_stack'], p['m_air'], cp_exh,
                                  p['eta_hrsg'], p['pinch_dT'], steam['q_boiler'])
            exg = calculate_exergy(gas, steam, p)
        except Exception as e:
            messagebox.showerror('Calculation Error',str(e)); return
        warns_out = validate_results(gas, steam, ad, hrsg)
        all_warns = warns_in + warns_out
        plot_hs_chart(self.fig_hs, steam); self.cv_hs.draw()
        plot_th_chart(self.fig_th, gas); self.cv_th.draw()
        self._fill_tables(gas, steam, ad, hrsg, exg, p, all_warns)
        self.nb.select(1)

    def _fill_tables(self, gas, steam, ad, hrsg, exg, p, warns):
        for w in self.tbl_f.winfo_children(): w.destroy()
        cv = tk.Canvas(self.tbl_f,bg=DARK_BG,highlightthickness=0)
        sb = ttk.Scrollbar(self.tbl_f,orient='vertical',command=cv.yview)
        inner = tk.Frame(cv,bg=DARK_BG)
        inner.bind('<Configure>',lambda e:cv.configure(scrollregion=cv.bbox('all')))
        cv.create_window((0,0),window=inner,anchor='nw'); cv.configure(yscrollcommand=sb.set)
        sb.pack(side='right',fill='y'); cv.pack(side='left',fill='both',expand=True)

        def card(parent,title,fg=ACCENT_CYAN):
            lf = tk.LabelFrame(parent,text=title,bg=CARD_BG,fg=fg,font=('Segoe UI',10,'bold'),padx=10,pady=8)
            lf.pack(fill='x',padx=10,pady=6); return lf

        def metric_grid(parent, items, fg=ACCENT_CYAN, cols=4):
            fr = tk.Frame(parent,bg=CARD_BG); fr.pack(fill='x')
            for i,(lbl,val) in enumerate(items):
                r,c = divmod(i,cols)
                f = tk.Frame(fr,bg='#1a1a2e',padx=8,pady=5); f.grid(row=r,column=c,padx=3,pady=2,sticky='nsew')
                fr.columnconfigure(c,weight=1)
                tk.Label(f,text=val,font=('Consolas',11,'bold'),fg=fg,bg='#1a1a2e').pack()
                tk.Label(f,text=lbl,font=('Segoe UI',8),fg='#94a3b8',bg='#1a1a2e').pack()

        # ── Warnings ──
        if warns:
            wf = card(inner,'Engineering Warnings',fg=ACCENT_RED)
            for sev,msg in warns:
                col = ACCENT_RED if sev=='danger' else ACCENT_YELLOW
                icon = 'X' if sev=='danger' else '!'
                tk.Label(wf,text=f'[{icon}] {msg}',fg=col,bg=CARD_BG,font=('Segoe UI',8),
                         anchor='w',wraplength=600,justify='left').pack(fill='x',pady=1)

        # ── Combined Performance ──
        W_GT = gas['w_t'] * p['m_air'] / 1000  # MW
        W_comp = gas['w_c'] * p['m_air'] / 1000
        W_net_gas = gas['w_net'] * p['m_air'] / 1000
        m_st = hrsg['m_steam']
        W_ST = steam['w_st'] * m_st / 1000
        W_pump = steam['w_fp'] * m_st / 1000
        W_net_steam = steam['w_net'] * m_st / 1000
        W_net_comb = W_net_gas + W_net_steam
        Q_in_total = gas['q_in'] * p['m_air'] / 1000
        eta_comb = W_net_comb / Q_in_total * 100 if Q_in_total > 0 else 0
        cf = card(inner,'Combined Cycle Performance')
        metric_grid(cf, [
            ('W_GT [MW]', f'{W_GT:.2f}'), ('W_Comp [MW]', f'{W_comp:.2f}'),
            ('W_net,gas [MW]', f'{W_net_gas:.2f}'), ('Gas eta [%]', f'{gas["eta"]:.1f}'),
            ('W_ST [MW]', f'{W_ST:.2f}'), ('W_Pump [MW]', f'{W_pump:.3f}'),
            ('W_net,steam [MW]', f'{W_net_steam:.2f}'), ('Steam eta [%]', f'{steam["eta"]:.1f}'),
            ('W_net,combined [MW]', f'{W_net_comb:.2f}'), ('eta_combined [%]', f'{eta_comb:.1f}'),
            ('BWR [%]', f'{gas["bwr"]:.1f}'), ('Q_in [MW]', f'{Q_in_total:.2f}'),
        ])

        # ── Gas Properties ──
        gf = card(inner,'Gas Cycle Properties',fg=ACCENT_BLUE)
        metric_grid(gf, [
            ('Cp_avg [kJ/kgK]', f'{gas["cp_avg"]:.4f}'), ('gamma_avg', f'{gas["gamma_avg"]:.4f}'),
            ('Fuel-Air Ratio f', f'{gas["f_ratio"]:.4f}'), ('m_fuel [kg/s]', f'{gas["m_fuel"]:.2f}'),
            ('T_exhaust [K]', f'{gas["T_exhaust"]:.0f}'), ('Combustion model', 'T3-constrained'),
        ], fg=ACCENT_BLUE)

        # ── HRSG ──
        hf = card(inner,'HRSG Heat Recovery',fg=ACCENT_ORANGE)
        metric_grid(hf, [
            ('Q_available [kW]', f'{hrsg["Q_available"]:.0f}'), ('Q_recovered [kW]', f'{hrsg["Q_recovered"]:.0f}'),
            ('m_steam [kg/s]', f'{hrsg["m_steam"]:.2f}'), ('T_stack [K]', f'{hrsg["T_stack_actual"]:.0f}'),
            ('Pinch OK', 'Yes' if hrsg['pinch_ok'] else 'NO'),
        ], fg=ACCENT_ORANGE)

        # ── AD Balance ──
        af = card(inner,'AD-HTC Mass & Energy Balance',fg=ACCENT_GREEN)
        status = 'SURPLUS' if ad['surplus'] >= 0 else 'DEFICIT'
        metric_grid(af, [
            ('Biomass [kg/s]', f'{ad["m_total"]:.2f}'), ('Moist-Rich->AD [kg/s]', f'{ad["m_rich"]:.2f}'),
            ('Moist-Lean->HTC [kg/s]', f'{ad["m_lean"]:.2f}'), ('Biogas [m3/s]', f'{ad["biogas_vol"]:.2f}'),
            ('m_biogas [kg/s]', f'{ad["m_biogas"]:.2f}'), ('E_biogas [kW]', f'{ad["E_biogas"]:.0f}'),
            ('Fuel Demand [kg/s]', f'{gas["m_fuel"]:.2f}'), ('Renewable [%]', f'{ad["renewable_frac"]:.0f}'),
            ('Supply vs Demand', status), ('HTC Energy [kW]', f'{ad["htc_energy"]:.0f}'),
        ], fg=ACCENT_GREEN)

        # ── Second-Law ──
        ef = card(inner,'Second-Law (Exergy) Analysis',fg=ACCENT_PURPLE)
        metric_grid(ef, [
            ('I_comp [kW]', f'{exg["I_comp"]:.0f}'), ('I_comb [kW]', f'{exg["I_comb"]:.0f}'),
            ('I_turb [kW]', f'{exg["I_turb"]:.0f}'), ('I_total [kW]', f'{exg["I_total"]:.0f}'),
            ('S_gen [kW/K]', f'{exg["S_gen_total"]:.2f}'), ('E_fuel [kW]', f'{exg["E_fuel"]:.0f}'),
            ('eta_II [%]', f'{exg["eta_II"]:.1f}'),
        ], fg=ACCENT_PURPLE)

        # ── State Tables ──
        self._tree(inner,'Gas Cycle State Points',gas['df'],
                   ['State','Location','T(K)','P(kPa)','h(kJ/kg)','s(kJ/kgK)'])
        self._tree(inner,'Steam Cycle State Points',steam['df'],
                   ['State','Location','T(K)','P(kPa)','h(kJ/kg)','s(kJ/kgK)','x'])

    def _tree(self, parent, title, df, cols):
        lf = tk.LabelFrame(parent,text=title,bg=CARD_BG,fg=ACCENT_CYAN,font=('Segoe UI',10,'bold'),padx=10,pady=8)
        lf.pack(fill='x',padx=10,pady=6)
        tr = ttk.Treeview(lf,columns=cols,show='headings',height=len(df))
        for c in cols:
            tr.heading(c,text=c); tr.column(c,width=100,anchor='center')
        tr.column('Location',width=200,anchor='w')
        for _,row in df.iterrows():
            vals = []
            for c in cols:
                v = row.get(c,'')
                vals.append(f'{v:.4f}' if isinstance(v,float) and 's(' in c else (f'{v:.1f}' if isinstance(v,float) else str(v)))
            tr.insert('','end',values=vals)
        tr.pack(fill='x')

if __name__ == '__main__':
    app = ADHTCApp()
    app.mainloop()
