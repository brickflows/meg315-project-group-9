"""
AD-HTC Fuel-Enhanced Power Gas Cycle — Thermodynamic Property Tables
=====================================================================
Property sources:
  • Gas  — Polynomial Cp(T) correlations for dry air (ideal gas)
  • Steam — IAPWS-IF97 via `iapws` library (primary), CoolProp (fallback),
            or built-in Antoine/correlation (final fallback)

Additional modules:
  • HRSG heat-recovery model
  • Exergy / second-law analysis helpers
  • Engineering validation checks
"""

import numpy as np

# ─────────────────────────────────────────────
#  Library detection: iapws (fast) → CoolProp
# ─────────────────────────────────────────────
try:
    from iapws import IAPWS97
    HAS_IAPWS = True
except ImportError:
    HAS_IAPWS = False

HAS_COOLPROP = False
if not HAS_IAPWS:
    try:
        import CoolProp.CoolProp as CP
        HAS_COOLPROP = True
    except ImportError:
        HAS_COOLPROP = False

STEAM_SOURCE = "IAPWS-IF97 (iapws)" if HAS_IAPWS else (
    "CoolProp" if HAS_COOLPROP else "Built-in correlations")


# ═══════════════════════════════════════════════
#  GAS TABLE API  (Dry air — ideal gas)
#  Source: polynomial fit to JANAF data
# ═══════════════════════════════════════════════

class GasTable:
    """Ideal-gas property functions for air using polynomial cp(T)."""

    R = 0.287  # kJ/(kg·K) — specific gas constant for air

    @staticmethod
    def cp(T):
        """Specific heat cp [kJ/(kg·K)] at temperature T [K]."""
        t = T / 1000.0
        return (1.0483
                - 0.3717 * t
                + 0.9483 * t**2
                - 0.6271 * t**3
                + 0.1507 * t**4)

    @staticmethod
    def cv(T):
        """Specific heat cv [kJ/(kg·K)] at temperature T [K]."""
        return GasTable.cp(T) - GasTable.R

    @staticmethod
    def gamma(T):
        """Ratio of specific heats γ = cp/cv at temperature T [K]."""
        c = GasTable.cp(T)
        return c / (c - GasTable.R)

    @staticmethod
    def cp_avg(T1, T2):
        """Average cp over [T1, T2] via numerical integration."""
        temps = np.linspace(T1, T2, 100)
        cp_vals = np.vectorize(GasTable.cp)(temps)
        return float(np.mean(cp_vals))

    @staticmethod
    def gamma_avg(T1, T2):
        """Average γ over [T1, T2]."""
        temps = np.linspace(T1, T2, 100)
        g_vals = np.vectorize(GasTable.gamma)(temps)
        return float(np.mean(g_vals))

    @staticmethod
    def h(T):
        """Enthalpy h [kJ/kg] relative to 0 K via trapezoidal integration."""
        n = 200
        T = max(T, 1.0)
        temps = np.linspace(1.0, T, n + 1)
        cp_vals = np.vectorize(GasTable.cp)(temps)
        return float(np.trapezoid(cp_vals, temps))

    @staticmethod
    def s(T, P):
        """
        Entropy s [kJ/(kg·K)] relative to T_ref=298.15 K, P_ref=101.325 kPa.
        s(T,P) = integral(cp/T)dT from T_ref to T  -  R·ln(P/P_ref)
        """
        T_ref, P_ref = 298.15, 101.325
        n = 200
        temps = np.linspace(T_ref, T, n + 1)
        integrand = np.vectorize(GasTable.cp)(temps) / temps
        s_thermal = float(np.trapezoid(integrand, temps))
        return s_thermal - GasTable.R * np.log(P / P_ref)

    @staticmethod
    def T_isentropic(T1, P1, P2):
        """Find T2 such that s(T2, P2) = s(T1, P1) via bisection."""
        s_target = GasTable.s(T1, P1)
        lo, hi = T1 * 0.4, T1 * 3.5
        for _ in range(80):
            mid = (lo + hi) / 2.0
            if GasTable.s(mid, P2) < s_target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0


# ═══════════════════════════════════════════════
#  STEAM TABLE API  (Water/Steam)
#  Source: IAPWS-IF97 → CoolProp → correlations
# ═══════════════════════════════════════════════

class SteamTable:
    """
    Steam/water property functions.
    Uses IAPWS97 (primary, fast import), CoolProp (fallback),
    or built-in correlations.
    """

    # ── Saturation properties ──

    @staticmethod
    def T_sat(P_kPa):
        """Saturation temperature [K] at pressure P [kPa]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, x=0).T
        elif HAS_COOLPROP:
            return CP.PropsSI('T', 'P', P_kPa * 1000, 'Q', 0, 'Water')
        else:
            return SteamTable._T_sat_correlation(P_kPa)

    @staticmethod
    def P_sat(T_K):
        """Saturation pressure [kPa] at temperature T [K]."""
        if HAS_IAPWS:
            return IAPWS97(T=T_K, x=0).P * 1000
        elif HAS_COOLPROP:
            return CP.PropsSI('P', 'T', T_K, 'Q', 0, 'Water') / 1000
        else:
            lo, hi = 0.5, 25000
            for _ in range(60):
                mid = (lo + hi) / 2
                if SteamTable._T_sat_correlation(mid) < T_K:
                    lo = mid
                else:
                    hi = mid
            return (lo + hi) / 2

    @staticmethod
    def hf(P_kPa):
        """Saturated liquid enthalpy [kJ/kg]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, x=0).h
        elif HAS_COOLPROP:
            return CP.PropsSI('H', 'P', P_kPa * 1000, 'Q', 0, 'Water') / 1000
        else:
            Tc = SteamTable._T_sat_correlation(P_kPa) - 273.15
            return 4.18 * Tc + 0.00088 * Tc**2

    @staticmethod
    def hg(P_kPa):
        """Saturated vapour enthalpy [kJ/kg]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, x=1).h
        elif HAS_COOLPROP:
            return CP.PropsSI('H', 'P', P_kPa * 1000, 'Q', 1, 'Water') / 1000
        else:
            Tc = SteamTable._T_sat_correlation(P_kPa) - 273.15
            return 2501.3 + 1.86 * Tc

    @staticmethod
    def hfg(P_kPa):
        """Latent heat of vaporisation [kJ/kg]."""
        return SteamTable.hg(P_kPa) - SteamTable.hf(P_kPa)

    @staticmethod
    def sf(P_kPa):
        """Saturated liquid entropy [kJ/(kg·K)]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, x=0).s
        elif HAS_COOLPROP:
            return CP.PropsSI('S', 'P', P_kPa * 1000, 'Q', 0, 'Water') / 1000
        else:
            Ts = SteamTable._T_sat_correlation(P_kPa)
            return 4.18 * np.log(Ts / 273.15)

    @staticmethod
    def sg(P_kPa):
        """Saturated vapour entropy [kJ/(kg·K)]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, x=1).s
        elif HAS_COOLPROP:
            return CP.PropsSI('S', 'P', P_kPa * 1000, 'Q', 1, 'Water') / 1000
        else:
            return SteamTable.sf(P_kPa) + SteamTable.hfg(P_kPa) / SteamTable._T_sat_correlation(P_kPa)

    @staticmethod
    def vf(P_kPa):
        """Saturated liquid specific volume [m^3/kg]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, x=0).v
        elif HAS_COOLPROP:
            return 1.0 / CP.PropsSI('D', 'P', P_kPa * 1000, 'Q', 0, 'Water')
        else:
            return 0.001

    # ── Superheated properties ──

    @staticmethod
    def h_super(P_kPa, T_K):
        """Superheated steam enthalpy [kJ/kg]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, T=T_K).h
        elif HAS_COOLPROP:
            return CP.PropsSI('H', 'P', P_kPa * 1000, 'T', T_K, 'Water') / 1000
        else:
            cp_s = 1.87 + 0.00036 * (T_K - 373.15)
            return SteamTable.hg(P_kPa) + cp_s * (T_K - SteamTable._T_sat_correlation(P_kPa))

    @staticmethod
    def s_super(P_kPa, T_K):
        """Superheated steam entropy [kJ/(kg·K)]."""
        if HAS_IAPWS:
            return IAPWS97(P=P_kPa / 1000, T=T_K).s
        elif HAS_COOLPROP:
            return CP.PropsSI('S', 'P', P_kPa * 1000, 'T', T_K, 'Water') / 1000
        else:
            Ts = SteamTable._T_sat_correlation(P_kPa)
            cp_s = 1.87 + 0.00036 * (T_K - 373.15)
            return SteamTable.sg(P_kPa) + cp_s * np.log(T_K / Ts)

    @staticmethod
    def T_from_s_super(P_kPa, s_target):
        """Temperature [K] for given entropy in superheated region."""
        Ts = SteamTable.T_sat(P_kPa)
        lo, hi = Ts + 0.1, 1200.0
        for _ in range(80):
            mid = (lo + hi) / 2
            if SteamTable.s_super(P_kPa, mid) < s_target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    @staticmethod
    def x_from_s(P_kPa, s_val):
        """Quality x from entropy in wet region."""
        sfv = SteamTable.sf(P_kPa)
        sgv = SteamTable.sg(P_kPa)
        if s_val <= sfv:
            return 0.0
        if s_val >= sgv:
            return 1.0
        return (s_val - sfv) / (sgv - sfv)

    @staticmethod
    def h_from_x(P_kPa, x):
        """Enthalpy from quality."""
        return SteamTable.hf(P_kPa) + x * SteamTable.hfg(P_kPa)

    @staticmethod
    def saturation_dome(n_points=60):
        """Generate saturation dome data for H-S chart overlay."""
        pressures = np.linspace(5, 18000, n_points)
        sf_arr, hf_arr, sg_arr, hg_arr = [], [], [], []
        for P in pressures:
            try:
                sf_arr.append(SteamTable.sf(P))
                hf_arr.append(SteamTable.hf(P))
                sg_arr.append(SteamTable.sg(P))
                hg_arr.append(SteamTable.hg(P))
            except Exception:
                continue
        return (np.array(sf_arr), np.array(hf_arr),
                np.array(sg_arr), np.array(hg_arr))

    # ── Fallback correlation ──
    @staticmethod
    def _T_sat_correlation(P_kPa):
        """Antoine-type saturation temperature correlation [K]."""
        P = P_kPa / 1000.0  # MPa
        n1, n2, n3 = 1167.0521452767, -724213.16703206, -17.073846940092
        n4, n5, n6 = 12020.82470247, -3232555.0322333, 14.91510861353
        n7, n8 = -4823.2657361591, 405113.40542057
        n9, n10 = -0.23855557567849, 650.17534844798
        beta = P ** 0.25
        E = beta**2 + n3 * beta + n6
        F = n1 * beta**2 + n4 * beta + n7
        G = n2 * beta**2 + n5 * beta + n8
        D = 2 * G / (-F - np.sqrt(F**2 - 4 * E * G))
        return (n10 + D - np.sqrt((n10 + D)**2 - 4 * (n9 + n10 * D))) / 2


# ═══════════════════════════════════════════════
#  HRSG MODEL
# ═══════════════════════════════════════════════

def calculate_hrsg(T_exhaust, T_stack, m_gas, cp_avg_gas, eta_hrsg,
                   pinch_dT, q_boiler_per_kg):
    """
    Heat Recovery Steam Generator model.

    Parameters
    ----------
    T_exhaust : float  — Gas turbine exhaust temperature [K]
    T_stack   : float  — Stack / exit temperature [K]
    m_gas     : float  — Gas mass flow rate [kg/s]
    cp_avg_gas: float  — Average Cp of exhaust gas [kJ/(kg·K)]
    eta_hrsg  : float  — HRSG effectiveness (0-1)
    pinch_dT  : float  — Minimum pinch temperature difference [K]
    q_boiler_per_kg : float — Steam cycle heat input per unit steam [kJ/kg]

    Returns
    -------
    dict with Q_available, Q_recovered, m_steam, T_stack_actual, pinch_ok
    """
    Q_available = m_gas * cp_avg_gas * (T_exhaust - T_stack)           # kW
    Q_recovered = eta_hrsg * Q_available                                # kW
    m_steam = Q_recovered / q_boiler_per_kg if q_boiler_per_kg > 0 else 0  # kg/s
    T_stack_actual = T_exhaust - Q_recovered / (m_gas * cp_avg_gas) if (m_gas * cp_avg_gas) > 0 else T_stack
    pinch_ok = (T_stack_actual - SteamTable.T_sat(100)) >= pinch_dT  # simplified check

    return {
        'Q_available': Q_available,
        'Q_recovered': Q_recovered,
        'm_steam': m_steam,
        'T_stack_actual': T_stack_actual,
        'pinch_ok': pinch_ok,
    }


# ═══════════════════════════════════════════════
#  EXERGY / SECOND-LAW ANALYSIS
# ═══════════════════════════════════════════════

def exergy_flow_gas(h, s, h0, s0, T0):
    """Specific flow exergy for gas [kJ/kg].  e = (h-h0) - T0*(s-s0)"""
    return (h - h0) - T0 * (s - s0)


def exergy_destruction_component(m, s_out, s_in, Q=0, T_boundary=None, T0=298.15):
    """
    Exergy destruction for a steady-state component [kW].
    I_dot = T0 * S_gen
    S_gen = m*(s_out - s_in) - Q/T_boundary  (for heat exchange)
    """
    S_gen = m * (s_out - s_in)
    if Q != 0 and T_boundary is not None and T_boundary > 0:
        S_gen -= Q / T_boundary
    return T0 * max(S_gen, 0)  # can't be negative physically


def second_law_efficiency(W_net, E_fuel):
    """Second-law (exergetic) efficiency."""
    return (W_net / E_fuel * 100) if E_fuel > 0 else 0


def fuel_exergy(m_fuel, LHV, phi=1.04):
    """
    Chemical exergy of fuel [kW].
    phi ~ 1.04 for methane/biogas (ratio of exergy to LHV).
    """
    return m_fuel * LHV * phi


# ═══════════════════════════════════════════════
#  ENGINEERING VALIDATION
# ═══════════════════════════════════════════════

def validate_inputs(params):
    """
    Return list of warning strings for problematic input values.
    Each item is (severity, message) where severity is 'warning' or 'danger'.
    """
    warnings = []

    TIT = params.get('TIT', 0)
    if TIT > 1600:
        warnings.append(('danger',
            f'TIT = {TIT:.0f} K exceeds typical turbine blade limit (~1600 K). '
            'Risk of blade failure without advanced cooling.'))
    elif TIT > 1500:
        warnings.append(('warning',
            f'TIT = {TIT:.0f} K is near upper limit. '
            'Requires advanced blade cooling technology.'))

    T_steam = params.get('T_steam', 0)
    if T_steam > 873:
        warnings.append(('danger',
            f'Steam temperature {T_steam:.0f} K exceeds typical material limit (~873 K / 600 C).'))

    P_cond = params.get('P_cond', 0)
    if P_cond < 3:
        warnings.append(('warning',
            f'Condenser pressure {P_cond:.1f} kPa is very low. '
            'Requires deep vacuum — verify feasibility.'))

    rp = params.get('rp', 0)
    if rp < 4:
        warnings.append(('warning',
            f'Pressure ratio {rp:.1f} is unusually low for a gas turbine.'))
    elif rp > 40:
        warnings.append(('warning',
            f'Pressure ratio {rp:.1f} is very high — multi-stage compression recommended.'))

    eta_cc = params.get('eta_cc', 1.0)
    if eta_cc > 1.0:
        warnings.append(('danger', 'Combustion efficiency > 100% is non-physical.'))

    return warnings


def validate_results(gas, steam, ad, hrsg):
    """Post-calculation validation warnings."""
    warnings = []

    if gas.get('w_net', 0) <= 0:
        warnings.append(('danger',
            f"Negative net gas cycle work ({gas['w_net']:.1f} kJ/kg). "
            'Compressor work exceeds turbine work.'))

    if steam.get('w_net', 0) <= 0:
        warnings.append(('danger',
            f"Negative net steam cycle work ({steam['w_net']:.1f} kJ/kg)."))

    m_biogas = ad.get('m_biogas', 0)
    m_fuel_req = gas.get('m_fuel', 0)
    if m_fuel_req > 0 and m_biogas < m_fuel_req:
        deficit = (1 - m_biogas / m_fuel_req) * 100
        warnings.append(('warning',
            f'Biogas supply ({m_biogas:.2f} kg/s) is insufficient for fuel demand '
            f'({m_fuel_req:.2f} kg/s). Deficit: {deficit:.0f}%.'))

    if not hrsg.get('pinch_ok', True):
        warnings.append(('warning',
            'HRSG pinch point temperature difference may be violated. '
            'Increase stack temperature or reduce HRSG effectiveness.'))

    return warnings

