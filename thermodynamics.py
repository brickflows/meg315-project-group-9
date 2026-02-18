"""
AD-HTC Fuel-Enhanced Power Gas Cycle — Steam & Gas Property Tables
===================================================================
Uses CoolProp for accurate steam/water properties and polynomial
correlations for ideal-gas (air) properties.
"""

import numpy as np

# ─────────────────────────────────────────────
#  Try iapws first (fast import), then CoolProp
# ─────────────────────────────────────────────
try:
    from iapws import IAPWS97
    HAS_IAPWS = True
except ImportError:
    HAS_IAPWS = False

# CoolProp is very slow to import; only use if iapws unavailable
HAS_COOLPROP = False
if not HAS_IAPWS:
    try:
        import CoolProp.CoolProp as CP
        HAS_COOLPROP = True
    except ImportError:
        HAS_COOLPROP = False


# ═══════════════════════════════════════════════
#  GAS TABLE API  (Dry air — ideal gas)
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
    def h(T):
        """Enthalpy h [kJ/kg] relative to 0 K via Simpson integration."""
        n = 200
        T = max(T, 1.0)
        temps = np.linspace(1.0, T, n + 1)
        cp_vals = np.vectorize(GasTable.cp)(temps)
        return float(np.trapezoid(cp_vals, temps))

    @staticmethod
    def s(T, P):
        """
        Entropy s [kJ/(kg·K)] relative to T_ref=298.15 K, P_ref=101.325 kPa.
        s(T,P) = ∫(cp/T)dT from T_ref to T  −  R·ln(P/P_ref)
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

    @staticmethod
    def gamma(T):
        c = GasTable.cp(T)
        return c / (c - GasTable.R)


# ═══════════════════════════════════════════════
#  STEAM TABLE API  (Water/Steam)
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
            st = IAPWS97(P=P_kPa / 1000, x=0)
            return st.T
        elif HAS_COOLPROP:
            return CP.PropsSI('T', 'P', P_kPa * 1000, 'Q', 0, 'Water')
        else:
            return SteamTable._T_sat_correlation(P_kPa)

    @staticmethod
    def P_sat(T_K):
        """Saturation pressure [kPa] at temperature T [K]."""
        if HAS_IAPWS:
            st = IAPWS97(T=T_K, x=0)
            return st.P * 1000
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
        """Saturated liquid specific volume [m³/kg]."""
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
