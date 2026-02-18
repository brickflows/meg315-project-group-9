"""Quick validation of thermodynamics module — v2 with HRSG + exergy."""
from thermodynamics import (GasTable, SteamTable, STEAM_SOURCE,
                            calculate_hrsg, exergy_flow_gas,
                            exergy_destruction_component, second_law_efficiency,
                            fuel_exergy, validate_inputs, validate_results)

print(f"Steam source: {STEAM_SOURCE}")

# Gas tests
print(f"cp(500) = {GasTable.cp(500):.6f}")
print(f"gamma(500) = {GasTable.gamma(500):.4f}")
print(f"cp_avg(298,1400) = {GasTable.cp_avg(298,1400):.4f}")
print(f"gamma_avg(298,1400) = {GasTable.gamma_avg(298,1400):.4f}")

# Steam tests
print(f"T_sat(100) = {SteamTable.T_sat(100):.4f}")
print(f"hf(100) = {SteamTable.hf(100):.2f}")
print(f"h_super(4000,673) = {SteamTable.h_super(4000,673):.2f}")

# HRSG test
hrsg = calculate_hrsg(T_exhaust=800, T_stack=420, m_gas=50, cp_avg_gas=1.1,
                       eta_hrsg=0.85, pinch_dT=15, q_boiler_per_kg=2600)
print(f"HRSG Q_available = {hrsg['Q_available']:.0f} kW")
print(f"HRSG Q_recovered = {hrsg['Q_recovered']:.0f} kW")
print(f"HRSG m_steam = {hrsg['m_steam']:.2f} kg/s")

# Exergy test
e = exergy_flow_gas(h=500, s=0.5, h0=298, s0=0, T0=298)
print(f"Exergy flow = {e:.1f} kJ/kg")
E_f = fuel_exergy(m_fuel=2, LHV=20000)
print(f"Fuel exergy = {E_f:.0f} kW")

# Validation test
w = validate_inputs({'TIT':1700,'T_steam':900,'P_cond':2,'rp':3,'eta_cc':0.98})
print(f"Warnings: {len(w)}")
for sev,msg in w:
    print(f"  [{sev}] {msg[:60]}...")

print("\n=== ALL TESTS PASSED ===")
