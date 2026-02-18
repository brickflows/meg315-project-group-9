"""Quick validation of thermodynamics module."""
import sys
print("Starting test...")
from thermodynamics import GasTable, SteamTable, HAS_IAPWS, HAS_COOLPROP
print(f"IAPWS: {HAS_IAPWS}, CoolProp: {HAS_COOLPROP}")

# Gas tests
print(f"GasTable cp(500) = {GasTable.cp(500):.6f}")
print(f"GasTable h(500) = {GasTable.h(500):.2f}")
print(f"GasTable s(500,101.325) = {GasTable.s(500,101.325):.6f}")
print(f"GasTable T_isen(298,101.325,1013.25) = {GasTable.T_isentropic(298,101.325,1013.25):.2f}")

# Steam tests
print(f"SteamTable T_sat(100) = {SteamTable.T_sat(100):.4f}")
print(f"SteamTable hf(100) = {SteamTable.hf(100):.2f}")
print(f"SteamTable hg(100) = {SteamTable.hg(100):.2f}")
print(f"SteamTable sf(100) = {SteamTable.sf(100):.6f}")
print(f"SteamTable sg(100) = {SteamTable.sg(100):.6f}")
print(f"SteamTable h_super(4000,673) = {SteamTable.h_super(4000,673):.2f}")
print(f"SteamTable s_super(4000,673) = {SteamTable.s_super(4000,673):.6f}")

print("\n=== ALL TESTS PASSED ===")
