# AD-HTC Fuel-Enhanced Power Gas Cycle Analyser
### MEG 315 — Assignment 7 | Energhx Research Group, University of Lagos

> An interactive thermodynamic simulation of a combined **Anaerobic Digestion (AD)** + **Hydrothermal Carbonization (HTC)** fuel-enhanced gas power cycle — with a coupled Rankine steam cycle and real-gas property models.

---

## 🔗 Live Demo & Repository

| | Link |
|---|---|
| 🌐 **Live Web App** | [https://meg315-assignment7.vercel.app](https://meg315-assignment7.vercel.app) |
| 📦 **GitHub Repository** | [https://github.com/brickflows/meg315-assignment7](https://github.com/brickflows/meg315-assignment7) |

---

## 📸 Schematic Preview

The interactive animated system schematic — showing all energy and mass flows across the AD-HTC combined cycle:

![AD-HTC Cycle Schematic](https://raw.githubusercontent.com/brickflows/meg315-assignment7/main/preview.png)

> **Live version**: Open [`index.html`](index.html) in your browser to see the fully animated schematic with flowing particle animations, spinning generator symbol, pulsing state-point indicators, and an animated background grid.

---

## 📋 Table of Contents

- [Overview](#overview)
- [System Architecture](#system-architecture)
- [Features](#features)
- [File Structure](#file-structure)
- [Thermodynamic Models](#thermodynamic-models)
- [Installation & Setup](#installation--setup)
- [Usage](#usage)
- [Validation](#validation)

---

## Overview

This project models an integrated **AD-HTC Fuel-Enhanced Gas Power Cycle** — a novel combined system where:

- **Moisture-rich biomass** is fed to an **Anaerobic Digester (AD)** to produce biogas
- **Moisture-lean biomass** is fed to a **Hydrothermal Carbonization (HTC) Reactor** to produce hydrochar and process heat
- The biogas fuels a **Brayton (Gas) Cycle** turbine
- The HTC process heat drives a **Rankine (Steam) Cycle** via a boiler
- Both cycles generate shaft work / electricity

The result is a renewable, waste-to-energy system suitable for distributed power generation and building energy supply.

---

## System Architecture

```
Biomass Feedstock
       │
       ▼
┌─────────────────────┐
│ Biomass Feedstock   │
│   Homogenizer       │
└──────────┬──────────┘
           │
     ┌─────┴──────┐
     │             │
     ▼             ▼
Moisture-lean  Moisture-rich
Biomass         Biomass
     │               │
     ▼               ▼
┌────────┐     ┌─────────┐
│  HTC   │     │   AD    │
│ Reactor│◄────│(Biogas) │
└────┬───┘     └────┬────┘
     │  ▲           │
     │  │           ▼
     ▼  │     ┌──────────────────┐
┌─────────┐   │ Enhanced Biogas  │
│  Boiler │──►│   Collector      │
└─────────┘   └──────┬───────────┘
     │               │         │
     ▼               ▼         ▼
  Rankine     Combustion   Biogas to
  Cycle        Chamber     Buildings
                   │
        ┌──────────┴──────────┐
        ▼                     ▼
   Compressor ◄── Shaft ──► Turbine
        ▲                     │
        │ Air            Exhaust + Generator
```

---

## Features

### 🌐 Web Application (`index.html`)
- **Animated schematic** — glowing particle dots travel along every pipe path in real-time
- **Engineering-accurate turbomachinery** — trapezoid Compressor (blue) and Turbine (red)
- **Spinning GEN symbol** at turbine exit (generator output)
- **Pulsing state-point indicators** ①②③④⑤ at each thermodynamic state
- **Animated background grid** for a premium visual feel
- **Real-time Brayton cycle** calculations from input parameters
- **Real-time Rankine cycle** calculations with saturation dome on H-s chart
- **State-point tables** — T, P, h, s, x at every state
- **AD-HTC mass & energy balance** — biogas production vs combustor fuel demand

### 🖥️ Desktop Application (`app.py`)
- **Tkinter GUI** — dual-panel layout with scrollable parameter inputs and tabbed charts
- **HRSG coupling** — gas turbine exhaust heats the steam cycle with pinch-point enforcement
- **Exergy / Second-Law analysis** — irreversibility per component, `η_II`, `Ṡ_gen`
- **Engineering validation** — warns on low back-work ratio, unrealistic efficiencies, etc.
- **Collapsible advanced sections** — HRSG, Biomass/AD, optional parameters

---

## File Structure

```
meg315-assignment7/
│
├── index.html              # Web application (open directly in browser)
├── styles.css              # Dark-theme CSS, glassmorphism, animations
├── script.js               # Thermodynamic engine + Chart.js visualisations
│
├── app.py                  # Python/Tkinter desktop GUI application
├── thermodynamics.py       # Core models: GasTable, SteamTable, HRSG, Exergy
├── validate_benchmark.py   # Benchmark validation script (standalone)
├── test_thermo.py          # Unit tests for thermodynamic functions
│
├── requirements.txt        # Python package dependencies
└── README.md               # This file
```

---

## Thermodynamic Models

### Gas Cycle (Brayton)

| Property | Method |
|----------|--------|
| `cp(T)` [kJ/kg·K] | 4th-order polynomial (Borgnakke & Sonntag) — valid 250–2000 K |
| `h(T)` [kJ/kg] | Simpson rule integration of `cp(T)` from 0 K |
| `s(T,P)` [kJ/kg·K] | Simpson integration + `−R·ln(P/P_ref)` term |
| Isentropic temperature | Bisection solver on entropy equality |
| Compressor/Turbine work | Isentropic efficiency `η` applied to ideal `Δh` |
| Combustor | Turbine Inlet Temperature (TIT) constrained; optional combustion efficiency `η_cc` |

### Steam Cycle (Rankine)

| Property | Method |
|----------|--------|
| `T_sat(P)` | IAPWS-IF97 4th-degree approximation |
| `hf`, `hg`, `hfg` | Correlation fits vs saturation pressure |
| Superheated `h(P,T)` | Linear `cp_steam` correction above saturation |
| Superheated `s(P,T)` | Log-law correction above saturation |
| Quality `x` | Two-phase entropy interpolation |
| Pump work | `w_p = v_f × ΔP / η_fp` |

### HRSG Coupling

- Gas turbine exhaust `→` HRSG `→` boiler feedwater preheating
- Pinch-point minimum `ΔT = 15 K` enforced
- Steam mass flow: `ṁ_steam = Q_recovered / Δh_boiler`

### Second Law (Exergy)

- Flow exergy: `ė = (h − h₀) − T₀(s − s₀)` [kJ/kg]
- Component irreversibility: `İ = T₀ · ΔṠ_gen` [kW]
- Overall second-law efficiency: `η_II = Ẇ_net / Ė_fuel`
- Entropy generation rate: `Ṡ_gen = İ_total / T₀` [kW/K]

### AD-HTC Biomass Balance

| Parameter | Expression |
|-----------|-----------|
| Biogas volumetric flow | `V̇ = ṁ_rich × y_AD` m³/s |
| Biogas mass flow | `ṁ_biogas = V̇ × ρ_biogas` (≈ 1.15 kg/m³) |
| Biogas energy available | `Ė_biogas = ṁ_biogas × LHV` kW |
| Combust fuel demand | `ṁ_fuel = (q_in × ṁ_air) / LHV` kg/s |
| Renewable fraction | `min(Ė_biogas / Ė_demand, 100%)` % |
| HTC thermal output | `Q̇_HTC ≈ ṁ_lean × 1.5 × (T_HTC − 298)` kJ/s |

---

## Installation & Setup

### Web App

No installation needed:

```bash
# Windows
start index.html

# macOS/Linux
open index.html
```

### Python Desktop App

```bash
# Clone
git clone https://github.com/brickflows/meg315-assignment7.git
cd meg315-assignment7

# Install
pip install -r requirements.txt

# Run
python app.py
```

### Run Validation

```bash
python validate_benchmark.py
```

---

## Usage

### Web App Inputs

| Section | Parameters |
|---------|-----------|
| **Gas Cycle (Brayton)** | T₁, P₁, rp, TIT, η_c, η_t, LHV, ṁ_air |
| **Steam Cycle (Rankine)** | P_boiler, T_steam, P_cond, η_st, η_fp |
| **Biomass / AD** | ṁ_biomass, moisture split, AD yield, T_HTC |

### Default Design Point

| Parameter | Value |
|-----------|-------|
| Ambient Temp T₁ | 298 K |
| Pressure Ratio rp | 10 |
| Turbine Inlet Temp | 1400 K |
| Compressor Eff. η_c | 0.86 |
| Turbine Eff. η_t | 0.89 |
| Boiler Pressure | 4.0 MPa |
| Steam Temp | 673 K |
| Condenser Pressure | 10 kPa |
| Biomass Feed | 5 kg/s |
| Moisture-rich Fraction | 60% |
| AD Biogas Yield | 0.4 m³/kg |

---

## Validation

The `validate_benchmark.py` script runs a full cycle against a published benchmark specification and reports:

- State-point temperatures (K) — States 1–5
- State-point enthalpies (kJ/kg)
- Net specific work `w_net` (kJ/kg)
- Thermal efficiency `η_Brayton` (%)
- Back-work ratio BWR (%)

Results are printed as `[value | reference | Δ | Δ%]` for each quantity.

---

*© 2025 Energhx Research Group — Faculty of Engineering, University of Lagos | MEG 315 Assignment 7*
