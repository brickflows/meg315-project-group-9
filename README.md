# AD-HTC Fuel-Enhanced Power Gas Cycle Analyser
### MEG 315 — Group 9 | Energhx Research Group, University of Lagos

> An interactive thermodynamic simulation of a combined **Anaerobic Digestion (AD)** + **Hydrothermal Carbonization (HTC)** fuel-enhanced gas power cycle — with a coupled Rankine steam cycle and real-gas property models.

---

## 🔗 Live Demo & Repository

| | Link |
|---|---|
| 🌐 **Live Web App** | [https://meg315-project-group-9.vercel.app](https://meg315-project-group-9.vercel.app) |
| 📦 **GitHub Repository** | [https://github.com/brickflows/meg315-project-group-9](https://github.com/brickflows/meg315-project-group-9) |

---

## 📸 Schematic Preview

The interactive animated system schematic — showing all energy and mass flows across the AD-HTC combined cycle:

![AD-HTC Cycle Schematic](https://raw.githubusercontent.com/brickflows/meg315-project-group-9/main/ass9-schematic-meg315.png)

> **Live version**: Open [`index.html`](index.html) in your browser or visit the [Live Demo](https://meg315-project-group-9.vercel.app) to see the fully animated schematic with flowing particle animations, spinning generator symbol, and pulsing state-point indicators.

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
- [Team](#team)

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
     │               │
     └───────────────┘
              │
              ▼
         Net Power Output
```

---

## Features

### 🌐 Web Application (`index.html`)
- **Interactive animated schematic** — live particle flow animations on every mass/energy stream
- **Engineering-accurate turbomachinery** — trapezoid Compressor (blue) and Turbine (red)
- **Spinning GEN symbol** at turbine exit (generator output)
- **Pulsing state-point indicators** ①②③④⑤ at each thermodynamic state
- **Animated background grid** for a premium visual feel
- **Real-time thermodynamic calculations** on parameter input
- **H-s diagram** — Rankine steam cycle (with saturation dome)
- **T-h diagram** — Brayton gas cycle
- **State-point tables** — full thermodynamic properties at each state
- **AD-HTC mass & energy balance** — biogas supply vs fuel demand

### 🖥️ Desktop Application (`app.py`)
- **Tkinter GUI** with scrollable parameter panel and tabbed charts
- **HRSG coupling** between gas and steam cycles with pinch-point enforcement
- **Exergy / Second-Law analysis** — irreversibility per component, `η_II`, `Ṡ_gen`
- **Engineering validation warnings**
- **Collapsible advanced sections** for optional parameters

---

## File Structure

```
meg315-project-group-9/
│
├── index.html                  # Web application (open directly in browser)
├── styles.css                  # Dark-theme CSS with glassmorphism
├── script.js                   # Thermodynamic engine + Chart.js rendering
│
├── app.py                      # Python/Tkinter desktop GUI
├── thermodynamics.py           # Core models: GasTable, SteamTable, HRSG, Exergy
├── validate_benchmark.py       # Benchmark validation script
│
├── ass9-schematic-meg315.png   # Schematic screenshot (this README preview)
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

---

## Thermodynamic Models

### Gas Cycle (Brayton)

| Property | Method |
|----------|--------|
| `cp(T)` | 4th-order polynomial — valid 250–2000 K |
| `h(T)` | Simpson integration of `cp(T)` |
| `s(T,P)` | Simpson integration + ideal-gas pressure correction |
| Isentropic T | Bisection solver on entropy equality |
| Work | Isentropic efficiency applied to ideal `Δh` |

### Steam Cycle (Rankine)

| Property | Method |
|----------|--------|
| `T_sat(P)` | IAPWS-IF97 approximation |
| `hf`, `hg`, `hfg` | Correlation fits vs saturation pressure |
| Superheated `h`, `s` | Linear `cp` correction above saturation |
| Quality `x` | Two-phase entropy interpolation |
| Pump work | `w_p = v_f × ΔP / η_fp` |

### HRSG Coupling
- Pinch-point constraint: `ΔT_min = 15 K`
- Steam mass flow: `ṁ_steam = Q_recovered / Δh_boiler`

### Second Law (Exergy)
- Flow exergy: `ė = (h − h₀) − T₀(s − s₀)`
- Irreversibility: `İ = T₀ · ΔṠ_gen`
- Second-law efficiency: `η_II = Ẇ_net / Ė_fuel`

### AD-HTC Balance

| Parameter | Model |
|-----------|-------|
| Biogas flow | `V̇ = ṁ_rich × y_AD` m³/s |
| Biogas energy | `Ė = ṁ_biogas × LHV` kW |
| HTC output | `Q̇ ≈ ṁ_lean × 1.5 × (T_HTC − 298)` kJ/s |
| Renewable % | `min(Ė_biogas / Ė_demand, 100%)` |

---

## Installation & Setup

### Web App

```bash
start index.html        # Windows
open index.html         # macOS/Linux
```

### Python Desktop App

```bash
git clone https://github.com/brickflows/meg315-project-group-9.git
cd meg315-project-group-9
pip install -r requirements.txt
python app.py
```

### Validation

```bash
python validate_benchmark.py
```

---

## Usage

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

```bash
python validate_benchmark.py
```

Verifies solver accuracy against reference benchmark values and prints absolute + percentage errors for each state point and performance metric.

---

## Team

**MEG 315 — Group 9**
Faculty of Engineering, University of Lagos
*Energhx Research Group*

---

*© 2025 Energhx Research Group — Faculty of Engineering, University of Lagos | MEG 315 Group 9*
