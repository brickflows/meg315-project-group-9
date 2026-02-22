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

![AD-HTC Cycle Schematic](https://raw.githubusercontent.com/brickflows/meg315-project-group-9/main/preview.png)

> **Live version**: Open [`index.html`](index.html) in your browser to see the fully animated schematic with flowing particle animations, spinning generator, and pulsing state-point indicators.

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
- **Correct turbomachinery symbols** — trapezoid Compressor (blue) and Turbine (red) shapes
- **Spinning generator** symbol at turbine exit
- **Pulsing state-point indicators** (①②③④⑤)
- **Real-time thermodynamic calculations** on parameter input
- **H-s diagram** — Rankine steam cycle (with saturation dome)
- **T-h diagram** — Brayton gas cycle
- **State-point tables** — full thermodynamic properties at each state
- **AD-HTC mass & energy balance** — biogas supply vs fuel demand

### 🖥️ Desktop Application (`app.py`)
- Full **Tkinter GUI** with scrollable parameter panel
- **HRSG coupling** between gas and steam cycles
- **Second-law (exergy) analysis** — irreversibility per component
- **Engineering validation warnings** — flags physically unreasonable inputs
- **Collapsible sections** for optional parameters

---

## File Structure

```
meg315-project-group-9/
│
├── index.html              # Web application (main entry point)
├── styles.css              # Dark-theme styling with glassmorphism
├── script.js               # Thermodynamic engine & Chart.js rendering
│
├── app.py                  # Python/Tkinter desktop GUI
├── thermodynamics.py       # Core property models (GasTable, SteamTable, HRSG, Exergy)
├── validate_benchmark.py   # Standalone benchmark validation script
│
├── requirements.txt        # Python dependencies
└── README.md               # This file
```

---

## Thermodynamic Models

### Gas Cycle (Brayton)

| Property | Method |
|----------|--------|
| `cp(T)` | 4th-order polynomial — valid 250–2000 K |
| `h(T)` | Simpson integration of `cp(T)` |
| `s(T,P)` | Simpson integration + ideal-gas pressure correction |
| Isentropic T | Bisection solver on `s(T,P) = const` |
| Work | Isentropic efficiency applied to ideal work |

### Steam Cycle (Rankine)

| Property | Method |
|----------|--------|
| `T_sat(P)` | IAPWS-IF97 approximation |
| `hf`, `hg`, `hfg` | Correlation fits vs saturation pressure |
| Superheated `h`, `s` | Linear `cp` correction above saturation |
| Quality `x` | Two-phase interpolation on entropy |

### HRSG Coupling
- Pinch-point constraint enforced (`ΔT_pinch ≥ 15 K`)
- Steam mass flow from: `Q_recovered = ṁ_steam × Δh_boiler`

### Second Law (Exergy)
- Flow exergy: `e = (h − h₀) − T₀(s − s₀)`
- Irreversibility: `I = T₀ × Ṡ_gen`
- Second-law efficiency: `η_II = Ẇ_net / Ė_fuel`

### AD-HTC Balance

| Parameter | Model |
|-----------|-------|
| Biogas production | `V̇ = ṁ_rich × AD_yield` |
| Biogas energy | `Ė = ṁ_biogas × LHV` |
| HTC thermal output | `Q̇ ≈ ṁ_lean × 1.5 × (T_HTC − 298)` kJ/s |
| Renewable fraction | `min(Ė_biogas / Ė_demand, 100%)` |

---

## Installation & Setup

### Web App (No install needed)

```bash
# Just open in browser
start index.html
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
| Pressure Ratio | 10 |
| Turbine Inlet Temp | 1400 K |
| Compressor Efficiency | 86% |
| Turbine Efficiency | 89% |
| Boiler Pressure | 4.0 MPa |
| Steam Temperature | 673 K |
| Condenser Pressure | 10 kPa |
| Biomass Feed | 5 kg/s |
| Moisture-rich Fraction | 60% |
| AD Biogas Yield | 0.4 m³/kg |

---

## Validation

```bash
python validate_benchmark.py
```

Checks solver accuracy against reference values for state-point temperatures, enthalpies, specific work, and thermal efficiency. Errors are printed as absolute and percentage deviations.

---

## Team

**MEG 315 — Group 9**
Faculty of Engineering, University of Lagos
*Energhx Research Group*

---

## License

Submitted as a course assignment for **MEG 315 (Engineering Thermodynamics)**, University of Lagos. All thermodynamic models and code are original work by Group 9.
