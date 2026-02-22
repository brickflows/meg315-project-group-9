# AD-HTC Fuel-Enhanced Power Gas Cycle Analyser
### MEG 315 — Group 9 | Energhx Research Group, University of Lagos

> An interactive thermodynamic simulation of a combined **Anaerobic Digestion (AD)** + **Hydrothermal Carbonization (HTC)** fuel-enhanced gas power cycle — with a coupled Rankine steam cycle and real-gas property models.

---

## 📋 Table of Contents

- [Overview](#overview)
- [System Architecture](#system-architecture)
- [Features](#features)
- [Live Demo (Web App)](#live-demo-web-app)
- [Desktop App (Python/Tkinter)](#desktop-app-pythontkinter)
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
└─────────┬───────────┘
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
│ Reactor│     │(Biogas) │
└────┬───┘     └────┬────┘
     │              │
     ▼              ▼
┌─────────┐   ┌──────────────────┐
│  Boiler │   │ Enhanced Biogas  │
│(Steam)  │   │   Collector      │
└────┬────┘   └──────┬───────────┘
     │               │         │
     ▼               ▼         ▼
  Rankine     Combustion   Biogas to
  Cycle        Chamber     Buildings
     │               │
     ▼               ▼
  Steam         Brayton Cycle
  Turbine    (Compressor → Turbine)
     │               │
     └───────────────┘
              │
              ▼
         Net Power Output
         + Exhaust Gases
```

---

## Features

### 🌐 Web Application (`index.html`)
- **Interactive animated schematic** — live particle flow animations showing every mass/energy stream
- **Trapezoid compressor & turbine** shapes (engineering-accurate turbomachinery symbols)
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
- **Navigation toolbar** on all charts

---

## Live Demo (Web App)

To run the web app locally:

1. Open `index.html` in any modern browser (Chrome, Edge, Firefox)
2. No server or build step required — it's pure HTML/CSS/JavaScript

**Or** view the deployed version on GitHub Pages *(if enabled on this repo)*.

---

## Desktop App (Python/Tkinter)

### Requirements

```
Python 3.10+
matplotlib
numpy
pandas
tkinter (built-in)
```

Install dependencies:

```bash
pip install -r requirements.txt
```

### Run

```bash
python app.py
```

---

## File Structure

```
meg315-project-group-9/
│
├── index.html              # Web application (main entry point)
├── styles.css              # Web app styling (dark theme, glassmorphism)
├── script.js               # Web app thermodynamic engine & chart rendering
│
├── app.py                  # Python/Tkinter desktop GUI
├── thermodynamics.py       # Core thermodynamic property models (GasTable, SteamTable, HRSG, Exergy)
├── validate_benchmark.py   # Standalone validation against benchmark specifications
│
├── requirements.txt        # Python dependencies
└── README.md               # This file
```

---

## Thermodynamic Models

### Gas Cycle (Brayton)

| Property | Method |
|----------|--------|
| Specific heat `cp(T)` | 4th-order polynomial (Borgnakke & Sonntag) — valid 250–2000 K |
| Enthalpy `h(T)` | Simpson integration of `cp(T)` |
| Entropy `s(T,P)` | Simpson integration + ideal-gas pressure correction |
| Isentropic temperature | Bisection solver on `s(T,P) = const` |
| Compressor & Turbine | Isentropic efficiency applied to ideal work |

### Steam Cycle (Rankine)

| Property | Method |
|----------|--------|
| Saturation temperature | IAPWS-IF97 approximation |
| `hf`, `hg`, `hfg` | Correlation fits vs saturation pressure |
| Superheated steam `h`, `s` | Linear `cp` correction above saturation |
| Quality `x` | Two-phase interpolation on entropy |

### HRSG Coupling

- Gas turbine exhaust recovers heat via a **Heat Recovery Steam Generator (HRSG)**
- Pinch-point constraint enforced (`ΔT_pinch ≥ 15 K` default)
- Steam mass flow computed from energy balance: `Q_recovered = m_steam × (h_boiler_exit − h_pump_exit)`

### Second Law (Exergy)

- Flow exergy per state: `e = (h − h₀) − T₀(s − s₀)`
- Irreversibility per component: `I = T₀ × Ṡ_gen`
- Second-law efficiency: `η_II = Ẇ_net / Ė_fuel`

### AD-HTC Biomass Balance

| Parameter | Model |
|-----------|-------|
| Biogas production | `V̇_biogas = ṁ_rich × AD_yield` |
| Biogas energy | `Ė_biogas = ṁ_biogas × LHV` |
| HTC thermal output | `Q̇_HTC ≈ ṁ_lean × 1.5 × (T_HTC − 298)` kJ/s |
| Renewable fraction | `min(Ė_biogas / Ė_demand, 100%)` |

---

## Installation & Setup

### Web App

No installation needed. Open `index.html` directly in a browser.

### Python Desktop App

```bash
# Clone the repo
git clone https://github.com/brickflows/meg315-project-group-9.git
cd meg315-project-group-9

# Install dependencies
pip install -r requirements.txt

# Launch the GUI
python app.py
```

### Validation Script

```bash
python validate_benchmark.py
```

Runs a known benchmark case and prints errors relative to reference values for each state point and performance metric.

---

## Usage

### Web App Parameter Guide

| Section | Key Inputs |
|---------|-----------|
| **Gas Cycle** | Ambient T & P, Pressure Ratio `rp`, Turbine Inlet Temp `TIT`, Efficiencies |
| **Steam Cycle** | Boiler pressure & temperature, Condenser pressure, Efficiencies |
| **Biomass / AD** | Biomass feed rate, Moisture-rich fraction, AD biogas yield, HTC reactor temp |

Click **Analyse** to compute all cycles, populate charts, and display state-point tables.

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

Run the benchmark validator to verify solver accuracy:

```bash
python validate_benchmark.py
```

**Expected outputs** are checked against:
- Compressor exit temperature (State 2)
- Turbine exhaust temperature (State 5)
- Net specific work `w_net`
- Thermal efficiency `η_Brayton`
- Steam cycle state enthalpies

Errors are reported as **absolute** and **percentage** deviations from reference values.

---

## Team

**MEG 315 — Group 9**
Faculty of Engineering, University of Lagos

> *Energhx Research Group*

---

## License

This project is submitted as a course assignment for **MEG 315 (Engineering Thermodynamics)** at the University of Lagos. All thermodynamic models and code are original work by Group 9.
