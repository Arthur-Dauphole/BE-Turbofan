# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a SUPAERO 2nd-year aerospace engineering assignment (BE 10 — Modèle Turbofan). The goal is to build a parametric thermodynamic cycle model of a dual-spool, dual-flux (turbofan) engine, modeled after the LEAP-1A at cruise conditions.

The model is to be implemented in Python (or Excel/Matlab). No code exists yet — this repo contains only the assignment PDF (`BE-Turbofan-2025.pdf`).

## Flight Conditions & Reference Parameters (LEAP-1A at cruise)

- M₀ = 0.8, altitude = 35 000 ft → P₀ = 227 hPa, T₀ = 217 K
- Target thrust: F = 21 000 N
- BPR (bypass ratio) λ = ṁs/ṁp = 11
- Fan compression ratio πf = 1.45
- OPR πc = Pt2/Pt2 = 40 (global compression ratio)
- HP compressor ratio πCHP = Pt3/Pt25 = 22
- Combustion temperature Tt4 = 1600 K
- Size constraint: diameter < 2 m (rmax), hub-to-tip ratio rmin/rmax ≈ 0.3

## Gas Properties

- Air (before combustion): r = 287 J·kg⁻¹·K⁻¹, γ = 1.4
- Burnt gases (after combustion): r* = 291.6 J·kg⁻¹·K⁻¹, γ* = 1.33
- Kerosene lower heating value: Pk = 42 800 kJ/kg

## Component Efficiencies & Loss Coefficients

| Parameter | Symbol | Value |
|---|---|---|
| Inlet pressure loss | ξe | 0.98 |
| Polytropic efficiency (compressors BP & HP) | ηc,BP / ηc,HP | 0.90 |
| Polytropic efficiency (fan) | ηf | 0.92 |
| Combustion efficiency | ηcomb | 0.99 |
| Combustion chamber pressure loss | ξcc | 0.95 |
| Mechanical shaft efficiency | ηm | 0.98 |
| Polytropic efficiency (HP turbine) | ηt,HP | 0.89 |
| Polytropic efficiency (BP turbine) | ηt,BP | 0.90 |
| Nozzle pressure loss | ξtuy | 0.98 |

## Engine Architecture & Station Numbering

Dual-spool, dual-flux turbofan with adapted nozzles. Station indices from the diagram:

- 0: freestream
- 2: fan inlet (after intake)
- 2.5: fan exit / LP compressor exit / HP compressor inlet
- 21: secondary (bypass) duct inlet
- 3: HP compressor exit
- 4: combustion chamber exit (turbine inlet)
- 4.5: HP turbine exit / LP turbine inlet
- 5: LP turbine exit
- 7/9: core nozzle exit
- 12/17/18/19: bypass duct and secondary nozzle stations

Key constraint: LP compressor ratio πcBP is derived from πf and πCHP to satisfy OPR = πf · πcBP · πCHP.

## Outputs the Model Must Compute

- Primary mass flow rate ṁp (and secondary ṁs = λ · ṁp)
- Engine outer radius rmax
- Fuel mass flow rate ṁk
- Specific fuel consumption: Cs = ṁk/F (kg·h⁻¹·daN⁻¹)
- Specific thrust: f^sp = F/ṁ (m·s⁻¹)
- Thermal efficiency: ηth = Ẇcy / Ẇchim
- Propulsive efficiency: ηpr = Ẇpr / Ẇcy
- Overall efficiency: η = ηth · ηpr

Where:
- Ẇpr = F · V₀
- Ẇcy = ½ṁp(1+λ)V₉² − V₀²) + ½ṁs(V₁₉² − V₀²)
- Ẇchim = ṁk · Pk

## Parametric Studies Required (hors séance report)

1. Turbofan vs turbojet efficiency comparison
2. Thermal efficiency vs OPR for several Tt4 values
3. Overall efficiency vs BPR for several fan compression ratio values
4. How to evolve OPR, BPR, Tt4, πf to improve LEAP performance within size constraint
5. Effect of reducing M₀ from 0.8 to 0.7 on engine size and performance (T-s diagram)

## Fan Inlet Mach Number Constraint

Fan tip Mach number should be around 0.6 to minimize compressibility effects at blade tips. This constrains the fan inlet area and hence rmax given the hub-to-tip ratio rmin/rmax ≈ 0.3.
