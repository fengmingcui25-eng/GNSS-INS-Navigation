# GNSS-INS-Navigation
C++ implementation of GNSS (SPP/RTK) positioning and INS Mechanization algorithms for high-precision navigation.
# GNSS-INS-Navigation

A professional C++ implementation of high-precision positioning and navigation algorithms, featuring GNSS (SPP/RTK) and Inertial Navigation System (INS) mechanization.

## Overview

This repository contains core algorithms for modern navigation systems:
1. **GNSS Engine**: From raw data parsing to centimeter-level RTK positioning.
2. **INS Mechanization**: Precise state updates (Position, Velocity, Attitude) using strapdown inertial navigation principles.

## Features

### 1. GNSS RTK Engine (`/GNSS_RTK_Engine`)
- **Data Parsing**: Support for RINEX 3.x Observation (O), Navigation (N), and Precise Orbit (SP3) files.
- **Positioning Modes**: 
  - **SPP**: Pseudorange-based Single Point Positioning.
  - **RTK**: Carrier-phase Real-Time Kinematic positioning for high-precision solutions.
- **Ambiguity Resolution**: Integrated **LAMBDA** algorithm for integer ambiguity fixing.
- **Satellite Geometry**: Computation of satellite positions and clock errors.

### 2. INS Mechanization (`/INS_Mechanization`)
- **Earth Model**: WGS-84 ellipsoidal model implementation.
- **Rotation & Quaternions**: Robust attitude representation using Quaternions and Rotation Matrices to avoid singularity.
- **State Update**: Mechanization equations for updating position, velocity, and attitude based on IMU inputs.

## Dependencies

- **Eigen 3**: A high-level C++ library for linear algebra, used extensively for matrix operations and coordinate transformations.
- **LAMBDA**: Integrated third-party module for GNSS ambiguity resolution.

## Directory Structure

```text
.
├── GNSS_RTK_Engine/      # GNSS processing core
│   ├── include/          # Header files (.h)
│   ├── src/              # Implementation files (.cpp)
│   └── main.cpp          # GNSS engine entry point
├── INS_Mechanization/    # Inertial navigation core
│   ├── include/          # Header files (.h)
│   ├── src/              # Implementation files (.cpp)
│   └── main.cpp          # INS simulation entry point
└── README.md             # Project documentation
