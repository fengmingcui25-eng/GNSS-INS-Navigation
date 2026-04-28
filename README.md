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
<img width="1527" height="353" alt="17" src="https://github.com/user-attachments/assets/68ee02c4-254b-48ce-801a-9fc13c6736a5" />
<img width="1567" height="354" alt="18" src="https://github.com/user-attachments/assets/7f4284f1-c745-4cbd-b31c-d1a35afe199b" />
<img width="1545" height="353" alt="19" src="https://github.com/user-attachments/assets/fa48ad0c-a377-47e2-abca-a3a40ddaecd4" />


### 2. INS Mechanization (`/INS_Mechanization`)
- **Earth Model**: WGS-84 ellipsoidal model implementation.
- **Rotation & Quaternions**: Robust attitude representation using Quaternions and Rotation Matrices to avoid singularity.
- **State Update**: Mechanization equations for updating position, velocity, and attitude based on IMU inputs.
<img width="416" height="208" alt="image" src="https://github.com/user-attachments/assets/5b309c89-3874-42f2-b616-adcf2a9f9c61" />
<img width="416" height="208" alt="image" src="https://github.com/user-attachments/assets/8946ada6-5765-410e-a715-3014fff37287" />


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
