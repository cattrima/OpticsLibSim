# 🌌 OpticsLibSim — Geant4 Optical Simulation Library

## 🎄 OpticsLibSim - Holiday Update - First Release 🎄

The first release of OpticsLibSim is expected in January 2026. This initial version will include two optical configurations.

**Author:** Caterina Trimarelli
**License:** MIT  
**Description:**  
OpticsLibSim is a modular **optical component library** for **Geant4**, designed to let users easily add and configure optical elements such as mirrors, corrector plates, and focal plane arrays (FPA).  
It is structured to be reusable as a **standalone library** or integrated into any Geant4-based telescope or detector simulation.
User Guide and first application "OpticsLibSim: A Modular and Lightweight Optical Library for Geant4" soon in arxiv. The initial version of the code was developed for the Terzina Cherenkov Telescope (first version October-December 2023), with collaboration and scientific guidance from Leonid Burmistrov.
https://arxiv.org/abs/2510.02844 
---

## 📁 Project Structure
<pre>
OpticsLibSim/
├── CMakeLists.txt
├── LICENSE
├── include/
│   └── OpticsLib/
│       ├── OpticalMirror.hh        # Generic mirror class (primary/secondary)       #follow the beginning for other modules
│       ├── CorrectorPlate.hh       # Optical corrector element    //currently not used 
│       ├── FocalPlaneArray.hh      # Focal plane / detector surface   
│       └── OpticalMaterial.hh      # Shared material definitions (refractive indices, etc.)
│       └── OpticalComponent.hh     # 
├── src/
│   ├── OpticalMirror.cc             
│   ├── CorrectorPlate.cc
│   ├── FocalPlaneArray.cc
├── config/
|   ├──optics_config.txt
|   |──data/
│	└── mirror_reflectivity.txt    # Example reflectivity file (energy [eV], reflectivity)
│   	└── corrector_lens_refractive_index.dat    # Example refraction for corrector lens quartz file (energy [eV], refraction)
│       └── corrector_lens_absorption.dat    # Example absorption for corrector lens quartz file (energy [eV], absorption)
└── demo/
    ├── CMakeLists.txt
    ├── DetectorConstruction.hh
    ├── DetectorConstruction.cc
    ├── PrimaryGeneratorAction.hh
    ├── PrimaryGeneratorAction.cc
    ├── SteppingAction.hh
    ├── SteppingAction.cc
    ├── FPA_SD.cc                   #not used useful for output
    ├── FPA_SD.hh		    #not used useful for output
    ├── main.cc                     # Example simulation using the library
    ├── init_vis.mac                # Visualization macro
    └── run.mac                     # Example run macro
</pre>

---

## 🪞 Core Classes Overview

### **1. OpticalMirror**
> Generic optical mirror used for **primary**, **secondary** reflectors.

**Key features:**
- Configurable reflectivity (constant or from file)
- Adjustable curvature and diameter
- Automatic creation of optical surface in Geant4

---

### **2. CorrectorPlate**
> Represents an optical corrector (e.g. Schmidt or Maksutov plate).

**Key features:**
- Configurable refractive index and thickness  
- Can use predefined material (BK7, fused silica, etc.)

---

### **3. FocalPlaneArray (FPA)**
> Represents a focal surface or detector plane.

**Key features:**
- Defined by active area and pixel size  
- Placeholder for photodetector behavior (e.g. PMT, SiPM)

---

### **4. OpticalMaterial**
> Centralized factory for shared materials and optical constants.


---

## ⚙️  Dependencies
  

OpticsLibSim depends on:

Geant4 (≥ 11.0)

CMake (≥ 3.16)

A C++17 compiler (e.g. GCC ≥ 9.0)

Optional: Qt (for visualization /vis/open Qt)

Make sure to source your Geant4 environment before building:
source /path/to/geant4-install/bin/geant4.sh


---


##  🏗️  CMake Configuration



Top-level CMakeLists.txt


---


##  🚀  Demo Guide


The demo folder provides a ready-to-run Geant4 example that uses the library.

---


## 🧪 Building & Running


mkdir build && cd build
cmake ..
make
./demo/OpticsDemo
If Qt is enabled, the Geant4 visualization window will appear showing your optical system.

---


##  🧾  License (MIT)
Each source file starts with:
// OpticsLibSim - Optical Simulation Library for Geant4
// Copyright (c) 2025 Caterina Trimarelli
// Licensed under the MIT License (see LICENSE file in the project root)
And the LICENSE file contains the full MIT text.