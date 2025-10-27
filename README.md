# 🌌 OpticsLibSim — Geant4 Optical Simulation Library

**Author:** Caterina  
**License:** MIT  
**Description:**  
OpticsLibSim is a modular **optical component library** for **Geant4**, designed to let users easily add and configure optical elements such as mirrors, corrector plates, and focal plane arrays (FPA).  
It is structured to be reusable as a **standalone library** or integrated into any Geant4-based telescope or detector simulation.

---

## 📁 Project Structure

OpticsLibSim/
├── CMakeLists.txt
├── LICENSE
│
├── include/
│ └── OpticsLib/
│ ├── OpticalMirror.hh # Generic mirror class (primary/secondary)
│ ├── CorrectorPlate.hh # Optical corrector element
│ ├── FocalPlaneArray.hh # Focal plane / detector surface
│ └── OpticalMaterial.hh # Shared material definitions (refractive indices, etc.)
│
├── src/
│ ├── OpticalMirror.cc
│ ├── CorrectorPlate.cc
│ ├── FocalPlaneArray.cc
│ └── OpticalMaterial.cc
│
├── data/
│ └── mirror_reflectivity.txt # Example reflectivity file (energy [eV], reflectivity)
│
└── demo/
├── CMakeLists.txt
├── DetectorConstruction.hh
├── DetectorConstruction.cc
├── main.cc # Example simulation using the library
├── init_vis.mac # Visualization macro
└── run.mac # Example run macro


---

## 🪞 Core Classes Overview

### **1. OpticalMirror**
> Generic optical mirror used for **primary**, **secondary**, or **tertiary** reflectors.

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

**Example:**
```cpp
auto glass = OpticalMaterial::Get("BK7");
auto aluminum = OpticalMaterial::Get("Aluminum");


## ⚙️ Dependencies
OpticsLibSim depends on:

Geant4 (≥ 11.0)

CMake (≥ 3.16)

A C++17 compiler (e.g. GCC ≥ 9.0)

Optional: Qt (for visualization /vis/open Qt)

Make sure to source your Geant4 environment before building:
source /path/to/geant4-install/bin/geant4.sh



## 🏗️ CMake Configuration

Top-level CMakeLists.txt

## 🚀 Demo Guide
The demo folder provides a ready-to-run Geant4 example that uses the library.


## 🧪 Building & Running

mkdir build && cd build
cmake ..
make
./OpticsDemo
If Qt is enabled, the Geant4 visualization window will appear showing your optical system.

## 🧾 License (MIT)
Each source file starts with:
// OpticsLibSim - Optical Simulation Library for Geant4
// Copyright (c) 2025 Caterina
// Licensed under the MIT License (see LICENSE file in the project root)
And the LICENSE file contains the full MIT text.