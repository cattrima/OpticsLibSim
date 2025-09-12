# OpticsLibSim
A dedicated optical library with a detailed modeling of the optical system and Point Spread Function (PSF) for Geant4 simulation

1. Install Geant4
Make sure you have Geant4 compiled with visualization and UI support.
For example:
cmake -DCMAKE_INSTALL_PREFIX=/path/to/geant4-install ..
make -jN
make install

Then source Geant4 each time before building/running (put in the bash)

source /path/to/geant4-install/bin/geant4.sh


2. Build this project

From the project root (OpticsLibSim)
mkdir build
cd build
cmake ..
make -jN

---> build/

2. Run the program

Interactive mode with visualization (no macro argument)
./G4OpticsLibSim

This will load the visualization macro (init_vis.mac)

Batch mode
./G4OpticsLibSim run.mac

This exectues the macro run.mac without opening the UI


Structure
OpticsLibSim/
├── CMakeLists.txt
├── include/
│   ├── DetectorConstruction.hh
│   ├── ActionInitialization.hh
│   └── PrimaryGeneratorAction.hh
└── src/
    ├── G4OpticsLibSim.cc     
    ├── DetectorConstruction.cc
    ├── ActionInitialization.cc
    ├── PrimaryGeneratorAction.cc
