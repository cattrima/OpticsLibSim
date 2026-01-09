// OpticsLibSim - Optical Simulation Library for Geant4
// Copyright (c) 2025 Caterina Trimarelli
// Licensed under the MIT License (see LICENSE file in the project root)
// SPDX-License-Identifier: MIT
#include "OpticsLib/FocalPlaneArray.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "OpticsLib/OpticalMaterial.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
using namespace OpticsLib;

FocalPlaneArray::FocalPlaneArray(const std::string& name,
                                 G4double width,
                                 G4double height,
                                 G4double thickness,
                                 G4Material* mat)
    : OpticalComponent(name),
      fWidth(width), fHeight(height), fThickness(thickness),
      fMaterial(mat ? mat : OpticalMaterial::GetAir())
{
}
G4VPhysicalVolume* FocalPlaneArray::Place(G4LogicalVolume* mother,
                                          const G4ThreeVector& pos,
                                          G4RotationMatrix* rot)
{
    auto solid = new G4Box(fName, fWidth/2., fHeight/2., fThickness/2.);
    logical_ = new G4LogicalVolume(solid, fMaterial, fName+"LV");

    // vis attributes
    auto visAttrib = new G4VisAttributes(G4Color(1.0,0.0,0.0)); // red
    visAttrib->SetForceSolid(true);
    logical_->SetVisAttributes(visAttrib);

    return new G4PVPlacement(rot, pos, logical_, logical_->GetName(), mother, false, 0, true);
}
