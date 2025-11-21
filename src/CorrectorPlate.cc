// SPDX-License-Identifier: MIT
#include "OpticsLib/CorrectorPlate.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "OpticsLib/OpticalMaterial.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh" 
using namespace CLHEP;
using namespace OpticsLib;

CorrectorPlate::CorrectorPlate(const std::string& name,
                               G4double width,
                               G4double height,
                               G4double thickness,
                               G4Material* mat)
    : OpticalComponent(name),
      fWidth(width), fHeight(height), fThickness(thickness),
      fMaterial(mat ? mat : OpticalMaterial::GetAir())

{

  
}
G4VPhysicalVolume* CorrectorPlate::Place(G4LogicalVolume* mother,
                                         const G4ThreeVector& pos,
                                         G4RotationMatrix* rot)
{
    solid_ = new G4Box(fName, fWidth/2., fHeight/2., fThickness/2.);
    logical_ = new G4LogicalVolume(solid_, fMaterial, fName+"LV");

    // vis attributes
    auto visAttrib = new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.3)); // cyan trasparente
    visAttrib->SetForceSolid(true);
    logical_->SetVisAttributes(visAttrib);

    return new G4PVPlacement(rot, pos, logical_, logical_->GetName(), mother, false, 0, true);
}
