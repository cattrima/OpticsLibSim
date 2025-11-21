// SPDX-License-Identifier: MIT
#include "OpticsLib/OpticalMirror.hh"

// Geant4 headers
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4RotationMatrix.hh"
#include <cmath>
#include <vector>
#include <sstream>
using namespace OpticsLib;
namespace OpticsLib {
    OpticsConfig global_cfg;  
    OpticsConfig* cfg = &global_cfg;  
}

// ----------------------------- Profiles ------------------------------------

G4double OpticalMirror::M1ProfileFunction(G4double x) {
    if (!cfg_) return -999.0;
    if (x == 0) x = 1e-20;

    G4double c       = cfg_->M1_c;
    G4double k       = cfg_->M1_k;
    G4double alpha_4 = cfg_->M1_alpha_4;
    G4double alpha_6 = cfg_->M1_alpha_6;

    G4double term1 = (c*pow(x,2)) / (1 + sqrt(1 - (1+k)*c*c*pow(x,2)));
    G4double term2 = alpha_4*pow(x,4);
    G4double term3 = alpha_6*pow(x,6);

    return -(term1 + term2 + term3);
}

G4double OpticalMirror::M2ProfileFunction(G4double x) {
    if (!cfg_) return -999.0;
    if (x == 0) x = 1e-21;

    G4double c       = cfg_->M2_c;
    G4double k       = cfg_->M2_k;
    G4double alpha_4 = cfg_->M2_alpha_4;
    G4double alpha_6 = cfg_->M2_alpha_6;
    G4double alpha_8 = cfg_->M2_alpha_8;
    G4double alpha_10= cfg_->M2_alpha_10;

    G4double term1 = (c*pow(x,2)) / (1 + sqrt(1-(1+k)*c*c*pow(x,2)));
    G4double term2 = alpha_4*pow(x,4);
    G4double term3 = alpha_6*pow(x,6);
    G4double term4 = alpha_8*pow(x,8);
    G4double term5 = alpha_10*pow(x,10);

    return -(term1 + term2 + term3 + term4 + term5);
}

// ----------------------------- Constructors --------------------------------
// ----------------------------- Constructors (config-based) -----------------

// Primary mirror constructor
OpticalMirror::OpticalMirror(const std::string& name,
                             G4Material* mirrorMaterial,
                             MirrorKind kind)
    : kind_(kind), fName(name), fMaterial(mirrorMaterial),
      reflectivity_(0.9), cfg_(cfg) // default reflectivity 0.9 se non presente
{
    if(cfg_) {
        fRadius          = cfg_->M1_radius;
        fThickness       = cfg_->M1_thickness;
        nn_ring_M1_      = cfg_->M1_nn_ring;
        nn_ring_start_M1_= cfg_->M1_nn_ring_start;
    }

    // fallback solid & logical
    solid_ = new G4Tubs(fName, 0., fRadius, fThickness/2., 0.*deg, CLHEP::twopi);
    logical_ = new G4LogicalVolume(solid_, fMaterial, fName + "LV");

    opticalSurface_ = new G4OpticalSurface(fName + "_OpticalSurface");
    opticalSurface_->SetType(dielectric_metal);
    opticalSurface_->SetFinish(polished);
    opticalSurface_->SetModel(unified);

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    G4double photonEnergy[] = {2.0*eV, 3.5*eV};
    G4double refl[] = {reflectivity_, reflectivity_};
    mpt->AddProperty("REFLECTIVITY", photonEnergy, refl, 2);
    opticalSurface_->SetMaterialPropertiesTable(mpt);

    new G4LogicalSkinSurface(fName + "Skin", logical_, opticalSurface_);

    auto visAttrib = new G4VisAttributes(G4Colour(0.78,0.78,0.80));
    visAttrib->SetForceSolid(true);
    logical_->SetVisAttributes(visAttrib);
}

// Secondary simple mirror constructor
OpticalMirror::OpticalMirror(const std::string& name,
                             G4Material* mirrorMaterial,
                             MirrorKind kind,
                             bool /*secondaryMarker*/)
    : kind_(kind), fName(name), fMaterial(mirrorMaterial),
      reflectivity_(0.9), cfg_(cfg)
{
    if(cfg_) {
        fRadius    = cfg_->M2_radius;
        fThickness = cfg_->M2_thickness;
        nn_ring_M2_= cfg_->M2_nn_ring;
    }

    solid_ = new G4Tubs(fName, 0., fRadius, fThickness/2., 0.*deg, CLHEP::twopi);
    logical_ = new G4LogicalVolume(solid_, fMaterial, fName + "LV");

    opticalSurface_ = new G4OpticalSurface(fName + "_OpticalSurface");
    opticalSurface_->SetType(dielectric_metal);
    opticalSurface_->SetFinish(polished);
    opticalSurface_->SetModel(unified);

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    G4double photonEnergy[] = {2.0*eV, 3.5*eV};
    G4double refl[] = {reflectivity_, reflectivity_};
    mpt->AddProperty("REFLECTIVITY", photonEnergy, refl, 2);
    opticalSurface_->SetMaterialPropertiesTable(mpt);

    new G4LogicalSkinSurface(fName + "Skin", logical_, opticalSurface_);

    auto visAttrib = new G4VisAttributes(G4Colour(0.78,0.78,0.80));
    visAttrib->SetForceSolid(true);
    logical_->SetVisAttributes(visAttrib);
}

// Secondary + Lens constructor
OpticalMirror::OpticalMirror(const std::string& name,
                             G4Material* lensMaterial,
                             G4Material* mirrorMaterial,
                             MirrorKind kind)
    : kind_(kind), fName(name), fMaterial(mirrorMaterial),
      fLensMaterial(lensMaterial), reflectivity_(0.9), cfg_(cfg)
{
    if(cfg_) {
        lensRadius_    = cfg_->Lens_radius;
        lensThickness_ = cfg_->Lens_thickness;
        fRadius        = cfg_->M2_radius;
        fThickness     = cfg_->M2_thickness;
        nn_ring_M2_    = cfg_->M2_nn_ring;
    }

    solid_ = new G4Tubs(fName, 0., fRadius, fThickness/2., 0.*deg, CLHEP::twopi);
    logical_ = new G4LogicalVolume(solid_, fMaterial, fName + "LV");

    opticalSurface_ = new G4OpticalSurface(fName + "_OpticalSurface");
    opticalSurface_->SetType(dielectric_metal);
    opticalSurface_->SetFinish(polished);
    opticalSurface_->SetModel(unified);

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    G4double photonEnergy[] = {2.0*eV, 3.5*eV};
    G4double refl[] = {reflectivity_, reflectivity_};
    mpt->AddProperty("REFLECTIVITY", photonEnergy, refl, 2);
    opticalSurface_->SetMaterialPropertiesTable(mpt);

    new G4LogicalSkinSurface(fName + "Skin", logical_, opticalSurface_);

    auto visAttrib = new G4VisAttributes(G4Colour(0.78,0.78,0.80));
    visAttrib->SetForceSolid(true);
    logical_->SetVisAttributes(visAttrib);
}

// --------------------------- Build() dispatcher -----------------------------

void OpticalMirror::Build(G4LogicalVolume* mother, const G4ThreeVector& pos, G4RotationMatrix* rot) {
    switch (kind_) {
        case PRIMARY_M1:
            BuildPrimaryM1(mother);
            break;
        case SECONDARY_M2_WITH_LENS:
            BuildSecondaryM2WithLens(mother);
            break;
        default:
            BuildPrimaryM1(mother);
            break;
    }
}

// Simple Place (single-cylindrical mirror)
G4VPhysicalVolume* OpticalMirror::Place(G4LogicalVolume* mother, const G4ThreeVector& pos, G4RotationMatrix* rot) {
    return new G4PVPlacement(rot, pos, logical_, fName + "_Phys", mother, false, 0, true);
}

// -------------------------- Helper (placeholder) ---------------------------
void OpticalMirror::get_cons_ring_par(G4double x0, G4double y0, G4double x2, G4double y2,
                                      G4double dR, G4double& Rmin_out, G4double& Rmax_out, G4double& h_out)
{
    Rmin_out = x0;
    Rmax_out = x2;
    h_out = std::abs(y2 - y0);
}
// -------------------------- Build Primary M1 -------------------------------
void OpticalMirror::BuildPrimaryM1(G4LogicalVolume* mother) {
    if(!cfg_) return; // sanity check

    G4RotationMatrix* Ra_optics = new G4RotationMatrix();
    Ra_optics->rotateY(180.0*deg);

    // Envelope
    G4Material* envelopeMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4double envelopeThickness = 20.08*mm;
    G4Tubs* M1_cyl = new G4Tubs(fName + "_Envelope", 110.0*mm, cfg_->M1_radius+1.*mm, envelopeThickness/2.0, 0, CLHEP::twopi);
    G4LogicalVolume* envelopeLogical = new G4LogicalVolume(M1_cyl, envelopeMat, fName + "_EnvelopeLV");
    new G4PVPlacement(Ra_optics,
                      G4ThreeVector(0,0,-cfg_->M1_to_FPA + envelopeThickness/2.0),
                      envelopeLogical, fName + "_EnvelopePV", mother, false, 0, true);

    // Rings
    int nn_ring = cfg_->M1_nn_ring;
    int start_ring = cfg_->M1_nn_ring_start;
    G4double dl_ring = cfg_->M1_radius / nn_ring;
    G4double ring_opening = 360.*deg;

    // Optical surface
    G4OpticalSurface* ringSurface = new G4OpticalSurface(fName + "_RingSurface");
    ringSurface->SetType(dielectric_metal);
    ringSurface->SetFinish(polished);
    ringSurface->SetModel(unified);

    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    G4double photonEnergy[] = {2.0*eV, 3.5*eV};
    G4double refl[] = {reflectivity_, reflectivity_};
    mpt->AddProperty("REFLECTIVITY", photonEnergy, refl, 2);
    ringSurface->SetMaterialPropertiesTable(mpt);

    for(int i = start_ring; i < nn_ring; ++i) {
        G4double Rmin = dl_ring * i;
        G4double Rmax = dl_ring * (i + 1);

        G4double y0 = M1ProfileFunction(Rmin);
        G4double y2 = M1ProfileFunction(Rmax);
        G4double htubs_ring = y2;
        G4double hcons_ring = y2 - y0;

        G4Tubs* lens_ring_solid = new G4Tubs(fName + "_lens_ring_" + std::to_string(i),
                                             Rmin, Rmax, htubs_ring/2.0, 0, ring_opening);

        G4double cons_Rmin, cons_Rmax, cons_h;
        get_cons_ring_par(Rmin, y0, Rmax, y2, 0.05, cons_Rmin, cons_Rmax, cons_h);

        G4Cons* cons_solid = new G4Cons(fName + "_cons_" + std::to_string(i),
                                        0.0, cons_Rmax, 0.0, cons_Rmin, cons_h/2.0, 0.0, ring_opening);

        G4RotationMatrix* rotMatrix = new G4RotationMatrix();
        G4ThreeVector transVector(0.0, 0.0, -htubs_ring/2.0 + hcons_ring/2.0);

        G4SubtractionSolid* ring_m_cons_solid = new G4SubtractionSolid(fName + "_ring_m_cons_" + std::to_string(i),
                                                                       lens_ring_solid, cons_solid, rotMatrix, transVector);

        G4LogicalVolume* ringLogical = new G4LogicalVolume(ring_m_cons_solid, fMaterial, fName + "_RingLV_" + std::to_string(i));
        G4ThreeVector ringPos(0,0,-htubs_ring/2.0 + envelopeThickness/2.0);

        if(if_M1_) {
            new G4PVPlacement(nullptr, ringPos, ringLogical, fName + "_RingPV_" + std::to_string(i), envelopeLogical, false, i, true);
        }

        new G4LogicalSkinSurface(fName + "_RingSkin_" + std::to_string(i), ringLogical, ringSurface);
    }
}

// --------------------- Build Secondary + Lens + Rings ----------------------
void OpticalMirror::BuildSecondaryM2WithLens(G4LogicalVolume* mother) {
    if(!cfg_) return;

    G4RotationMatrix* Ra_optics = new G4RotationMatrix();
    Ra_optics->rotateY(180.0*deg);

    // Lens
    G4VSolid* sphere_solid = new G4Sphere(fName + "_sphere", 0.0, cfg_->Lens_radius, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
    G4Tubs* secondary_cyl = new G4Tubs(fName + "_secondary_cyl", 0, 354.*mm, cfg_->Lens_radius, 0, CLHEP::twopi);


    G4SubtractionSolid* final_secondary = new G4SubtractionSolid(fName + "_final_lens", sphere_solid, secondary_cyl, nullptr, G4ThreeVector(0,0,-36.8*mm)); //the m2 thickness/2+ lens thickness/2

    G4Tubs* secondary_monster = new G4Tubs("secondary_monster",102*mm, 205*mm, 359.168*mm,0,CLHEP::twopi);


    G4SubtractionSolid* final_lens = new G4SubtractionSolid("final_secondary_monster", final_secondary,secondary_monster,0,G4ThreeVector( 0, 0, 0));

    G4Material* lensMat = fLensMaterial ? fLensMaterial : G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4LogicalVolume* lensLogical = new G4LogicalVolume(final_lens, lensMat, fName + "_LensLV");


    G4VisAttributes* visAttr = new G4VisAttributes(G4Colour(0.2,0.6,1.0,0.3));
    visAttr->SetForceSolid(true);  // forza il rendering della geometria interna
    lensLogical->SetVisAttributes(visAttr);

    new G4PVPlacement(Ra_optics, G4ThreeVector(0,0,cfg_->M2_to_FocalPlane+cfg_->Lens_radius), lensLogical, fName + "_LensPV", mother, false, 0, true);

    // M2 Rings
    int nn_ring = cfg_->M2_nn_ring;
    G4double dl_ring = cfg_->M2_radius / nn_ring;

    for(int i=0; i<nn_ring; ++i) {
        G4double Rmin = dl_ring * i;
        G4double Rmax = dl_ring * (i + 1);
        G4double y0 = M2ProfileFunction(Rmin);
        G4double y2 = M2ProfileFunction(Rmax);
        G4double hcons = y2 - y0;
        G4Cons* M2_cons_solid = new G4Cons(fName + "_M2_cons_" + std::to_string(i),
                                           0.0, Rmax, 0.0, Rmin, hcons/2.0, 0.0, CLHEP::twopi);
        G4LogicalVolume* M2_ring_logical = new G4LogicalVolume(M2_cons_solid,
                                                               G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"),
                                                               fName + "_M2_ringLV_" + std::to_string(i));
        G4ThreeVector ringPos(0,0,-y2 + hcons/2.0 + 325.35*mm);
        new G4PVPlacement(nullptr, ringPos, M2_ring_logical, fName + "_M2_ringPV_" + std::to_string(i), lensLogical, false, i, true);
    }
}
