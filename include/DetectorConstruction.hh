//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file Pnnlsurf/include/DetectorConstruction.hh
/// \brief Definition of the Pnnlsurf::DetectorConstruction class

#ifndef PnnlsurfDetectorConstruction_h
#define PnnlsurfDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GenericMessenger;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    G4Material* CreateCustomDensityMaterial(const G4String& baseName, const G4String& newName, G4double newDensity);
    void SetMaterial(const G4String& name);
    void SetDetector(const G4String& name);
    void SetSampleHeight(G4double height);
    void SetSampleMass(G4double mass);

    void PrintSampleDensity();
    void CalculateAndSetDensity();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4double GetSampleHeight() const { return Sample_height; }
    G4double GetSampleMass() const { return Sample_mass; }
    G4double GetSampleDensity() const { return Sample_density; }

  private:

    void DefineCommands();
    G4GenericMessenger* fMessenger = nullptr;
    G4LogicalVolume* Sample_log = nullptr;  // Your main logical volume
    G4PVPlacement* Sample_phys = nullptr;
    G4double Sample_height = 10.0 * mm;
    G4double Sample_mass = 1000*g;
    G4double Sample_density = 0;
    G4Material* Sample_mat = nullptr;
    G4String DetectorName = "Maeve";

  protected:

    G4LogicalVolume* fScoringVolume = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
