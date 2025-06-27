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
/// \file Pnnlsurf/src/DetectorConstruction.cc
/// \brief Implementation of the Pnnlsurf::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"

#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyhedron.hh"

#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();


  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;


  // materials
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material * Beaker_mat = nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
  G4Material * Ge_mat = nist->FindOrBuildMaterial("G4_Ge");
  G4Material * Cap_mat = nist->FindOrBuildMaterial("G4_C");
  G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  Sample_mat = nist->FindOrBuildMaterial("G4_Ni");
  //
  // World
  //
  G4double world_sizeXY = 2*m;
  G4double world_sizeZ = 2*m;
  auto solidWorld =
    new G4Box("World",  // its name
              0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size
  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        world_mat,  // its material
                                        "World");  // its name
  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     checkOverlaps);  // overlaps checking
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());



   //Ge Detector Cap
  G4double numZPlanesCap = 7;
  G4double rInnerCap[] = { 0,       0*mm,       34.25*mm, 46*mm,  46*mm,    0,      0};
  G4double rOuterCap[] = { 34.25*mm,35.74*mm,   35.75*mm, 47.5*mm,47.5*mm,  47.5*mm, 47.5*mm};
  G4double zPlaneCap[] = {  0,      0.759*mm,   0.76*mm,  4*mm,   159.4*mm, 159.5*mm,161*mm};

  G4Polycone * Cap_sol = new G4Polycone("Cap_sol",0,2.*CLHEP::pi,numZPlanesCap,zPlaneCap,rInnerCap,rOuterCap);
  G4LogicalVolume* Cap_log = new G4LogicalVolume(Cap_sol, Cap_mat, "Cap_log");
  G4VPhysicalVolume* Cap_phys= new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, -4*mm),Cap_log, "Cap_phys",  logicWorld, false, 0,checkOverlaps);
  Cap_log->SetVisAttributes(new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.3)));

  // inner cavity shape (make sure it's entirely inside the cap!)
  G4double numZPlanesVac = 3;
  G4double rInnerVac[] = {0, 0, 0};
  G4double rOuterVac[] = {34.25*mm, 46*mm,46*mm};
  G4double zVac[]      = {0.76*mm,4*mm,159.4*mm,159.5*mm};

  G4Polycone* VacCavity_sol = new G4Polycone("VacCavity_sol", 0, 2*M_PI, numZPlanesVac, zVac, rInnerVac, rOuterVac);
  G4LogicalVolume* VacCavity_log = new G4LogicalVolume(VacCavity_sol, vacuum, "VacCavity_log");
  G4VPhysicalVolume* VacCavity_phys= new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, -4*mm),VacCavity_log, "VacCavity_phys",  logicWorld, false, 0,checkOverlaps);
  VacCavity_log->SetVisAttributes(G4VisAttributes::GetInvisible());


// Ge Detector configuration
struct DetectorConfig {
  G4double GeLen;
  G4double GeDia;
  G4double HoleLen;
  G4double HoleDia;
  G4double FaceCurve;
};
DetectorConfig config;

if (DetectorName == "Mordred") {
  config = {68.7 * mm, 68.5 * mm, 62.7 * mm, 10.3 * mm, 8.0 * mm};
} else if (DetectorName == "Maeve") {
  config = {69.8 * mm, 77.6 * mm, 56.1 * mm, 9.4 * mm, 8.0 * mm};
} else if (DetectorName == "Morgan") {
  config = {95.3 * mm, 73.0 * mm, 86.3 * mm, 10.8 * mm, 8.0 * mm};
} else {
  G4Exception("Invalid detector name", "", FatalException, "Unknown DetectorName.");
}

G4double GeLen     = config.GeLen;
G4double GeDia     = config.GeDia;
G4double HoleLen   = config.HoleLen;
G4double HoleDia   = config.HoleDia;
G4double FaceCurve = config.FaceCurve;

// Build angle list and transition point
std::vector<G4double> angles_deg = {0, 22.5, 45, 67.5, 90};

G4double z = GeLen - HoleLen;
G4double phi_deg_transition = 90.0;

if (z <= FaceCurve) {
  G4double arg = (FaceCurve - z) / FaceCurve;
  arg = std::clamp(arg, -1.0, 1.0);  // safe domain for acos
  G4double phi_rad = std::acos(arg);
  phi_deg_transition = phi_rad * 180.0 / M_PI;
}

// Always insert rInner transition
angles_deg.push_back(phi_deg_transition);
angles_deg.push_back(phi_deg_transition + 0.001);
std::sort(angles_deg.begin(), angles_deg.end());

// Build polycone points
std::vector<G4double> rInner;
std::vector<G4double> rOuter;
std::vector<G4double> zPlane;

for (G4double angle_deg : angles_deg) {
  G4double phi = angle_deg * M_PI / 180.0;

  if (angle_deg <= 90.0) {
    rOuter.push_back(GeDia / 2.0 - FaceCurve + FaceCurve * std::sin(phi));
  }
  else {
    rOuter.push_back(GeDia / 2.0);
  }

  zPlane.push_back(FaceCurve - FaceCurve * std::cos(phi));
  rInner.push_back(angle_deg <= phi_deg_transition ? 0.0 : HoleDia / 2.0);
}

// Final flat end
rOuter.push_back(GeDia / 2.0);
zPlane.push_back(GeLen);
rInner.push_back(HoleDia / 2.0);

G4int numZPlanes = rOuter.size();

// Create Ge detector solid and volume
G4Polycone* Ge_sol = new G4Polycone("Ge_sol", 0, 2. * CLHEP::pi,
                                    numZPlanes, zPlane.data(),
                                    rInner.data(), rOuter.data());

G4LogicalVolume* Ge_log = new G4LogicalVolume(Ge_sol, Ge_mat, "Ge_log");
G4VPhysicalVolume* Ge_phys = new G4PVPlacement(
    nullptr, G4ThreeVector(0, 0, +4.0 * mm), Ge_log, "Ge_phys",
    VacCavity_log, false, 0, checkOverlaps);
Ge_log->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0, 0)));

// Optional: Debug printout
G4cout << "============ " << G4endl;
G4cout << "Detector : " << DetectorName << G4endl;
G4cout << "rInner[] = { ";
for (size_t i = 0; i < rInner.size(); ++i) {
  G4cout << rInner[i] / mm << " mm";
  if (i < rInner.size() - 1) G4cout << ", ";
}
G4cout << " }" << G4endl;

G4cout << "rOuter[] = { ";
for (size_t i = 0; i < rOuter.size(); ++i) {
  G4cout << rOuter[i] / mm << " mm";
  if (i < rOuter.size() - 1) G4cout << ", ";
}
G4cout << " }" << G4endl;

G4cout << "zPlane[] = { ";
for (size_t i = 0; i < zPlane.size(); ++i) {
  G4cout << zPlane[i] / mm << " mm";
  if (i < zPlane.size() - 1) G4cout << ", ";
}
G4cout << " }" << G4endl;
G4cout << "Volume : " << Ge_sol->GetCubicVolume() / cm3 << " cm³" << G4endl;
G4cout << "============ " << G4endl;

fScoringVolume = Ge_log;

  //Ge Beaker
  G4double numZPlanesBeaker = 6;
  G4double rInnerBeaker[] = {75.5*mm,  75.5*mm, 0,        0,        48.5*mm,  48.5*mm};
  G4double rOuterBeaker[] = {77.5*mm,  77.5*mm, 77.5*mm,  77.5*mm,  77.5*mm,  77.5*mm};
  G4double zPlaneBeaker[] = {0,        66*mm,   66.01*mm, 67.99*mm, 68*mm,    165*mm};

  G4Polycone* Beaker_sol = new G4Polycone("Beaker_sol",0,2.*CLHEP::pi,numZPlanesBeaker,zPlaneBeaker,rInnerBeaker,rOuterBeaker);
  G4Tubs* Beaker_sol2 = new G4Tubs("Beaker2_sol", 50.49*mm, 75.51*mm, 48.51*mm, 0.*deg, 360.*deg);
  G4SubtractionSolid* Beaker_sol3 = new G4SubtractionSolid("Beaker_sol3", Beaker_sol, Beaker_sol2, 0, G4ThreeVector(0,0,114.5*mm));
  G4LogicalVolume* Beaker_log = new G4LogicalVolume(Beaker_sol3, Beaker_mat, "Beaker_log");
  G4VPhysicalVolume* Beaker_phys= new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, -72.5*mm),
                                                    Beaker_log,"Beaker_phys",  logicWorld, false, 0,checkOverlaps);
  Beaker_log->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0.5,0.9)));
  G4Polyhedron* test = Beaker_sol3->GetPolyhedron();



  G4Tubs* Sample_sol = new G4Tubs("Sample_sol", 50.5*mm, 75.5 * mm, Sample_height/2, 0., 360. * deg);
  Sample_log = new G4LogicalVolume(Sample_sol, Sample_mat, "Sample_log");
  Sample_phys = new G4PVPlacement(0, G4ThreeVector(0,0,165*mm-72.5*mm-2*mm-Sample_height/2), Sample_log, "Sample_phys",
              logicWorld,false, 0, checkOverlaps);
  Sample_log->SetVisAttributes(new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.9)));





  //
  // always return the physical World
  //
  return physWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& name)
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* newMaterial = nullptr;
    // Map of predefined materials
    static const std::map<G4String, G4String> materialMap = {
        {"Ti", "G4_Ti"},
        {"Nb", "G4_Nb"},
        {"Cr", "G4_Cr"}
    };

    if (name == "Sand") {
        // Define custom "Sand" material
        G4double density = 2.6 * g/cm3; // 2.6 normal, this gravel so assume 2 ?
        newMaterial = new G4Material("Sand", density, 14);

        // Add elements by mass fraction
        newMaterial->AddElement(nist->FindOrBuildElement("O"), 0.526);
        newMaterial->AddElement(nist->FindOrBuildElement("Si"), 0.217);
        newMaterial->AddElement(nist->FindOrBuildElement("Al"), 0.100);
        newMaterial->AddElement(nist->FindOrBuildElement("Fe"), 0.048);
        newMaterial->AddElement(nist->FindOrBuildElement("Ca"), 0.040);
        newMaterial->AddElement(nist->FindOrBuildElement("K"), 0.024);
        newMaterial->AddElement(nist->FindOrBuildElement("Na"), 0.020);
        newMaterial->AddElement(nist->FindOrBuildElement("Mg"), 0.018);
        newMaterial->AddElement(nist->FindOrBuildElement("P"), 0.002);
        newMaterial->AddElement(nist->FindOrBuildElement("Ba"), 0.001);
        newMaterial->AddElement(nist->FindOrBuildElement("Ti"), 0.001);
        newMaterial->AddElement(nist->FindOrBuildElement("Sr"), 0.001);
        newMaterial->AddElement(nist->FindOrBuildElement("Mn"), 0.001);
        newMaterial->AddElement(nist->FindOrBuildElement("Zn"), 0.001); // Representing trace elements

        // Use newMaterial as needed
    } else {
        auto it = materialMap.find(name);
        if (it != materialMap.end()) {
            newMaterial = nist->FindOrBuildMaterial(it->second);
            // Use newMaterial as needed
        } else {
            G4cerr << "ERROR: Material '" << name << "' is not supported." << G4endl;
            return;
        }
    }

    Sample_log->SetMaterial(newMaterial);
    Sample_mat = newMaterial;
    CalculateAndSetDensity();
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << Sample_log->GetMaterial() << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSampleHeight(G4double height) {
    if (height <= 0.0) {
        G4cerr << "ERROR: Sample height must be positive!" << G4endl;
        return;
    }
    Sample_height = height;

    // Recreate the solid with new height
    G4Tubs* solid = dynamic_cast<G4Tubs*>(Sample_log->GetSolid());
    if (!solid) {
        G4cerr << "ERROR: Solid is not a G4Tubs!" << G4endl;
        return;
    }
    // Update half-height
    solid->SetZHalfLength(Sample_height/2);

     // Update placement Z to keep the same base
    G4double newZ = 165*mm - 72.5*mm - 2*mm - height / 2.0;
    Sample_phys->SetTranslation(G4ThreeVector(0, 0, newZ));
    CalculateAndSetDensity();

    // Notify run manager
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSampleMass(G4double mass) {
    if (mass <= 0.0) {
        G4cerr << "ERROR: Sample mass must be positive!" << G4endl;
        return;
    }
    Sample_mass = mass;
    CalculateAndSetDensity();

    // Notify run manager
    //G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintSampleDensity() {

    G4Material* mat = Sample_log->GetMaterial();
    G4cout << "Sample Material: " << mat->GetName()
           << " | Density: " << mat->GetDensity() / (g/cm3) << " g/cm³" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CalculateAndSetDensity() {
  G4double volume = Sample_log->GetSolid()->GetCubicVolume();  // in mm³
  G4cout << "Sample geometric volume: " << volume / cm3 << " cm³" << G4endl;
  Sample_density = Sample_mass/volume;
  PrintSampleDensity();
  G4cout << "Desired Sample density: " << Sample_density / (g/cm3) << " g/cm³" << G4endl;

  G4Material* newMaterial = CreateCustomDensityMaterial(Sample_mat->GetName(), Sample_mat->GetName(), Sample_density / (g/cm3) * (g/cm3) );
  Sample_log->SetMaterial(newMaterial);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4Material* DetectorConstruction::CreateCustomDensityMaterial(const G4String& baseName,
                                                               const G4String& newName,
                                                               G4double newDensity) {
    auto* baseMat = G4NistManager::Instance()->FindOrBuildMaterial(baseName);
    if (!baseMat) return nullptr;
    G4cout << baseMat << G4endl;

    auto* newMat = new G4Material(newName, newDensity, baseMat->GetNumberOfElements());
    for (size_t i = 0; i < baseMat->GetNumberOfElements(); ++i) {
        auto* elem = const_cast<G4Element*>(baseMat->GetElement(i));
        newMat->AddElement(elem, baseMat->GetFractionVector()[i]);
    }
    return newMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetector(const G4String& name)
{
  if (name == "Mordred" || name == "Morgan" || name == "Maeve") {
    DetectorName = name;
    G4cout << "Detector changed to: " << DetectorName << G4endl;

    // Reinitialize geometry to apply change
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  } else {
    G4cerr << "ERROR: Unknown detector name '" << name << "'" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/detector/", "Detector control");

  // Material command
  auto& materialCmd = fMessenger->DeclareMethod(
    "setMaterial", &DetectorConstruction::SetMaterial, "Set the material of the sample.");
  materialCmd.SetParameterName("material", false);  // Mandatory parameter
  materialCmd.SetGuidance("Choose material: Ti, Nb, Cr or Sand.");

  // Material command
  auto& detectorCmd = fMessenger->DeclareMethod(
      "setDetector", &DetectorConstruction::SetDetector, "Set the detector.");
  detectorCmd.SetParameterName("Detector", false);  // Mandatory parameter
  detectorCmd.SetGuidance("Choose Mordred Morgan Maeve");

  auto& heightCmd = fMessenger->DeclareMethodWithUnit(
      "setSampleHeight", "mm", &DetectorConstruction::SetSampleHeight,
      "Set the half-height (Z extent) of the sample volume.");
  heightCmd.SetParameterName("height", false);  // required
  heightCmd.SetRange("height > 0.");

  auto& massCmd = fMessenger->DeclareMethodWithUnit(
      "setSampleMass", "kg", &DetectorConstruction::SetSampleMass,
      "Set the mass of the sample volume.");
  massCmd.SetParameterName("mass", false);  // required
  massCmd.SetRange("mass > 0.");


  fMessenger->DeclareMethod("printSampleDensity", &DetectorConstruction::PrintSampleDensity,
                          "Print the current sample material and density.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
