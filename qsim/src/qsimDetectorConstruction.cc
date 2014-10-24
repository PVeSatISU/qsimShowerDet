
#include "qsimDetectorConstruction.hh"

#include "qsimDetector.hh"
#include "qsimScintDetector.hh"
#include "G4SDManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4GenericTrap.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4TwoVector.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::qsimDetectorConstruction()
{
  det_x = det_y = det_z = 275*cm;
  quartz_x = 1.75*cm; 
  quartz_y = 7.*cm;  //2.5 
//Half cm
    quartz_z = 0.5*cm;
//One cm
//  quartz_z = 0.5*cm;

	quartz_zPos = -.0*cm;//-1.1*cm; //-.9*cm; //-.6*cm;

  cone_rmin1 = 2.1*cm;
  cone_rmax1 = cone_rmin1+.05*cm;
  cone_rmin2 = 2.5*cm;  // normally 2.5*cm;
  cone_rmax2 = cone_rmin2+.05*cm;
  cone_z = quartz_y+.5*cm;    //3
  cone_sphi = 0.;
  cone_fphi = 2*3.1415;

  rin = cone_rmin2;  // normally 2.5*cm;
  rout = rin+.05*cm;
  lngth = 1.9*cm;  // PMT dist. = 2*lngth +1cm  (10.4 == 4.5, 6.8 == 2.9)


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::~qsimDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* qsimDetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Air
// 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

// Quartz
// 
  G4Element* Si = new G4Element("Silicon", "Si", z=14 , a=28*g/mole);

  G4Material* Quartz = new G4Material("Quartz", density= 2.203*g/cm3, nelements=2);
  Quartz->AddElement(Si, 1);
  Quartz->AddElement(O, 2);

//Tungsten
//
  G4Element* W = new G4Element("Tungsten", "W", z=74 , a=183.84*g/mole);

  G4Material* Tungsten = new G4Material("Tungsten", density= 19.25*g/cm3, nelements=1);
  Tungsten->AddElement(W, 1);

// Mirror
// 
  G4Element* Al = new G4Element("Aluminum", "Al", z=13 , a=27*g/mole);

  G4Element* Pb = new G4Element("Lead", "Pb", z=82 , a=207.2*g/mole);

  G4Material* Alu_Mat = new G4Material("Alu_Mat", 2.7*g/cm3, nelements=1);
  Alu_Mat->AddElement(Al, 1);

  G4Material* Pb_Mat = new G4Material("Pb_Mat", 11.34*g/cm3, nelements=1);
  Pb_Mat->AddElement(Pb, 1);

	//G4Material* Pb_Mat=Air; // To remove lead bricks, uncomment.
	
  G4Material* Mirror = new G4Material("Mirror", density= 2.7*g/cm3, nelements=1);
  Mirror->AddElement(Al, 1);


// Let us make cathode from a special metal (reflectivity 0, efficiency of photoelectrons 25%)
  G4Material* CATH = new G4Material("CATH", density= 2.7*g/cm3, nelements=1);
  CATH->AddElement(Al, 1);


//
// ------------ Generate & Add Material Properties Table ------------
//



const G4int nEntries = 190;

	G4double PhotonEnergy[nEntries] =
		{  2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,
		2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,
		2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,
		3,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,
		3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,
		3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,
		3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,
		3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,
		4,4.02,4.04,4.06,4.08,4.1,4.12 ,4.14,4.16,4.18,  //Glass cuts off above 4.135eV, 87 entries
		4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,
		4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,
		4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,
		4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98, //  Cut off -> 4.96eV ~ 250nm
		5,5.02,5.04   ,   5.06,5.08,5.1,5.12,5.14,5.16,5.18,   // 5.04eV = 246 nm is the 30% cutoff, 133 entries
		5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,
		5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,	
		5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,
		5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,
		6,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18   };  // 200 nm

	G4double RefractiveIndex1[nEntries];
	G4double Absorption1[nEntries];
	G4double RefractiveIndex2[nEntries];
	G4double RefractiveIndex3[nEntries];
	G4double Reflectivity4[nEntries];
	G4double Efficiency4[nEntries];
	G4double Reflectivity3[nEntries];

 	for (int i = 0; i < nEntries; i++) {
		PhotonEnergy[i] = PhotonEnergy[i]*eV;
		RefractiveIndex1[i]= 1.455 -(.005836*PhotonEnergy[i])+(.003374*PhotonEnergy[i]*PhotonEnergy[i]);

//Aluminum
//		Reflectivity3[i] = 0; //.6;
	
//Aluminum Real
   if (PhotonEnergy[i] < 4.135) Reflectivity3[i] = .75;  // regularly .75, .7 below  .56/.53/.46 tunes to 50 PEs
		else if (PhotonEnergy[i] >= 4.135 && PhotonEnergy[i] < 6.203) Reflectivity3[i] = .7;
   else Reflectivity3[i] = .6;		// .6
		
//ALZAK		
//		if (PhotonEnergy[i] < 3.26*eV) {
//			Reflectivity3[i]=.93; }
//		else { Reflectivity3[i] = 0;}

// No Mirror
//		Reflectivity3[i] = 0;
		
//		Absorption1[i] = 50.*cm;  //Uniform
   
        Absorption1[i] = (exp(4.325)*exp(1.191*PhotonEnergy[i])*exp(-.213*PhotonEnergy[i]*PhotonEnergy[i])*exp(-.04086*PhotonEnergy[i]*PhotonEnergy[i]*PhotonEnergy[i]))*m;

       if (Absorption1[i] > 25*m) {Absorption1[i] = 25*m;}

		RefractiveIndex2[i]=1;
		RefractiveIndex3[i]=0;
		Reflectivity4[i]=0;
		Efficiency4[i]=.25;

		
	}
	
//QUARTZ
	
  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  
  Quartz->SetMaterialPropertiesTable(myMPT1);

//
// Air
//

/*  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };  */

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  
  Air->SetMaterialPropertiesTable(myMPT2);

//
// Mirror (refractive index = 0) 
//


/*  G4double RefractiveIndex3[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };   */


 

  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex3, nEntries);
  // myMPT3->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption3, nEntries);  
  //  myMPT3->AddProperty("REFLECTIVITY", PhotonEnergy, Reflectivity3, nEntries);
  //  myMPT3->AddProperty("EFFICIENCY",    PhotonEnergy, Efficiency3, nEntries);  

  Mirror->SetMaterialPropertiesTable(myMPT3);

//
// CATH
//

/*
 G4double Reflectivity4[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };

  G4double Efficiency4[nEntries] =
           {0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25, 0.25 };  */



 
  G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();
  myMPT4->AddProperty("REFLECTIVITY",       PhotonEnergy, Reflectivity4,nEntries);
  myMPT4->AddProperty("EFFICIENCY",    PhotonEnergy, Efficiency4, nEntries);  

  CATH->SetMaterialPropertiesTable(myMPT4);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

//
//	------------- Volumes --------------

// The detector
//
  G4Box* det_box = new G4Box("World",det_x,det_y,det_z);

  G4LogicalVolume* det_log
    = new G4LogicalVolume(det_box,Air,"World",0,0,0);

  G4VPhysicalVolume* det_phys
    = new G4PVPlacement(0,G4ThreeVector(),det_log,"World",0,false,0);


/*
// Al plate
	G4Box* Al_box = new G4Box("Al_box", 20*cm, 20*cm, .0176*cm);
	
  G4LogicalVolume* Al_log
    = new G4LogicalVolume(Al_box,Alu_Mat,"Al_log",0,0,0);
	
	G4VPhysicalVolume* Al_phys = new G4PVPlacement(0, G4ThreeVector(0,0,-9.0*cm), Al_log, "Aluminum",
									det_log, false, 0);
*/

// The quartz

  G4double Qthickness = 5*mm;
  G4double Wthickness = 2.4*mm;
  G4double Mthickness = 0.01*mm;

  G4Trap* tungsten_box_1 = new G4Trap("tungsten_box_1",245.195*mm/2,210.493*mm/2,Wthickness/2,Wthickness/2,79.497*mm);

  G4Box* tungsten_box_2 = new G4Box("tungsten_box_2",245.195*mm,Wthickness/2,Qthickness/2);

  G4ThreeVector zTrans_1(0, 0.0*mm,-79.497*mm - Qthickness/2); 

  G4RotationMatrix* rotW_1 = new G4RotationMatrix;
     rotW_1->rotateX(0*deg);

  G4UnionSolid* union_tungsten_box_1_tungsten_box_2 =
      new G4UnionSolid("union_tungsten_box_1_tungsten_box_2", tungsten_box_1, tungsten_box_2, rotW_1, zTrans_1);

 G4LogicalVolume* union_tungsten_box_1_tungsten_box_2_log
    = new G4LogicalVolume(union_tungsten_box_1_tungsten_box_2,Tungsten,"union_tungsten_box_1_tungsten_box_2_log",0,0,0);

 //G4VPhysicalVolume* union_tungsten_box_1_tungsten_box_2_phys
    //= new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,-Wthickness/2-Qthickness/2,0.0*mm),union_tungsten_box_1_tungsten_box_2_log,"Tungsten_1",
                        //det_log,false,0);  // normally zero vector

  G4Trap* quartz_box_1 = new G4Trap("quartz_box_1",245.195*mm/2,210.493*mm/2,Qthickness/2,Qthickness/2,79.497*mm);

  //G4Trap* quartz_box_1_single_piece = new G4Trap("quartz_box_1_single_piece", 210.493*mm, Qthickness, 2*79.497*mm + Qthickness, 2*79.497*mm);

  G4int npoints = 8;
    std::vector<G4TwoVector> points_quartz_1(npoints);
    points_quartz_1[0] = G4TwoVector(0*mm,  0*mm);
    points_quartz_1[1] = G4TwoVector(0*mm,  Qthickness);
    points_quartz_1[2] = G4TwoVector(Qthickness, 0*mm);
    points_quartz_1[3] = G4TwoVector(Qthickness, 0*mm);
    points_quartz_1[4] = G4TwoVector(0*mm,  0*mm);
    points_quartz_1[5] = G4TwoVector(0*mm,  Qthickness);
    points_quartz_1[6] = G4TwoVector(Qthickness, 0*mm);
    points_quartz_1[7] = G4TwoVector(Qthickness, 0*mm);

  G4GenericTrap* quartz_box_2 = new G4GenericTrap("quartz_box_2",123.19*mm,points_quartz_1);
    
  G4RotationMatrix* rotQ_1 = new G4RotationMatrix;
     rotQ_1->rotateX(0*deg);
     rotQ_1->rotateY(-90*deg);

 G4ThreeVector zTrans_2(0, -Qthickness/2, -79.497*mm); 

 G4UnionSolid* union_quartz_box_1_quartz_box_2 =
      new G4UnionSolid("union_quartz_box_1_quartz_box_2", quartz_box_1, quartz_box_2, rotQ_1, zTrans_2);

 G4LogicalVolume* union_quartz_box_1_quartz_box_2_log
    = new G4LogicalVolume(union_quartz_box_1_quartz_box_2,Quartz,"union_quartz_box_1_quartz_box_2_log",0,0,0);

 qsimScintDetector* quartzSD = new qsimScintDetector("QuartzSD", 10);
 SDman->AddNewDetector(quartzSD);
 union_quartz_box_1_quartz_box_2_log->SetSensitiveDetector(quartzSD);

 G4RotationMatrix* rotQ_2 = new G4RotationMatrix;
    rotQ_2->rotateX(0*deg);
    //rotQ_2->rotateY(-90*deg);

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_1
    = new G4PVPlacement(rotQ_2,G4ThreeVector(0.0*cm,0,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_1",
                        det_log,false,0);  // normally zero vector

 G4LogicalVolume* tungsten_box_3_log
    = new G4LogicalVolume(tungsten_box_1,Tungsten,"tungsten_box_3_log",0,0,0);

G4VPhysicalVolume* union_tungsten_box_1_tungsten_box_2_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,-Wthickness/2 - Qthickness/2 - Mthickness,0.0*mm),tungsten_box_3_log,"Tungsten_1",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* tungsten_box_3_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,Qthickness/2 + Wthickness/2 + Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_3_phys",
                        det_log,false,0);  // normally zero vector

    std::vector<G4TwoVector> points_quartz_2(npoints);
    points_quartz_2[0] = G4TwoVector(-Qthickness/2,  0*mm);
    points_quartz_2[1] = G4TwoVector(0*mm,  Qthickness/2);
    points_quartz_2[2] = G4TwoVector(Qthickness/2, 0*mm);
    points_quartz_2[3] = G4TwoVector(Qthickness/2, 0*mm);
    points_quartz_2[4] = G4TwoVector(-Qthickness/2,  0*mm);
    points_quartz_2[5] = G4TwoVector(0*mm,  Qthickness/2);
    points_quartz_2[6] = G4TwoVector(Qthickness/2, 0*mm);
    points_quartz_2[7] = G4TwoVector(Qthickness/2, 0*mm);

 G4GenericTrap* quartz_box_3 = new G4GenericTrap("quartz_box_3",123.19*mm,points_quartz_2);

 G4RotationMatrix* rotQ_3 = new G4RotationMatrix;
  rotQ_3->rotateX(90*deg);  
  rotQ_3->rotateY(-90*deg);
    

 G4ThreeVector zTrans_3(0, 0.0*mm, -79.497*mm); 

 G4UnionSolid* union_quartz_box_1_quartz_box_3 =
      new G4UnionSolid("union_quartz_box_1_quartz_box_3", quartz_box_1, quartz_box_3, rotQ_3, zTrans_3);

 G4LogicalVolume* union_quartz_box_1_quartz_box_3_log
    = new G4LogicalVolume(union_quartz_box_1_quartz_box_3,Quartz,"union_quartz_box_1_quartz_box_3_log",0,0,0);

 union_quartz_box_1_quartz_box_3_log->SetSensitiveDetector(quartzSD);

 G4RotationMatrix* rotQ_4 = new G4RotationMatrix;
    //rotQ_4->rotateY(-90*deg);
    rotQ_4->rotateX(0*deg);

 G4RotationMatrix* rotQ_5 = new G4RotationMatrix;
    //rotQ_5->rotateY(90*deg);
    rotQ_5->rotateZ(180*deg);

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_2
    = new G4PVPlacement(rotQ_5,G4ThreeVector(0.0*cm,Qthickness + Wthickness + 2*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_2",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* tungsten_box_4_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,3*Qthickness/2 + 3*Wthickness/2 + 3*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_4_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_3
    = new G4PVPlacement(rotQ_4,G4ThreeVector(0.0*cm,2*Qthickness + 2*Wthickness + 4*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_3",
                        det_log,false,0);  // normally zero vector
 
 G4VPhysicalVolume* tungsten_box_5_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,5*Qthickness/2 + 5*Wthickness/2 + 5*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_5_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_4
    = new G4PVPlacement(rotQ_5,G4ThreeVector(0.0*cm,3*Qthickness + 3*Wthickness + 6*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_4",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* tungsten_box_6_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,7*Qthickness/2 + 7*Wthickness/2 + 7*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_6_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_5
    = new G4PVPlacement(rotQ_4,G4ThreeVector(0.0*cm,4*Qthickness + 4*Wthickness + 8*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_5",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* tungsten_box_7_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,9*Qthickness/2 + 9*Wthickness/2 + 9*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_7_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_6
    = new G4PVPlacement(rotQ_5,G4ThreeVector(0.0*cm,5*Qthickness + 5*Wthickness + 10*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_6",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* tungsten_box_8_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,11*Qthickness/2 + 11*Wthickness/2 + 11*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_8_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_7
    = new G4PVPlacement(rotQ_4,G4ThreeVector(0.0*cm,6*Qthickness + 6*Wthickness + 12*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_7",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* tungsten_box_9_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,13*Qthickness/2 + 13*Wthickness/2 + 13*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_9_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_8
    = new G4PVPlacement(rotQ_5,G4ThreeVector(0.0*cm,7*Qthickness + 7*Wthickness + 14*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_8",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* tungsten_box_10_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,15*Qthickness/2 + 15*Wthickness/2 + 15*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_10_phys",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_9
    = new G4PVPlacement(rotQ_4,G4ThreeVector(0.0*cm,8*Qthickness + 8*Wthickness + 16*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_9",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* tungsten_box_11_phys
    = new G4PVPlacement(rotW_1,G4ThreeVector(0.0*cm,17*Qthickness/2 + 17*Wthickness/2 + 17*Mthickness,0.0*mm),tungsten_box_3_log,"tungsten_box_11_phys",
                        det_log,false,0);  // normally zero vector

G4VPhysicalVolume* union_quartz_box_1_quartz_box_2_phys_10
    = new G4PVPlacement(rotQ_5,G4ThreeVector(0.0*cm,9*Qthickness + 9*Wthickness + 18*Mthickness,0.0*mm),union_quartz_box_1_quartz_box_2_log,"union_quartz_box_1_quartz_box_2_phys_10",
                        det_log,false,0);  // normally zero vector

//Mirrors

 G4Trap* mirror_box_1 = new G4Trap("mirror_box_1",70.0*mm/2,89.54*mm/2,0.5*mm/2,0.5*mm/2,58.3*mm/2);

 G4ThreeVector zTrans_4(-49.48*mm,33.30*mm,-304.80*mm);

 G4RotationMatrix* rotM_1 = new G4RotationMatrix;
     rotM_1->rotateY(32.08*deg);
     rotM_1->rotateZ(90*deg);

 G4LogicalVolume* mirror_box_1_log
    = new G4LogicalVolume(mirror_box_1,Mirror,"mirror_box_1_log",0,0,0);

 G4VPhysicalVolume* mirror_box_1_phys
    = new G4PVPlacement(rotM_1,zTrans_4,mirror_box_1_log,"mirror_box_1_phys",
                        det_log,false,0);  // normally zero vector

 G4Trap* mirror_box_1_1 = new G4Trap("mirror_box_1_1",89.54*mm/2,128.0*mm/2,0.5*mm/2,0.5*mm/2,100.91*mm/2);

 G4ThreeVector zTrans_4_1(-78.89*mm,33.30*mm,-231.61*mm);

 G4RotationMatrix* rotM_1_1 = new G4RotationMatrix;
     rotM_1_1->rotateY(16.03*deg);
     rotM_1_1->rotateZ(90*deg);

 G4LogicalVolume* mirror_box_1_1_log
    = new G4LogicalVolume(mirror_box_1_1,Mirror,"mirror_box_1_1_log",0,0,0);

 G4VPhysicalVolume* mirror_box_1_1_phys
    = new G4PVPlacement(rotM_1_1,zTrans_4_1,mirror_box_1_1_log,"mirror_box_1_1_phys",
                        det_log,false,0);  // normally zero vector


 G4Trap* mirror_box_2 = new G4Trap("mirror_box_2",71.6*mm/2,128/2*mm,0.5*mm/2,0.5*mm/2,107.81*mm/2);

 G4ThreeVector zTrans_5(-107.71*mm,33.30*mm,-131.31*mm);

 G4RotationMatrix* rotM_2 = new G4RotationMatrix;
     rotM_2->rotateY(16.03*deg+180*deg);
     rotM_2->rotateZ(90*deg);

 G4LogicalVolume* mirror_box_2_log
    = new G4LogicalVolume(mirror_box_2,Mirror,"mirror_box_2_log",0,0,0);

 G4VPhysicalVolume* mirror_box_2_phys
    = new G4PVPlacement(rotM_2,zTrans_5,mirror_box_2_log,"mirror_box_2_phys",
                        det_log,false,0);  // normally zero vector

/////////////////////////////////////

 G4ThreeVector zTrans_6(49.48*mm,33.30*mm,-304.80*mm);

 G4RotationMatrix* rotM_3 = new G4RotationMatrix;
     rotM_3->rotateY(-32.08*deg);
     rotM_3->rotateZ(90*deg);

 G4VPhysicalVolume* mirror_box_3_phys
    = new G4PVPlacement(rotM_3,zTrans_6,mirror_box_1_log,"mirror_box_3_phys",
                        det_log,false,0);  // normally zero vector

G4ThreeVector zTrans_6_1(78.90*mm,33.30*mm,-231.61*mm);

 G4RotationMatrix* rotM_3_1 = new G4RotationMatrix;
     rotM_3_1->rotateY(-16.03*deg);
     rotM_3_1->rotateZ(90*deg);

 G4VPhysicalVolume* mirror_box_3_1_phys
    = new G4PVPlacement(rotM_3_1,zTrans_6_1,mirror_box_1_1_log,"mirror_box_3_1_phys",
                        det_log,false,0);  // normally zero vector

 G4ThreeVector zTrans_7(107.71*mm,33.30*mm,-131.31*mm);

 G4RotationMatrix* rotM_2_1 = new G4RotationMatrix;
     rotM_2_1->rotateY(-16.03*deg+180*deg);
     rotM_2_1->rotateZ(90*deg);

 G4VPhysicalVolume* mirror_box_4_phys
    = new G4PVPlacement(rotM_2_1,zTrans_7,mirror_box_2_log,"mirror_box_4_phys",
                        det_log,false,0);  // normally zero vector

/////////////////////////////////////////

 G4Trap* mirror_box_3 = new G4Trap("mirror_box_3",185.66*mm/2,245.20*mm/2,0.5*mm/2,0.5*mm/2,107.39*mm/2);

 G4ThreeVector zTrans_8(0.0*mm,-16.60*mm,-131.31*mm);

 G4RotationMatrix* rotM_4 = new G4RotationMatrix;
     rotM_4->rotateX(15.22*deg);
     rotM_4->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_3_log
    = new G4LogicalVolume(mirror_box_3,Mirror,"mirror_box_3_log",0,0,0);

 G4VPhysicalVolume* mirror_box_5_phys
    = new G4PVPlacement(rotM_4,zTrans_8,mirror_box_3_log,"mirror_box_5_phys",
                        det_log,false,0);  // normally zero vector

 G4ThreeVector zTrans_9(0.0*mm,83.20*mm,-131.31*mm);

 G4RotationMatrix* rotM_5 = new G4RotationMatrix;
     rotM_5->rotateX(-15.22*deg);
     rotM_5->rotateZ(0.0*deg);

 G4VPhysicalVolume* mirror_box_6_phys
    = new G4PVPlacement(rotM_5,zTrans_9,mirror_box_3_log,"mirror_box_5_phys",
                        det_log,false,0);  // normally zero vector

 G4Trap* mirror_box_4 = new G4Trap("mirror_box_4",129.92*mm/2,185.66*mm/2,0.5*mm/2,0.5*mm/2,98.87*mm/2);

 G4ThreeVector zTrans_10(0.0*mm,-21.08*mm,-231.61*mm);

 G4RotationMatrix* rotM_6 = new G4RotationMatrix;
     rotM_6->rotateX(-10.32*deg);
     rotM_6->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_4_log
    = new G4LogicalVolume(mirror_box_4,Mirror,"mirror_box_4_log",0,0,0);

 G4VPhysicalVolume* mirror_box_7_phys
    = new G4PVPlacement(rotM_6,zTrans_10,mirror_box_4_log,"mirror_box_7_phys",
                        det_log,false,0);  // normally zero vector

 G4ThreeVector zTrans_11(0.0*mm,87.69*mm,-231.61*mm);

 G4RotationMatrix* rotM_7 = new G4RotationMatrix;
     rotM_7->rotateX(10.32*deg);
     rotM_7->rotateZ(0.0*deg);

 G4VPhysicalVolume* mirror_box_8_phys
    = new G4PVPlacement(rotM_7,zTrans_11,mirror_box_4_log,"mirror_box_8_phys",
                        det_log,false,0);  // normally zero vector

 G4Trap* mirror_box_4_1 = new G4Trap("mirror_box_4_1",68*mm/2,129.92*mm/2,0.5*mm/2,0.5*mm/2,50.36*mm/2);

 G4ThreeVector zTrans_10_1(0.0*mm,-6.58*mm,-304.80*mm);

 G4RotationMatrix* rotM_6_1 = new G4RotationMatrix;
     rotM_6_1->rotateX(-11.19*deg);
     rotM_6_1->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_4_1_log
    = new G4LogicalVolume(mirror_box_4_1,Mirror,"mirror_box_4_1_log",0,0,0);

 G4VPhysicalVolume* mirror_box_7_1_phys
    = new G4PVPlacement(rotM_6_1,zTrans_10_1,mirror_box_4_1_log,"mirror_box_7_1_phys",
                        det_log,false,0);  // normally zero vector

 G4ThreeVector zTrans_11_1(0.0*mm,73.19*mm,-304.80*mm);

 G4RotationMatrix* rotM_7_1 = new G4RotationMatrix;
     rotM_7_1->rotateX(11.19*deg);
     rotM_7_1->rotateZ(0.0*deg);

 G4VPhysicalVolume* mirror_box_8_1_phys
    = new G4PVPlacement(rotM_7_1,zTrans_11_1,mirror_box_4_1_log,"mirror_box_8_1_phys",
                        det_log,false,0);  // normally zero vector


/////////////////////////////////////

 G4Box* mirror_box_5 = new G4Box("mirror_box_5",0.5*mm/2,71.6*mm/2,159.94*mm/2);

 G4ThreeVector zTrans_12(113.92*mm+0.5*mm/2,33.0*mm,0.0*mm);

 G4RotationMatrix* rotM_8 = new G4RotationMatrix;
     rotM_8->rotateY(6.43*deg);
     rotM_8->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_5_log
    = new G4LogicalVolume(mirror_box_5,Mirror,"mirror_box_5_log",0,0,0);

 G4VPhysicalVolume* mirror_box_9_phys
    = new G4PVPlacement(rotM_8,zTrans_12,mirror_box_5_log,"mirror_box_9_phys",
                        det_log,false,0);  // normally zero vector

 G4ThreeVector zTrans_13(-113.92*mm-0.5*mm/2,33.0*mm,0.0*mm);

 G4RotationMatrix* rotM_9 = new G4RotationMatrix;
     rotM_9->rotateY(-6.43*deg);
     rotM_9->rotateZ(0.0*deg);

 G4VPhysicalVolume* mirror_box_10_phys
    = new G4PVPlacement(rotM_9,zTrans_13,mirror_box_5_log,"mirror_box_10_phys",
                        det_log,false,0);  // normally zero vector

 G4Box* mirror_box_6 = new G4Box("mirror_box_6",210.49*mm/2,71.6*mm/2,0.5*mm/2);

 G4ThreeVector zTrans_14(0.0*mm,33.3*mm,79.5*mm+0.5*mm/2);

 G4RotationMatrix* rotM_10 = new G4RotationMatrix;
     rotM_10->rotateY(0.0*deg);
     rotM_10->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_6_log
    = new G4LogicalVolume(mirror_box_6,Mirror,"mirror_box_6_log",0,0,0);

 G4VPhysicalVolume* mirror_box_11_phys
    = new G4PVPlacement(rotM_10,zTrans_14,mirror_box_6_log,"mirror_box_11_phys",
                        det_log,false,0);  // normally zero vector

 G4Trap* mirror_box_7 = new G4Trap("mirror_box_7",245.2*mm/2,210.49*mm/2,0.5*mm/2,0.5*mm/2,158.99*mm/2);

 G4ThreeVector zTrans_15(0.0*mm,69.10*mm+0.5*mm/2,0.0*mm);

 G4RotationMatrix* rotM_11 = new G4RotationMatrix;
     rotM_11->rotateY(0.0*deg);
     rotM_11->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_7_log
    = new G4LogicalVolume(mirror_box_7,Mirror,"mirror_box_7_log",0,0,0);

 G4VPhysicalVolume* mirror_box_12_phys
    = new G4PVPlacement(rotM_11,zTrans_15,mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector


/*
 G4Trap* mirror_box_7 = new G4Trap("mirror_box_7",245.195*mm/2,210.493*mm/2,Mthickness/2,Mthickness/2,79.497*mm);

 G4RotationMatrix* rotM_11 = new G4RotationMatrix;
     rotM_11->rotateY(0.0*deg);
     rotM_11->rotateZ(0.0*deg);

 G4LogicalVolume* mirror_box_7_log
    = new G4LogicalVolume(mirror_box_7,Mirror,"mirror_box_7_log",0,0,0);

 G4VPhysicalVolume* mirror_box_12_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,-Qthickness/2 - Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_13_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,Qthickness/2 + Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_14_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,Qthickness/2 + Wthickness + 3*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_15_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,3*Qthickness/2 + Wthickness + 5*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_16_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,3*Qthickness/2 + 2*Wthickness + 7*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_17_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,5*Qthickness/2 + 2*Wthickness + 9*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_18_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,5*Qthickness/2 + 3*Wthickness + 11*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_19_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,7*Qthickness/2 + 3*Wthickness + 13*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_20_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,7*Qthickness/2 + 4*Wthickness + 15*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_21_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,9*Qthickness/2 + 4*Wthickness + 17*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_22_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,9*Qthickness/2 + 5*Wthickness + 19*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_23_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,11*Qthickness/2 + 5*Wthickness + 21*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_24_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,11*Qthickness/2 + 6*Wthickness + 23*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_25_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,13*Qthickness/2 + 6*Wthickness + 25*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_26_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,13*Qthickness/2 + 7*Wthickness + 27*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_27_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,15*Qthickness/2 + 7*Wthickness + 29*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_28_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,15*Qthickness/2 + 8*Wthickness + 31*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_29_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,17*Qthickness/2 + 8*Wthickness + 33*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_30_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,17*Qthickness/2 + 9*Wthickness + 35*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector

 G4VPhysicalVolume* mirror_box_31_phys
    = new G4PVPlacement(rotM_11,G4ThreeVector(0.0*mm,19*Qthickness/2 + 9*Wthickness + 37*Mthickness/2,0.0*mm),mirror_box_7_log,"mirror_box_12_phys",
                        det_log,false,0);  // normally zero vector
*/
// Front Plate
//
  G4Box* front_plate_box = new G4Box("front_plate", 8.89*cm,19.863*cm,0.9525*cm);

  //G4LogicalVolume* front_plate_log
    //= new G4LogicalVolume(front_plate_box,Alu_Mat,"front_plate",0,0,0);

  //G4VPhysicalVolume* front_plate_phys
    //= new G4PVPlacement(0,G4ThreeVector(15.625*cm,0*cm,0*cm),front_plate_log,"front_plate",det_log,false,0);

  G4Tubs* plate_hole = new G4Tubs("plate_hole",0*cm,3.1813*cm,0.9525*cm,0*deg,360*deg);

  G4SubtractionSolid* FrontPlateWithHole
    = new G4SubtractionSolid("FrontPlateWithHole", front_plate_box, plate_hole);

  G4LogicalVolume* FrontPlateWithHole_log
    = new G4LogicalVolume(FrontPlateWithHole,Alu_Mat,"front_plate",0,0,0);

  G4RotationMatrix* rotPlate = new G4RotationMatrix;
    rotPlate->rotateY(90*deg);

  //G4VPhysicalVolume* FrontPlateWithHole_phys
    //= new G4PVPlacement(rotPlate,G4ThreeVector(15.625*cm,0*cm,0*cm),FrontPlateWithHole_log,"front_plate",det_log,false,0);

  //G4LogicalVolume* plate_hole_log
    //= new G4LogicalVolume(plate_hole,Air,"plate_hole",0,0,0);

  //G4RotationMatrix* rotPlate = new G4RotationMatrix;
    //rotPlate->rotateY(90*deg);

  //G4VPhysicalVolume* plate_hole_phys
    //= new G4PVPlacement(rotPlate,G4ThreeVector(15.625*cm,0*cm,0*cm),plate_hole_log,"plate_hole",det_log,false,0);


//quartz holder 

  G4Box* quartz_holder_bar = new G4Box("quartz_holder",7.0*cm,0.30*cm,0.30*cm);

  G4Box* quartz_holder_bar_cut = new G4Box("quartz_holder",7.0*cm,0.15*cm,0.15*cm);

  G4SubtractionSolid* QuartzHolderLeft
    = new G4SubtractionSolid("QuartzHolder", quartz_holder_bar, quartz_holder_bar_cut,0,G4ThreeVector(0*cm,-0.15*cm,0.15*cm));

  G4LogicalVolume* quartz_holder_log_left
    = new G4LogicalVolume(QuartzHolderLeft,Alu_Mat,"quartz_holder",0,0,0);

  G4SubtractionSolid* QuartzHolderRight
    = new G4SubtractionSolid("QuartzHolder", quartz_holder_bar, quartz_holder_bar_cut,0,G4ThreeVector(0*cm,0.15*cm,0.15*cm));

  G4LogicalVolume* quartz_holder_log_right
    = new G4LogicalVolume(QuartzHolderRight,Alu_Mat,"quartz_holder",0,0,0);

  //G4LogicalVolume* quartz_holder_log
    //= new G4LogicalVolume(quartz_holder_bar,Alu_Mat,"quartz_holder",0,0,0);

  //G4VPhysicalVolume* quartz_holder_phys_left
    //= new G4PVPlacement(0,G4ThreeVector(7.5*cm,1.78*cm,-0.5*cm),quartz_holder_log_left,"quartz_holder_left",det_log,false,0);

  //G4VPhysicalVolume* quartz_holder_phys_right
    //= new G4PVPlacement(0,G4ThreeVector(7.5*cm,-1.78*cm,-0.5*cm),quartz_holder_log_right,"quartz_holder_right",det_log,false,0);

	//////////////////////Air
	
/*	  G4Box* ruler_box = new G4Box("Ruler",1.5*cm,1*cm,1*cm);
	 
	 G4LogicalVolume* ruler_log
	 = new G4LogicalVolume(ruler_box,Air,"Ruler",0,0,0);
	 
	 G4VPhysicalVolume* ruler_phys
	 = new G4PVPlacement(0,G4ThreeVector(quartz_y+1.4*cm,0.5*cm,0*cm),ruler_log,"Ruler",
	 det_log,false,0);  
*/	 
	 
	
// The small mirror on the quartz
//
/*
  G4Box* mirr_Box = new G4Box("QuMirror", .05*cm, quartz_x, quartz_z*1.4142);
  
  G4LogicalVolume* mirr_log = new G4LogicalVolume(mirr_Box, Mirror, "QuMirror",0,0,0);
  
    G4RotationMatrix* rotM = new G4RotationMatrix;

	//rotM->rotateZ(M_PI*rad);
    rotM->rotateY(-M_PI/4.*rad);
  
  
  G4VPhysicalVolume* mirr_Phys = new G4PVPlacement(rotM,G4ThreeVector(-1*quartz_y+.05*cm,0,quartz_zPos+.06*cm), mirr_log, "QuMirror",
													det_log,false, 0);  // normally z= .06
*/

// Trapezoid tube plates

G4Box* frontPlate_1 = new G4Box("frontPlate_1", 0.025*cm, 2.19*cm, 1.01*cm);
   
	G4LogicalVolume* frontPlate_1_log = new G4LogicalVolume(frontPlate_1, Mirror, "FrontPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotF_1 = new G4RotationMatrix;
	rotF_1->rotateY(0*deg);//-M_PI*.4462*rad
	 //G4VPhysicalVolume* frontPlate_1_phys
		//= new G4PVPlacement(rotF_1,G4ThreeVector(-7.41*cm,0.0*cm,0.0*cm),frontPlate_1_log,"FrontPlate_1_phys",
							//det_log,false,0);  // normally zero vector

// Top

   G4Box* topPlate_1 = new G4Box("topPlate_1", 6.85*cm, 2.19*cm, 0.025*cm);
   
	G4LogicalVolume* topPlate_1_log = new G4LogicalVolume(topPlate_1, Mirror, "TopPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotT_1 = new G4RotationMatrix;
	rotT_1->rotateY(0*deg);//-M_PI*.4462*rad
	 //G4VPhysicalVolume* topPlate_1_phys
		//= new G4PVPlacement(rotT_1,G4ThreeVector(-0.56*cm,0.0*cm,-0.93*cm),topPlate_1_log,"TopPlate_1_phys",
							//det_log,false,0);  // normally zero vector

   G4Box* topPlate_2 = new G4Box("topPlate_2", 0.62*cm, 2.19*cm, 0.025*cm);
   
	G4LogicalVolume* topPlate_2_log = new G4LogicalVolume(topPlate_2, Mirror, "TopPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotT_2 = new G4RotationMatrix;
	rotT_2->rotateY(-255*deg);//-M_PI*.4462*rad
	 //G4VPhysicalVolume* topPlate_2_phys
		//= new G4PVPlacement(rotT_2,G4ThreeVector(6.45*cm,0.0*cm,-1.52*cm),topPlate_2_log,"TopPlate_2_phys",
							//det_log,false,0);  // normally zero vector

G4Box* topPlate_3 = new G4Box("topPlate_3", 1.92*cm, 2.19*cm, 0.025*cm);
   
	G4LogicalVolume* topPlate_3_log = new G4LogicalVolume(topPlate_3, Mirror, "TopPlate_3_log",0,0,0);
	
	G4RotationMatrix* rotT_3 = new G4RotationMatrix;
	rotT_3->rotateY(-45*deg);//-M_PI*.4462*rad
	 //G4VPhysicalVolume* topPlate_3_phys
		//= new G4PVPlacement(rotT_3,G4ThreeVector(7.97*cm,0.0*cm,-3.48*cm),topPlate_3_log,"TopPlate_3_phys",
							//det_log,false,0);  // normally zero vector

// Bottom 
   G4Box* botPlate_1 = new G4Box("botPlate_1", 8.07*cm, 2.19*cm, 0.025*cm);
   
	G4LogicalVolume* botPlate_1_log = new G4LogicalVolume(botPlate_1, Mirror, "botPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotB_1 = new G4RotationMatrix;
	rotB_1->rotateY(0*deg);//M_PI*.4794*rad
        rotB_1->rotateY(0*deg); 	

	 //G4VPhysicalVolume* botPlate_1_phys
		//= new G4PVPlacement(rotB_1,G4ThreeVector(0.67*cm, 0.0*cm,0.93*cm),botPlate_1_log,"botPlate_1_phys",
			//det_log,false,0);  // normally zero vector

   G4Box* botPlate_2 = new G4Box("botPlate_2", 2.31*cm, 2.19*cm, 0.025*cm);
   
	G4LogicalVolume* botPlate_2_log = new G4LogicalVolume(botPlate_2, Mirror, "botPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotB_2 = new G4RotationMatrix;
	rotB_2->rotateY(-45*deg);//M_PI*.4794*rad 	

	 //G4VPhysicalVolume* botPlate_2_phys
		//= new G4PVPlacement(rotB_2,G4ThreeVector(10.36*cm, 0.0*cm,-0.54*cm),botPlate_2_log,"botPlate_2_phys",
			//det_log,false,0);  // normally zero vector


												
    G4Trap* RPlate_1 = new G4Trap("RPlate_1", 0.05*cm, 2.02*cm, 16.14*cm, 13.70*cm);
	  
	G4LogicalVolume* RPlate_1_log = new G4LogicalVolume(RPlate_1, Mirror, "RPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotR_1 = new G4RotationMatrix;
	rotR_1->rotateY(0*deg);//-M_PI*.485*rad
	rotR_1->rotateX(90*deg);//-M_PI*.0127*rad

	 //G4VPhysicalVolume* RPlate_1_phys
		//= new G4PVPlacement(rotR_1,G4ThreeVector(0.0*cm,2.19*cm,0.0*cm),RPlate_1_log,"RPlate_1_phys",
							//det_log,false,0);  // normally zero vector

    G4Trap* RPlate_2 = new G4Trap("RPlate_2", 0.05*cm, 3.77*cm, 4.61*cm, 3.84*cm);
	  
	G4LogicalVolume* RPlate_2_log = new G4LogicalVolume(RPlate_2, Mirror, "RPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotR_2 = new G4RotationMatrix;
        rotR_2->rotateX(90*deg);
	rotR_2->rotateY(180*deg);//-M_PI*.485*rad
	rotR_2->rotateZ(45*deg);//-M_PI*.0127*rad
        
	 //G4VPhysicalVolume* RPlate_2_phys
		//= new G4PVPlacement(rotR_2,G4ThreeVector(9.05*cm,2.19*cm,-2.14*cm),RPlate_2_log,"RPlate_2_phys",
							//det_log,false,0);  // normally zero vector


    G4int nCVtx = 8;
    std::vector<G4TwoVector> cvtx(nCVtx);
    cvtx[0] = G4TwoVector(  0.0*cm,   0.0*mm);
    cvtx[1] = G4TwoVector(-0.51*cm,  1.12*cm);
    cvtx[2] = G4TwoVector( 3.17*cm,   0.0*cm);
    cvtx[3] = G4TwoVector( 3.17*cm,   0.0*cm);
    cvtx[4] = G4TwoVector(  0.0*cm,   0.0*mm);
    cvtx[5] = G4TwoVector(-0.51*cm,  1.12*cm);
    cvtx[6] = G4TwoVector( 3.17*cm,   0.0*cm);
    cvtx[7] = G4TwoVector( 3.17*cm,   0.0*cm);
     
    G4GenericTrap* RPlate_3 = new G4GenericTrap("RPlate_3",0.025*cm,cvtx);
    
    G4LogicalVolume* RPlate_3_log = new G4LogicalVolume(RPlate_3, Mirror, "RPlate_3_log",0,0,0);
    
    G4RotationMatrix* rotR_3 = new G4RotationMatrix;
        rotR_3->rotateX(90*deg);
	rotR_3->rotateY(0*deg);//-M_PI*.485*rad
	rotR_3->rotateZ(39.61*deg);//-M_PI*.0127*rad
        
	 //G4VPhysicalVolume* RPlate_3_phys
		//= new G4PVPlacement(rotR_3,G4ThreeVector(6.3*cm,2.19*cm,-1.02*cm),RPlate_3_log,"RPlate_3_phys",
							//det_log,false,0);  // normally zero vector

    G4Trap* LPlate_1 = new G4Trap("LPlate_1", 0.05*cm, 2.02*cm, 16.14*cm, 13.70*cm);
	  
	G4LogicalVolume* LPlate_1_log = new G4LogicalVolume(LPlate_1, Mirror, "LPlate_log",0,0,0);
	
	G4RotationMatrix* rotL_1 = new G4RotationMatrix;
	rotL_1->rotateY(0*deg);//-M_PI*.485*rad
	rotL_1->rotateX(90*deg);//M_PI*.0127*rad
	
	 //G4VPhysicalVolume* LPlate_1_phys
		//= new G4PVPlacement(rotL_1,G4ThreeVector(0.0*cm,-2.19*cm,0.0*cm),LPlate_1_log,"RPlate_phys",
							//det_log,false,0);  // normally zero vector

    G4Trap* LPlate_2 = new G4Trap("LPlate_2", 0.05*cm, 3.77*cm, 4.61*cm, 3.84*cm);
	  
	G4LogicalVolume* LPlate_2_log = new G4LogicalVolume(LPlate_2, Mirror, "LPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotL_2 = new G4RotationMatrix;
        rotL_2->rotateX(90*deg);
	rotL_2->rotateY(180*deg);//-M_PI*.485*rad
	rotL_2->rotateZ(45*deg);//-M_PI*.0127*rad
        
	 //G4VPhysicalVolume* LPlate_2_phys
		//= new G4PVPlacement(rotL_2,G4ThreeVector(9.05*cm,-2.19*cm,-2.14*cm),LPlate_2_log,"LPlate_2_phys",
							//det_log,false,0);  // normally zero vector

    G4GenericTrap* LPlate_3 = new G4GenericTrap("LPlate_3",0.025*cm,cvtx);
    
    G4LogicalVolume* LPlate_3_log = new G4LogicalVolume(LPlate_3, Mirror, "LPlate_3_log",0,0,0);
    
    G4RotationMatrix* rotL_3 = new G4RotationMatrix;
        rotL_3->rotateX(90*deg);
	rotL_3->rotateY(0*deg);//-M_PI*.485*rad
	rotL_3->rotateZ(39.61*deg);//-M_PI*.0127*rad
        
	 //G4VPhysicalVolume* LPlate_3_phys
		//= new G4PVPlacement(rotL_3,G4ThreeVector(6.3*cm,-2.19*cm,-1.02*cm),LPlate_3_log,"LPlate_3_phys",
							//det_log,false,0);  // normally zero vector


// The cone mirror
//
	
  G4Cons* cmirror_cone = new G4Cons("CMirror",cone_rmin1,cone_rmax1,
  cone_rmin2,cone_rmax2,cone_z,cone_sphi,cone_fphi);

  G4LogicalVolume* cmirror_log
    = new G4LogicalVolume(cmirror_cone,Mirror,"CMirror",0,0,0);

  // Rotation

  G4VPhysicalVolume* cmirror_phys;

    G4double rtphi = 0.0*deg;
    G4RotationMatrix rm;
    rm.rotateY(rtphi);
    G4double tmpx = 0.5*cm;

//    cmirror_phys
//    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(tmpx,0., 0.)),  
//                        cmirror_log,"CMirror",
//                        det_log,false,0);


// The tube mirror
//	//	//	//	//	//	//	//	//	//	

    G4double rin = 2.5*cm;
    G4double rout = 2.55*cm;
    G4double lngth = 2.3*cm;  // reg. tube
//    G4double lngth = 3.5*cm;  // long tube	
    G4double anini = 0*deg;
    G4double anspan = 360*deg;
/*


G4Cons* mirror_tube = new G4Cons("TMirror",cone_rmin2,cone_rmax2,
  2.5*cm,2.55*cm,lngth,cone_sphi,cone_fphi);

  G4LogicalVolume* tmirror_log
    = new G4LogicalVolume(mirror_tube,Mirror,"TMirror",0,0,0);


  // Rotation

  G4VPhysicalVolume* tmirror_phys;

 
    G4double tmp = lngth+cone_z+tmpx;

    tmirror_phys      // Only for tube at full distance (7.1cm)
   // = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(11.8*cm,0.,0.)),  // Long tube
    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(12.4*cm,0.,0.)),  // Long tube
						tmirror_log,"TMirror",
                        det_log,false,0);
*/
 
 /*   
  // Rotation

  G4VPhysicalVolume* cmirror_phys;

 

    rm.rotateY(rtphi);



*/
 //	//	//	//	//	//	//	//	//	//	

// The photomultiplier
//	quartz window

    G4double prin = 0;
    G4double prout = 76.2*mm/2;
    G4double plngth = 1.5*mm;
    
/*
  G4Tubs* quartz_window = new G4Tubs("QuartzWin",prin,prout,3*mm,anini,anspan);

  G4LogicalVolume* QuartzWin_log
  = new G4LogicalVolume(quartz_window,Quartz,"QuartzWin",0,0,0);

  G4VPhysicalVolume* QuartzWin_phys;

  QuartzWin_phys = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector( 15.150*cm, 0., 0.)), QuartzWin_log,"QuartzWin", det_log, false, 0);
*/

    
  //pmt
  
  G4Tubs* pmt = new G4Tubs("PMT",prin,3.0*cm,plngth,anini,anspan);

  G4LogicalVolume* pmt_log
    = new G4LogicalVolume(pmt,Air,"PMT",0,0,0);

  // Make PMT Sensitive
	
  
  G4String DetSDname = "tracker1";

  qsimDetector* trackerSD = new qsimDetector(DetSDname, 1);
  
  SDman->AddNewDetector(trackerSD);
  pmt_log->SetSensitiveDetector(trackerSD);

  // Rotation

  G4VPhysicalVolume* pmt_phys;

//  G4double sep = 0.5*cm;
//    G4double ptmp = tmp+lngth+sep;

    pmt_phys
        //= new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(10.1*cm,0.,0.)),  // Original sim. position
	      //= new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(10.7*cm,0.,0.)),  // Cosmic test sim. position
		  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(0.0*mm,33.30*mm,-329.5*mm)),  // Final detector-pmt length position
                        pmt_log,"PMT",
                        det_log,false,0);

 // metal cathode
/*
    G4double cin = 0;
    G4double cout = 2.6*cm;
    G4double clngth = 0.1*mm;


  G4Tubs* cath = new G4Tubs("CATH",cin,cout,clngth,anin,anspan);

  G4LogicalVolume* cath_log
    = new G4LogicalVolume(cath,CATH,"CATH",0,0,0);

  qsimDetector* cathSD = new qsimDetector("cath", 2);
  
  SDman->AddNewDetector(cathSD);
  cath_log->SetSensitiveDetector(cathSD);

  // Rotation

  G4VPhysicalVolume* cath_phys;

    G4double ctmp = 10.1*cm+plngth;

    cath_phys
    //    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(ctmp,0.,0.)),cath_log,"CATH",det_log,false,0);  // Original sim. position
	//	  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(10.85*cm,0.,0.)),cath_log,"CATH",det_log,false,0);  // Cosmic test position (3.7cm quartz-pmt)
		  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(14.25*cm,0.,0.)),cath_log,"CATH",det_log,false,0);  // Cosmic test position (7.1cm quartz-pmt)
*/	
	
 // Coincidence volumes **** NOTE: Upper scint is below the quartz (First coincidence w/ e-)
 
	G4Box* upperScint = new G4Box("upperScint",20.0*cm,3.5*cm,0.25*cm);
	G4LogicalVolume* uScint_log = new G4LogicalVolume(upperScint,Air,"upperScint",0,0,0);

	// Make sensitive
			//G4String 
	DetSDname = "tracker3";

	qsimScintDetector* upScintSD = new qsimScintDetector(DetSDname, 1);
  
	SDman->AddNewDetector(upScintSD);
	uScint_log->SetSensitiveDetector(upScintSD);
	
	G4double scintAngle = 0.0*deg;

	G4RotationMatrix* scintRoll = new G4RotationMatrix;
	scintRoll->rotateY(scintAngle);
		
	G4double upScint_pos=(-1*quartz_z)+(-2.5*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
		
	G4PVPlacement* uScint_phys;
	//uScint_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(upScint_pos-1.*cm,0.0*cm,upScint_pos-1.*cm),
	//						uScint_log,"upperScint",det_log,false,0);

        //uScint_phys 
		//= new G4PVPlacement(scintRoll,G4ThreeVector(0.0*cm,0.0*cm,-55*cm),
							//uScint_log,"upperScint",det_log,false,0);
								

 /////////////
 
	G4Box* lowScint = new G4Box("lowScint",20.0*cm,3.5*cm,0.25*cm);
	G4LogicalVolume* lScint_log = new G4LogicalVolume(lowScint,Air,"lowScint",0,0,0);

    // Make sensitive
/*			//G4String 
	DetSDname = "/tracker3";

	qsimTrackerSD* loScintSD = new qsimTrackerSD(DetSDname);
  
	SDman->AddNewDetector(loScintSD);
	lScint_log->SetSensitiveDetector(loScintSD);
*/	
	
		DetSDname = "tracker2";
	 
	 qsimScintDetector* loScintSD = new qsimScintDetector(DetSDname, 2);
	 
	 SDman->AddNewDetector(loScintSD);
	 lScint_log->SetSensitiveDetector(loScintSD);
	 
	
	G4double loScint_pos=upScint_pos+70.71*cm;

//(-1*quartz_z)+(41.25*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
		
	G4PVPlacement* lScint_phys;
	//lScint_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(loScint_pos,0.0*cm,loScint_pos),
	//						lScint_log,"lowerScint",det_log,false,0);

        //lScint_phys 
		//= new G4PVPlacement(scintRoll,G4ThreeVector(0.0*cm,0.0*cm,55*cm),
	     						//lScint_log,"lowerScint",det_log,false,0);


 /////////////
 
	G4Box* Pb_blox = new G4Box("Pb_blox",20.0*cm,10.0*cm,13.0*cm);
		// Really 10x5x10 cm half-lengths, expanded to ensure nothing
		//   can hit the scint. w/o the lead.
	G4LogicalVolume* Pb_log = new G4LogicalVolume(Pb_blox,Pb_Mat,"Lead",0,0,0);

	G4double Pb_pos=loScint_pos-11.25*cm; //(-1*quartz_z)+(30.0*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
		
	G4PVPlacement* Pb_phys;
	//Pb_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(Pb_pos,0.0*cm,Pb_pos),
	//						Pb_log,"Pb",det_log,false,0);

        //Pb_phys 
		//= new G4PVPlacement(scintRoll,G4ThreeVector(0.0*cm,0.0*cm,41*cm),
							//Pb_log,"Pb",det_log,false,0);
  


//	------------- Surfaces --------------
//
// Quartz
//
  G4OpticalSurface* OpQuartzSurface = new G4OpticalSurface("QuartzSurface");
  OpQuartzSurface->SetType(dielectric_dielectric);
  OpQuartzSurface->SetFinish(ground);
  OpQuartzSurface->SetModel(unified);

  //  G4LogicalBorderSurface* QuartzSurface = 
  //                                 new G4LogicalBorderSurface("QuartzSurface",
  //                                 quartz_phys,det_phys,OpQuartzSurface);

  //  if(QuartzSurface->GetVolume1() == quartz_phys) G4cout << "Equal" << G4endl;
  //  if(QuartzSurface->GetVolume2() == det_phys  ) G4cout << "Equal" << G4endl;

// Mirrors and cathode

  G4OpticalSurface* MOpSurface = new G4OpticalSurface("MirrorOpSurface");
  G4OpticalSurface* CTHOpSurface = new G4OpticalSurface("CathodeOpSurface");

  MOpSurface -> SetType(dielectric_metal);
  MOpSurface -> SetFinish(ground);
  MOpSurface -> SetModel(glisur);

  //  G4double polish = 0.8;

/*  G4double Reflectivity3[nEntries] =
            { 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90 };  */


  G4MaterialPropertiesTable* MOpSurfaceProperty = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* COpSurfaceProperty = new G4MaterialPropertiesTable();

  MOpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity3,nEntries);

  MOpSurface -> SetMaterialPropertiesTable(MOpSurfaceProperty);

  COpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity4,nEntries);
  COpSurfaceProperty -> AddProperty("EFFICIENCY",PhotonEnergy,Efficiency4,nEntries);

  CTHOpSurface -> SetMaterialPropertiesTable(COpSurfaceProperty);
 

  //G4LogicalSkinSurface* TSurface = new
               //G4LogicalSkinSurface("TMirrorOpS",tmirror_log,MOpSurface);

	G4LogicalSkinSurface* FrontSurface_1 = new
	G4LogicalSkinSurface("FrontMirrorOpS_1",frontPlate_1_log,MOpSurface);

	G4LogicalSkinSurface* TopSurface_1 = new
	G4LogicalSkinSurface("TopMirrorOpS_1",topPlate_1_log,MOpSurface);

	G4LogicalSkinSurface* TopSurface_2 = new
	G4LogicalSkinSurface("TopMirrorOpS_2",topPlate_2_log,MOpSurface);

	G4LogicalSkinSurface* TopSurface_3 = new
	G4LogicalSkinSurface("TopMirrorOpS_3",topPlate_3_log,MOpSurface);
	
	G4LogicalSkinSurface* BotSurface_1 = new
	G4LogicalSkinSurface("BotMirrorOpS_1",botPlate_1_log,MOpSurface);

	G4LogicalSkinSurface* BotSurface_2 = new
	G4LogicalSkinSurface("BotMirrorOpS_2",botPlate_2_log,MOpSurface);

	G4LogicalSkinSurface* LSurface_1 = new
	G4LogicalSkinSurface("LMirrorOpS_1",LPlate_1_log,MOpSurface);

        G4LogicalSkinSurface* LSurface_2 = new
	G4LogicalSkinSurface("LMirrorOpS_2",LPlate_2_log,MOpSurface);

        G4LogicalSkinSurface* LSurface_3 = new
	G4LogicalSkinSurface("LMirrorOpS_3",LPlate_3_log,MOpSurface);	

	G4LogicalSkinSurface* RSurface_1 = new
	G4LogicalSkinSurface("RMirrorOpS_1",RPlate_1_log,MOpSurface);

        G4LogicalSkinSurface* RSurface_2 = new
	G4LogicalSkinSurface("RMirrorOpS_2",RPlate_2_log,MOpSurface);

        G4LogicalSkinSurface* RSurface_3 = new
	G4LogicalSkinSurface("RMirrorOpS_3",RPlate_3_log,MOpSurface);
	
        G4LogicalSkinSurface* mirror_1_OpS = new
               G4LogicalSkinSurface("mirror_1_OpS",mirror_box_1_log,MOpSurface);

        G4LogicalSkinSurface* mirror_1_1_OpS = new
               G4LogicalSkinSurface("mirror_1_1_OpS",mirror_box_1_1_log,MOpSurface);

        G4LogicalSkinSurface* mirror_2_OpS = new
               G4LogicalSkinSurface("mirror_2_OpS",mirror_box_2_log,MOpSurface);

        G4LogicalSkinSurface* mirror_3_OpS = new
               G4LogicalSkinSurface("mirror_3_OpS",mirror_box_3_log,MOpSurface);

        G4LogicalSkinSurface* mirror_4_OpS = new
               G4LogicalSkinSurface("mirror_4_OpS",mirror_box_4_log,MOpSurface);

        G4LogicalSkinSurface* mirror_4_1_OpS = new
               G4LogicalSkinSurface("mirror_4_1_OpS",mirror_box_4_1_log,MOpSurface);

        G4LogicalSkinSurface* mirror_5_OpS = new
               G4LogicalSkinSurface("mirror_5_OpS",mirror_box_5_log,MOpSurface);

        G4LogicalSkinSurface* mirror_6_OpS = new
               G4LogicalSkinSurface("mirror_6_OpS",mirror_box_6_log,MOpSurface);

        G4LogicalSkinSurface* mirror_7_OpS = new
               G4LogicalSkinSurface("mirror_7_OpS",mirror_box_7_log,MOpSurface);

  //G4LogicalSkinSurface* CSurface = new
               //G4LogicalSkinSurface("CMirrorOpS",cmirror_log,MOpSurface);

  //G4LogicalSkinSurface* CathSurface = new
               //G4LogicalSkinSurface("CathOpS",cath_log,CTHOpSurface);


  //G4LogicalSkinSurface* QuartWinSurface = new 
               //G4LogicalSkinSurface("QuartzWinOpS", QuartzWin_log, OpQuartzSurface);




//
// Generate & Add Material Properties Table attached to the optical surfaces
//
  const G4int num = 2;
  G4double Ephoton[num] = {2.038*eV, 4.144*eV};

  //OpticalQuartzSurface 
  G4double RefractiveIndex[num] = {1.46, 1.46};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OpQuartzSurface->SetMaterialPropertiesTable(myST1);




//always return the physical World
  return det_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
