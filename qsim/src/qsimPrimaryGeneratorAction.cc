#include "qsimPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "qsimIO.hh"
#include "qsimEvent.hh"
#include "qsimtypes.hh"
#include "globals.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

//
#include <fstream>
#include <iostream>
//

//
static std::ifstream myfile;
//

qsimPrimaryGeneratorAction::qsimPrimaryGeneratorAction() {
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);


  fDefaultEvent = new qsimEvent();

  fXmin =  -7.0*cm;
  fXmax =   7.0*cm;

  fYmin =  -1.75*cm;
  fYmax =   1.75*cm;

  fZmin =  -65.0*cm;
  fZmax =  -47.0*cm;

  fEmin = 220*MeV;
  fEmax = 100000*MeV;

  fThetaMin = 40.0*deg;
  fThetaMax = 50.0*deg;
  fPhiMin = -5.0*deg;
  fPhiMax = 5.0*deg;
}

qsimPrimaryGeneratorAction::~qsimPrimaryGeneratorAction() {

  delete fParticleGun;
  delete fDefaultEvent;
}


void qsimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    /*  Generate event, set IO data */

    // Use default, static single generator
    // Update this just in case things changed
    // from the command user interface
    fDefaultEvent->Reset();

    // Set data //////////////////////////////////
    // Magic happens here

    double xPos = CLHEP::RandFlat::shoot( fXmin, fXmax );
    double yPos = CLHEP::RandFlat::shoot( fYmin, fYmax );
    double zPos = CLHEP::RandFlat::shoot( fZmin, fZmax );
    double theta = CLHEP::RandFlat::shoot( fThetaMin, fThetaMax );
    double phi = CLHEP::RandFlat::shoot( fPhiMin, fPhiMax );
    double E = CLHEP::RandFlat::shoot( fEmin, fEmax );
    //
    /*double Ene = 0, Prob = 0, thisEne = 0, thisProb = 0;
    myfile.open ("/home/bulacarl/solid/qsim_A/build/MuonDistribution.dat");
    if (myfile.is_open()) {
      while (!myfile.eof()) {
        myfile >> Ene >> Prob;
        //G4cout << Ene << "\t" << Prob << G4endl;
        if (Ene <= E*MeV) {
          thisProb = Prob;
          thisEne = Ene;
        }
      }
    }
    myfile.close();*/
    //G4cout << E << "\t" << thisEne << "\t" << thisProb << G4endl; 
    //
    double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
    
    assert( E > 0.0 );
    assert( E > mass );

    double p = sqrt( E*E - mass*mass );

    double pX = sin(theta)*cos(phi)*p;
    double pY = sin(theta)*sin(phi)*p;
    double pZ = cos(theta)*p;

    //if (CLHEP::RandFlat::shoot( 0.0, 1.0 ) <= thisProb) {
    fDefaultEvent->ProduceNewParticle(
	    G4ThreeVector(xPos, yPos, zPos),
	    G4ThreeVector(pX, pY, pZ ),
	    fParticleGun->GetParticleDefinition()->GetParticleName() );

    /////////////////////////////////////////////////////////////
    // Register and create event
    // 
    double kinE = sqrt(fDefaultEvent->fPartMom[0].mag()*fDefaultEvent->fPartMom[0].mag()
	    + fDefaultEvent->fPartType[0]->GetPDGMass()*fDefaultEvent->fPartType[0]->GetPDGMass() )
	-  fDefaultEvent->fPartType[0]->GetPDGMass();

      fParticleGun->SetParticleDefinition(fDefaultEvent->fPartType[0]);
      fParticleGun->SetParticleMomentumDirection(fDefaultEvent->fPartMom[0].unit());
      fParticleGun->SetParticleEnergy( kinE  );
      fParticleGun->SetParticlePosition( fDefaultEvent->fPartPos[0] );
      fIO->SetEventData(fDefaultEvent);
      fParticleGun->GeneratePrimaryVertex(anEvent);
      //thisProb = 0;      
    //}
}

G4ParticleGun* qsimPrimaryGeneratorAction::GetParticleGun() {
  return fParticleGun;
} 

