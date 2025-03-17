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
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Evaporation.cc
//
//      Authors:        V.Ivanchenko, N.Chalyi
//
//      Creation date: 25 april 2024
//
//      Modifications: 
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeList.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"

#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"

#include "G4Fragment.hh"
#include "G4NucleiProperties.hh"
#include "G4LorentzVector.hh"

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4ParticleInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4GammaNuclearXS.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4EvaporationChannel.hh"
#include "G4EvaporationProbability.hh"
#include "G4CoulombBarrier.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4KalbachCrossSection.hh"
#include "G4ChatterjeeCrossSection.hh"
#include "G4ChargeExchangeXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4DecayPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Log.hh"

int main(int argc, char** argv)
{
  
  // Control on input
  if(argc < 2) {
    G4cout << "Input parameters are not specified! Exit" << G4endl;
    return 1;
  }
  G4int verbose = 2;
  if(0 < verbose) {
    G4cout << "====================================================" << G4endl;
    G4cout << "======         ChargeExchangeXS Test        ========" << G4endl;
    G4cout << "====================================================" << G4endl;
  }
  
  G4int Z = 6;
  G4int A = 12;
  G4double de = 1*GeV;
  G4double emax = 100.*GeV + de*0.5;
  G4double emin = de*0.5;
  G4String partname;
  G4String m;

  // convert string to Z 
  G4String sz = "";
  if(argc >= 2) {
    sz = argv[1];
    std::istringstream is(sz);
    is >> Z;
  }
  // convert string to A
  if(argc >= 3) {
    sz = argv[2];
    std::istringstream is(sz);
    is >> A;
  }
  //
  auto nist = G4NistManager::Instance();
  auto elm = nist->FindOrBuildElement(Z);
  if(argc >= 4) {
    sz = argv[3];
    std::istringstream is(sz);
    is >> partname; 
  }
  /*if(argc >= 5) {
    m = argv[4];
    std::istringstream is(m);
    is >> m; 
  }*/

  auto fParticleList = new G4DecayPhysics(1);
  fParticleList->ConstructParticle();
  //G4Material* mat = G4Material::GetMaterial(m);

  G4int nbins = (G4int) ((emax-emin)/de);
  
  //if(0 < verbose) {
  //	G4cout << "### Z=" << Z << " A=" << A  << " nbins=" << nbins << " Emax(MeV)=" << emax << G4endl;
  //	}
  // particles
  //Нужно. Добавить нужные и реорганизовать pi+-0 K+- LS eta eta_prime F2 omega electron
  
  // histograms name
  Histo histo;
  G4String hname;
  hname = "hist_" + partname + "_z" + std::to_string(Z) + "_a" + std::to_string(A);
  histo.Add1D("0","Total cross section",nbins,0,emax,1);
  histo.Add1D("1","Charge Exchange cross section",nbins,0,emax,1);
  histo.Add1D("2","Ration of Charge Exchage to total cross section",nbins,0,emax,1);
  
  histo.SetFileName(hname);
  histo.Book();

  if(0 < verbose) {
    G4cout << "Histograms are booked output file <" << hname << "> " << G4endl;
  }

  const G4ParticleDefinition* part;
  G4double plmass[8] = {0.0};
  /*part[0] = G4Gamma::Gamma();
  part[1] = G4Electron::Electron();
  part[2] = G4Neutron::Neutron();
  part[3] = G4Proton::Proton();
  part[4] = G4Deuteron::Deuteron();
  part[5] = G4Triton::Triton();
  part[6] = G4He3::He3();
  part[7] = G4Alpha::Alpha();*/
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();
  G4cout << "### Z=" << Z << " A=" << A  << " nbins=" << nbins << " Emax(MeV)=" << emax << G4endl;
  //------------------------
  G4ThreeVector dir(0.0,0.0,1.0);
  part = partTable->FindParticle(partname);
  auto dp = new G4DynamicParticle(part, dir, emax);
  
  // x-sections
  /*G4VCrossSectionDataSet* xs[8];
  xs[0] = new G4ChargeExchangeXS();
  xs[0]->BuildPhysicsTable(*(part[0]));
  xs[1] = nullptr;
  xs[2] = new BGGPionInelasticXS();
  xs[2]->BuildPhysicsTable(*(part[2]));
  pmass[2] = neutron_mass_c2;
  for(G4int i=3; i<8; ++i) {
    pmass[i] = part[i]->GetPDGMass();
    xs[i] = new G4ParticleInelasticXS(part[i]);
    xs[i]->BuildPhysicsTable(*(part[i]));
  }*/
  
  G4VCrossSectionDataSet* xs1 = new G4ChargeExchangeXS();
  G4VCrossSectionDataSet* xs2 = new G4BGGPionInelasticXS(part);

  G4double ArrXS1[nbins];
  G4double ArrXS2[nbins];
  
  G4double e = emin+de*0.5;
  for (G4int i = 0; i < nbins; i++)
  {
    //G4cout << e << G4endl;
    G4bool Ap1 = xs1->IsElementApplicable(dp, Z);
    G4bool Ap2 = xs2->IsElementApplicable(dp, Z);
    if (Ap1)
    {
      ArrXS1[i]=xs1->GetElementCrossSection(dp, Z);
    }
    if (Ap2)
    {
      ArrXS2[i]=xs2->GetElementCrossSection(dp, Z);
    }
    G4cout << e << std::setw(35) << ArrXS1[i] << std::setw(35) << ArrXS2[i] <<G4endl;
    
    histo.Fill(0, e, ArrXS1[i]/barn);
    histo.Fill(1, e, ArrXS2[i]/barn);
    histo.Fill(2, e, ArrXS1[i]/ArrXS2[i]);
    
    e += de;
    dp->SetKineticEnergy(e);
  }
  
  // fragment
  /*if (A <= 0) A = lrint(aeff[Z]);
  G4double mass = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double massEx = mass + eexc;
  G4LorentzVector lv(0.0,0.0,0.0,massEx);
  G4Fragment* frag = new G4Fragment(A, Z, lv);

  // handler
  G4ExcitationHandler handler;
  handler.Initialise();
  G4VEvaporation* evap = handler.GetEvaporation();
  G4CoulombBarrier* bCoulomb[8] = {nullptr};
  G4double CB[8] = {0.0};

  G4int prec = G4cout.precision(6);

  for (G4int ii=0; ii<8; ++ii) {
    G4VEvaporationChannel* ch = evap->GetChannel(ii);
    G4EvaporationChannel* cevap = nullptr;
    G4EvaporationProbability* pevap = nullptr;
    if (ii > 1) {
      cevap = static_cast<G4EvaporationChannel*>(ch);
      pevap = cevap->GetEvaporationProbability();
    }
    G4int theZ = 0;
    G4int theA = 0;
    G4int idx  = 0;
    if (ii == 2) {
      theA = 1;
    } else if (ii == 3) {
      theZ = 1;
      theA = 1;
      idx = 1;
    } else if (ii == 4) {
      theZ = 1;
      theA = 2;
      idx = 2;
    } else if (ii == 5) {
      theZ = 1;
      theA = 3;
      idx = 3;
    } else if (ii == 6) {
      theZ = 2;
      theA = 3;
      idx = 4;
    } else if (ii == 7) {
      theZ = 2;
      theA = 4;
      idx = 5;
    }
   
    G4double evapmass = G4NucleiProperties::GetNuclearMass(theA, theZ);
    G4double evapmass2 = evapmass*evapmass;
    G4int resA = A - theA;
    G4int resZ = Z - theZ;
    G4double rmass = G4NucleiProperties::GetNuclearMass(resA, resZ);
    G4cout << "&& rmass=" << rmass << "; rmass-nist= " << rmass - nist->GetIsotopeMass(resZ,resA) << G4endl;
    G4LorentzVector Rlv(0.0,0.0,0.0,rmass+eexc);
    G4Fragment* Rfrag = new G4Fragment(resA, resZ, Rlv);
    if (resA < theA || resA < resZ || resZ < 0 || (resA == theA && resZ < theZ)
       || ((resA > 1) && (resA == resZ || resZ == 0))) {
      if (0 < verbose) {
        G4cout << "WARN: Z=" << Z << " A=" << A << " channel=" << ii
               << " not allowed resZ=" << resZ << " resA=" << resA
               << " Rmass(GeV)=" << rmass/GeV << G4endl;
      }
      continue;
    }
    G4double resA13 = nist->GetZ13(resA);
    G4double muu = G4KalbachCrossSection::ComputePowerParameter(resA, idx); // ?
    if (ii > 1) {
      bCoulomb[ii] = new G4CoulombBarrier(theA, theZ);
      CB[ii] = bCoulomb[ii]->GetCoulombBarrier(resA, resZ, 0.0);
    }
    // total probability - provides initialisation of the channel
    G4double prob1 = ch->GetEmissionProbability(frag);
    if (1 < verbose) {
    	G4cout << "=================================================" << G4endl;
    	G4cout << "### Channel " << ii << " Z=" << Z << " A=" << A << " resZ=" << resZ << " resA=" << resA << " Eexc(MeV)=" << eexc << " CB(MeV)=" << CB[ii] << " prob=" << prob1 << G4endl;
    	G4cout << "=================================================" << G4endl;
    	}
    G4double e, xshad{0.0}, xskalb{0.0}, xschat{0.0};
    for (G4int i=0; i<nbins; ++i) {
      e = emin + i*de;
      if (ii > 1) {
        xshad = xs[ii]->ComputeIsoCrossSection(e, G4Log(e), part[ii], resZ, resA)/CLHEP::millibarn;
        xskalb = G4KalbachCrossSection::ComputeCrossSection(e, CB[ii], resA13, muu, idx, theZ, theA, resA);
        xschat = G4ChatterjeeCrossSection::ComputeCrossSection(e, CB[ii], resA13, muu, idx, theZ, resA);
      } else if (0 == ii) {
        dp->SetKineticEnergy(e);
        xshad = xs[ii]->ComputeCrossSection(dp, elm)/CLHEP::millibarn;
      }
      G4cout << i << ".  E(MeV)=" << e << " XS(mb)=" << xshad << " XSKalbach(mb)=" << xskalb << " XSChatterjee(mb)=" << xschat << G4endl;
      histo.Fill(16 + ii, e, xshad);
      histo.Fill(24 + ii, e, xschat);
      histo.Fill(ii, e, xskalb);
    }
    if (prob1 <= 0.0) { continue; }
    for (G4int i=0; i<=ibins; ++i) {
      e = emin + i*de;
      G4double prob = ch->ComputeProbability(frag, e);
      if (prob > 0.0) { histo.Fill(8 + ii, e, prob); }
    }
  }
  */
  if (verbose > 0) { G4cout << "###### Save histograms" << G4endl; }
  histo.Save();
 
  if (verbose > 0) {
    G4cout << "###### End of run # " << G4endl;
  }
}
