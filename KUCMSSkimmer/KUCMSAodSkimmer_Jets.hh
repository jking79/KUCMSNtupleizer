//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSVSkimmer.hh"
#include "KUCMSHelperFunctions.hh"

//#define DEBUG true
#define DEBUG false

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//  do any processing and calulations for objects and save values to output varibles 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// void KUCMSAodSkimmer::processTemplate(){
// 	
// 		Clear out branch vector varibles &/or initilize other out branch vars
//	------------------------------------------------
// 		Do any calculations / cuts
//  ------------------------------------------------
//		Fill out branch varibles
//
//}//<<>>void KUCMSAodSkimmer::processTemplate()

void KUCMSAodSkimmer::processJets(){

  // intilize
  selJets.clearBranches(); // <<<<<<<   must do

  // calc
  //if( DEBUG ) std::cout << "Finding jets" << std::endl;
  //--------- jets --------------------------
  isSelJet.clear(); 

  uInt nSelJets = 0;
  uInt nJets = Jet_energy->size();
  std::vector<float> seljetpt;
  std::vector<float> seljeteta;
  std::vector<float> seljetphi;
  std::vector<float> seljetmass;
  if( DEBUG ) std::cout << " - Looping over for " << nJets << " jets" << std::endl;
  for( uInt it = 0; it < nJets; it++ ){

    // pull values ---------------------------------------------------
    auto energy = (*Jet_energy)[it];
    auto mass = (*Jet_mass)[it];
    //if( DEBUG ) std::cout << " - Finding Jet Quality" << std::endl;
    auto quality = getJetQuality(it);
    auto pt = (*Jet_pt)[it];
    auto eta = (*Jet_eta)[it];
    auto phi = (*Jet_phi)[it];
    auto rhids = (*Jet_drRhIds)[it];
    if( DEBUG ) std::cout << " -- jettDrId size: " << rhids.size() << std::endl;		

    auto area = (*Jet_area)[it];
    auto chEmEF = (*Jet_chEmEF)[it];
    auto chHEF = (*Jet_chHEF)[it];
    auto chHM = (*Jet_chHM)[it];
    auto muEF = (*Jet_muEF)[it];
    auto neEmEF = (*Jet_neEmEF)[it];
    auto neHEF = (*Jet_neHEF)[it];
    auto neHM = (*Jet_neHM)[it];

    int qrkllpId = 0; //*Jet_genLlpId)[it];
    int gjllpId = 0; //*Jet_genLlpId)[it];
    float jgdpt = 0; //*Jet_genDptMatch)[it];
    float jgdr = 0; //*Jet_genDrMatch)[it];
    float jge = 0; //*Jet_genEnergy)[it];
    float jgeta = 0; //*Jet_genEta)[it];
    float jgimpang = 0; //*Jet_genImpactAngle)[it];
    float llpdp = 0; //*Jet_genLlpDp)[it];
    float llpdr = 0; //*Jet_genLlpDr)[it];
    float jgpt = 0; //*Jet_genPt)[it];
    float jgtof = 0; //*Jet_genTOF)[it];
    float jgt = 0; //(*Jet_genTime)[it];
    float jgllpt = 0; //*Jet_genTimeLLP)[it];
    if( hasGenInfoFlag ){
      //qrkllpId = int((*Jet_genQrkLlpId)[it]);
      gjllpId = int((*Jet_genJetLlpId)[it]);
      jgdpt = (*Jet_genDptMatch)[it];
      jgdr = (*Jet_genDrMatch)[it];
      jge = (*Jet_genEnergy)[it];
      jgeta = (*Jet_genEta)[it];
      jgimpang = (*Jet_genImpactAngle)[it];
      //llpdp = (*Jet_genQrkLlpDp)[it];
      //llpdr = (*Jet_genQrkLlpDr)[it];
      jgpt = (*Jet_genPt)[it];
      jgtof = (*Jet_genTOF)[it];
      jgt = (*Jet_genTime)[it];
      jgllpt = (*Jet_genTimeLLP)[it];
    }//if( doGenInfo )

    // get cacluated values-------------------------------------------
    if( DEBUG ) std::cout << " - Finding Jet Time." << std::endl;
    auto rhenergies = getRhGrpEnergies( rhids );
    auto rhtimes = getRhGrpTimes( rhids );
    auto timedist = getDistStats( rhtimes, rhenergies );
    auto time = timedist[6];

    if( DEBUG ) std::cout << " - Jet Obj selection." << std::endl;
    // jet object selection ------------------------------------------
    auto overMinPt = pt > 30; 
    auto underMaxEta = std::abs(eta) <= 2.4;
    auto isMinQuality = quality > 1; // 2 = "tight" 3 = "tighter"

    bool isNotPhoJet = not phoJetVeto[it];

	bool hemEligible1 = overMinPt && isMinQuality;
    bool hemEligible2 = pt > 50 && quality > 2;
    bool isInHemRegion = inHEMRegion( eta, phi );
    hemBits["jet1"] = isInHemRegion && hemEligible1;
    hemBits["jet2"] = isInHemRegion && hemEligible2;

    //if( hemEligible && inHEMRegion( eta, phi ) ) hasHemObj = true;

    auto jetSelected = underMaxEta && isMinQuality && overMinPt && isNotPhoJet;
    if( not jetSelected ){ isSelJet.push_back(false); continue; }
	isSelJet.push_back(true);
    nSelJets++;

    seljetpt.push_back(pt);
    seljeteta.push_back(eta);
    seljetphi.push_back(phi);
    seljetmass.push_back(mass);

    // fill vectors
    selJets.fillBranch( "selJetGenLlpId", gjllpId );
    selJets.fillBranch( "selJetQrkLlpId", qrkllpId );
    selJets.fillBranch( "selJetQuality", quality );
    selJets.fillBranch( "selJetPt", pt);
    selJets.fillBranch( "selJetMass", mass);
    selJets.fillBranch( "selJetEnergy", energy);
    selJets.fillBranch( "selJetEta", eta);
    selJets.fillBranch( "selJetPhi", phi);
    selJets.fillBranch( "selJetTime", time);

    selJets.fillBranch( "selJetArea", area ); //*Jet_area)[it]; 
    selJets.fillBranch( "selJetChEmEF", chEmEF ); //*Jet_chEmEF)[it]; 
    selJets.fillBranch( "selJetchHEF", chHEF ); //*Jet_chHEF)[it]; 
    selJets.fillBranch( "selJetChHM", chHM ); //*Jet_chHM)[it]; 
    selJets.fillBranch( "selJetMuEF", muEF ); //*Jet_muEF)[it]; 
    selJets.fillBranch( "selJetNeEmEF", neEmEF ); //*Jet_neEmEF)[it]; 
    selJets.fillBranch( "selJetNeHEF", neHEF ); //*Jet_neHEF)[it]; 
    selJets.fillBranch( "selJetNeHM", neHM ); //*Jet_neHM)[it]; 

    selJets.fillBranch( "selGenJetDpt", jgdpt ); //*Jet_genDptMatch)[it]; 
    selJets.fillBranch( "selGenJetdr", jgdr ); //*Jet_genDrMatch)[it]; 
    selJets.fillBranch( "selGenJetEnergy", jge ); //*Jet_genEnergy)[it]; 
    selJets.fillBranch( "selGenJeteta", jgeta ); //*Jet_genEta)[it]; 
    selJets.fillBranch( "selGenJetImpAng", jgimpang ); //*Jet_genImpactAngle)[it]; 
    selJets.fillBranch( "selJetLlpDp", llpdp ); //*Jet_genLlpDp)[it]; 
    selJets.fillBranch( "selJetLlpDr", llpdr ); //*Jet_genLlpDr)[it]; 
    selJets.fillBranch( "selGenJetPt", jgpt ); //*Jet_genPt)[it]; 
    selJets.fillBranch( "selGenJetTof", jgtof ); //*Jet_genTOF)[it]; 
    selJets.fillBranch( "selGenJetTime", jgt ); // (*Jet_genTime)[it]; 
    selJets.fillBranch( "selGenJetLlpTime", jgllpt ); //*Jet_genTimeLLP)[it]; 

  }//<<>>for( int it = 0; it < nJets; it++ )
  if( DEBUG ) std::cout << " - Finished Jet loop." << std::endl;

  // fill other
  geVects.set( "selJetPt", seljetpt );
  geVects.set( "selJetEta", seljeteta );
  geVects.set( "selJetPhi", seljetphi );
  geVects.set( "selJetMass", seljetmass );
  geCnts.set( "nSelJets", nSelJets );
  selJets.fillBranch( "nJets", nJets);
  selJets.fillBranch( "nSelJets", nSelJets );

}//<<>>void KUCMSAodSkimmer::processJets()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setJetsBranches( TTree* fOutTree ){

    //selJets.makeBranch( "JetHt", &JetHt );
    selJets.makeBranch( "nJets", UINT );
    selJets.makeBranch( "nSelJets", UINT );
    selJets.makeBranch( "selJetGenLlpId", VINT ); 
    selJets.makeBranch( "selJetQrkLlpId", VINT );
    selJets.makeBranch( "selJetQuality", VINT ); 
    selJets.makeBranch( "selJetPt", VFLOAT ); 
    selJets.makeBranch( "selJetEnergy", VFLOAT );
    selJets.makeBranch( "selJetEta", VFLOAT ); 
    selJets.makeBranch( "selJetPhi", VFLOAT );
    selJets.makeBranch( "selJetTime", VFLOAT ); 
    selJets.makeBranch( "selJetMass", VFLOAT );

    selJets.makeBranch( "selJetArea", VFLOAT ); //*Jet_area)[it]; 
    selJets.makeBranch( "selJetChEmEF", VFLOAT ); //*Jet_chEmEF)[it]; 
    selJets.makeBranch( "selJetchHEF", VFLOAT ); //*Jet_chHEF)[it]; 
    selJets.makeBranch( "selJetChHM", VFLOAT ); //*Jet_chHM)[it]; 
    selJets.makeBranch( "selJetMuEF", VFLOAT ); //*Jet_muEF)[it]; 
    selJets.makeBranch( "selJetNeEmEF", VFLOAT ); //*Jet_neEmEF)[it]; 
    selJets.makeBranch( "selJetNeHEF", VFLOAT ); //*Jet_neHEF)[it]; 
    selJets.makeBranch( "selJetNeHM", VFLOAT ); //*Jet_neHM)[it]; 

    selJets.makeBranch( "selGenJetDpt", VFLOAT ); //*Jet_genDptMatch)[it]; 
    selJets.makeBranch( "selGenJetdr", VFLOAT ); //*Jet_genDrMatch)[it]; 
    selJets.makeBranch( "selGenJetEnergy", VFLOAT ); //*Jet_genEnergy)[it]; 
    selJets.makeBranch( "selGenJeteta", VFLOAT ); //*Jet_genEta)[it]; 
    selJets.makeBranch( "selGenJetImpAng", VFLOAT ); //*Jet_genImpactAngle)[it]; 
    selJets.makeBranch( "selJetLlpDp", VFLOAT ); //*Jet_genLlpDp)[it]; 
    selJets.makeBranch( "selJetLlpDr", VFLOAT ); //*Jet_genLlpDr)[it]; 
    selJets.makeBranch( "selGenJetPt", VFLOAT ); //*Jet_genPt)[it]; 
    selJets.makeBranch( "selGenJetTof", VFLOAT ); //*Jet_genTOF)[it]; 
    selJets.makeBranch( "selGenJetTime", VFLOAT ); // (*Jet_genTime)[it]; 
    selJets.makeBranch( "selGenJetLlpTime", VFLOAT ); //*Jet_genTimeLLP)[it]; 

    // add new jet branches above
    selJets.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

//------------------------------------------------------------------------------------------------------------
// object quality ids
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getJetQuality( int it ){

  const auto eta  = std::abs((*Jet_eta)[it]);
  const auto NHF  = (*Jet_neHEF)[it];
  const auto NEMF = (*Jet_neEmEF)[it];
  const auto CHF  = (*Jet_chHEF)[it];
  const auto CEMF = (*Jet_chEmEF)[it];
  const auto NHM  = (*Jet_neHM)[it];
  const auto CHM  = (*Jet_chHM)[it];
  const auto SHM  = NHM + CHM;
  const auto MUF  = (*Jet_muEF)[it];

  int tighter = 3;
  int tightLepVeto = 0;
  int tight = 2;
  int loose = 1;

  bool nhfup  = NHF  <= 0.90;   // delpho
  bool nhflw  = NHF  >= 0.2;

  bool nemfup1 = NEMF <= 0.90; // delpho
  bool nemfup2 = NEMF <= 0.99;
  bool nemf80 = NEMF <= 0.80;
  bool nemflw = NEMF >= 0.01;
  bool nemf10 = NEMF >= 0.1;

  bool shm1  = SHM  >= 1;
  bool muf8  = MUF  <= 0.80;
  bool chf0  = CHF  >= 0;       // delpho
  bool chf10  = CHF  >= 0.10;
  bool chm0  = CHM  >= 0;       // delpho
  bool cemf8 = CEMF >= 0.80;
  bool nhm2  = NHM  >= 1;
  bool nhm10 = NHM  >= 10;

  bool eta1 = eta <= 2.6;
  bool eta2 = eta <= 2.7;
  bool eta3 = eta <= 3.0;

  if (eta1){
    if      (nhfup && nemfup1 && shm1 && muf8 && chf0 && chm0 && cemf8) return tightLepVeto;
    else if (nhfup && nemf80 && shm1 && chf10 && chm0) return tighter;
    else if (nhfup && nemfup1 && shm1 && chf0 && chm0) return tight;
    else    return loose;
  } else if (!eta1 && eta2 ){ //<<>>if (eta1)
    if      (nhfup && nemfup2 && chm0 && muf8 && cemf8) return tightLepVeto;
    else if (nhfup && nemf80 && chm0) return tighter;
    else if (nhfup && nemfup2 && chm0) return tight;
    else    return loose;
  } else if (!eta2 && eta3){ //<<>>if (eta1) : else if
    if      (nemf10 && nemf80 && nhm2) return tighter;
    else if (nemflw && nemfup2 && nhm2) return tight;
    else    return loose;
  } else { //<<>>if (eta1) : else if : else if
    if      (nhflw && nemfup1 && nhm10) return tight;
    else    return loose;
  }//<<>>if (eta1) : else if : else if : else

  return -1; // should not happen

}//<<>>int KUCMSAodSkimmer::getJetQuality( int iter )


