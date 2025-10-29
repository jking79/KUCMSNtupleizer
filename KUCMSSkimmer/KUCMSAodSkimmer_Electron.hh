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

void KUCMSAodSkimmer::processElectrons(){

  //if( DEBUG ) std::cout << "Finding electrons" << std::endl;
  //-------- electrons --------------------------------------

  // intilize
  selElectrons.clearBranches(); // <<<<<<<   must do

  uInt nElectrons = Electron_energy->size();
  uInt nLooseEle = 0;
  uInt nSelIsoElectrons = 0;
  float nEpDrMatch = 0;
  float nEpSeedIdMatch = 0;
  uInt nEpsidsolo = 0;
  float elePhoIsoMinDr(10.0); 
  //if( DEBUG ) std::cout << " -- Looping electrons: " << std::endl;
  for( uInt itr = 0; itr < nElectrons; itr++ ){

    if( ( hasGenInfoFlag && ( (*Electron_genIdx)[itr] > -1 ) ) || not hasGenInfoFlag ){
      hist1d[10]->Fill((*Electron_pt)[itr]);
      if( (*Electron_isLVeto)[itr] ) hist1d[11]->Fill((*Electron_pt)[itr]);
      if( (*Electron_isVeto)[itr] ) hist1d[12]->Fill((*Electron_pt)[itr]);
      if( (*Electron_isLoose)[itr] ) hist1d[13]->Fill((*Electron_pt)[itr]);
      if( (*Electron_isMedium)[itr] ) hist1d[14]->Fill((*Electron_pt)[itr]);
    }//<<>>if( (*Electron_hasGenMatch)[itr] )

    if( DEBUG ) std::cout << " ---- Processing  electron: " << itr << " of: " << nElectrons << std::endl;
    float eta = (*Electron_eta)[itr];
    float phi = (*Electron_phi)[itr];
    auto scIdx = (*Electron_scIndex)[itr];
    if( DEBUG ) std::cout << " -- doing sc sstuff " << std::endl;
    auto esid = getLeadRhID( (*SuperCluster_rhIds)[scIdx] );
    //std::cout << " -- esid : " << esid << std::endl;
    auto eletime = (*Electron_seedTOFTime)[itr];
    if( std::abs(eletime) > 30 ) continue;

    uInt nPhotons = Photon_excluded->size();
    float elePhoIsoDr(10.0);
    bool epSeedIdMatch = false;
    bool epsidsolo = true;
    //if( DEBUG ) std::cout << " -- doing ele pho iso sstuff " << std::endl;
    for( uInt pit2 = 0; pit2 < nPhotons; pit2++ ){
      float eta2 = (*Photon_eta)[pit2];
      float phi2 = (*Photon_phi)[pit2];
      auto scIdx2 = (*Photon_scIndex)[pit2];
      auto dphi12 = dPhi( phi, phi2 );
      float elephodr = hypo( eta-eta2, dphi12 );
      //std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << elephodr << std::endl;
      if( elephodr < elePhoIsoDr ) elePhoIsoDr = elephodr;
      //std::cout << " -- scIdx : " << scIdx2 << std::endl; 
      auto rhIds = (*SuperCluster_rhIds)[scIdx2];
      //std::cout << " -- scRh1 : " << rhIds[0] << std::endl;
      auto psid = getLeadRhID(rhIds);
      //std::cout << " -- psid : " << psid << std::endl;
      if( esid == psid ){ if( epSeedIdMatch == true ) epsidsolo = false ; epSeedIdMatch = true; }
    }//for( uInt pit2 = 0; pit2 < nPhotons; pit2++ )
    //std::cout << " ---- next electrons: " << std::endl;
    elePhoIsoMinDr = elePhoIsoDr;
    bool epDrMatch = ( elePhoIsoDr < 0.2  ) ? true : false;
    bool eleminpt = (*Electron_pt)[itr] > 10;
    if( not epDrMatch & eleminpt ) nSelIsoElectrons++;

    if( epDrMatch ) nEpDrMatch++;
    if( epSeedIdMatch ) nEpSeedIdMatch++;
    if( epsidsolo ) nEpsidsolo++;
    if( DEBUG ) std::cout << " ---- next electrons: " << std::endl;
    //std::cout << " ---- filling electrons: " << std::endl;
    selElectrons.fillBranch( "elePhoMinDr", elePhoIsoMinDr );
    selElectrons.fillBranch( "eleHasSVMatch", (*Electron_hasSVMatch)[itr] );
    selElectrons.fillBranch( "eleIsLoose", (*Electron_isLoose)[itr] );
    if( (*Electron_isLoose)[itr] ) nLooseEle++;

	bool hemEligible1 = (*Electron_isLoose)[itr] && (*Electron_pt)[itr] > 30;
    bool hemEligible2 = (*Electron_isMedium)[itr] && (*Electron_pt)[itr] > 50;
    bool isInHemRegion = inHEMRegion( eta, phi );
    hemBits["el1"] = isInHemRegion && hemEligible1;
    hemBits["el2"] = isInHemRegion && hemEligible2;

    //if( hemEligible && inHEMRegion( eta, phi ) ) hasHemObj = true;	

  }//<<>>for( int itr = 0; itr < nElectrons; itr++ )
  if( DEBUG ) std::cout << " -- Finishd looping electrons " << std::endl;

  geCnts.set("nSelEle",nSelIsoElectrons);// passes min pt, isolated from photons, && good time criteria

  selElectrons.fillBranch( "nElectrons", nElectrons );
  if( nElectrons == 0 ) nElectrons = 1;
  //std::cout << " -- Finish nEle looping electrons " << std::endl;
  selElectrons.fillBranch( "epDrMatch", nEpDrMatch/float(nElectrons) );
  selElectrons.fillBranch( "epSeedIdMatch", nEpSeedIdMatch/float(nElectrons) );
  selElectrons.fillBranch( "epSeedIdUnique", nEpsidsolo/float(nElectrons) );
  //std::cout << " -- Finish New looping electrons " << std::endl;
  selElectrons.fillBranch( "nSelElectrons", Electron_nSelElectrons );
  selElectrons.fillBranch( "nEleSVMatched", Electron_nSVMatched );
  selElectrons.fillBranch( "nSelIsoElectrons", nLooseEle );

}//<<>>void KUCMSAodSkimmer::processElectrons

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setElectronBranches( TTree* fOutTree ){

  selElectrons.makeBranch( "nElectrons", UINT );
  selElectrons.makeBranch( "nSelIsoElectrons", UINT );
  selElectrons.makeBranch( "epDrMatch", VFLOAT );
  selElectrons.makeBranch( "epSeedIdMatch", VFLOAT );
  selElectrons.makeBranch( "epSeedIdUnique", VFLOAT );
  selElectrons.makeBranch( "elePhoMinDr", VFLOAT );
  selElectrons.makeBranch( "eleIsLoose", VBOOL );
  selElectrons.makeBranch( "nSelElectrons", UINT );
  selElectrons.makeBranch( "eleHasSVMatch", VBOOL );
  selElectrons.makeBranch( "nEleSVMatched", INT );

  selElectrons.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

