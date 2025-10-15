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

void KUCMSAodSkimmer::processMuons(){

  if( DEBUG ) std::cout << "Finding muons" << std::endl;
  //-------- muons --------------------------------------
    
  selMuons.clearBranches(); // <<<<<<<   must do

  uInt nMuons = Muon_energy->size();
  uInt nLooseMuons = 0;
  if( DEBUG ) std::cout << " -- Looping muons: " << std::endl;
  for( uInt itr = 0; itr < nMuons; itr++ ){
	
    //float energy = (*Muons_energy)[itr];
    selMuons.fillBranch( "muonIsLoose", (*Muon_isLoose)[itr] );
    selMuons.fillBranch( "muonIsMedium", (*Muon_isMedium)[itr] );

    if( (*Muon_isLoose)[itr] ) nLooseMuons++;

  }//<<>>for( uInt itr = 0; itr < nMuons; itr++ )
  selMuons.fillBranch( "nMuons", nMuons );
  selMuons.fillBranch( "nLooseMuons", nLooseMuons );
  selMuons.fillBranch( "nMuonSVMatched", Muon_nSVMatched );

}//<<>>void KUCMSAodSkimmer::processMuons

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setMuonsBranches( TTree* fOutTree ){

    selMuons.makeBranch( "nLooseMuons", UINT );
    selMuons.makeBranch( "nMuons", UINT );
    selMuons.makeBranch( "muonIsLoose", VBOOL );
    selMuons.makeBranch( "muonIsMedium", VBOOL );
    //selMuons.makeBranch( "muonHasSVMatch", VBOOL );
    selMuons.makeBranch( "nMuonSVMatched", INT );

    selMuons.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

