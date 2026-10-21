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

// Track-derived branch filling for the SV-enabled workflow.

void KUCMSAodSkimmer::processTracks(){

  int nTracks = Track_nIsoTracks;
  selTracks.fillBranch( "Trk_nIsoTrack", nTracks );

}//<<>>void KUCMSAodSkimmer::processTracks()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setTrackBranches( TTree* fOutTree ){

  std::cout << " - Making Branches for Tracks." << std::endl;
  selTracks.makeBranch( "Trk_nIsoTrack", INT );

  selTracks.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setTrackBranches( TTree* fOutTree )

