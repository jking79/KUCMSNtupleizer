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

void KUCMSAodSkimmer::processEvntVars(){
  
  if( DEBUG ) std::cout << "Finding PV Vars" << std::endl;
  selEvtVars.clearBranches(); // <<<<<<<   must do
  
  // calc
  selEvtVars.fillBranch( "PVx", PV_x );
  selEvtVars.fillBranch( "PVy", PV_y );
  selEvtVars.fillBranch( "PVz", PV_z ); 
  
  //fill

  if( DEBUG ) std::cout << "Finding Event wts Vars" << std::endl;
  selEvtVars.fillBranch( "dsKey", dataSetKey );
  float evtGenWgt = 1;
  if( useEvtWgt && doGenInfo ) evtGenWgt = Evt_genWgt;
  selEvtVars.fillBranch( "evtGenWgt", evtGenWgt );
  //sumEvtGenWgt += Evt_genWgt;
  selEvtVars.fillBranch( "evtXSection", xsctn );
  
  float fillWgt = 1;
  if( doGenInfo ) fillWgt = ( ( xsctn * 1000 ) * evtGenWgt  ) / configWgts["sumEvtWgt"];

  selEvtVars.fillBranch( "evtFillWgt", fillWgt );

  bool BadChargedCandidateFilter = true;
  bool BadPFMuonDzFilter = true;
  bool BadPFMuonFilter = true;
  bool EcalDeadCellTriggerPrimitiveFilter = true;
  bool HBHENoiseFilter = true;
  bool HBHENoiseIsoFilter = true;
  bool ecalBadCalibFilter = true;
  bool eeBadScFilter = true;
  bool globalSuperTightHalo2016Filter = true;
  bool goodVertices = true;
  bool hfNoisyHitsFilter = true;

  if( not doGenInfo ){

	BadChargedCandidateFilter = Flag_BadChargedCandidateFilter;
	BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
	BadPFMuonFilter = Flag_BadPFMuonFilter;
	EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
	HBHENoiseFilter = Flag_HBHENoiseFilter;
	HBHENoiseIsoFilter = Flag_HBHENoiseIsoFilter;
	ecalBadCalibFilter = Flag_ecalBadCalibFilter;
	eeBadScFilter = Flag_eeBadScFilter;
	globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
	goodVertices = Flag_goodVertices;
	hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;

  }//<<>>if( not doGenInfo )

  selEvtVars.fillBranch( "Flag_BadChargedCandidateFilter", BadChargedCandidateFilter );//not suggested
  selEvtVars.fillBranch( "Flag_BadPFMuonDzFilter", BadPFMuonDzFilter );//suggested
  selEvtVars.fillBranch( "Flag_BadPFMuonFilter", BadPFMuonFilter );//suggested
  selEvtVars.fillBranch( "Flag_EcalDeadCellTriggerPrimitiveFilter", EcalDeadCellTriggerPrimitiveFilter );//suggested
  selEvtVars.fillBranch( "Flag_HBHENoiseFilter", HBHENoiseFilter );//suggested
  selEvtVars.fillBranch( "Flag_HBHENoiseIsoFilter", HBHENoiseIsoFilter );//suggested
  selEvtVars.fillBranch( "Flag_ecalBadCalibFilter", ecalBadCalibFilter );//suggested
  selEvtVars.fillBranch( "Flag_eeBadScFilter", eeBadScFilter );//suggested
  selEvtVars.fillBranch( "Flag_globalSuperTightHalo2016Filter", globalSuperTightHalo2016Filter );//suggested
  selEvtVars.fillBranch( "Flag_goodVertices", goodVertices );//suggested
  selEvtVars.fillBranch( "Flag_hfNoisyHitsFilter", hfNoisyHitsFilter );//optional

  bool PFMET120_PFMHT120_IDTight = true;
  bool PFMETNoMu120_PFMHTNoMu120_IDTight = true;
  bool PFMET120_PFMHT120_IDTight_PFHT60 = true;
  bool PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = true;

  if( not doGenInfo ){

	PFMET120_PFMHT120_IDTight = HLT_PFMET120_PFMHT120_IDTight_v;
	PFMETNoMu120_PFMHTNoMu120_IDTight = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
	PFMET120_PFMHT120_IDTight_PFHT60 = HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
	PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;

  }//<<>>if( not doGenInfo )

  selEvtVars.fillBranch( "Trigger_PFMET120_PFMHT120_IDTight", PFMET120_PFMHT120_IDTight );
  selEvtVars.fillBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight", PFMETNoMu120_PFMHTNoMu120_IDTight );
  selEvtVars.fillBranch( "Trigger_PFMET120_PFMHT120_IDTight_PFHT60", PFMET120_PFMHT120_IDTight_PFHT60 );
  selEvtVars.fillBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 );

}//<<>>void KUCMSAodSkimmer::processEvntVars()

void KUCMSAodSkimmer::processMet(){

  //intilize
  selMet.clearBranches(); // <<<<<<<   must do

  //calc
  //auto met = std::sqrt(sq2(Met_Cpx)+sq2(Met_Cpy));

  geVars.set("metPx", Met_Cpx );
  geVars.set("metPy", Met_Cpy );
  geVars.set("cmet", Met_CPt );

  //fill
  selMet.fillBranch( "cmet", Met_CPt );
  selMet.fillBranch( "cmetPx", Met_Cpx );
  selMet.fillBranch( "cmetPy", Met_Cpy );
  selMet.fillBranch( "met", Met_pt );
  selMet.fillBranch( "metPx", Met_px );
  selMet.fillBranch( "metPy", Met_py );

}//<<>>void KUCMSAodSkimmer::processMet()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setEvtVarMetBranches( TTree* fOutTree ){

  //fOutTree->Branch( "RunNumber", &RunNumber );
  selEvtVars.makeBranch( "dsKey", "DataSetKey", STR, "Key for source data set of event" );

  selEvtVars.makeBranch( "evtGenWgt", FLOAT );
  selEvtVars.makeBranch( "evtXSection", FLOAT );
  selEvtVars.makeBranch( "evtFillWgt", FLOAT );

  selEvtVars.makeBranch( "PVx", FLOAT );
  selEvtVars.makeBranch( "PVy", FLOAT );
  selEvtVars.makeBranch( "PVz", FLOAT );
  
  selEvtVars.makeBranch( "Flag_BadChargedCandidateFilter", BOOL );
  selEvtVars.makeBranch( "Flag_BadPFMuonDzFilter", BOOL );
  selEvtVars.makeBranch( "Flag_BadPFMuonFilter", BOOL );
  selEvtVars.makeBranch( "Flag_EcalDeadCellTriggerPrimitiveFilter", BOOL );
  selEvtVars.makeBranch( "Flag_HBHENoiseFilter", BOOL );
  selEvtVars.makeBranch( "Flag_HBHENoiseIsoFilter", BOOL );
  selEvtVars.makeBranch( "Flag_ecalBadCalibFilter", BOOL );
  selEvtVars.makeBranch( "Flag_eeBadScFilter", BOOL );
  selEvtVars.makeBranch( "Flag_globalSuperTightHalo2016Filter", BOOL );
  selEvtVars.makeBranch( "Flag_goodVertices", BOOL );
  selEvtVars.makeBranch( "Flag_hfNoisyHitsFilter", BOOL );

  selEvtVars.makeBranch( "Trigger_PFMET120_PFMHT120_IDTight", BOOL );
  selEvtVars.makeBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight", BOOL );
  selEvtVars.makeBranch( "Trigger_PFMET120_PFMHT120_IDTight_PFHT60", BOOL );
  selEvtVars.makeBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", BOOL );

  selEvtVars.attachBranches( fOutTree );

  selMet.makeBranch( "cmet", "selCMet", FLOAT, "Magnitude of event Met corrected for OOT photons" );
  selMet.makeBranch( "cmetPx", "selCMetPx", FLOAT, "Magnitude of event MetPx corrected for OOT photons" );
  selMet.makeBranch( "cmetPy", "selCMetPy", FLOAT, "Magnitude of event MetPy corrected for OOT photons" );
  selMet.makeBranch( "met", "selMet", FLOAT, "Magnitude of event Met" );
  selMet.makeBranch( "metPx", "selMetPx", FLOAT, "Magnitude of event MetPx" );
  selMet.makeBranch( "metPy", "selMetPy", FLOAT, "Magnitude of event MetPy" );

  selMet.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

