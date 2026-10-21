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
#include <cmath>
#include <cstdlib>

//#define DEBUG true
#define DEBUG false

// Event-level and MET quantities.  These values are filled after the object
// modules so event-level flags can include information from the full event.

void KUCMSAodSkimmer::processEvntVars(){

  if( DEBUG ) std::cout << "Finding PV Vars" << std::endl;
  selEvtVars.clearBranches(); // <<<<<<<   must do

  // calc
  selEvtVars.fillBranch( "PVx", PV_x );
  selEvtVars.fillBranch( "PVy", PV_y );
  selEvtVars.fillBranch( "PVz", PV_z );

  //fill
  bool useSMSWgts = false;
  if( ( dataSetKey.find("gogo") != std::string::npos ) || ( dataSetKey.find("sqsq") != std::string::npos ) ){ useSMSWgts = true; }
  if( useEvtGenWgtFlag == false ) useSMSWgts = true;

  if( DEBUG ) std::cout << "Finding Event wts Vars" << std::endl;
  selEvtVars.fillBranch( "dsKey", dataSetKey );
  float evtGenWgt = 1;
  if( hasGenInfoFlag ) evtGenWgt = Evt_genWgt;
  if( useSMSWgts ) evtGenWgt = 1;
  selEvtVars.fillBranch( "evtGenWgt", evtGenWgt );
  selEvtVars.fillBranch( "evtXSection", xsctn );

  float fillWgt = 1;
  //if( !isLocal ) std::cout << " -- ECT : isLocal " << isLocal << " useSMSWgts : " << useSMSWgts << " use: " << useEvtGenWgtFlag << std::endl;
  // mctype != 1 (not data): both FullSim (0) and FastSim (2) MC get normalized the
  // same way, against whichever tier-specific EventCount.txt entry dataSetKey points at.
  // Only look up EventCount.txt for MC: data dataSetKeys (e.g. "MET18A") are never
  // added to EventCount.txt (skim_submit.py only tracks/checks it for non-data
  // samples), so calling lookupConfigData() unconditionally here would hard-exit
  // every real data job on a "missing" key that was never meant to exist.
  if( mctype != 1 ){
    // Not isLocal (real condor/grid production): configData[dataSetKey] via operator[]
    // would silently default-construct {0.0, 0.0} for a missing key, producing an
    // inf/nan evtFillWgt instead of failing the job where the actual problem is.
    float configSumEvtWgt = isLocal ? configWgts["sumEvtWgt"]
    : useSMSWgts ? lookupConfigData( dataSetKey ).first
    : lookupConfigData( dataSetKey ).second;
    if( !std::isfinite(configSumEvtWgt) || configSumEvtWgt == 0.f ){
      std::cerr << "ERROR: invalid sum of event weights for dataset '" << dataSetKey
      << "': " << configSumEvtWgt << std::endl;
      std::exit(EXIT_FAILURE);
    }//<<>>if( !std::isfinite(configSumEvtWgt) || configSumEvtWgt == 0.f )
    fillWgt = ( ( xsctn * 1000 ) * evtGenWgt  ) / configSumEvtWgt;
    if( !std::isfinite(fillWgt) ){
      std::cerr << "ERROR: non-finite event fill weight for dataset '" << dataSetKey
      << "' (xsec=" << xsctn << ", evtGenWgt=" << evtGenWgt
      << ", sumEvtWgt=" << configSumEvtWgt << ")." << std::endl;
      std::exit(EXIT_FAILURE);
    }//<<>>if( !std::isfinite(fillWgt) )
  }//<<>>if( mctype != 1 )
  //if( !isLocal ) std::cout << " -- ECT : fillWgt " << fillWgt << " " << xsctn << " " << evtGenWgt << " " << configSumEvtWgt << std::endl;

  selEvtVars.fillBranch( "evtFillWgt", fillWgt );
  //std::cout << " -- ECT : fillWgt " << fillWgt << " " << xsctn << " " << evtGenWgt << " " << configSumEvtWgt << std::endl;

  bool BadChargedCandidateFilter = Flag_BadChargedCandidateFilter;
  bool BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
  bool BadPFMuonFilter = Flag_BadPFMuonFilter;
  bool EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool HBHENoiseFilter = Flag_HBHENoiseFilter;
  bool HBHENoiseIsoFilter = Flag_HBHENoiseIsoFilter;
  bool ecalBadCalibFilter = Flag_ecalBadCalibFilter;
  bool eeBadScFilter = Flag_eeBadScFilter;
  bool globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
  bool goodVertices = Flag_goodVertices;
  bool hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;


  bool metFlags = BadChargedCandidateFilter && BadPFMuonFilter && EcalDeadCellTriggerPrimitiveFilter;
  metFlags = metFlags && HBHENoiseFilter && HBHENoiseIsoFilter && ecalBadCalibFilter && eeBadScFilter && goodVertices;

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
  selEvtVars.fillBranch( "Flag_MetFilters", metFlags );//optional

  bool PFMET120_PFMHT120_IDTight = false;
  bool PFMETNoMu120_PFMHTNoMu120_IDTight = false;
  bool PFMET120_PFMHT120_IDTight_PFHT60 = false;
  bool PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = false;

  if( doHTLPathsBase ){

    //Ele35_WPTight_Gsf = HLT_Ele35_WPTight_Gsf_v;
    //Mu55 = HLT_Mu55_v;
    //IsoMu27 = HLT_IsoMu27_v;
    //IsoMu20 = HLT_IsoMu20_v;

    PFMET120_PFMHT120_IDTight = HLT_PFMET120_PFMHT120_IDTight_v;
    PFMETNoMu120_PFMHTNoMu120_IDTight = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
    PFMET120_PFMHT120_IDTight_PFHT60 = HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
    PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;

  }//<<>>if( not doGenInfo )

  selEvtVars.fillBranch( "Trigger_PFMET120_PFMHT120_IDTight", PFMET120_PFMHT120_IDTight );
  selEvtVars.fillBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight", PFMETNoMu120_PFMHTNoMu120_IDTight );
  selEvtVars.fillBranch( "Trigger_PFMET120_PFMHT120_IDTight_PFHT60", PFMET120_PFMHT120_IDTight_PFHT60 );
  selEvtVars.fillBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 );

  //selEvtVars.fillBranch( "Trigger_Ele35_WPTight_Gsf", Ele35_WPTight_Gsf);
  //selEvtVars.fillBranch( "Trigger_Mu55", Mu55 );
  //selEvtVars.fillBranch( "Trigger_IsoMu27", IsoMu27 );
  //selEvtVars.fillBranch( "Trigger_IsoMu20", IsoMu20 );

  bool hemRun = ( ( Evt_run > 319076 ) && ( Evt_run < 326635 ) ) ? true : false;
  bool eLHem( hemBits("el1hvl") );
  bool eMHem( hemBits("el2hvm") );
  bool jLHem( hemBits("jet1hvl") );
  bool jMHem( hemBits("jet2hvm") );
  bool pLHem( hemBits("pho1hvl") );
  bool pMHem( hemBits("pho2hvm") );
  bool mLHem; //( hemBits("mu1hvl") );
  mLHem = doSVs ? hemBits("mu1hvl") : false;
  bool mMHem; //( hemBits("mu2hvm") );
  mMHem = doSVs ? hemBits("mu2hvm") : false;
  bool fullHemVeto( ( eLHem || jLHem || pMHem || mLHem ) && hemRun && ( mctype == 1 ) );

  selEvtVars.fillBranch( "Flag_hemRun", hemRun );
  selEvtVars.fillBranch( "Flag_eLHemVeto", eLHem );
  selEvtVars.fillBranch( "Flag_eMHemVeto", eMHem );
  selEvtVars.fillBranch( "Flag_jLHemVeto", jLHem );
  selEvtVars.fillBranch( "Flag_jMHemVeto", jMHem );
  selEvtVars.fillBranch( "Flag_pLHemVeto", pLHem );
  selEvtVars.fillBranch( "Flag_pMHemVeto", pMHem );
  selEvtVars.fillBranch( "Flag_mLHemVeto", mLHem );
  selEvtVars.fillBranch( "Flag_mMHemVeto", mMHem );
  selEvtVars.fillBranch( "Flag_hemVeto", fullHemVeto );

}//<<>>void KUCMSAodSkimmer::processEvntVars()

void KUCMSAodSkimmer::processMet(){

  //intilize
  selMet.clearBranches(); // <<<<<<<   must do

  // Rebuild the GED/OOT overlap decision from the photons retained in the
  // ntuple.  The input Photon_excluded branch can reflect the last raw OOT
  // candidate rather than an overlap with any retained OOT photon.
  if( Photon_pt == nullptr || Photon_eta == nullptr || Photon_phi == nullptr || Photon_isOot == nullptr ){
    std::cerr << "ERROR: missing photon input needed to rebuild corrected MET." << std::endl;
    std::exit(EXIT_FAILURE);
  }//<<>>if( Photon_pt == nullptr || Photon_eta == nullptr || Photon_phi == nullptr || Photon_isOot == nullptr )

  const std::size_t nPhotons = Photon_pt->size();
  if( Photon_eta->size() != nPhotons || Photon_phi->size() != nPhotons || Photon_isOot->size() != nPhotons ){
    std::cerr << "ERROR: inconsistent photon input sizes while rebuilding corrected MET." << std::endl;
    std::exit(EXIT_FAILURE);
  }//<<>>if( Photon_eta->size() != nPhotons || Photon_phi->size() != nPhotons || Photon_isOot->size() != nPhotons )

  correctedPhotonExcluded.assign( nPhotons, false );
  for( std::size_t gedIdx = 0; gedIdx < nPhotons; gedIdx++ ){
    if( (*Photon_isOot)[gedIdx] ) continue;
    for( std::size_t ootIdx = 0; ootIdx < nPhotons; ootIdx++ ){
      if( not (*Photon_isOot)[ootIdx] ) continue;
      if( dR1( (*Photon_eta)[gedIdx], (*Photon_phi)[gedIdx],
      (*Photon_eta)[ootIdx], (*Photon_phi)[ootIdx] ) < 0.3 ){
        correctedPhotonExcluded[gedIdx] = true;
        break;
      }//<<>>if( dR1( (*Photon_eta)[gedIdx], (*Photon_phi)[gedIdx], (*Photon_eta)[ootIdx], (*Photon_phi)[ootIdx] ) < 0.3 )
    }//<<>>for( std::size_t ootIdx = 0; ootIdx < nPhotons; ootIdx++ )
  }//<<>>for( std::size_t gedIdx = 0; gedIdx < nPhotons; gedIdx++ )

  double correctedMetPx = Met_px;
  double correctedMetPy = Met_py;
  double correctedMetSumEt = Met_sumEt;
  for( std::size_t phoIdx = 0; phoIdx < nPhotons; phoIdx++ ){
    const double pt = (*Photon_pt)[phoIdx];
    const double phi = (*Photon_phi)[phoIdx];
    if( (*Photon_isOot)[phoIdx] ){
      correctedMetPx -= pt * std::cos(phi);
      correctedMetPy -= pt * std::sin(phi);
      correctedMetSumEt += pt;
    } else if( correctedPhotonExcluded[phoIdx] ){
      correctedMetPx += pt * std::cos(phi);
      correctedMetPy += pt * std::sin(phi);
      correctedMetSumEt -= pt;
    }//<<>>else if( correctedPhotonExcluded[phoIdx] )
  }//<<>>for( std::size_t phoIdx = 0; phoIdx < nPhotons; phoIdx++ )
  const double correctedMet = std::hypot( correctedMetPx, correctedMetPy );

  if( DEBUG ){
    std::cout << "Rebuilt corrected MET: pt=" << correctedMet
    << " px=" << correctedMetPx << " py=" << correctedMetPy
    << " sumEt=" << correctedMetSumEt << std::endl;
  }//<<>>if( DEBUG )

  geVars.set("metPx", correctedMetPx );
  geVars.set("metPy", correctedMetPy );
  geVars.set("metBPx", Met_px );
  geVars.set("metBPy", Met_py );
  geVars.set("cmet", correctedMet );

  //fill
  selMet.fillBranch( "cmet", static_cast<float>(correctedMet) );
  selMet.fillBranch( "cmetPx", static_cast<float>(correctedMetPx) );
  selMet.fillBranch( "cmetPy", static_cast<float>(correctedMetPy) );
  selMet.fillBranch( "met", Met_pt );
  selMet.fillBranch( "metPx", Met_px );
  selMet.fillBranch( "metPy", Met_py );

}//<<>>void KUCMSAodSkimmer::processMet()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setEvtVarMetBranches( TTree* fOutTree ){

  std::cout << " - Making Branches for EvtVarMet." << std::endl;
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
  selEvtVars.makeBranch( "Flag_MetFilters", BOOL );//optional

  selEvtVars.makeBranch( "Flag_hemRun", BOOL );
  selEvtVars.makeBranch( "Flag_hemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_eLHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_eMHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_jLHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_jMHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_pLHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_pMHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_mLHemVeto", BOOL );
  selEvtVars.makeBranch( "Flag_mMHemVeto", BOOL );

  selEvtVars.makeBranch( "Trigger_PFMET120_PFMHT120_IDTight", BOOL );
  selEvtVars.makeBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight", BOOL );
  selEvtVars.makeBranch( "Trigger_PFMET120_PFMHT120_IDTight_PFHT60", BOOL );
  selEvtVars.makeBranch( "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", BOOL );
  selEvtVars.makeBranch( "Trigger_Ele35_WPTight_Gsf", BOOL );
  selEvtVars.makeBranch( "Trigger_Photon20", BOOL );
  selEvtVars.makeBranch( "Trigger_Mu55", BOOL );
  selEvtVars.makeBranch( "Trigger_Mu12", BOOL );
  selEvtVars.makeBranch( "Trigger_IsoMu27", BOOL );
  selEvtVars.makeBranch( "Trigger_IsoMu20", BOOL );
  selEvtVars.makeBranch( "Trigger_Ele27_WPTight_Gsf", BOOL );

  selEvtVars.attachBranches( fOutTree );

  selMet.makeBranch( "cmet", "selCMet", FLOAT, "Magnitude of event Met corrected for OOT photons" );
  selMet.makeBranch( "cmetPx", "selCMetPx", FLOAT, "Magnitude of event MetPx corrected for OOT photons" );
  selMet.makeBranch( "cmetPy", "selCMetPy", FLOAT, "Magnitude of event MetPy corrected for OOT photons" );
  selMet.makeBranch( "met", "selMet", FLOAT, "Magnitude of event Met" );
  selMet.makeBranch( "metPx", "selMetPx", FLOAT, "Magnitude of event MetPx" );
  selMet.makeBranch( "metPy", "selMetPy", FLOAT, "Magnitude of event MetPy" );

  selMet.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setEvtVarMetBranches( TTree* fOutTree )
