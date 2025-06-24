//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 28 14:33:01 2024 by ROOT version 6.26/11
// from TTree llpgtree/KUCMSNtuple
// found on file: gmsb_AODSIM_KUCMSNtuplizer_justin_test_fullgen_oottrue.root
//////////////////////////////////////////////////////////

#ifndef llpgtree_h
#define llpgtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

//#define sbDEBUG true
#define sbDEBUG false

class llpgtree {

public :

   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<int>     *DisplacedElectron_charge;
   Float_t         DisplacedElectron_costGen;
   Float_t         DisplacedElectron_costSC;
   std::vector<float>   *DisplacedElectron_dRGen;
   std::vector<float>   *DisplacedElectron_dRSC;
   std::vector<float>   *DisplacedElectron_dxy;
   std::vector<float>   *DisplacedElectron_energy;
   std::vector<float>   *DisplacedElectron_eta;
   std::vector<float>   *DisplacedElectron_etaECAL;
   std::vector<int>     *DisplacedElectron_genSigIndex;
   std::vector<int>     *DisplacedElectron_indexGen;
   std::vector<int>     *DisplacedElectron_indexSC;
   std::vector<int>     *DisplacedElectron_indexTrack;
   std::vector<float>   *DisplacedElectron_ip2D;
   std::vector<float>   *DisplacedElectron_ip2Dbs;
   std::vector<float>   *DisplacedElectron_ip3D;
   std::vector<float>   *DisplacedElectron_ip3Dbs;
   std::vector<bool>    *DisplacedElectron_isBarrel;
   std::vector<bool>    *DisplacedElectron_isConversion;
   std::vector<bool>    *DisplacedElectron_isEndCap;
   std::vector<bool>    *DisplacedElectron_isGeneralTrack;
   std::vector<bool>    *DisplacedElectron_isGsfTrack;
   std::vector<bool>    *DisplacedElectron_isHeavy;
   std::vector<bool>    *DisplacedElectron_isLight;
   std::vector<bool>    *DisplacedElectron_isOOTSC;
   std::vector<bool>    *DisplacedElectron_isSC;
   std::vector<bool>    *DisplacedElectron_isSignal;
   std::vector<bool>    *DisplacedElectron_isTau;
   std::vector<bool>    *DisplacedElectron_isUnmatched;
   Int_t           DisplacedElectron_nGenSig;
   std::vector<int>     *DisplacedElectron_nHits;
   Int_t           DisplacedElectron_nOOTSCs;
   Int_t           DisplacedElectron_nSCs;
   Int_t           DisplacedElectron_nTotal;
   std::vector<float>   *DisplacedElectron_p;
   std::vector<int>     *DisplacedElectron_pdgId;
   std::vector<float>   *DisplacedElectron_pfIso0p4;
   std::vector<float>   *DisplacedElectron_phi;
   std::vector<float>   *DisplacedElectron_phiECAL;
   std::vector<float>   *DisplacedElectron_pt;
   std::vector<float>   *DisplacedElectron_scIso0p4;
   std::vector<float>   *DisplacedElectron_sip2D;
   std::vector<float>   *DisplacedElectron_sip2Dbs;
   std::vector<float>   *DisplacedElectron_sip3D;
   std::vector<float>   *DisplacedElectron_sip3Dbs;
   std::vector<float>   *DisplacedElectron_trackIso0p4;
   std::vector<float>   *DisplacedElectron_x;
   std::vector<float>   *DisplacedElectron_xECAL;
   std::vector<float>   *DisplacedElectron_y;
   std::vector<float>   *DisplacedElectron_yECAL;
   std::vector<float>   *DisplacedElectron_z;
   std::vector<float>   *DisplacedElectron_zECAL;
   std::vector<float>   *ECALRecHit_energy;
   std::vector<float>   *ECALRecHit_timeError;
   std::vector<bool>    *ECALRecHit_hasGS1;
   std::vector<bool>    *ECALRecHit_hasGS6;
   std::vector<unsigned int> *ECALRecHit_ID;
   std::vector<float>   *ECALRecHit_swCross;
   std::vector<float>   *ECALRecHit_0TOF;
   std::vector<float>   *ECALRecHit_pvTOF;
   std::vector<float>   *ECALRecHit_time;
   std::vector<float>   *ECALRecHit_amplitude;
   std::vector<float>   *ECALRecHit_ampres;
   std::vector<float>   *ECALRecHit_eta;
   std::vector<bool>    *ECALRecHit_isDead;
   std::vector<bool>    *ECALRecHit_isDiWrd;
   std::vector<bool>    *ECALRecHit_isOOT;
   std::vector<bool>    *ECALRecHit_isOther;
   std::vector<bool>    *ECALRecHit_isPoor;
   std::vector<bool>    *ECALRecHit_isRecov;
   std::vector<bool>    *ECALRecHit_isWrd;
   std::vector<float>   *ECALRecHit_phi;
   std::vector<float>   *ECALRecHit_precentUsed;
   std::vector<float>   *ECALRecHit_rhx;
   std::vector<float>   *ECALRecHit_rhy;
   std::vector<float>   *ECALRecHit_rhz;
   std::vector<float>   *SuperCluster_covEtaEta;
   std::vector<float>   *SuperCluster_covEtaPhi;
   std::vector<float>   *SuperCluster_covPhiPhi;
   std::vector<float>   *SuperCluster_energyRaw;
   std::vector<float>   *SuperCluster_etaWidth;
   std::vector<bool>    *SuperCluster_excluded;
   std::vector<bool>    *SuperCluster_seedIsEB;
   std::vector<bool>    *SuperCluster_isScEtaEB;
   std::vector<bool>    *SuperCluster_isScEtaEE;
   std::vector<bool>    *SuperCluster_isOot;
   std::vector<bool>    *SuperCluster_original;
   std::vector<float>   *SuperCluster_phiWidth;
   std::vector<float>   *SuperCluster_salp;
   std::vector<float>   *SuperCluster_smaj;
   std::vector<float>   *SuperCluster_smin;
   std::vector<int>     *SuperCluster_seediEtaOriX;
   std::vector<int>     *SuperCluster_seediPhiOriY;
   std::vector<unsigned int> *SuperCluster_nBasicClusters;
   std::vector<int>     *SuperCluster_ElectronIndx;
   std::vector<float>   *SuperCluster_energy;
   std::vector<float>   *SuperCluster_eta;
   std::vector<float>   *SuperCluster_clcx;
   std::vector<float>   *SuperCluster_clcy;
   std::vector<float>   *SuperCluster_clcz;
   Int_t           SuperCluster_nSuperCluster;
   std::vector<unsigned int> *SuperCluster_nRHXtals;
   std::vector<int>     *SuperCluster_ObjectPdgId;
   std::vector<float>   *SuperCluster_phi;
   std::vector<int>     *SuperCluster_PhotonIndx;
   std::vector<std::vector<float> > *SuperCluster_rhFracs;
   std::vector<std::vector<unsigned int> > *SuperCluster_rhIds;
   std::vector<unsigned int> *SuperCluster_XtalSeedID;
   std::vector<unsigned int> *SuperCluster_nHFXtals;
   std::vector<float>   *SuperCluster_x_calo;
   std::vector<float>   *SuperCluster_y_calo;
   std::vector<float>   *SuperCluster_z_calo;
   std::vector<unsigned int> *ECALTrackDetID_detId;
   std::vector<bool>    *ECALTrackDetID_isECAL;
   std::vector<bool>    *ECALTrackDetID_isHCAL;
   std::vector<unsigned int> *ECALTrackDetID_trackIndex;
   std::vector<int>     *ECALTrack_charge;
   std::vector<float>   *ECALTrack_chi2;
   std::vector<float>   *ECALTrack_dxy;
   std::vector<float>   *ECALTrack_eta;
   std::vector<float>   *ECALTrack_etaError;
   std::vector<float>   *ECALTrack_genDeltaR;
   std::vector<int>     *ECALTrack_genIndex;
   std::vector<float>   *ECALTrack_lambda;
   Int_t           ECALTrack_nGeneral;
   Int_t           ECALTrack_nGsf;
   std::vector<int>     *ECALTrack_nLostHits;
   Int_t           ECALTrack_nTracks;
   std::vector<int>     *ECALTrack_nValidHits;
   std::vector<float>   *ECALTrack_ndof;
   std::vector<float>   *ECALTrack_normalizedChi2;
   std::vector<float>   *ECALTrack_p;
   std::vector<int>     *ECALTrack_pdgId;
   std::vector<float>   *ECALTrack_phi;
   std::vector<float>   *ECALTrack_phiError;
   std::vector<float>   *ECALTrack_pt;
   std::vector<float>   *ECALTrack_ptError;
   std::vector<float>   *ECALTrack_px;
   std::vector<float>   *ECALTrack_py;
   std::vector<float>   *ECALTrack_pz;
   std::vector<float>   *ECALTrack_qoverp;
   std::vector<int>     *ECALTrack_qualityMask;
   std::vector<float>   *ECALTrack_theta;
   std::vector<float>   *ECALTrack_vx;
   std::vector<float>   *ECALTrack_vy;
   std::vector<float>   *ECALTrack_vz;
   std::vector<float>   *TrackSCMatch_deltaR;
   std::vector<float>   *TrackSCMatch_scEnergyToTrackPRatio;
   std::vector<float>   *Electron_energy;
   std::vector<float>   *Electron_eta;
   std::vector<int>     *Electron_genIdx;
   std::vector<int>     *Electron_genSigWZId;
   std::vector<int>     *Electron_genSigXMomId;
   std::vector<float>   *Electron_phi;
   std::vector<float>   *Electron_pt;
   std::vector<float>   *Electron_px;
   std::vector<float>   *Electron_py;
   std::vector<float>   *Electron_pz;
   std::vector<float>   *Electron_seedTOFTime;
   std::vector<float>   *Electron_trackz;
   std::vector<bool>   *Electron_isLoose;
   Int_t           Electron_nSelElectrons;
   Int_t           Electron_nElectrons;
   std::vector<int>     *Electron_scIndex;

   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_hfNoisyHitsFilter;
   UInt_t          Evt_luminosityBlock;
   UInt_t          Evt_run;
   UInt_t          Evt_event;
   Bool_t          Trigger_hltEG20EtFilterLooseHoverE;
   Bool_t          Trigger_hltEG20EtL1TripleEGFilter;
   Bool_t          Trigger_hltEG20HEFilterLooseHoverE;
   Bool_t          Trigger_hltEG20HEL1TripleEGFilter;
   Bool_t          Trigger_hltEG20L1SingleEGLowETOrEtFilter;
   Bool_t          Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter;
   Bool_t          Trigger_hltEG30EBHE10R9Id50b80eHEFilter;
   Bool_t          Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter;
   Bool_t          Trigger_hltEG30EBR9Id50b80eR9IdFilter;
   Bool_t          Trigger_hltEG30EIso15HE30EcalIsoLastFilter;
   Bool_t          Trigger_hltEG30EtFilterLooseHoverE;
   Bool_t          Trigger_hltEG30EtL1TripleEGFilter;
   Bool_t          Trigger_hltEG30HE30HEFilter;
   Bool_t          Trigger_hltEG30HEFilterLooseHoverE;
   Bool_t          Trigger_hltEG30HEL1TripleEGFilter;
   Bool_t          Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter;
   Bool_t          Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter;
   Bool_t          Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter;
   Bool_t          Trigger_hltEG30LHE12R9Id50b80eHEFilter;
   Bool_t          Trigger_hltEG30LR9Id50b80eR9IdFilter;
   Bool_t          Trigger_hltEG30PVHE10R9Id50b80eHEFilter;
   Bool_t          Trigger_hltEG30PVR9Id50b80eR9IdFilter;
   Bool_t          Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter;
   Bool_t          Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter;
   Bool_t          Trigger_hltEG60EtFilter;
   Bool_t          Trigger_hltEG60HEFilter;
   Bool_t          Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter;
   Bool_t          Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter;
   Bool_t          Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter;
   Bool_t          Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter;
   Bool_t          Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter;
   Bool_t          Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter;
   Bool_t          Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter;
   Bool_t          Trigger_hltHT130Jet30;
   Bool_t          Trigger_hltHT175Jet10;
   Bool_t          Trigger_hltHT200Jet30;
   Bool_t          Trigger_hltL1sSingleAndDoubleEG;
   Bool_t          Trigger_hltL1sSingleAndDoubleEGNonIsoOr;
   Bool_t          Trigger_hltL1sSingleAndDoubleEGor;
   Bool_t          Trigger_hltL1sSingleEG15;
   Bool_t          Trigger_hltL1sSingleEG18;
   Bool_t          Trigger_hltL1sSingleEG24;
   Bool_t          Trigger_hltL1sSingleEG26;
   Bool_t          Trigger_hltL1sSingleEG34to45;
   Bool_t          Trigger_hltL1sSingleEG34to50;
   Bool_t          Trigger_hltL1sSingleEG40to50;
   Bool_t          Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS;
   Bool_t          Trigger_hltL1sSingleEGor;
   Bool_t          Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG;
   Bool_t          Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter;
   Bool_t          Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt;
   Bool_t          Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt;
   Bool_t          Trigger_hltPFHT180Jet30;
   Bool_t          Trigger_hltPFHT250Jet30;
   Bool_t          Trigger_hltPFHT350Jet15;
   Bool_t          Trigger_hltPFMET100;
   Bool_t          Trigger_hltPFMET110;
   Bool_t          Trigger_hltPFMET120;
   Bool_t          Trigger_hltPFMET130;
   Bool_t          Trigger_hltPFMET140;
   Bool_t          Trigger_hltPFMET200;
   Bool_t          Trigger_hltPFMET250;
   Bool_t          Trigger_hltPFMET300;
   Bool_t          Trigger_hltPFMET50;
   Bool_t          Trigger_hltPFMET70;
   Bool_t          Trigger_hltPFMET90;
   Bool_t          Trigger_hltPFMETNoMu100;
   Bool_t          Trigger_hltPFMETNoMu110;
   Bool_t          Trigger_hltPFMETNoMu120;
   Bool_t          Trigger_hltPFMETNoMu130;
   Bool_t          Trigger_hltPFMETNoMu140;
   Bool_t          Trigger_hltPFMETNoMu60;
   Bool_t          Trigger_hltPFMHTNoMuTightID100;
   Bool_t          Trigger_hltPFMHTNoMuTightID110;
   Bool_t          Trigger_hltPFMHTNoMuTightID120;
   Bool_t          Trigger_hltPFMHTNoMuTightID130;
   Bool_t          Trigger_hltPFMHTNoMuTightID140;
   Bool_t          Trigger_hltPFMHTNoMuTightID70;
   Bool_t          Trigger_hltPFMHTNoMuTightID90;
   Bool_t          Trigger_hltPFMHTTightID120;
   Bool_t          Trigger_hltPFMHTTightID130;
   Bool_t          Trigger_hltPFMHTTightID140;

   UInt_t          PV_npvs;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;

   std::vector<int>     *Gen_charge;
   std::vector<float>   *Gen_energy;
   std::vector<float>   *Gen_eta;
   std::vector<float>   *Gen_mass;
   std::vector<float>   *Gen_momDisplacment;
   std::vector<int>     *Gen_motherIdx;
   UInt_t          Gen_nTotal;
   std::vector<unsigned int> *Gen_pdgId;
   std::vector<float>   *Gen_phi;
   std::vector<float>   *Gen_pt;
   std::vector<float>   *Gen_px;
   std::vector<float>   *Gen_py;
   std::vector<float>   *Gen_pz;
   std::vector<int>     *Gen_status;
   std::vector<int>     *Gen_susEvtType;
   std::vector<int>     *Gen_susId;
   std::vector<float>   *Gen_vx;
   std::vector<float>   *Gen_vy;
   std::vector<float>   *Gen_vz;
   Float_t         Evt_genWgt;

   std::vector<float>   *Muon_energy;
   std::vector<float>   *Muon_eta;
   std::vector<float>   *Muon_phi;
   std::vector<float>   *Muon_pt;
   std::vector<float>   *Muon_px;
   std::vector<float>   *Muon_py;
   std::vector<float>   *Muon_pz;
   std::vector<bool>    *Muon_isLoose;
   //std::vector<bool>    *Muon_isMedium;
   Int_t           Muon_nMuons;
   Int_t           Muon_nSelMuons;

   std::vector<float>   *Jet_area;
   std::vector<float>   *Jet_chEmEF;
   std::vector<float>   *Jet_chHEF;
   std::vector<float>   *Jet_chHM;
   std::vector<std::vector<unsigned int> > *Jet_drRhIds;
   std::vector<float>   *Jet_energy;
   std::vector<float>   *Jet_eta;
   std::vector<float>   *Jet_genJetLlpDr;
   std::vector<float>   *Jet_genDptMatch;
   std::vector<float>   *Jet_genDrMatch;
   std::vector<float>   *Jet_genEnergy;
   std::vector<float>   *Jet_genEta;
   std::vector<float>   *Jet_genImpactAngle;
   std::vector<int>     *Jet_genJetLlpId;
   std::vector<float>   *Jet_genJetLlpRp;
   std::vector<float>   *Jet_genPhi;
   std::vector<float>   *Jet_genPt;
   std::vector<float>   *Jet_genQrkLlpDr;
   std::vector<int>     *Jet_genQrkLlpId;
   std::vector<float>   *Jet_genQrkLlpRp;
   std::vector<float>   *Jet_genTOF;
   std::vector<float>   *Jet_genTime;
   std::vector<float>   *Jet_genTimeLLP;
   std::vector<float>   *Jet_mass;
   std::vector<float>   *Jet_muEF;
   std::vector<float>   *Jet_neEmEF;
   std::vector<float>   *Jet_neHEF;
   std::vector<float>   *Jet_neHM;
   std::vector<std::vector<unsigned int> > *Jet_egIndxs;
   std::vector<float>   *Jet_phi;
   std::vector<float>   *Jet_pt;
   std::vector<int>     *Jet_nConstituents;

   Float_t         Met_CPt;
   Float_t         Met_Cpx;
   Float_t         Met_Cpy;
   Float_t         Met_CsumEt;
   Float_t         Met_eta;
   Float_t         Met_phi;
   Float_t         Met_pt;
   Float_t         Met_px;
   Float_t         Met_py;
   Float_t         Met_sumEt;

   std::vector<float>   *Photon_ecalRHSumEtConeDR04;
   std::vector<float>   *Photon_energy;
   std::vector<float>   *Photon_energyErr;
   std::vector<float>   *Photon_eta;
   std::vector<bool>    *Photon_excluded;
   std::vector<int>     *Photon_genIdx;
   std::vector<int>     *Photon_genSigXMomId;
   std::vector<float>   *Photon_hadOverEM;
   std::vector<float>   *Photon_hadTowOverEM;
   std::vector<float>   *Photon_hcalTowerSumEtBcConeDR04;
   std::vector<float>   *Photon_hcalTowerSumEtConeDR04;
   std::vector<bool>    *Photon_isOot;
   std::vector<float>   *Photon_nTrkHollowConeDR04;
   std::vector<float>   *Photon_nTrkSolidConeDR04;
   std::vector<float>   *Photon_phi;
   std::vector<float>   *Photon_pt;
   std::vector<float>   *Photon_px;
   std::vector<float>   *Photon_py;
   std::vector<float>   *Photon_pz;
   std::vector<float>   *Photon_r9;
   std::vector<float>   *Photon_s4;
   std::vector<float>   *Photon_seedTOFTime;
   std::vector<float>   *Photon_SigmaIEtaIEta;
   std::vector<float>   *Photon_trkSumPtHollowConeDR03;
   std::vector<float>   *Photon_trkSumPtHollowConeDR04;
   std::vector<float>   *Photon_trkSumPtSolidConeDR04;
   std::vector<float>   *Photon_ecalPFClusterIso;
   std::vector<bool>    *Photon_electronVeto;
   std::vector<std::vector<unsigned int> > *Photon_gloResRhId;
   std::vector<bool>    *hasConversionTracks;
   std::vector<bool>    *Photon_pixelSeed;
   std::vector<float>   *Photon_hcalPFClusterIso;
   std::vector<float>   *Photon_Hoe_PUcorr;
   std::vector<std::vector<unsigned int> > *Photon_locResRhId;
   Int_t           Photon_nPhoton;
   std::vector<float>   *Photon_pfChargedIso;
   std::vector<float>   *Photon_pfChargedIsoPFPV;
   std::vector<float>   *Photon_pfChargedIsoWorstVtx;
   std::vector<float>   *Photon_pfPhoIso03;
   std::vector<float>   *pfRelIso03_all_quadratic;
   std::vector<float>   *pfRelIso03_chg_quadratic;
   std::vector<int>     *Photon_scIndex;
   std::vector<float>   *Photon_sieie;
   std::vector<float>   *Photon_sieip;
   std::vector<float>   *Photon_sipip;

   std::vector<unsigned int> *Vertex_nTracks;
   std::vector<float>   *Vertex_mass;
   std::vector<bool>    *Vertex_passLooseElectronID;
   std::vector<bool>    *Vertex_passLooseMuonID;

   // List of branches
   TBranch        *b_DisplacedElectron_charge;   //!
   TBranch        *b_DisplacedElectron_costGen;   //!
   TBranch        *b_DisplacedElectron_costSC;   //!
   TBranch        *b_DisplacedElectron_dRGen;   //!
   TBranch        *b_DisplacedElectron_dRSC;   //!
   TBranch        *b_DisplacedElectron_dxy;   //!
   TBranch        *b_DisplacedElectron_energy;   //!
   TBranch        *b_DisplacedElectron_eta;   //!
   TBranch        *b_DisplacedElectron_etaECAL;   //!
   TBranch        *b_DisplacedElectron_genSigIndex;   //!
   TBranch        *b_DisplacedElectron_indexGen;   //!
   TBranch        *b_DisplacedElectron_indexSC;   //!
   TBranch        *b_DisplacedElectron_indexTrack;   //!
   TBranch        *b_DisplacedElectron_ip2D;   //!
   TBranch        *b_DisplacedElectron_ip2Dbs;   //!
   TBranch        *b_DisplacedElectron_ip3D;   //!
   TBranch        *b_DisplacedElectron_ip3Dbs;   //!
   TBranch        *b_DisplacedElectron_isBarrel;   //!
   TBranch        *b_DisplacedElectron_isConversion;   //!
   TBranch        *b_DisplacedElectron_isEndCap;   //!
   TBranch        *b_DisplacedElectron_isGeneralTrack;   //!
   TBranch        *b_DisplacedElectron_isGsfTrack;   //!
   TBranch        *b_DisplacedElectron_isHeavy;   //!
   TBranch        *b_DisplacedElectron_isLight;   //!
   TBranch        *b_DisplacedElectron_isOOTSC;   //!
   TBranch        *b_DisplacedElectron_isSC;   //!
   TBranch        *b_DisplacedElectron_isSignal;   //!
   TBranch        *b_DisplacedElectron_isTau;   //!
   TBranch        *b_DisplacedElectron_isUnmatched;   //!
   TBranch        *b_DisplacedElectron_nGenSig;   //!
   TBranch        *b_DisplacedElectron_nHits;   //!
   TBranch        *b_DisplacedElectron_nOOTSCs;   //!
   TBranch        *b_DisplacedElectron_nSCs;   //!
   TBranch        *b_DisplacedElectron_nTotal;   //!
   TBranch        *b_DisplacedElectron_p;   //!
   TBranch        *b_DisplacedElectron_pdgId;   //!
   TBranch        *b_DisplacedElectron_pfIso0p4;   //!
   TBranch        *b_DisplacedElectron_phi;   //!
   TBranch        *b_DisplacedElectron_phiECAL;   //!
   TBranch        *b_DisplacedElectron_pt;   //!
   TBranch        *b_DisplacedElectron_scIso0p4;   //!
   TBranch        *b_DisplacedElectron_sip2D;   //!
   TBranch        *b_DisplacedElectron_sip2Dbs;   //!
   TBranch        *b_DisplacedElectron_sip3D;   //!
   TBranch        *b_DisplacedElectron_sip3Dbs;   //!
   TBranch        *b_DisplacedElectron_trackIso0p4;   //!
   TBranch        *b_DisplacedElectron_x;   //!
   TBranch        *b_DisplacedElectron_xECAL;   //!
   TBranch        *b_DisplacedElectron_y;   //!
   TBranch        *b_DisplacedElectron_yECAL;   //!
   TBranch        *b_DisplacedElectron_z;   //!
   TBranch        *b_DisplacedElectron_zECAL;   //!

   TBranch        *b_ECALRecHit_energy;   //!
   TBranch        *b_ECALRecHit_timeError;
   TBranch        *b_ECALRecHit_hasGS1;
   TBranch        *b_ECALRecHit_hasGS6;
   TBranch        *b_ECALRecHit_ID;   //!
   TBranch        *b_ECALRecHit_swCross;   //!
   TBranch        *b_ECALRecHit_0TOF;   //!
   TBranch        *b_ECALRecHit_pvTOF;   //!
   TBranch        *b_ECALRecHit_time;   //!
   TBranch        *b_ECALRecHit_amplitude;   //!
   TBranch        *b_ECALRecHit_ampres;   //!
   TBranch        *b_ECALRecHit_eta;   //!
   TBranch        *b_ECALRecHit_isDead;   //!
   TBranch        *b_ECALRecHit_isDiWrd;   //!
   TBranch        *b_ECALRecHit_isOOT;   //!
   TBranch        *b_ECALRecHit_isOther;   //!
   TBranch        *b_ECALRecHit_isPoor;   //!
   TBranch        *b_ECALRecHit_isRecov;   //!
   TBranch        *b_ECALRecHit_isWrd;   //!
   TBranch        *b_ECALRecHit_phi;   //!
   TBranch        *b_ECALRecHit_precentUsed;   //!
   TBranch        *b_ECALRecHit_rhx;   //!
   TBranch        *b_ECALRecHit_rhy;   //!
   TBranch        *b_ECALRecHit_rhz;   //!
   TBranch        *b_SuperCluster_covEtaEta;   //!
   TBranch        *b_SuperCluster_covEtaPhi;   //!
   TBranch        *b_SuperCluster_covPhiPhi;   //!
   TBranch        *b_SuperCluster_energyRaw;   //!
   TBranch        *b_SuperCluster_etaWidth;   //!
   TBranch        *b_SuperCluster_excluded;   //!
   TBranch        *b_SuperCluster_seedIsEB;   //!
   TBranch        *b_SuperCluster_isScEtaEB;   //!
   TBranch        *b_SuperCluster_isScEtaEE;   //!
   TBranch        *b_SuperCluster_isOot;   //!
   TBranch        *b_SuperCluster_original;   //!
   TBranch        *b_SuperCluster_phiWidth;   //!
   TBranch        *b_SuperCluster_salp;   //!
   TBranch        *b_SuperCluster_smaj;   //!
   TBranch        *b_SuperCluster_smin;   //!
   TBranch        *b_SuperCluster_seediEtaOriX;   //!
   TBranch        *b_SuperCluster_seediPhiOriY;   //!
   TBranch        *b_SuperCluster_nBasicClusters;   //!
   TBranch        *b_SuperCluster_ElectronIndx;   //!
   TBranch        *b_SuperCluster_energy;   //!
   TBranch        *b_SuperCluster_eta;   //!
   TBranch        *b_SuperCluster_clcx;   //!
   TBranch        *b_SuperCluster_clcy;   //!
   TBranch        *b_SuperCluster_clcz;   //!
   TBranch        *b_SuperCluster_nSuperCluster;   //!
   TBranch        *b_SuperCluster_nRHXtals;   //!
   TBranch        *b_SuperCluster_ObjectPdgId;   //!
   TBranch        *b_SuperCluster_phi;   //!
   TBranch        *b_SuperCluster_PhotonIndx;   //!
   TBranch        *b_SuperCluster_rhFracs;   //!
   TBranch        *b_SuperCluster_rhIds;   //!
   TBranch        *b_SuperCluster_XtalSeedID;   //!
   TBranch        *b_SuperCluster_nHFXtals;   //!
   TBranch        *b_SuperCluster_x_calo;   //!
   TBranch        *b_SuperCluster_y_calo;   //!
   TBranch        *b_SuperCluster_z_calo;   //!
   TBranch        *b_ECALTrackDetID_detId;   //!
   TBranch        *b_ECALTrackDetID_isECAL;   //!
   TBranch        *b_ECALTrackDetID_isHCAL;   //!
   TBranch        *b_ECALTrackDetID_trackIndex;   //!
   TBranch        *b_ECALTrack_charge;   //!
   TBranch        *b_ECALTrack_chi2;   //!
   TBranch        *b_ECALTrack_dxy;   //!
   TBranch        *b_ECALTrack_eta;   //!
   TBranch        *b_ECALTrack_etaError;   //!
   TBranch        *b_ECALTrack_genDeltaR;   //!
   TBranch        *b_ECALTrack_genIndex;   //!
   TBranch        *b_ECALTrack_lambda;   //!
   TBranch        *b_ECALTrack_nGeneral;   //!
   TBranch        *b_ECALTrack_nGsf;   //!
   TBranch        *b_ECALTrack_nLostHits;   //!
   TBranch        *b_ECALTrack_nTracks;   //!
   TBranch        *b_ECALTrack_nValidHits;   //!
   TBranch        *b_ECALTrack_ndof;   //!
   TBranch        *b_ECALTrack_normalizedChi2;   //!
   TBranch        *b_ECALTrack_p;   //!
   TBranch        *b_ECALTrack_pdgId;   //!
   TBranch        *b_ECALTrack_phi;   //!
   TBranch        *b_ECALTrack_phiError;   //!
   TBranch        *b_ECALTrack_pt;   //!
   TBranch        *b_ECALTrack_ptError;   //!
   TBranch        *b_ECALTrack_px;   //!
   TBranch        *b_ECALTrack_py;   //!
   TBranch        *b_ECALTrack_pz;   //!
   TBranch        *b_ECALTrack_qoverp;   //!
   TBranch        *b_ECALTrack_qualityMask;   //!
   TBranch        *b_ECALTrack_theta;   //!
   TBranch        *b_ECALTrack_vx;   //!
   TBranch        *b_ECALTrack_vy;   //!
   TBranch        *b_ECALTrack_vz;   //!
   TBranch        *b_TrackSCMatch_deltaR;   //!
   TBranch        *b_TrackSCMatch_scEnergyToTrackPRatio;   //!
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_genIdx;   //!
   TBranch        *b_Electron_genSigWZId;   //!
   TBranch        *b_Electron_genSigXMomId;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_px;   //!
   TBranch        *b_Electron_py;   //!
   TBranch        *b_Electron_pz;   //!
   TBranch        *b_Electron_seedTOFTime;   //!
   TBranch        *b_Electron_trackz;   //!
   TBranch        *b_Electron_nElectrons;   //!
   TBranch        *b_Electron_scIndex;   //!
   TBranch        *b_Electron_isLoose;   //!
   TBranch        *b_Electron_nSelElectrons;   //!

   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_Evt_luminosityBlock;   //!
   TBranch        *b_Evt_run;   //!
   TBranch        *b_Evt_event;   //!
   TBranch        *b_Trigger_hltEG20EtFilterLooseHoverE;   //!
   TBranch        *b_Trigger_hltEG20EtL1TripleEGFilter;   //!
   TBranch        *b_Trigger_hltEG20HEFilterLooseHoverE;   //!
   TBranch        *b_Trigger_hltEG20HEL1TripleEGFilter;   //!
   TBranch        *b_Trigger_hltEG20L1SingleEGLowETOrEtFilter;   //!
   TBranch        *b_Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter;   //!
   TBranch        *b_Trigger_hltEG30EBHE10R9Id50b80eHEFilter;   //!
   TBranch        *b_Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter;   //!
   TBranch        *b_Trigger_hltEG30EBR9Id50b80eR9IdFilter;   //!
   TBranch        *b_Trigger_hltEG30EIso15HE30EcalIsoLastFilter;   //!
   TBranch        *b_Trigger_hltEG30EtFilterLooseHoverE;   //!
   TBranch        *b_Trigger_hltEG30EtL1TripleEGFilter;   //!
   TBranch        *b_Trigger_hltEG30HE30HEFilter;   //!
   TBranch        *b_Trigger_hltEG30HEFilterLooseHoverE;   //!
   TBranch        *b_Trigger_hltEG30HEL1TripleEGFilter;   //!
   TBranch        *b_Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter;   //!
   TBranch        *b_Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter;   //!
   TBranch        *b_Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter;   //!
   TBranch        *b_Trigger_hltEG30LHE12R9Id50b80eHEFilter;   //!
   TBranch        *b_Trigger_hltEG30LR9Id50b80eR9IdFilter;   //!
   TBranch        *b_Trigger_hltEG30PVHE10R9Id50b80eHEFilter;   //!
   TBranch        *b_Trigger_hltEG30PVR9Id50b80eR9IdFilter;   //!
   TBranch        *b_Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter;   //!
   TBranch        *b_Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter;   //!
   TBranch        *b_Trigger_hltEG60EtFilter;   //!
   TBranch        *b_Trigger_hltEG60HEFilter;   //!
   TBranch        *b_Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter;   //!
   TBranch        *b_Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter;   //!
   TBranch        *b_Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter;   //!
   TBranch        *b_Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter;   //!
   TBranch        *b_Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter;   //!
   TBranch        *b_Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter;   //!
   TBranch        *b_Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter;   //!
   TBranch        *b_Trigger_hltHT130Jet30;   //!
   TBranch        *b_Trigger_hltHT175Jet10;   //!
   TBranch        *b_Trigger_hltHT200Jet30;   //!
   TBranch        *b_Trigger_hltL1sSingleAndDoubleEG;   //!
   TBranch        *b_Trigger_hltL1sSingleAndDoubleEGNonIsoOr;   //!
   TBranch        *b_Trigger_hltL1sSingleAndDoubleEGor;   //!
   TBranch        *b_Trigger_hltL1sSingleEG15;   //!
   TBranch        *b_Trigger_hltL1sSingleEG18;   //!
   TBranch        *b_Trigger_hltL1sSingleEG24;   //!
   TBranch        *b_Trigger_hltL1sSingleEG26;   //!
   TBranch        *b_Trigger_hltL1sSingleEG34to45;   //!
   TBranch        *b_Trigger_hltL1sSingleEG34to50;   //!
   TBranch        *b_Trigger_hltL1sSingleEG40to50;   //!
   TBranch        *b_Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS;   //!
   TBranch        *b_Trigger_hltL1sSingleEGor;   //!
   TBranch        *b_Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG;   //!
   TBranch        *b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter;   //!
   TBranch        *b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt;   //!
   TBranch        *b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt;   //!
   TBranch        *b_Trigger_hltPFHT180Jet30;   //!
   TBranch        *b_Trigger_hltPFHT250Jet30;   //!
   TBranch        *b_Trigger_hltPFHT350Jet15;   //!
   TBranch        *b_Trigger_hltPFMET100;   //!
   TBranch        *b_Trigger_hltPFMET110;   //!
   TBranch        *b_Trigger_hltPFMET120;   //!
   TBranch        *b_Trigger_hltPFMET130;   //!
   TBranch        *b_Trigger_hltPFMET140;   //!
   TBranch        *b_Trigger_hltPFMET200;   //!
   TBranch        *b_Trigger_hltPFMET250;   //!
   TBranch        *b_Trigger_hltPFMET300;   //!
   TBranch        *b_Trigger_hltPFMET50;   //!
   TBranch        *b_Trigger_hltPFMET70;   //!
   TBranch        *b_Trigger_hltPFMET90;   //!
   TBranch        *b_Trigger_hltPFMETNoMu100;   //!
   TBranch        *b_Trigger_hltPFMETNoMu110;   //!
   TBranch        *b_Trigger_hltPFMETNoMu120;   //!
   TBranch        *b_Trigger_hltPFMETNoMu130;   //!
   TBranch        *b_Trigger_hltPFMETNoMu140;   //!
   TBranch        *b_Trigger_hltPFMETNoMu60;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID100;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID110;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID120;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID130;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID140;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID70;   //!
   TBranch        *b_Trigger_hltPFMHTNoMuTightID90;   //!
   TBranch        *b_Trigger_hltPFMHTTightID120;   //!
   TBranch        *b_Trigger_hltPFMHTTightID130;   //!
   TBranch        *b_Trigger_hltPFMHTTightID140;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_Gen_charge;   //!
   TBranch        *b_Gen_energy;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_mass;   //!
   TBranch        *b_Gen_momDisplacment;
   TBranch        *b_Gen_motherIdx;   //!
   TBranch        *b_Gen_nTotal;   //!
   TBranch        *b_Gen_pdgId;   //!
   TBranch        *b_Gen_phi;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_px;   //!
   TBranch        *b_Gen_py;   //!
   TBranch        *b_Gen_pz;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_susEvtType;   //!
   TBranch        *b_Gen_susId;   //!
   TBranch        *b_Gen_vx;   //!
   TBranch        *b_Gen_vy;   //!
   TBranch        *b_Gen_vz;   //!
   TBranch        *b_Evt_genWgt;   //!

   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_px;   //!
   TBranch        *b_Muon_py;   //!
   TBranch        *b_Muon_pz;   //!
   TBranch        *b_Muon_isLoose;   //!
   //TBranch        *b_Muon_isMedium;   //!
   TBranch        *b_Muon_nMuons;   //!
   TBranch        *b_Muon_nSelMuons;   //!

   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_chHM;   //!
   TBranch        *b_Jet_drRhIds;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_genJetLlpDr;   //!
   TBranch        *b_Jet_genDptMatch;   //!
   TBranch        *b_Jet_genDrMatch;   //!
   TBranch        *b_Jet_genEnergy;   //!
   TBranch        *b_Jet_genEta;   //!
   TBranch        *b_Jet_genImpactAngle;   //!
   TBranch        *b_Jet_genJetLlpId;   //!
   TBranch        *b_Jet_genJetLlpRp;   //!
   TBranch        *b_Jet_genPhi;   //!
   TBranch        *b_Jet_genPt;   //!
   TBranch        *b_Jet_genQrkLlpDr;   //!
   TBranch        *b_Jet_genQrkLlpId;   //!
   TBranch        *b_Jet_genQrkLlpRp;   //!
   TBranch        *b_Jet_genTOF;   //!
   TBranch        *b_Jet_genTime;   //!
   TBranch        *b_Jet_genTimeLLP;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_neHM;   //!
   TBranch        *b_Jet_egIndxs;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Met_CPt;   //!
   TBranch        *b_Met_Cpx;   //!
   TBranch        *b_Met_Cpy;   //!
   TBranch        *b_Met_CsumEt;   //!
   TBranch        *b_Met_eta;   //!
   TBranch        *b_Met_phi;   //!
   TBranch        *b_Met_pt;   //!
   TBranch        *b_Met_px;   //!
   TBranch        *b_Met_py;   //!
   TBranch        *b_Met_sumEt;   //!
   TBranch        *b_Photon_ecalRHSumEtConeDR04;   //!
   TBranch        *b_Photon_energy;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_excluded;   //!
   TBranch        *b_Photon_genIdx;   //!
   TBranch        *b_Photon_genSigXMomId;   //!
   TBranch        *b_Photon_hadOverEM;   //!
   TBranch        *b_Photon_hadTowOverEM;   //!
   TBranch        *b_Photon_hcalTowerSumEtBcConeDR04;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_Photon_isOot;   //!
   TBranch        *b_Photon_nTrkHollowConeDR04;   //!
   TBranch        *b_Photon_nTrkSolidConeDR04;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_px;   //!
   TBranch        *b_Photon_py;   //!
   TBranch        *b_Photon_pz;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_s4;   //!
   TBranch        *b_Photon_seedTOFTime;   //!
   TBranch        *b_Photon_SigmaIEtaIEta;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_ecalPFClusterIso;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_gloResRhId;   //!
   TBranch        *b_hasConversionTracks;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_hcalPFClusterIso;   //!
   TBranch        *b_Photon_Hoe_PUcorr;   //!
   TBranch        *b_Photon_locResRhId;   //!
   TBranch        *b_Photon_nPhoton;   //!
   TBranch        *b_Photon_pfChargedIso;   //!
   TBranch        *b_Photon_pfChargedIsoPFPV;   //!
   TBranch        *b_Photon_pfChargedIsoWorstVtx;   //!
   TBranch        *b_Photon_pfPhoIso03;   //!
   TBranch        *b_pfRelIso03_all_quadratic;   //!
   TBranch        *b_pfRelIso03_chg_quadratic;   //!
   TBranch        *b_Photon_scIndex;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_sieip;   //!
   TBranch        *b_Photon_sipip;   //!

   TBranch        *b_Vertex_nTracks;   //!
   TBranch        *b_Vertex_mass;   //!
   TBranch        *b_Vertex_passLooseElectronID;   //!
   TBranch        *b_Vertex_passLooseMuonID;   //!

   //llpgtree( TChain *tree=0 );
   //virtual ~llpgtree();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   void Init( TChain *tree, bool doGenInfo );
   void getBranches( Long64_t entry, bool doGenInfo );
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);

};//<<>>class llpgtree

/*
llpgtree::llpgtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gmsb_AODSIM_KUCMSNtuplizer_justin_test_fullgen_oottrue.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gmsb_AODSIM_KUCMSNtuplizer_justin_test_fullgen_oottrue.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("gmsb_AODSIM_KUCMSNtuplizer_justin_test_fullgen_oottrue.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

llpgtree::~llpgtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t llpgtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t llpgtree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
*/

void llpgtree::Init( TChain *tree, bool doGenInfo ){

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   Evt_genWgt = 1;

   // Set object pointer
   DisplacedElectron_charge = 0;
   DisplacedElectron_dRGen = 0;
   DisplacedElectron_dRSC = 0;
   DisplacedElectron_dxy = 0;
   DisplacedElectron_energy = 0;
   DisplacedElectron_eta = 0;
   DisplacedElectron_etaECAL = 0;
   DisplacedElectron_genSigIndex = 0;
   DisplacedElectron_indexGen = 0;
   DisplacedElectron_indexSC = 0;
   DisplacedElectron_indexTrack = 0;
   DisplacedElectron_ip2D = 0;
   DisplacedElectron_ip2Dbs = 0;
   DisplacedElectron_ip3D = 0;
   DisplacedElectron_ip3Dbs = 0;
   DisplacedElectron_isBarrel = 0;
   DisplacedElectron_isConversion = 0;
   DisplacedElectron_isEndCap = 0;
   DisplacedElectron_isGeneralTrack = 0;
   DisplacedElectron_isGsfTrack = 0;
   DisplacedElectron_isHeavy = 0;
   DisplacedElectron_isLight = 0;
   DisplacedElectron_isOOTSC = 0;
   DisplacedElectron_isSC = 0;
   DisplacedElectron_isSignal = 0;
   DisplacedElectron_isTau = 0;
   DisplacedElectron_isUnmatched = 0;
   DisplacedElectron_nHits = 0;
   DisplacedElectron_p = 0;
   DisplacedElectron_pdgId = 0;
   DisplacedElectron_pfIso0p4 = 0;
   DisplacedElectron_phi = 0;
   DisplacedElectron_phiECAL = 0;
   DisplacedElectron_pt = 0;
   DisplacedElectron_scIso0p4 = 0;
   DisplacedElectron_sip2D = 0;
   DisplacedElectron_sip2Dbs = 0;
   DisplacedElectron_sip3D = 0;
   DisplacedElectron_sip3Dbs = 0;
   DisplacedElectron_trackIso0p4 = 0;
   DisplacedElectron_x = 0;
   DisplacedElectron_xECAL = 0;
   DisplacedElectron_y = 0;
   DisplacedElectron_yECAL = 0;
   DisplacedElectron_z = 0;
   DisplacedElectron_zECAL = 0;
   ECALRecHit_energy = 0;
   ECALRecHit_timeError = 0;
   ECALRecHit_hasGS1 = 0;
   ECALRecHit_hasGS6 = 0;
   ECALRecHit_ID = 0;
   ECALRecHit_swCross = 0;
   ECALRecHit_0TOF = 0;
   ECALRecHit_pvTOF = 0;
   ECALRecHit_time = 0;
   ECALRecHit_amplitude = 0;
   ECALRecHit_ampres = 0;
   ECALRecHit_eta = 0;
   ECALRecHit_isDead = 0;
   ECALRecHit_isDiWrd = 0;
   ECALRecHit_isOOT = 0;
   ECALRecHit_isOther = 0;
   ECALRecHit_isPoor = 0;
   ECALRecHit_isRecov = 0;
   ECALRecHit_isWrd = 0;
   ECALRecHit_phi = 0;
   ECALRecHit_precentUsed = 0;
   ECALRecHit_rhx = 0;
   ECALRecHit_rhy = 0;
   ECALRecHit_rhz = 0;
   SuperCluster_covEtaEta = 0;
   SuperCluster_covEtaPhi = 0;
   SuperCluster_covPhiPhi = 0;
   SuperCluster_energyRaw = 0;
   SuperCluster_etaWidth = 0;
   SuperCluster_excluded = 0;
   SuperCluster_seedIsEB = 0;
   SuperCluster_isScEtaEB = 0;
   SuperCluster_isScEtaEE = 0;
   SuperCluster_isOot = 0;
   SuperCluster_original = 0;
   SuperCluster_phiWidth = 0;
   SuperCluster_salp = 0;
   SuperCluster_smaj = 0;
   SuperCluster_smin = 0;
   SuperCluster_seediEtaOriX = 0;
   SuperCluster_seediPhiOriY = 0;
   SuperCluster_nBasicClusters = 0;
   SuperCluster_ElectronIndx = 0;
   SuperCluster_energy = 0;
   SuperCluster_eta = 0;
   SuperCluster_clcx = 0;
   SuperCluster_clcy = 0;
   SuperCluster_clcz = 0;
   SuperCluster_nRHXtals = 0;
   SuperCluster_ObjectPdgId = 0;
   SuperCluster_phi = 0;
   SuperCluster_PhotonIndx = 0;
   SuperCluster_rhFracs = 0;
   SuperCluster_rhIds = 0;
   SuperCluster_XtalSeedID = 0;
   SuperCluster_nHFXtals = 0;
   SuperCluster_x_calo = 0;
   SuperCluster_y_calo = 0;
   SuperCluster_z_calo = 0;
   ECALTrackDetID_detId = 0;
   ECALTrackDetID_isECAL = 0;
   ECALTrackDetID_isHCAL = 0;
   ECALTrackDetID_trackIndex = 0;
   ECALTrack_charge = 0;
   ECALTrack_chi2 = 0;
   ECALTrack_dxy = 0;
   ECALTrack_eta = 0;
   ECALTrack_etaError = 0;
   ECALTrack_genDeltaR = 0;
   ECALTrack_genIndex = 0;
   ECALTrack_lambda = 0;
   ECALTrack_nLostHits = 0;
   ECALTrack_nValidHits = 0;
   ECALTrack_ndof = 0;
   ECALTrack_normalizedChi2 = 0;
   ECALTrack_p = 0;
   ECALTrack_pdgId = 0;
   ECALTrack_phi = 0;
   ECALTrack_phiError = 0;
   ECALTrack_pt = 0;
   ECALTrack_ptError = 0;
   ECALTrack_px = 0;
   ECALTrack_py = 0;
   ECALTrack_pz = 0;
   ECALTrack_qoverp = 0;
   ECALTrack_qualityMask = 0;
   ECALTrack_theta = 0;
   ECALTrack_vx = 0;
   ECALTrack_vy = 0;
   ECALTrack_vz = 0;
   TrackSCMatch_deltaR = 0;
   TrackSCMatch_scEnergyToTrackPRatio = 0;
   Electron_energy = 0;
   Electron_eta = 0;
   Electron_genIdx = 0;
   Electron_genSigWZId = 0;
   Electron_genSigXMomId = 0;
   Electron_phi = 0;
   Electron_pt = 0;
   Electron_px = 0;
   Electron_py = 0;
   Electron_pz = 0;
   Electron_seedTOFTime = 0;
   Electron_trackz = 0;
   Electron_scIndex = 0;
   Electron_isLoose = 0;

   Gen_charge = 0;
   Gen_energy = 0;
   Gen_eta = 0;
   Gen_mass = 0;
   Gen_momDisplacment = 0;
   Gen_motherIdx = 0;
   Gen_pdgId = 0;
   Gen_phi = 0;
   Gen_pt = 0;
   Gen_px = 0;
   Gen_py = 0;
   Gen_pz = 0;
   Gen_status = 0;
   Gen_susEvtType = 0;
   Gen_susId = 0;
   Gen_vx = 0;
   Gen_vy = 0;
   Gen_vz = 0;

   Muon_energy = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_pt = 0;
   Muon_px = 0;
   Muon_py = 0;
   Muon_pz = 0;
   Muon_isLoose = 0;
   //Muon_isMedium = 0;

   Jet_area = 0;
   Jet_chEmEF = 0;
   Jet_chHEF = 0;
   Jet_chHM = 0;
   Jet_drRhIds = 0;
   Jet_energy = 0;
   Jet_eta = 0;
   Jet_genJetLlpDr = 0;
   Jet_genDptMatch = 0;
   Jet_genDrMatch = 0;
   Jet_genEnergy = 0;
   Jet_genEta = 0;
   Jet_genImpactAngle = 0;
   Jet_genJetLlpId = 0;
   Jet_genJetLlpRp = 0;
   Jet_genPhi = 0;
   Jet_genPt = 0;
   Jet_genQrkLlpDr = 0;
   Jet_genQrkLlpId = 0;
   Jet_genQrkLlpRp = 0;
   Jet_genTOF = 0;
   Jet_genTime = 0;
   Jet_genTimeLLP = 0;
   Jet_mass = 0;
   Jet_muEF = 0;
   Jet_neEmEF = 0;
   Jet_neHEF = 0;
   Jet_neHM = 0;
   Jet_egIndxs = 0;
   Jet_phi = 0;
   Jet_pt = 0;
   Jet_nConstituents = 0;
   Photon_ecalRHSumEtConeDR04 = 0;
   Photon_energy = 0;
   Photon_energyErr = 0;
   Photon_eta = 0;
   Photon_excluded = 0;
   Photon_genIdx = 0;
   Photon_genSigXMomId = 0;
   Photon_hadOverEM = 0;
   Photon_hadTowOverEM = 0;
   Photon_hcalTowerSumEtBcConeDR04 = 0;
   Photon_hcalTowerSumEtConeDR04 = 0;
   Photon_isOot = 0;
   Photon_nTrkHollowConeDR04 = 0;
   Photon_nTrkSolidConeDR04 = 0;
   Photon_phi = 0;
   Photon_pt = 0;
   Photon_px = 0;
   Photon_py = 0;
   Photon_pz = 0;
   Photon_r9 = 0;
   Photon_s4 = 0;
   Photon_seedTOFTime = 0;
   Photon_SigmaIEtaIEta = 0;
   Photon_trkSumPtHollowConeDR03 = 0;
   Photon_trkSumPtHollowConeDR04 = 0;
   Photon_trkSumPtSolidConeDR04 = 0;
   Photon_ecalPFClusterIso = 0;
   Photon_electronVeto = 0;
   Photon_gloResRhId = 0;
   hasConversionTracks = 0;
   Photon_pixelSeed = 0;
   Photon_hcalPFClusterIso = 0;
   Photon_Hoe_PUcorr = 0;
   Photon_locResRhId = 0;
   Photon_pfChargedIso = 0;
   Photon_pfChargedIsoPFPV = 0;
   Photon_pfChargedIsoWorstVtx = 0;
   Photon_pfPhoIso03 = 0;
   pfRelIso03_all_quadratic = 0;
   pfRelIso03_chg_quadratic = 0;
   Photon_scIndex = 0;
   Photon_sieie = 0;
   Photon_sieip = 0;
   Photon_sipip = 0;

   Vertex_nTracks = 0;
   Vertex_mass = 0;
   Vertex_passLooseElectronID = 0;
   Vertex_passLooseMuonID = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);

/*
   fChain->SetBranchAddress("DisplacedElectron_charge", &DisplacedElectron_charge, &b_DisplacedElectron_charge);
   fChain->SetBranchAddress("DisplacedElectron_costGen", &DisplacedElectron_costGen, &b_DisplacedElectron_costGen);
   fChain->SetBranchAddress("DisplacedElectron_costSC", &DisplacedElectron_costSC, &b_DisplacedElectron_costSC);
   fChain->SetBranchAddress("DisplacedElectron_dRGen", &DisplacedElectron_dRGen, &b_DisplacedElectron_dRGen);
   fChain->SetBranchAddress("DisplacedElectron_dRSC", &DisplacedElectron_dRSC, &b_DisplacedElectron_dRSC);
   fChain->SetBranchAddress("DisplacedElectron_dxy", &DisplacedElectron_dxy, &b_DisplacedElectron_dxy);
   fChain->SetBranchAddress("DisplacedElectron_energy", &DisplacedElectron_energy, &b_DisplacedElectron_energy);
   fChain->SetBranchAddress("DisplacedElectron_eta", &DisplacedElectron_eta, &b_DisplacedElectron_eta);
   fChain->SetBranchAddress("DisplacedElectron_etaECAL", &DisplacedElectron_etaECAL, &b_DisplacedElectron_etaECAL);
   fChain->SetBranchAddress("DisplacedElectron_genSigIndex", &DisplacedElectron_genSigIndex, &b_DisplacedElectron_genSigIndex);
   fChain->SetBranchAddress("DisplacedElectron_indexGen", &DisplacedElectron_indexGen, &b_DisplacedElectron_indexGen);
   fChain->SetBranchAddress("DisplacedElectron_indexSC", &DisplacedElectron_indexSC, &b_DisplacedElectron_indexSC);
   fChain->SetBranchAddress("DisplacedElectron_indexTrack", &DisplacedElectron_indexTrack, &b_DisplacedElectron_indexTrack);
   fChain->SetBranchAddress("DisplacedElectron_ip2D", &DisplacedElectron_ip2D, &b_DisplacedElectron_ip2D);
   fChain->SetBranchAddress("DisplacedElectron_ip2Dbs", &DisplacedElectron_ip2Dbs, &b_DisplacedElectron_ip2Dbs);
   fChain->SetBranchAddress("DisplacedElectron_ip3D", &DisplacedElectron_ip3D, &b_DisplacedElectron_ip3D);
   fChain->SetBranchAddress("DisplacedElectron_ip3Dbs", &DisplacedElectron_ip3Dbs, &b_DisplacedElectron_ip3Dbs);
   fChain->SetBranchAddress("DisplacedElectron_isBarrel", &DisplacedElectron_isBarrel, &b_DisplacedElectron_isBarrel);
   fChain->SetBranchAddress("DisplacedElectron_isConversion", &DisplacedElectron_isConversion, &b_DisplacedElectron_isConversion);
   fChain->SetBranchAddress("DisplacedElectron_isEndCap", &DisplacedElectron_isEndCap, &b_DisplacedElectron_isEndCap);
   fChain->SetBranchAddress("DisplacedElectron_isGeneralTrack", &DisplacedElectron_isGeneralTrack, &b_DisplacedElectron_isGeneralTrack);
   fChain->SetBranchAddress("DisplacedElectron_isGsfTrack", &DisplacedElectron_isGsfTrack, &b_DisplacedElectron_isGsfTrack);
   fChain->SetBranchAddress("DisplacedElectron_isHeavy", &DisplacedElectron_isHeavy, &b_DisplacedElectron_isHeavy);
   fChain->SetBranchAddress("DisplacedElectron_isLight", &DisplacedElectron_isLight, &b_DisplacedElectron_isLight);
   fChain->SetBranchAddress("DisplacedElectron_isOOTSC", &DisplacedElectron_isOOTSC, &b_DisplacedElectron_isOOTSC);
   fChain->SetBranchAddress("DisplacedElectron_isSC", &DisplacedElectron_isSC, &b_DisplacedElectron_isSC);
   fChain->SetBranchAddress("DisplacedElectron_isSignal", &DisplacedElectron_isSignal, &b_DisplacedElectron_isSignal);
   fChain->SetBranchAddress("DisplacedElectron_isTau", &DisplacedElectron_isTau, &b_DisplacedElectron_isTau);
   fChain->SetBranchAddress("DisplacedElectron_isUnmatched", &DisplacedElectron_isUnmatched, &b_DisplacedElectron_isUnmatched);

   fChain->SetBranchAddress("DisplacedElectron_nGenSig", &DisplacedElectron_nGenSig, &b_DisplacedElectron_nGenSig);
   fChain->SetBranchAddress("DisplacedElectron_pdgId", &DisplacedElectron_pdgId, &b_DisplacedElectron_pdgId);

   fChain->SetBranchAddress("DisplacedElectron_nHits", &DisplacedElectron_nHits, &b_DisplacedElectron_nHits);
   fChain->SetBranchAddress("DisplacedElectron_nOOTSCs", &DisplacedElectron_nOOTSCs, &b_DisplacedElectron_nOOTSCs);
   fChain->SetBranchAddress("DisplacedElectron_nSCs", &DisplacedElectron_nSCs, &b_DisplacedElectron_nSCs);
   fChain->SetBranchAddress("DisplacedElectron_nTotal", &DisplacedElectron_nTotal, &b_DisplacedElectron_nTotal);
   fChain->SetBranchAddress("DisplacedElectron_p", &DisplacedElectron_p, &b_DisplacedElectron_p);

   fChain->SetBranchAddress("DisplacedElectron_pfIso0p4", &DisplacedElectron_pfIso0p4, &b_DisplacedElectron_pfIso0p4);
   fChain->SetBranchAddress("DisplacedElectron_phi", &DisplacedElectron_phi, &b_DisplacedElectron_phi);
   fChain->SetBranchAddress("DisplacedElectron_phiECAL", &DisplacedElectron_phiECAL, &b_DisplacedElectron_phiECAL);
   fChain->SetBranchAddress("DisplacedElectron_pt", &DisplacedElectron_pt, &b_DisplacedElectron_pt);
   fChain->SetBranchAddress("DisplacedElectron_scIso0p4", &DisplacedElectron_scIso0p4, &b_DisplacedElectron_scIso0p4);
   fChain->SetBranchAddress("DisplacedElectron_sip2D", &DisplacedElectron_sip2D, &b_DisplacedElectron_sip2D);
   fChain->SetBranchAddress("DisplacedElectron_sip2Dbs", &DisplacedElectron_sip2Dbs, &b_DisplacedElectron_sip2Dbs);
   fChain->SetBranchAddress("DisplacedElectron_sip3D", &DisplacedElectron_sip3D, &b_DisplacedElectron_sip3D);
   fChain->SetBranchAddress("DisplacedElectron_sip3Dbs", &DisplacedElectron_sip3Dbs, &b_DisplacedElectron_sip3Dbs);
   fChain->SetBranchAddress("DisplacedElectron_trackIso0p4", &DisplacedElectron_trackIso0p4, &b_DisplacedElectron_trackIso0p4);
   fChain->SetBranchAddress("DisplacedElectron_x", &DisplacedElectron_x, &b_DisplacedElectron_x);
   fChain->SetBranchAddress("DisplacedElectron_xECAL", &DisplacedElectron_xECAL, &b_DisplacedElectron_xECAL);
   fChain->SetBranchAddress("DisplacedElectron_y", &DisplacedElectron_y, &b_DisplacedElectron_y);
   fChain->SetBranchAddress("DisplacedElectron_yECAL", &DisplacedElectron_yECAL, &b_DisplacedElectron_yECAL);
   fChain->SetBranchAddress("DisplacedElectron_z", &DisplacedElectron_z, &b_DisplacedElectron_z);
   fChain->SetBranchAddress("DisplacedElectron_zECAL", &DisplacedElectron_zECAL, &b_DisplacedElectron_zECAL);
*/

   fChain->SetBranchAddress("ECALRecHit_energy", &ECALRecHit_energy, &b_ECALRecHit_energy);
   //fChain->SetBranchAddress("ECALRecHit_timeError", &ECALRecHit_timeError, &b_ECALRecHit_timeError);
   //fChain->SetBranchAddress("ECALRecHit_hasGS1", &ECALRecHit_hasGS1, &b_ECALRecHit_hasGS1);
   //fChain->SetBranchAddress("ECALRecHit_hasGS6", &ECALRecHit_hasGS6, &b_ECALRecHit_hasGS6);
   fChain->SetBranchAddress("ECALRecHit_ID", &ECALRecHit_ID, &b_ECALRecHit_ID);
   fChain->SetBranchAddress("ECALRecHit_swCross", &ECALRecHit_swCross, &b_ECALRecHit_swCross);
   fChain->SetBranchAddress("ECALRecHit_0TOF", &ECALRecHit_0TOF, &b_ECALRecHit_0TOF);
   fChain->SetBranchAddress("ECALRecHit_pvTOF", &ECALRecHit_pvTOF, &b_ECALRecHit_pvTOF);
   fChain->SetBranchAddress("ECALRecHit_time", &ECALRecHit_time, &b_ECALRecHit_time);
   fChain->SetBranchAddress("ECALRecHit_amplitude", &ECALRecHit_amplitude, &b_ECALRecHit_amplitude);
   fChain->SetBranchAddress("ECALRecHit_ampres", &ECALRecHit_ampres, &b_ECALRecHit_ampres);
   fChain->SetBranchAddress("ECALRecHit_eta", &ECALRecHit_eta, &b_ECALRecHit_eta);
   fChain->SetBranchAddress("ECALRecHit_isDead", &ECALRecHit_isDead, &b_ECALRecHit_isDead);
   fChain->SetBranchAddress("ECALRecHit_isDiWrd", &ECALRecHit_isDiWrd, &b_ECALRecHit_isDiWrd);
   fChain->SetBranchAddress("ECALRecHit_isOOT", &ECALRecHit_isOOT, &b_ECALRecHit_isOOT);
   fChain->SetBranchAddress("ECALRecHit_isOther", &ECALRecHit_isOther, &b_ECALRecHit_isOther);
   fChain->SetBranchAddress("ECALRecHit_isPoor", &ECALRecHit_isPoor, &b_ECALRecHit_isPoor);
   fChain->SetBranchAddress("ECALRecHit_isRecov", &ECALRecHit_isRecov, &b_ECALRecHit_isRecov);
   fChain->SetBranchAddress("ECALRecHit_isWrd", &ECALRecHit_isWrd, &b_ECALRecHit_isWrd);
   fChain->SetBranchAddress("ECALRecHit_phi", &ECALRecHit_phi, &b_ECALRecHit_phi);
   fChain->SetBranchAddress("ECALRecHit_precentUsed", &ECALRecHit_precentUsed, &b_ECALRecHit_precentUsed);
   fChain->SetBranchAddress("ECALRecHit_rhx", &ECALRecHit_rhx, &b_ECALRecHit_rhx);
   fChain->SetBranchAddress("ECALRecHit_rhy", &ECALRecHit_rhy, &b_ECALRecHit_rhy);
   fChain->SetBranchAddress("ECALRecHit_rhz", &ECALRecHit_rhz, &b_ECALRecHit_rhz);

   fChain->SetBranchAddress("SuperCluster_covEtaEta", &SuperCluster_covEtaEta, &b_SuperCluster_covEtaEta);
   fChain->SetBranchAddress("SuperCluster_covEtaPhi", &SuperCluster_covEtaPhi, &b_SuperCluster_covEtaPhi);
   fChain->SetBranchAddress("SuperCluster_covPhiPhi", &SuperCluster_covPhiPhi, &b_SuperCluster_covPhiPhi);
   fChain->SetBranchAddress("SuperCluster_energyRaw", &SuperCluster_energyRaw, &b_SuperCluster_energyRaw);
   fChain->SetBranchAddress("SuperCluster_etaWidth", &SuperCluster_etaWidth, &b_SuperCluster_etaWidth);
   fChain->SetBranchAddress("SuperCluster_excluded", &SuperCluster_excluded, &b_SuperCluster_excluded);
   fChain->SetBranchAddress("SuperCluster_seedIsEB", &SuperCluster_seedIsEB, &b_SuperCluster_seedIsEB);
   fChain->SetBranchAddress("SuperCluster_isScEtaEB", &SuperCluster_isScEtaEB, &b_SuperCluster_isScEtaEB);
   fChain->SetBranchAddress("SuperCluster_isScEtaEE", &SuperCluster_isScEtaEE, &b_SuperCluster_isScEtaEE);
   fChain->SetBranchAddress("SuperCluster_isOot", &SuperCluster_isOot, &b_SuperCluster_isOot);
   fChain->SetBranchAddress("SuperCluster_original", &SuperCluster_original, &b_SuperCluster_original);
   fChain->SetBranchAddress("SuperCluster_phiWidth", &SuperCluster_phiWidth, &b_SuperCluster_phiWidth);
   fChain->SetBranchAddress("SuperCluster_salp", &SuperCluster_salp, &b_SuperCluster_salp);
   fChain->SetBranchAddress("SuperCluster_smaj", &SuperCluster_smaj, &b_SuperCluster_smaj);
   fChain->SetBranchAddress("SuperCluster_smin", &SuperCluster_smin, &b_SuperCluster_smin);
   fChain->SetBranchAddress("SuperCluster_seediEtaOriX", &SuperCluster_seediEtaOriX, &b_SuperCluster_seediEtaOriX);
   fChain->SetBranchAddress("SuperCluster_seediPhiOriY", &SuperCluster_seediPhiOriY, &b_SuperCluster_seediPhiOriY);
   fChain->SetBranchAddress("SuperCluster_nBasicClusters", &SuperCluster_nBasicClusters, &b_SuperCluster_nBasicClusters);
   fChain->SetBranchAddress("SuperCluster_ElectronIndx", &SuperCluster_ElectronIndx, &b_SuperCluster_ElectronIndx);
   fChain->SetBranchAddress("SuperCluster_energy", &SuperCluster_energy, &b_SuperCluster_energy);
   fChain->SetBranchAddress("SuperCluster_eta", &SuperCluster_eta, &b_SuperCluster_eta);
   fChain->SetBranchAddress("SuperCluster_clcx", &SuperCluster_clcx, &b_SuperCluster_clcx);
   fChain->SetBranchAddress("SuperCluster_clcy", &SuperCluster_clcy, &b_SuperCluster_clcy);
   fChain->SetBranchAddress("SuperCluster_clcz", &SuperCluster_clcz, &b_SuperCluster_clcz);
   fChain->SetBranchAddress("SuperCluster_nSuperCluster", &SuperCluster_nSuperCluster, &b_SuperCluster_nSuperCluster);
   fChain->SetBranchAddress("SuperCluster_nRHXtals", &SuperCluster_nRHXtals, &b_SuperCluster_nRHXtals);
   if( doGenInfo ) fChain->SetBranchAddress("SuperCluster_ObjectPdgId", &SuperCluster_ObjectPdgId, &b_SuperCluster_ObjectPdgId);  
   fChain->SetBranchAddress("SuperCluster_phi", &SuperCluster_phi, &b_SuperCluster_phi);
   fChain->SetBranchAddress("SuperCluster_PhotonIndx", &SuperCluster_PhotonIndx, &b_SuperCluster_PhotonIndx);
   fChain->SetBranchAddress("SuperCluster_rhFracs", &SuperCluster_rhFracs, &b_SuperCluster_rhFracs);
   fChain->SetBranchAddress("SuperCluster_rhIds", &SuperCluster_rhIds, &b_SuperCluster_rhIds);
   fChain->SetBranchAddress("SuperCluster_XtalSeedID", &SuperCluster_XtalSeedID, &b_SuperCluster_XtalSeedID);
   fChain->SetBranchAddress("SuperCluster_nHFXtals", &SuperCluster_nHFXtals, &b_SuperCluster_nHFXtals);
   fChain->SetBranchAddress("SuperCluster_x_calo", &SuperCluster_x_calo, &b_SuperCluster_x_calo);
   fChain->SetBranchAddress("SuperCluster_y_calo", &SuperCluster_y_calo, &b_SuperCluster_y_calo);
   fChain->SetBranchAddress("SuperCluster_z_calo", &SuperCluster_z_calo, &b_SuperCluster_z_calo);

   fChain->SetBranchAddress("ECALTrackDetID_detId", &ECALTrackDetID_detId, &b_ECALTrackDetID_detId);
   fChain->SetBranchAddress("ECALTrackDetID_isECAL", &ECALTrackDetID_isECAL, &b_ECALTrackDetID_isECAL);
   fChain->SetBranchAddress("ECALTrackDetID_isHCAL", &ECALTrackDetID_isHCAL, &b_ECALTrackDetID_isHCAL);
   fChain->SetBranchAddress("ECALTrackDetID_trackIndex", &ECALTrackDetID_trackIndex, &b_ECALTrackDetID_trackIndex);
   fChain->SetBranchAddress("ECALTrack_charge", &ECALTrack_charge, &b_ECALTrack_charge);
   fChain->SetBranchAddress("ECALTrack_chi2", &ECALTrack_chi2, &b_ECALTrack_chi2);
   fChain->SetBranchAddress("ECALTrack_dxy", &ECALTrack_dxy, &b_ECALTrack_dxy);
   fChain->SetBranchAddress("ECALTrack_eta", &ECALTrack_eta, &b_ECALTrack_eta);
   fChain->SetBranchAddress("ECALTrack_etaError", &ECALTrack_etaError, &b_ECALTrack_etaError);
   if( doGenInfo ){
   fChain->SetBranchAddress("ECALTrack_genDeltaR", &ECALTrack_genDeltaR, &b_ECALTrack_genDeltaR);
   fChain->SetBranchAddress("ECALTrack_genIndex", &ECALTrack_genIndex, &b_ECALTrack_genIndex);
   fChain->SetBranchAddress("ECALTrack_pdgId", &ECALTrack_pdgId, &b_ECALTrack_pdgId);
   }//<<>>if( doGenInfo )
   fChain->SetBranchAddress("ECALTrack_lambda", &ECALTrack_lambda, &b_ECALTrack_lambda);
   fChain->SetBranchAddress("ECALTrack_nGeneral", &ECALTrack_nGeneral, &b_ECALTrack_nGeneral);
   fChain->SetBranchAddress("ECALTrack_nGsf", &ECALTrack_nGsf, &b_ECALTrack_nGsf);
   fChain->SetBranchAddress("ECALTrack_nLostHits", &ECALTrack_nLostHits, &b_ECALTrack_nLostHits);
   fChain->SetBranchAddress("ECALTrack_nTracks", &ECALTrack_nTracks, &b_ECALTrack_nTracks);
   fChain->SetBranchAddress("ECALTrack_nValidHits", &ECALTrack_nValidHits, &b_ECALTrack_nValidHits);
   fChain->SetBranchAddress("ECALTrack_ndof", &ECALTrack_ndof, &b_ECALTrack_ndof);
   fChain->SetBranchAddress("ECALTrack_normalizedChi2", &ECALTrack_normalizedChi2, &b_ECALTrack_normalizedChi2);
   fChain->SetBranchAddress("ECALTrack_p", &ECALTrack_p, &b_ECALTrack_p);
   fChain->SetBranchAddress("ECALTrack_phi", &ECALTrack_phi, &b_ECALTrack_phi);
   fChain->SetBranchAddress("ECALTrack_phiError", &ECALTrack_phiError, &b_ECALTrack_phiError);
   fChain->SetBranchAddress("ECALTrack_pt", &ECALTrack_pt, &b_ECALTrack_pt);
   fChain->SetBranchAddress("ECALTrack_ptError", &ECALTrack_ptError, &b_ECALTrack_ptError);
   fChain->SetBranchAddress("ECALTrack_px", &ECALTrack_px, &b_ECALTrack_px);
   fChain->SetBranchAddress("ECALTrack_py", &ECALTrack_py, &b_ECALTrack_py);
   fChain->SetBranchAddress("ECALTrack_pz", &ECALTrack_pz, &b_ECALTrack_pz);
   fChain->SetBranchAddress("ECALTrack_qoverp", &ECALTrack_qoverp, &b_ECALTrack_qoverp);
   fChain->SetBranchAddress("ECALTrack_qualityMask", &ECALTrack_qualityMask, &b_ECALTrack_qualityMask);
   fChain->SetBranchAddress("ECALTrack_theta", &ECALTrack_theta, &b_ECALTrack_theta);
   fChain->SetBranchAddress("ECALTrack_vx", &ECALTrack_vx, &b_ECALTrack_vx);
   fChain->SetBranchAddress("ECALTrack_vy", &ECALTrack_vy, &b_ECALTrack_vy);
   fChain->SetBranchAddress("ECALTrack_vz", &ECALTrack_vz, &b_ECALTrack_vz);
   fChain->SetBranchAddress("TrackSCMatch_deltaR", &TrackSCMatch_deltaR, &b_TrackSCMatch_deltaR);
   fChain->SetBranchAddress("TrackSCMatch_scEnergyToTrackPRatio", &TrackSCMatch_scEnergyToTrackPRatio, &b_TrackSCMatch_scEnergyToTrackPRatio);
   fChain->SetBranchAddress("Electron_energy", &Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   if( doGenInfo ){
   fChain->SetBranchAddress("Electron_genIdx", &Electron_genIdx, &b_Electron_genIdx);
   fChain->SetBranchAddress("Electron_genSigWZId", &Electron_genSigWZId, &b_Electron_genSigWZId);
   fChain->SetBranchAddress("Electron_genSigXMomId", &Electron_genSigXMomId, &b_Electron_genSigXMomId);
   }//<<>>if( doGenInfo )
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_px", &Electron_px, &b_Electron_px);
   fChain->SetBranchAddress("Electron_py", &Electron_py, &b_Electron_py);
   fChain->SetBranchAddress("Electron_pz", &Electron_pz, &b_Electron_pz);
   fChain->SetBranchAddress("Electron_seedTOFTime", &Electron_seedTOFTime, &b_Electron_seedTOFTime);
   fChain->SetBranchAddress("Electron_trackz", &Electron_trackz, &b_Electron_trackz);
   fChain->SetBranchAddress("Electron_nElectrons", &Electron_nElectrons, &b_Electron_nElectrons);
   fChain->SetBranchAddress("Electron_scIndex", &Electron_scIndex, &b_Electron_scIndex);
   //fChain->SetBranchAddress("Electron_nSelElectrons", &Electron_nSelElectrons, &b_Electron_nSelElectrons);
   //fChain->SetBranchAddress("Electron_isLoose", &Electron_isLoose, &b_Electron_isLoose);

   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);

   fChain->SetBranchAddress("Evt_luminosityBlock", &Evt_luminosityBlock, &b_Evt_luminosityBlock);
   fChain->SetBranchAddress("Evt_run", &Evt_run, &b_Evt_run);
   fChain->SetBranchAddress("Evt_event", &Evt_event, &b_Evt_event);

/*
   fChain->SetBranchAddress("Trigger_hltEG20EtFilterLooseHoverE", &Trigger_hltEG20EtFilterLooseHoverE, &b_Trigger_hltEG20EtFilterLooseHoverE);
   fChain->SetBranchAddress("Trigger_hltEG20EtL1TripleEGFilter", &Trigger_hltEG20EtL1TripleEGFilter, &b_Trigger_hltEG20EtL1TripleEGFilter);
   fChain->SetBranchAddress("Trigger_hltEG20HEFilterLooseHoverE", &Trigger_hltEG20HEFilterLooseHoverE, &b_Trigger_hltEG20HEFilterLooseHoverE);
   fChain->SetBranchAddress("Trigger_hltEG20HEL1TripleEGFilter", &Trigger_hltEG20HEL1TripleEGFilter, &b_Trigger_hltEG20HEL1TripleEGFilter);
   fChain->SetBranchAddress("Trigger_hltEG20L1SingleEGLowETOrEtFilter", &Trigger_hltEG20L1SingleEGLowETOrEtFilter, &b_Trigger_hltEG20L1SingleEGLowETOrEtFilter);
   fChain->SetBranchAddress("Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter", &Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter, &b_Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter);
   fChain->SetBranchAddress("Trigger_hltEG30EBHE10R9Id50b80eHEFilter", &Trigger_hltEG30EBHE10R9Id50b80eHEFilter, &b_Trigger_hltEG30EBHE10R9Id50b80eHEFilter);
   fChain->SetBranchAddress("Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter", &Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter, &b_Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter);
   fChain->SetBranchAddress("Trigger_hltEG30EBR9Id50b80eR9IdFilter", &Trigger_hltEG30EBR9Id50b80eR9IdFilter, &b_Trigger_hltEG30EBR9Id50b80eR9IdFilter);
   fChain->SetBranchAddress("Trigger_hltEG30EIso15HE30EcalIsoLastFilter", &Trigger_hltEG30EIso15HE30EcalIsoLastFilter, &b_Trigger_hltEG30EIso15HE30EcalIsoLastFilter);
   fChain->SetBranchAddress("Trigger_hltEG30EtFilterLooseHoverE", &Trigger_hltEG30EtFilterLooseHoverE, &b_Trigger_hltEG30EtFilterLooseHoverE);
   fChain->SetBranchAddress("Trigger_hltEG30EtL1TripleEGFilter", &Trigger_hltEG30EtL1TripleEGFilter, &b_Trigger_hltEG30EtL1TripleEGFilter);
   fChain->SetBranchAddress("Trigger_hltEG30HE30HEFilter", &Trigger_hltEG30HE30HEFilter, &b_Trigger_hltEG30HE30HEFilter);
   fChain->SetBranchAddress("Trigger_hltEG30HEFilterLooseHoverE", &Trigger_hltEG30HEFilterLooseHoverE, &b_Trigger_hltEG30HEFilterLooseHoverE);
   fChain->SetBranchAddress("Trigger_hltEG30HEL1TripleEGFilter", &Trigger_hltEG30HEL1TripleEGFilter, &b_Trigger_hltEG30HEL1TripleEGFilter);
   fChain->SetBranchAddress("Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter", &Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter, &b_Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter);
   fChain->SetBranchAddress("Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter", &Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter, &b_Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter);
   fChain->SetBranchAddress("Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter", &Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter, &b_Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter);
   fChain->SetBranchAddress("Trigger_hltEG30LHE12R9Id50b80eHEFilter", &Trigger_hltEG30LHE12R9Id50b80eHEFilter, &b_Trigger_hltEG30LHE12R9Id50b80eHEFilter);
   fChain->SetBranchAddress("Trigger_hltEG30LR9Id50b80eR9IdFilter", &Trigger_hltEG30LR9Id50b80eR9IdFilter, &b_Trigger_hltEG30LR9Id50b80eR9IdFilter);
   fChain->SetBranchAddress("Trigger_hltEG30PVHE10R9Id50b80eHEFilter", &Trigger_hltEG30PVHE10R9Id50b80eHEFilter, &b_Trigger_hltEG30PVHE10R9Id50b80eHEFilter);
   fChain->SetBranchAddress("Trigger_hltEG30PVR9Id50b80eR9IdFilter", &Trigger_hltEG30PVR9Id50b80eR9IdFilter, &b_Trigger_hltEG30PVR9Id50b80eR9IdFilter);
   fChain->SetBranchAddress("Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter", &Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter, &b_Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter);
   fChain->SetBranchAddress("Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter", &Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter, &b_Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter);
   fChain->SetBranchAddress("Trigger_hltEG60EtFilter", &Trigger_hltEG60EtFilter, &b_Trigger_hltEG60EtFilter);
   fChain->SetBranchAddress("Trigger_hltEG60HEFilter", &Trigger_hltEG60HEFilter, &b_Trigger_hltEG60HEFilter);
   fChain->SetBranchAddress("Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter", &Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter, &b_Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter);
   fChain->SetBranchAddress("Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter", &Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter, &b_Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter);
   fChain->SetBranchAddress("Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter", &Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter, &b_Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter);
   fChain->SetBranchAddress("Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter", &Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter, &b_Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter);
   fChain->SetBranchAddress("Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter", &Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter, &b_Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter);
   fChain->SetBranchAddress("Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter", &Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter, &b_Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter);
   fChain->SetBranchAddress("Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter", &Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter, &b_Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter);
   fChain->SetBranchAddress("Trigger_hltHT130Jet30", &Trigger_hltHT130Jet30, &b_Trigger_hltHT130Jet30);
   fChain->SetBranchAddress("Trigger_hltHT175Jet10", &Trigger_hltHT175Jet10, &b_Trigger_hltHT175Jet10);
   fChain->SetBranchAddress("Trigger_hltHT200Jet30", &Trigger_hltHT200Jet30, &b_Trigger_hltHT200Jet30);
   fChain->SetBranchAddress("Trigger_hltL1sSingleAndDoubleEG", &Trigger_hltL1sSingleAndDoubleEG, &b_Trigger_hltL1sSingleAndDoubleEG);
   fChain->SetBranchAddress("Trigger_hltL1sSingleAndDoubleEGNonIsoOr", &Trigger_hltL1sSingleAndDoubleEGNonIsoOr, &b_Trigger_hltL1sSingleAndDoubleEGNonIsoOr);
   fChain->SetBranchAddress("Trigger_hltL1sSingleAndDoubleEGor", &Trigger_hltL1sSingleAndDoubleEGor, &b_Trigger_hltL1sSingleAndDoubleEGor);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG15", &Trigger_hltL1sSingleEG15, &b_Trigger_hltL1sSingleEG15);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG18", &Trigger_hltL1sSingleEG18, &b_Trigger_hltL1sSingleEG18);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG24", &Trigger_hltL1sSingleEG24, &b_Trigger_hltL1sSingleEG24);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG26", &Trigger_hltL1sSingleEG26, &b_Trigger_hltL1sSingleEG26);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG34to45", &Trigger_hltL1sSingleEG34to45, &b_Trigger_hltL1sSingleEG34to45);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG34to50", &Trigger_hltL1sSingleEG34to50, &b_Trigger_hltL1sSingleEG34to50);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEG40to50", &Trigger_hltL1sSingleEG40to50, &b_Trigger_hltL1sSingleEG40to50);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS", &Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS, &b_Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS);
   fChain->SetBranchAddress("Trigger_hltL1sSingleEGor", &Trigger_hltL1sSingleEGor, &b_Trigger_hltL1sSingleEGor);
   fChain->SetBranchAddress("Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG", &Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG, &b_Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG);
   fChain->SetBranchAddress("Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter", &Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter, &b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter);
   fChain->SetBranchAddress("Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt", &Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt, &b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt);
   fChain->SetBranchAddress("Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt", &Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt, &b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt);
   fChain->SetBranchAddress("Trigger_hltPFHT180Jet30", &Trigger_hltPFHT180Jet30, &b_Trigger_hltPFHT180Jet30);
   fChain->SetBranchAddress("Trigger_hltPFHT250Jet30", &Trigger_hltPFHT250Jet30, &b_Trigger_hltPFHT250Jet30);
   fChain->SetBranchAddress("Trigger_hltPFHT350Jet15", &Trigger_hltPFHT350Jet15, &b_Trigger_hltPFHT350Jet15);
   fChain->SetBranchAddress("Trigger_hltPFMET100", &Trigger_hltPFMET100, &b_Trigger_hltPFMET100);
   fChain->SetBranchAddress("Trigger_hltPFMET110", &Trigger_hltPFMET110, &b_Trigger_hltPFMET110);
   fChain->SetBranchAddress("Trigger_hltPFMET120", &Trigger_hltPFMET120, &b_Trigger_hltPFMET120);
   fChain->SetBranchAddress("Trigger_hltPFMET130", &Trigger_hltPFMET130, &b_Trigger_hltPFMET130);
   fChain->SetBranchAddress("Trigger_hltPFMET140", &Trigger_hltPFMET140, &b_Trigger_hltPFMET140);
   fChain->SetBranchAddress("Trigger_hltPFMET200", &Trigger_hltPFMET200, &b_Trigger_hltPFMET200);
   fChain->SetBranchAddress("Trigger_hltPFMET250", &Trigger_hltPFMET250, &b_Trigger_hltPFMET250);
   fChain->SetBranchAddress("Trigger_hltPFMET300", &Trigger_hltPFMET300, &b_Trigger_hltPFMET300);
   fChain->SetBranchAddress("Trigger_hltPFMET50", &Trigger_hltPFMET50, &b_Trigger_hltPFMET50);
   fChain->SetBranchAddress("Trigger_hltPFMET70", &Trigger_hltPFMET70, &b_Trigger_hltPFMET70);
   fChain->SetBranchAddress("Trigger_hltPFMET90", &Trigger_hltPFMET90, &b_Trigger_hltPFMET90);
   fChain->SetBranchAddress("Trigger_hltPFMETNoMu100", &Trigger_hltPFMETNoMu100, &b_Trigger_hltPFMETNoMu100);
   fChain->SetBranchAddress("Trigger_hltPFMETNoMu110", &Trigger_hltPFMETNoMu110, &b_Trigger_hltPFMETNoMu110);
   fChain->SetBranchAddress("Trigger_hltPFMETNoMu120", &Trigger_hltPFMETNoMu120, &b_Trigger_hltPFMETNoMu120);
   fChain->SetBranchAddress("Trigger_hltPFMETNoMu130", &Trigger_hltPFMETNoMu130, &b_Trigger_hltPFMETNoMu130);
   fChain->SetBranchAddress("Trigger_hltPFMETNoMu140", &Trigger_hltPFMETNoMu140, &b_Trigger_hltPFMETNoMu140);
   fChain->SetBranchAddress("Trigger_hltPFMETNoMu60", &Trigger_hltPFMETNoMu60, &b_Trigger_hltPFMETNoMu60);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID100", &Trigger_hltPFMHTNoMuTightID100, &b_Trigger_hltPFMHTNoMuTightID100);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID110", &Trigger_hltPFMHTNoMuTightID110, &b_Trigger_hltPFMHTNoMuTightID110);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID120", &Trigger_hltPFMHTNoMuTightID120, &b_Trigger_hltPFMHTNoMuTightID120);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID130", &Trigger_hltPFMHTNoMuTightID130, &b_Trigger_hltPFMHTNoMuTightID130);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID140", &Trigger_hltPFMHTNoMuTightID140, &b_Trigger_hltPFMHTNoMuTightID140);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID70", &Trigger_hltPFMHTNoMuTightID70, &b_Trigger_hltPFMHTNoMuTightID70);
   fChain->SetBranchAddress("Trigger_hltPFMHTNoMuTightID90", &Trigger_hltPFMHTNoMuTightID90, &b_Trigger_hltPFMHTNoMuTightID90);
   fChain->SetBranchAddress("Trigger_hltPFMHTTightID120", &Trigger_hltPFMHTTightID120, &b_Trigger_hltPFMHTTightID120);
   fChain->SetBranchAddress("Trigger_hltPFMHTTightID130", &Trigger_hltPFMHTTightID130, &b_Trigger_hltPFMHTTightID130);
   fChain->SetBranchAddress("Trigger_hltPFMHTTightID140", &Trigger_hltPFMHTTightID140, &b_Trigger_hltPFMHTTightID140);
*/

   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);

   if( doGenInfo ){
   fChain->SetBranchAddress("Gen_charge", &Gen_charge, &b_Gen_charge);
   fChain->SetBranchAddress("Gen_energy", &Gen_energy, &b_Gen_energy);
   fChain->SetBranchAddress("Gen_eta", &Gen_eta, &b_Gen_eta);
   fChain->SetBranchAddress("Gen_mass", &Gen_mass, &b_Gen_mass);
   fChain->SetBranchAddress("Gen_momDisplacment", &Gen_momDisplacment, &b_Gen_momDisplacment);
   fChain->SetBranchAddress("Gen_motherIdx", &Gen_motherIdx, &b_Gen_motherIdx);
   fChain->SetBranchAddress("Gen_nTotal", &Gen_nTotal, &b_Gen_nTotal);
   fChain->SetBranchAddress("Gen_pdgId", &Gen_pdgId, &b_Gen_pdgId);
   fChain->SetBranchAddress("Gen_phi", &Gen_phi, &b_Gen_phi);
   fChain->SetBranchAddress("Gen_pt", &Gen_pt, &b_Gen_pt);
   fChain->SetBranchAddress("Gen_px", &Gen_px, &b_Gen_px);
   fChain->SetBranchAddress("Gen_py", &Gen_py, &b_Gen_py);
   fChain->SetBranchAddress("Gen_pz", &Gen_pz, &b_Gen_pz);
   fChain->SetBranchAddress("Gen_status", &Gen_status, &b_Gen_status);
   fChain->SetBranchAddress("Gen_susEvtType", &Gen_susEvtType, &b_Gen_susEvtType);
   fChain->SetBranchAddress("Gen_susId", &Gen_susId, &b_Gen_susId);
   fChain->SetBranchAddress("Gen_vx", &Gen_vx, &b_Gen_vx);
   fChain->SetBranchAddress("Gen_vy", &Gen_vy, &b_Gen_vy);
   fChain->SetBranchAddress("Gen_vz", &Gen_vz, &b_Gen_vz);
   fChain->SetBranchAddress("Evt_genWgt", &Evt_genWgt, &b_Evt_genWgt);
   }//<<>>if( doGenInfo )

   fChain->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_px", &Muon_px, &b_Muon_px);
   fChain->SetBranchAddress("Muon_py", &Muon_py, &b_Muon_py);
   fChain->SetBranchAddress("Muon_pz", &Muon_pz, &b_Muon_pz);
   fChain->SetBranchAddress("Muon_isLoose", &Muon_isLoose, &b_Muon_isLoose);
   //fChain->SetBranchAddress("Muon_isMedium", &Muon_isMedium, &b_Muon_isMedium);
   fChain->SetBranchAddress("Muon_nMuons", &Muon_nMuons, &b_Muon_nMuons);
   fChain->SetBranchAddress("Muon_nSelMuons", &Muon_nSelMuons, &b_Muon_nSelMuons);

   fChain->SetBranchAddress("Jet_area", &Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", &Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_chHM", &Jet_chHM, &b_Jet_chHM);
   fChain->SetBranchAddress("Jet_drRhIds", &Jet_drRhIds, &b_Jet_drRhIds);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   if( doGenInfo ){
   fChain->SetBranchAddress("Jet_genJetLlpDr", &Jet_genJetLlpDr, &b_Jet_genJetLlpDr);
   fChain->SetBranchAddress("Jet_genDptMatch", &Jet_genDptMatch, &b_Jet_genDptMatch);
   fChain->SetBranchAddress("Jet_genDrMatch", &Jet_genDrMatch, &b_Jet_genDrMatch);
   fChain->SetBranchAddress("Jet_genEnergy", &Jet_genEnergy, &b_Jet_genEnergy);
   fChain->SetBranchAddress("Jet_genEta", &Jet_genEta, &b_Jet_genEta);
   fChain->SetBranchAddress("Jet_genImpactAngle", &Jet_genImpactAngle, &b_Jet_genImpactAngle);
   fChain->SetBranchAddress("Jet_genJetLlpId", &Jet_genJetLlpId, &b_Jet_genJetLlpId);
   fChain->SetBranchAddress("Jet_genJetLlpRp", &Jet_genJetLlpRp, &b_Jet_genJetLlpRp);
   fChain->SetBranchAddress("Jet_genPhi", &Jet_genPhi, &b_Jet_genPhi);
   fChain->SetBranchAddress("Jet_genPt", &Jet_genPt, &b_Jet_genPt);
   fChain->SetBranchAddress("Jet_genQrkLlpDr", &Jet_genQrkLlpDr, &b_Jet_genQrkLlpDr);
   fChain->SetBranchAddress("Jet_genQrkLlpId", &Jet_genQrkLlpId, &b_Jet_genQrkLlpId);
   fChain->SetBranchAddress("Jet_genQrkLlpRp", &Jet_genQrkLlpRp, &b_Jet_genQrkLlpRp);
   fChain->SetBranchAddress("Jet_genTOF", &Jet_genTOF, &b_Jet_genTOF);
   fChain->SetBranchAddress("Jet_genTime", &Jet_genTime, &b_Jet_genTime);
   fChain->SetBranchAddress("Jet_genTimeLLP", &Jet_genTimeLLP, &b_Jet_genTimeLLP);
   }//<<>>if( doGenInfo )
   fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", &Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_neEmEF", &Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", &Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_neHM", &Jet_neHM, &b_Jet_neHM);
   fChain->SetBranchAddress("Jet_egIndxs", &Jet_egIndxs, &b_Jet_egIndxs);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_nConstituents", &Jet_nConstituents, &b_Jet_nConstituents);

   fChain->SetBranchAddress("Met_CPt", &Met_CPt, &b_Met_CPt);
   fChain->SetBranchAddress("Met_Cpx", &Met_Cpx, &b_Met_Cpx);
   fChain->SetBranchAddress("Met_Cpy", &Met_Cpy, &b_Met_Cpy);
   fChain->SetBranchAddress("Met_CsumEt", &Met_CsumEt, &b_Met_CsumEt);
   fChain->SetBranchAddress("Met_eta", &Met_eta, &b_Met_eta);
   fChain->SetBranchAddress("Met_phi", &Met_phi, &b_Met_phi);
   fChain->SetBranchAddress("Met_pt", &Met_pt, &b_Met_pt);
   fChain->SetBranchAddress("Met_px", &Met_px, &b_Met_px);
   fChain->SetBranchAddress("Met_py", &Met_py, &b_Met_py);
   fChain->SetBranchAddress("Met_sumEt", &Met_sumEt, &b_Met_sumEt);

   fChain->SetBranchAddress("Photon_ecalRHSumEtConeDR04", &Photon_ecalRHSumEtConeDR04, &b_Photon_ecalRHSumEtConeDR04);
   fChain->SetBranchAddress("Photon_energy", &Photon_energy, &b_Photon_energy);
   fChain->SetBranchAddress("Photon_energyErr", &Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_excluded", &Photon_excluded, &b_Photon_excluded);
   if( doGenInfo ){
   fChain->SetBranchAddress("Photon_genIdx", &Photon_genIdx, &b_Photon_genIdx);
   fChain->SetBranchAddress("Photon_genSigXMomId", &Photon_genSigXMomId, &b_Photon_genSigXMomId);
   }//<<>>if( doGenInfo )
   fChain->SetBranchAddress("Photon_hadOverEM", &Photon_hadOverEM, &b_Photon_hadOverEM);
   fChain->SetBranchAddress("Photon_hadTowOverEM", &Photon_hadTowOverEM, &b_Photon_hadTowOverEM);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtBcConeDR04", &Photon_hcalTowerSumEtBcConeDR04, &b_Photon_hcalTowerSumEtBcConeDR04);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", &Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_isOot", &Photon_isOot, &b_Photon_isOot);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", &Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", &Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
   fChain->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_px", &Photon_px, &b_Photon_px);
   fChain->SetBranchAddress("Photon_py", &Photon_py, &b_Photon_py);
   fChain->SetBranchAddress("Photon_pz", &Photon_pz, &b_Photon_pz);
   fChain->SetBranchAddress("Photon_r9", &Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_s4", &Photon_s4, &b_Photon_s4);
   fChain->SetBranchAddress("Photon_seedTOFTime", &Photon_seedTOFTime, &b_Photon_seedTOFTime);
   fChain->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", &Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", &Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", &Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_ecalPFClusterIso", &Photon_ecalPFClusterIso, &b_Photon_ecalPFClusterIso);
   fChain->SetBranchAddress("Photon_electronVeto", &Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_gloResRhId", &Photon_gloResRhId, &b_Photon_gloResRhId);
   fChain->SetBranchAddress("hasConversionTracks", &hasConversionTracks, &b_hasConversionTracks);
   fChain->SetBranchAddress("Photon_pixelSeed", &Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", &Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_Hoe_PUcorr", &Photon_Hoe_PUcorr, &b_Photon_Hoe_PUcorr);
   fChain->SetBranchAddress("Photon_locResRhId", &Photon_locResRhId, &b_Photon_locResRhId);
   fChain->SetBranchAddress("Photon_nPhoton", &Photon_nPhoton, &b_Photon_nPhoton);
   fChain->SetBranchAddress("Photon_pfChargedIso", &Photon_pfChargedIso, &b_Photon_pfChargedIso);
   fChain->SetBranchAddress("Photon_pfChargedIsoPFPV", &Photon_pfChargedIsoPFPV, &b_Photon_pfChargedIsoPFPV);
   fChain->SetBranchAddress("Photon_pfChargedIsoWorstVtx", &Photon_pfChargedIsoWorstVtx, &b_Photon_pfChargedIsoWorstVtx);
   fChain->SetBranchAddress("Photon_pfPhoIso03", &Photon_pfPhoIso03, &b_Photon_pfPhoIso03);
   fChain->SetBranchAddress("pfRelIso03_all_quadratic", &pfRelIso03_all_quadratic, &b_pfRelIso03_all_quadratic);
   fChain->SetBranchAddress("pfRelIso03_chg_quadratic", &pfRelIso03_chg_quadratic, &b_pfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("Photon_scIndex", &Photon_scIndex, &b_Photon_scIndex);
   fChain->SetBranchAddress("Photon_sieie", &Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_sieip", &Photon_sieip, &b_Photon_sieip);
   fChain->SetBranchAddress("Photon_sipip", &Photon_sipip, &b_Photon_sipip);
   //Notify();

   fChain->SetBranchAddress("Vertex_nTracks", &Vertex_nTracks, &b_Vertex_nTracks);
   fChain->SetBranchAddress("Vertex_mass", &Vertex_mass, &b_Vertex_mass);
   fChain->SetBranchAddress("Vertex_passLooseElectronID", &Vertex_passLooseElectronID, &b_Vertex_passLooseElectronID);
   fChain->SetBranchAddress("Vertex_passLooseMuonID", &Vertex_passLooseMuonID, &b_Vertex_passLooseMuonID);

}//<<>>void llpgtree::Init( TTree *tree, bool doGenInfo )

void llpgtree::getBranches( Long64_t entry, bool doGenInfo ){

   // List of branches
/*
   b_DisplacedElectron_charge->GetEntry(entry);     //!
   b_DisplacedElectron_costSC->GetEntry(entry);     //!
   b_DisplacedElectron_dRSC->GetEntry(entry);     //!
   b_DisplacedElectron_dxy->GetEntry(entry);     //!
   b_DisplacedElectron_energy->GetEntry(entry);     //!
   b_DisplacedElectron_eta->GetEntry(entry);     //!
   b_DisplacedElectron_etaECAL->GetEntry(entry);     //!

   if( doGenInfo ){
   b_DisplacedElectron_costGen->GetEntry(entry);     //!
   b_DisplacedElectron_dRGen->GetEntry(entry);     //!
   b_DisplacedElectron_genSigIndex->GetEntry(entry);     //!
   b_DisplacedElectron_indexGen->GetEntry(entry);     //!
   b_DisplacedElectron_nGenSig->GetEntry(entry);     //!
   }//<<>>if( doGenInfo )

   b_DisplacedElectron_indexSC->GetEntry(entry);     //!
   b_DisplacedElectron_indexTrack->GetEntry(entry);     //!
   b_DisplacedElectron_ip2D->GetEntry(entry);     //!
   b_DisplacedElectron_ip2Dbs->GetEntry(entry);     //!
   b_DisplacedElectron_ip3D->GetEntry(entry);     //!
   b_DisplacedElectron_ip3Dbs->GetEntry(entry);     //!
   b_DisplacedElectron_isBarrel->GetEntry(entry);     //!
   b_DisplacedElectron_isConversion->GetEntry(entry);     //!
   b_DisplacedElectron_isEndCap->GetEntry(entry);     //!
   b_DisplacedElectron_isGeneralTrack->GetEntry(entry);     //!
   b_DisplacedElectron_isGsfTrack->GetEntry(entry);     //!
   b_DisplacedElectron_isHeavy->GetEntry(entry);     //!
   b_DisplacedElectron_isLight->GetEntry(entry);     //!
   b_DisplacedElectron_isOOTSC->GetEntry(entry);     //!
   b_DisplacedElectron_isSC->GetEntry(entry);     //!
   b_DisplacedElectron_isSignal->GetEntry(entry);     //!
   b_DisplacedElectron_isTau->GetEntry(entry);     //!
   b_DisplacedElectron_isUnmatched->GetEntry(entry);     //!
   b_DisplacedElectron_nHits->GetEntry(entry);     //!
   b_DisplacedElectron_nOOTSCs->GetEntry(entry);     //!
   b_DisplacedElectron_nSCs->GetEntry(entry);     //!
   b_DisplacedElectron_nTotal->GetEntry(entry);     //!
   b_DisplacedElectron_p->GetEntry(entry);     //!
   b_DisplacedElectron_pdgId->GetEntry(entry);     //!
   b_DisplacedElectron_pfIso0p4->GetEntry(entry);     //!
   b_DisplacedElectron_phi->GetEntry(entry);     //!
   b_DisplacedElectron_phiECAL->GetEntry(entry);     //!
   b_DisplacedElectron_pt->GetEntry(entry);     //!
   b_DisplacedElectron_scIso0p4->GetEntry(entry);     //!
   b_DisplacedElectron_sip2D->GetEntry(entry);     //!
   b_DisplacedElectron_sip2Dbs->GetEntry(entry);     //!
   b_DisplacedElectron_sip3D->GetEntry(entry);     //!
   b_DisplacedElectron_sip3Dbs->GetEntry(entry);     //!
   b_DisplacedElectron_trackIso0p4->GetEntry(entry);     //!
   b_DisplacedElectron_x->GetEntry(entry);     //!
   b_DisplacedElectron_xECAL->GetEntry(entry);     //!
   b_DisplacedElectron_y->GetEntry(entry);     //!
   b_DisplacedElectron_yECAL->GetEntry(entry);     //!
   b_DisplacedElectron_z->GetEntry(entry);     //!
   b_DisplacedElectron_zECAL->GetEntry(entry);     //!
*/

   if( sbDEBUG ) std::cout << "Getting Branches RH" << std::endl;
   b_ECALRecHit_energy->GetEntry(entry);     //!
   //b_ECALRecHit_timeError->GetEntry(entry);     //!
   //b_ECALRecHit_hasGS1->GetEntry(entry);     //!
   //b_ECALRecHit_hasGS6->GetEntry(entry);     //!
   b_ECALRecHit_ID->GetEntry(entry);     //!
   b_ECALRecHit_swCross->GetEntry(entry);     //!
   b_ECALRecHit_0TOF->GetEntry(entry);     //!
   b_ECALRecHit_pvTOF->GetEntry(entry);     //!
   b_ECALRecHit_time->GetEntry(entry);     //!
   b_ECALRecHit_amplitude->GetEntry(entry);     //!
   b_ECALRecHit_ampres->GetEntry(entry);     //!
   b_ECALRecHit_eta->GetEntry(entry);     //!
   b_ECALRecHit_isDead->GetEntry(entry);     //!
   b_ECALRecHit_isDiWrd->GetEntry(entry);     //!
   b_ECALRecHit_isOOT->GetEntry(entry);     //!
   b_ECALRecHit_isOther->GetEntry(entry);     //!
   b_ECALRecHit_isPoor->GetEntry(entry);     //!
   b_ECALRecHit_isRecov->GetEntry(entry);     //!
   b_ECALRecHit_isWrd->GetEntry(entry);     //!
   b_ECALRecHit_phi->GetEntry(entry);     //!
   b_ECALRecHit_precentUsed->GetEntry(entry);     //!
   b_ECALRecHit_rhx->GetEntry(entry);     //!
   b_ECALRecHit_rhy->GetEntry(entry);     //!
   b_ECALRecHit_rhz->GetEntry(entry);     //!

   if( sbDEBUG ) std::cout << "Getting Branches SC" << std::endl;
   b_SuperCluster_covEtaEta->GetEntry(entry);     //!
   b_SuperCluster_covEtaPhi->GetEntry(entry);     //!
   b_SuperCluster_covPhiPhi->GetEntry(entry);     //!
   b_SuperCluster_energyRaw->GetEntry(entry);     //!
   b_SuperCluster_etaWidth->GetEntry(entry);     //!
   b_SuperCluster_excluded->GetEntry(entry);     //!
   b_SuperCluster_seedIsEB->GetEntry(entry);     //!
   b_SuperCluster_isScEtaEB->GetEntry(entry);     //!
   b_SuperCluster_isScEtaEE->GetEntry(entry);     //!
   b_SuperCluster_isOot->GetEntry(entry);     //!
   b_SuperCluster_original->GetEntry(entry);     //!
   b_SuperCluster_phiWidth->GetEntry(entry);     //!
   b_SuperCluster_salp->GetEntry(entry);     //!
   b_SuperCluster_smaj->GetEntry(entry);     //!
   b_SuperCluster_smin->GetEntry(entry);     //!
   b_SuperCluster_seediEtaOriX->GetEntry(entry);     //!
   b_SuperCluster_seediPhiOriY->GetEntry(entry);     //!
   b_SuperCluster_nBasicClusters->GetEntry(entry);     //!
   b_SuperCluster_ElectronIndx->GetEntry(entry);     //!
   b_SuperCluster_energy->GetEntry(entry);     //!
   b_SuperCluster_eta->GetEntry(entry);     //!
   b_SuperCluster_clcx->GetEntry(entry);     //!
   b_SuperCluster_clcy->GetEntry(entry);     //!
   b_SuperCluster_clcz->GetEntry(entry);     //!
   b_SuperCluster_nSuperCluster->GetEntry(entry);     //!
   b_SuperCluster_nRHXtals->GetEntry(entry);     //!
   b_SuperCluster_phi->GetEntry(entry);     //!
   b_SuperCluster_PhotonIndx->GetEntry(entry);     //!
   b_SuperCluster_rhFracs->GetEntry(entry);     //!
   b_SuperCluster_rhIds->GetEntry(entry);     //!
   b_SuperCluster_XtalSeedID->GetEntry(entry);     //!
   b_SuperCluster_nHFXtals->GetEntry(entry);     //!
   b_SuperCluster_x_calo->GetEntry(entry);     //!
   b_SuperCluster_y_calo->GetEntry(entry);     //!
   b_SuperCluster_z_calo->GetEntry(entry);     //!
   if( doGenInfo ) b_SuperCluster_ObjectPdgId->GetEntry(entry);     //!

/*
   if( sbDEBUG ) std::cout << "Getting Branches TR" << std::endl;
   b_ECALTrackDetID_detId->GetEntry(entry);     //!
   b_ECALTrackDetID_isECAL->GetEntry(entry);     //!
   b_ECALTrackDetID_isHCAL->GetEntry(entry);     //!
   b_ECALTrackDetID_trackIndex->GetEntry(entry);     //!
   b_ECALTrack_charge->GetEntry(entry);     //!
   b_ECALTrack_chi2->GetEntry(entry);     //!
   b_ECALTrack_dxy->GetEntry(entry);     //!
   b_ECALTrack_eta->GetEntry(entry);     //!
   b_ECALTrack_etaError->GetEntry(entry);     //!
   b_ECALTrack_lambda->GetEntry(entry);     //!
   b_ECALTrack_nGeneral->GetEntry(entry);     //!
   b_ECALTrack_nGsf->GetEntry(entry);     //!
   b_ECALTrack_nLostHits->GetEntry(entry);     //!
   b_ECALTrack_nTracks->GetEntry(entry);     //!
   b_ECALTrack_nValidHits->GetEntry(entry);     //!
   b_ECALTrack_ndof->GetEntry(entry);     //!
   b_ECALTrack_normalizedChi2->GetEntry(entry);     //!
   b_ECALTrack_p->GetEntry(entry);     //!
   b_ECALTrack_pdgId->GetEntry(entry);     //!
   b_ECALTrack_phi->GetEntry(entry);     //!
   b_ECALTrack_phiError->GetEntry(entry);     //!
   b_ECALTrack_pt->GetEntry(entry);     //!
   b_ECALTrack_ptError->GetEntry(entry);     //!
   b_ECALTrack_px->GetEntry(entry);     //!
   b_ECALTrack_py->GetEntry(entry);     //!
   b_ECALTrack_pz->GetEntry(entry);     //!
   b_ECALTrack_qoverp->GetEntry(entry);     //!
   b_ECALTrack_qualityMask->GetEntry(entry);     //!
   b_ECALTrack_theta->GetEntry(entry);     //!
   b_ECALTrack_vx->GetEntry(entry);     //!
   b_ECALTrack_vy->GetEntry(entry);     //!
   b_ECALTrack_vz->GetEntry(entry);     //!
   b_TrackSCMatch_deltaR->GetEntry(entry);     //!
   b_TrackSCMatch_scEnergyToTrackPRatio->GetEntry(entry);     //!

   if( doGenInfo ){
   if( sbDEBUG ) std::cout << "Getting Branches GenTR" << std::endl;
   b_ECALTrack_genDeltaR->GetEntry(entry);     //!
   b_ECALTrack_genIndex->GetEntry(entry);     //!
   }//<<>>if( doGenInfo )
*/

   if( sbDEBUG ) std::cout << "Getting Branches Ele" << std::endl;
   b_Electron_energy->GetEntry(entry);     //!
   b_Electron_eta->GetEntry(entry);     //!
   b_Electron_phi->GetEntry(entry);     //!
   b_Electron_pt->GetEntry(entry);     //!
   b_Electron_px->GetEntry(entry);     //!
   b_Electron_py->GetEntry(entry);     //!
   b_Electron_pz->GetEntry(entry);     //!
   b_Electron_seedTOFTime->GetEntry(entry);     //!
   b_Electron_trackz->GetEntry(entry);     //!
   b_Electron_nElectrons->GetEntry(entry);     //!
   b_Electron_scIndex->GetEntry(entry);     //!
   //b_Electron_isLoose->GetEntry(entry);   //!
   //b_Electron_nSelElectrons->GetEntry(entry);   //!

   if( doGenInfo ){
   if( sbDEBUG ) std::cout << "Getting Branches GenEle" << std::endl;
   b_Electron_genIdx->GetEntry(entry);     //!
   b_Electron_genSigWZId->GetEntry(entry);     //!
   b_Electron_genSigXMomId->GetEntry(entry);     //!
   }//<<>>if( doGenInfo )

   if( sbDEBUG ) std::cout << "Getting Branches FLG" << std::endl;
   b_Flag_BadChargedCandidateFilter->GetEntry(entry);     //!
   b_Flag_BadPFMuonDzFilter->GetEntry(entry);     //!
   b_Flag_BadPFMuonFilter->GetEntry(entry);     //!
   b_Flag_EcalDeadCellTriggerPrimitiveFilter->GetEntry(entry);     //!
   b_Flag_HBHENoiseFilter->GetEntry(entry);     //!
   b_Flag_HBHENoiseIsoFilter->GetEntry(entry);     //!
   b_Flag_ecalBadCalibFilter->GetEntry(entry);     //!
   b_Flag_eeBadScFilter->GetEntry(entry);     //!
   b_Flag_globalSuperTightHalo2016Filter->GetEntry(entry);     //!
   b_Flag_goodVertices->GetEntry(entry);     //!
   b_Flag_hfNoisyHitsFilter->GetEntry(entry);     //!

   if( sbDEBUG ) std::cout << "Getting Branches EVT" << std::endl;
   b_Evt_luminosityBlock->GetEntry(entry);     //!
   b_Evt_run->GetEntry(entry);     //!
   b_Evt_event->GetEntry(entry);     //!

/*
   if( sbDEBUG ) std::cout << "Getting Branches TRG" << std::endl;
   b_Trigger_hltEG20EtFilterLooseHoverE->GetEntry(entry);     //!
   b_Trigger_hltEG20EtL1TripleEGFilter->GetEntry(entry);     //!
   b_Trigger_hltEG20HEFilterLooseHoverE->GetEntry(entry);     //!
   b_Trigger_hltEG20HEL1TripleEGFilter->GetEntry(entry);     //!
   b_Trigger_hltEG20L1SingleEGLowETOrEtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30EBHE10R9Id50b80eHEFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30EBR9Id50b80eR9IdFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30EIso15HE30EcalIsoLastFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30EtFilterLooseHoverE->GetEntry(entry);     //!
   b_Trigger_hltEG30EtL1TripleEGFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30HE30HEFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30HEFilterLooseHoverE->GetEntry(entry);     //!
   b_Trigger_hltEG30HEL1TripleEGFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30LHE12R9Id50b80eHEFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30LR9Id50b80eR9IdFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30PVHE10R9Id50b80eHEFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30PVR9Id50b80eR9IdFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter->GetEntry(entry);     //!
   b_Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60EtFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60HEFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter->GetEntry(entry);     //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter->GetEntry(entry);     //!
   b_Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter->GetEntry(entry);     //!
   b_Trigger_hltHT130Jet30->GetEntry(entry);     //!
   b_Trigger_hltHT175Jet10->GetEntry(entry);     //!
   b_Trigger_hltHT200Jet30->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleAndDoubleEG->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleAndDoubleEGNonIsoOr->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleAndDoubleEGor->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG15->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG18->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG24->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG26->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG34to45->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG34to50->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEG40to50->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS->GetEntry(entry);     //!
   b_Trigger_hltL1sSingleEGor->GetEntry(entry);     //!
   b_Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG->GetEntry(entry);     //!
   b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter->GetEntry(entry);     //!
   b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt->GetEntry(entry);     //!
   b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt->GetEntry(entry);     //!
   b_Trigger_hltPFHT180Jet30->GetEntry(entry);     //!
   b_Trigger_hltPFHT250Jet30->GetEntry(entry);     //!
   b_Trigger_hltPFHT350Jet15->GetEntry(entry);     //!
   b_Trigger_hltPFMET100->GetEntry(entry);     //!
   b_Trigger_hltPFMET110->GetEntry(entry);     //!
   b_Trigger_hltPFMET120->GetEntry(entry);     //!
   b_Trigger_hltPFMET130->GetEntry(entry);     //!
   b_Trigger_hltPFMET140->GetEntry(entry);     //!
   b_Trigger_hltPFMET200->GetEntry(entry);     //!
   b_Trigger_hltPFMET250->GetEntry(entry);     //!
   b_Trigger_hltPFMET300->GetEntry(entry);     //!
   b_Trigger_hltPFMET50->GetEntry(entry);     //!
   b_Trigger_hltPFMET70->GetEntry(entry);     //!
   b_Trigger_hltPFMET90->GetEntry(entry);     //!
   b_Trigger_hltPFMETNoMu100->GetEntry(entry);     //!
   b_Trigger_hltPFMETNoMu110->GetEntry(entry);     //!
   b_Trigger_hltPFMETNoMu120->GetEntry(entry);     //!
   b_Trigger_hltPFMETNoMu130->GetEntry(entry);     //!
   b_Trigger_hltPFMETNoMu140->GetEntry(entry);     //!
   b_Trigger_hltPFMETNoMu60->GetEntry(entry);     //!   
   b_Trigger_hltPFMHTNoMuTightID100->GetEntry(entry);     //!
   b_Trigger_hltPFMHTNoMuTightID110->GetEntry(entry);     //!
   b_Trigger_hltPFMHTNoMuTightID120->GetEntry(entry);     //!
   b_Trigger_hltPFMHTNoMuTightID130->GetEntry(entry);     //!
   b_Trigger_hltPFMHTNoMuTightID140->GetEntry(entry);     //!
   b_Trigger_hltPFMHTNoMuTightID70->GetEntry(entry);     //!
   b_Trigger_hltPFMHTNoMuTightID90->GetEntry(entry);     //!
   b_Trigger_hltPFMHTTightID120->GetEntry(entry);     //!
   b_Trigger_hltPFMHTTightID130->GetEntry(entry);     //!
   b_Trigger_hltPFMHTTightID140->GetEntry(entry);     //!
*/

   if( sbDEBUG ) std::cout << "Getting Branches PV" << std::endl;
   b_PV_npvs->GetEntry(entry);     //!
   b_PV_x->GetEntry(entry);     //!
   b_PV_y->GetEntry(entry);     //!
   b_PV_z->GetEntry(entry);     //!

   if( doGenInfo ){
   if( sbDEBUG ) std::cout << "Getting Branches Gen" << std::endl;
   b_Gen_charge->GetEntry(entry);     //!
   b_Gen_energy->GetEntry(entry);     //!
   b_Gen_eta->GetEntry(entry);     //!
   b_Gen_mass->GetEntry(entry);     //!
   b_Gen_momDisplacment->GetEntry(entry);     //!
   b_Gen_motherIdx->GetEntry(entry);     //!
   b_Gen_nTotal->GetEntry(entry);     //!
   b_Gen_pdgId->GetEntry(entry);     //!
   b_Gen_phi->GetEntry(entry);     //!
   b_Gen_pt->GetEntry(entry);     //!
   b_Gen_px->GetEntry(entry);     //!
   b_Gen_py->GetEntry(entry);     //!
   b_Gen_pz->GetEntry(entry);     //!
   b_Gen_status->GetEntry(entry);     //!
   b_Gen_susEvtType->GetEntry(entry);     //!
   b_Gen_susId->GetEntry(entry);     //!
   b_Gen_vx->GetEntry(entry);     //!
   b_Gen_vy->GetEntry(entry);     //!
   b_Gen_vz->GetEntry(entry);     //!
   b_Evt_genWgt->GetEntry(entry);     //!
   }//<<>>if( doGenInfo )

   if( sbDEBUG ) std::cout << "Getting Branches Muons" << std::endl;
   b_Muon_energy->GetEntry(entry);   //!
   b_Muon_eta->GetEntry(entry);   //!
   b_Muon_phi->GetEntry(entry);   //!
   b_Muon_pt->GetEntry(entry);   //!
   b_Muon_px->GetEntry(entry);   //!
   b_Muon_py->GetEntry(entry);   //!
   b_Muon_pz->GetEntry(entry);   //!
   b_Muon_isLoose->GetEntry(entry);   //!
   //b_Muon_isMedium->GetEntry(entry);   //!
   b_Muon_nMuons->GetEntry(entry);   //!
   b_Muon_nSelMuons->GetEntry(entry);   //!

   if( sbDEBUG ) std::cout << "Getting Branches Jet" << std::endl;
   b_Jet_area->GetEntry(entry);     //!
   b_Jet_chEmEF->GetEntry(entry);     //!
   b_Jet_chHEF->GetEntry(entry);     //!
   b_Jet_chHM->GetEntry(entry);     //!
   b_Jet_drRhIds->GetEntry(entry);     //!
   b_Jet_energy->GetEntry(entry);     //!
   b_Jet_eta->GetEntry(entry);     //!
   b_Jet_mass->GetEntry(entry);     //!
   b_Jet_muEF->GetEntry(entry);     //!
   b_Jet_neEmEF->GetEntry(entry);     //!
   b_Jet_neHEF->GetEntry(entry);     //!
   b_Jet_neHM->GetEntry(entry);     //!
   b_Jet_egIndxs->GetEntry(entry);     //!
   b_Jet_phi->GetEntry(entry);     //!
   b_Jet_pt->GetEntry(entry);     //!
   b_Jet_nConstituents->GetEntry(entry);     //!

   if( doGenInfo ){
   if( sbDEBUG ) std::cout << "Getting Branches GenJet" << std::endl;
   b_Jet_genJetLlpDr->GetEntry(entry);     //!
   b_Jet_genDptMatch->GetEntry(entry);     //!
   b_Jet_genDrMatch->GetEntry(entry);     //!
   b_Jet_genEnergy->GetEntry(entry);     //!
   b_Jet_genEta->GetEntry(entry);     //!
   b_Jet_genImpactAngle->GetEntry(entry);     //!
   b_Jet_genJetLlpId->GetEntry(entry);     //!
   b_Jet_genJetLlpRp->GetEntry(entry);     //!
   b_Jet_genPhi->GetEntry(entry);     //!
   b_Jet_genPt->GetEntry(entry);     //!
   b_Jet_genQrkLlpDr->GetEntry(entry);     //!
   b_Jet_genQrkLlpId->GetEntry(entry);     //!
   b_Jet_genQrkLlpRp->GetEntry(entry);     //!
   b_Jet_genTOF->GetEntry(entry);     //!
   b_Jet_genTime->GetEntry(entry);     //!
   b_Jet_genTimeLLP->GetEntry(entry);     //!
   }//<<>>if( doGenInfo )

   if( sbDEBUG ) std::cout << "Getting Branches Met" << std::endl;
   b_Met_CPt->GetEntry(entry);     //!
   b_Met_Cpx->GetEntry(entry);     //!
   b_Met_Cpy->GetEntry(entry);     //!
   b_Met_CsumEt->GetEntry(entry);     //!
   b_Met_eta->GetEntry(entry);     //!
   b_Met_phi->GetEntry(entry);     //!
   b_Met_pt->GetEntry(entry);     //!
   b_Met_px->GetEntry(entry);     //!
   b_Met_py->GetEntry(entry);     //!
   b_Met_sumEt->GetEntry(entry);     //!

   if( sbDEBUG ) std::cout << "Getting Branches Pho" << std::endl;
   b_Photon_ecalRHSumEtConeDR04->GetEntry(entry);     //!
   b_Photon_energy->GetEntry(entry);     //!
   b_Photon_energyErr->GetEntry(entry);     //!
   b_Photon_eta->GetEntry(entry);     //!
   b_Photon_excluded->GetEntry(entry);     //!
   b_Photon_hadOverEM->GetEntry(entry);     //!
   b_Photon_hadTowOverEM->GetEntry(entry);     //!
   b_Photon_hcalTowerSumEtBcConeDR04->GetEntry(entry);     //!
   b_Photon_hcalTowerSumEtConeDR04->GetEntry(entry);     //!
   b_Photon_isOot->GetEntry(entry);     //!
   b_Photon_nTrkHollowConeDR04->GetEntry(entry);     //!
   b_Photon_nTrkSolidConeDR04->GetEntry(entry);     //!
   b_Photon_phi->GetEntry(entry);     //!
   b_Photon_pt->GetEntry(entry);     //!
   b_Photon_px->GetEntry(entry);     //!
   b_Photon_py->GetEntry(entry);     //!
   b_Photon_pz->GetEntry(entry);     //!
   b_Photon_r9->GetEntry(entry);     //!
   b_Photon_s4->GetEntry(entry);     //!
   b_Photon_seedTOFTime->GetEntry(entry);     //!
   b_Photon_SigmaIEtaIEta->GetEntry(entry);     //!
   b_Photon_trkSumPtHollowConeDR03->GetEntry(entry);     //!
   b_Photon_trkSumPtHollowConeDR04->GetEntry(entry);     //!
   b_Photon_trkSumPtSolidConeDR04->GetEntry(entry);     //!
   b_Photon_ecalPFClusterIso->GetEntry(entry);     //!
   b_Photon_electronVeto->GetEntry(entry);     //!
   b_Photon_gloResRhId->GetEntry(entry);     //!
   b_hasConversionTracks->GetEntry(entry);     //!
   b_Photon_pixelSeed->GetEntry(entry);     //!
   b_Photon_hcalPFClusterIso->GetEntry(entry);     //!
   b_Photon_Hoe_PUcorr->GetEntry(entry);     //!
   b_Photon_locResRhId->GetEntry(entry);     //!
   b_Photon_nPhoton->GetEntry(entry);     //!
   b_Photon_pfChargedIso->GetEntry(entry);     //!
   b_Photon_pfChargedIsoPFPV->GetEntry(entry);     //!
   b_Photon_pfChargedIsoWorstVtx->GetEntry(entry);     //!
   b_Photon_pfPhoIso03->GetEntry(entry);     //!
   b_pfRelIso03_all_quadratic->GetEntry(entry);     //!
   b_pfRelIso03_chg_quadratic->GetEntry(entry);     //!
   b_Photon_scIndex->GetEntry(entry);     //!
   b_Photon_sieie->GetEntry(entry);     //!
   b_Photon_sieip->GetEntry(entry);     //!
   b_Photon_sipip->GetEntry(entry);     //!

   if( doGenInfo ){
   if( sbDEBUG ) std::cout << "Getting Branches GenPho" << std::endl;
   b_Photon_genIdx->GetEntry(entry);     //!
   b_Photon_genSigXMomId->GetEntry(entry);     //!
   }//<<>>if( doGenInfo )

   b_Vertex_nTracks->GetEntry(entry);   //!
   b_Vertex_mass->GetEntry(entry);   //!
   b_Vertex_passLooseElectronID->GetEntry(entry);   //!
   b_Vertex_passLooseMuonID->GetEntry(entry);   //!

}//<<>>void root_base::getBranches(Long64_t entry)

/*
Bool_t llpgtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void llpgtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t llpgtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/

#endif //<<>> define llpgtree_h
