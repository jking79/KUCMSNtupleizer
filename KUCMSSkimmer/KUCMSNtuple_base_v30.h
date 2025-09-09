//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 12 11:59:16 2025 by ROOT version 6.26/11
// from TTree llpgtree/KUCMSNtuple
// found on file: kucmsntuple.root
//////////////////////////////////////////////////////////

#ifndef kuntuple_v30_h
#define kuntuple_v30_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

//#define sbDEBUG true
#define sbDEBUG false

class llpgtree {

public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<float>   *GenVertex_dxy;
   std::vector<float>   *GenVertex_eta;
   std::vector<bool>    *GenVertex_isElectron;
   std::vector<bool>    *GenVertex_isHadronic;
   std::vector<bool>    *GenVertex_isMuon;
   std::vector<float>   *GenVertex_mass;
   UInt_t          GenVertex_nElectron;
   UInt_t          GenVertex_nHadronic;
   UInt_t          GenVertex_nMuon;
   UInt_t          GenVertex_nTotal;
   std::vector<unsigned int> *GenVertex_nTracks;
   std::vector<float>   *GenVertex_p;
   std::vector<bool>    *GenVertex_passSelection;
   std::vector<bool>    *GenVertex_passSelectionAndCuts;
   std::vector<float>   *GenVertex_phi;
   std::vector<float>   *GenVertex_pt;
   std::vector<float>   *GenVertex_px;
   std::vector<float>   *GenVertex_py;
   std::vector<float>   *GenVertex_pz;
   std::vector<float>   *GenVertex_x;
   std::vector<float>   *GenVertex_y;
   std::vector<float>   *GenVertex_z;
   std::vector<float>   *VertexTrack_SCDR;
   std::vector<float>   *VertexTrack_energySC;
   std::vector<bool>    *VertexTrack_isSignalElectron;
   std::vector<bool>    *VertexTrack_isSignalMuon;
   std::vector<bool>    *VertexTrack_isSignalTrack;
   std::vector<float>   *VertexTrack_ratioPToEnergySC;
   std::vector<float>   *VertexTrack_shift3DAfterTrackRemoval;
   std::vector<float>   *VertexTrack_shiftDzAfterTrackRemoval;
   std::vector<float>   *VertexTrack_trackCompatibility;
   std::vector<float>   *VertexTrack_trackCosTheta;
   std::vector<float>   *VertexTrack_trackCosThetaAtCM;
   std::vector<unsigned int> *VertexTrack_trackIndex;
   std::vector<unsigned int> *VertexTrack_vertexIndex;
   std::vector<float>   *Vertex_chi2;
   std::vector<float>   *Vertex_cosTheta;
   std::vector<float>   *Vertex_cxx;
   std::vector<float>   *Vertex_cxy;
   std::vector<float>   *Vertex_cxz;
   std::vector<float>   *Vertex_cyy;
   std::vector<float>   *Vertex_cyz;
   std::vector<float>   *Vertex_czz;
   std::vector<float>   *Vertex_decayAngle;
   std::vector<float>   *Vertex_dxy;
   std::vector<float>   *Vertex_dxyError;
   std::vector<float>   *Vertex_eta;
   std::vector<int>     *Vertex_genVertexIndex;
   std::vector<bool>    *Vertex_isBronze;
   std::vector<bool>    *Vertex_isGold;
   std::vector<bool>    *Vertex_isSilver;
   std::vector<bool>    *Vertex_isUnique;
   std::vector<float>   *Vertex_mass;
   std::vector<float>   *Vertex_matchRatio;
   std::vector<float>   *Vertex_min3D;
   UInt_t          Vertex_nTotal;
   std::vector<unsigned int> *Vertex_nTracks;
   std::vector<float>   *Vertex_ndof;
   std::vector<int>     *Vertex_nearestGenVertexIndex;
   std::vector<float>   *Vertex_normalizedChi2;
   std::vector<float>   *Vertex_p;
   std::vector<bool>    *Vertex_passLooseElectronID;
   std::vector<bool>    *Vertex_passLooseID;
   std::vector<bool>    *Vertex_passLooseMuonID;
   std::vector<float>   *Vertex_phi;
   std::vector<float>   *Vertex_pt;
   std::vector<float>   *Vertex_px;
   std::vector<float>   *Vertex_py;
   std::vector<float>   *Vertex_pz;
   std::vector<float>   *Vertex_scMatchRatio;
   std::vector<int>     *Vertex_sumCharge;
   std::vector<float>   *Vertex_x;
   std::vector<float>   *Vertex_y;
   std::vector<float>   *Vertex_z;
   std::vector<float>   *ECALRecHit_energy;
   std::vector<unsigned int> *ECALRecHit_ID;
   std::vector<float>   *ECALRecHit_swCross;
   std::vector<float>   *ECALRecHit_0TOF;
   std::vector<float>   *ECALRecHit_pvTOF;
   std::vector<float>   *ECALRecHit_time;
   std::vector<float>   *ECALRecHit_timeError;
   std::vector<float>   *ECALRecHit_amplitude;
   std::vector<float>   *ECALRecHit_ampres;
   std::vector<float>   *ECALRecHit_eta;
   std::vector<bool>    *ECALRecHit_hasGS1;
   std::vector<bool>    *ECALRecHit_hasGS6;
   std::vector<bool>    *ECALRecHit_isDead;
   std::vector<bool>    *ECALRecHit_isDiWrd;
   std::vector<bool>    *ECALRecHit_isOOT;
   std::vector<bool>    *ECALRecHit_isOther;
   std::vector<bool>    *ECALRecHit_isPoor;
   std::vector<bool>    *ECALRecHit_isRecov;
   std::vector<bool>    *ECALRecHit_isTimeValid;
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
   std::vector<float>   *SuperCluster_SeedTime;
   std::vector<std::vector<int> > *SuperCluster_TrackIndx;
   std::vector<float>   *SuperCluster_x_calo;
   std::vector<float>   *SuperCluster_y_calo;
   std::vector<float>   *SuperCluster_z_calo;
   std::vector<float>   *Electron_DetaSCTV;
   std::vector<float>   *Electron_DphiSCTV;
   std::vector<float>   *Electron_energy;
   std::vector<float>   *Electron_eta;
   std::vector<int>     *Electron_genIdx;
   std::vector<int>     *Electron_genSigWZId;
   std::vector<int>     *Electron_genSigXMomId;
   std::vector<float>   *Electron_HOE;
   std::vector<bool>    *Electron_isLVeto;
   std::vector<bool>    *Electron_isLoose;
   std::vector<bool>    *Electron_isMedium;
   std::vector<bool>    *Electron_isPrompt;
   std::vector<bool>    *Electron_isVeto;
   std::vector<float>   *Electron_phi;
   std::vector<float>   *Electron_pt;
   std::vector<float>   *Electron_px;
   std::vector<float>   *Electron_py;
   std::vector<float>   *Electron_pz;
   std::vector<float>   *Electron_seedTOFTime;
   std::vector<float>   *Electron_Sieie;
   std::vector<float>   *Electron_trackz;
   std::vector<bool>    *Electron_Vdetaseed;
   std::vector<bool>    *Electron_Vdpscvc;
   std::vector<bool>    *Electron_Vhoec;
   std::vector<bool>    *Electron_Vmistrack;
   std::vector<bool>    *Electron_Vooemoop;
   std::vector<bool>    *Electron_VrelIsoWithEA;
   std::vector<bool>    *Electron_Vsieiec;
   std::vector<bool>    *Electron_hasGenMatch;
   std::vector<bool>    *Electron_isLWZP;
   Int_t           Electron_nElectrons;
   Int_t           Electron_nGenMatch;
   Int_t           Electron_nPrompt;
   Int_t           Electron_nSVMatched;
   Int_t           Electron_nSelElectrons;
   std::vector<bool>    *Electron_passConversionVeto;
   std::vector<bool>    *Electron_passPrompt;
   std::vector<int>     *Electron_scIndex;
   std::vector<bool>    *Electron_hasSVMatch;
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
   std::vector<int>     *Gen_GMomPdgId;
   std::vector<bool>    *Gen_isLWZQ;
   std::vector<bool>    *Gen_isLWZX;
   std::vector<float>   *Gen_mass;
   std::vector<float>   *Gen_matchDe;
   std::vector<float>   *Gen_matchDr;
   std::vector<float>   *Gen_momDisplacment;
   std::vector<int>     *Gen_motherIdx;
   std::vector<int>     *Gen_MomPdgId;
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
   Bool_t          Evt_hasLWZQ;
   Bool_t          Evt_hasLWZX;
   Int_t           Evt_nLWZX;
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
   std::vector<float>   *Muon_energy;
   std::vector<float>   *Muon_eta;
   std::vector<int>     *Muon_genIdx;
   std::vector<int>     *Muon_genSigWZId;
   std::vector<int>     *Muon_genSigXMomId;
   std::vector<bool>    *Muon_isPromptMuon;
   std::vector<float>   *Muon_phi;
   std::vector<float>   *Muon_pt;
   std::vector<float>   *Muon_px;
   std::vector<float>   *Muon_py;
   std::vector<float>   *Muon_pz;
   std::vector<bool>    *Muon_hasGenMatch;
   std::vector<bool>    *Muon_isLWZP;
   std::vector<bool>    *Muon_isLoose;
   std::vector<bool>    *Muon_isMedium;
   Int_t           Muon_nGenMatch;
   Int_t           Muon_nMuons;
   Int_t           Muon_nPrompt;
   Int_t           Muon_nSelMuons;
   Int_t           Muon_nSVMatched;
   std::vector<bool>    *Muon_passPrompt;
   std::vector<bool>    *Muon_hasSVMatch;
   std::vector<bool>    *Muon_isTight;
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
   std::vector<bool>    *hasConversionTracks;
   std::vector<bool>    *Photon_pixelSeed;
   std::vector<float>   *Photon_hcalPFClusterIso;
   std::vector<float>   *Photon_Hoe_PUcorr;
   Int_t           Photon_nPhoton;
   std::vector<bool>    *Photon_isSelect;
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
   std::vector<float>   *Track_absIso;
   std::vector<int>     *Track_charge;
   std::vector<float>   *Track_chi2;
   std::vector<float>   *Track_dxy;
   std::vector<std::vector<unsigned int> > *Track_ecalDetIds;
   std::vector<float>   *Track_eta;
   std::vector<float>   *Track_etaError;
   std::vector<std::vector<unsigned int> > *Track_hcalDetIds;
   std::vector<float>   *Track_ip2D;
   std::vector<float>   *Track_ip3D;
   std::vector<bool>    *Track_isGenTrack;
   std::vector<float>   *Track_lambda;
   std::vector<float>   *Track_miniRelIso;
   std::vector<int>     *Track_pdgId;
   std::vector<float>   *Track_mrelisolep;
   Int_t           Track_nIsoTracks;
   std::vector<int>     *Track_nLostHits;
   std::vector<int>     *Track_nValidHits;
   std::vector<float>   *Track_ndof;
   std::vector<float>   *Track_normalizedChi2;
   std::vector<float>   *Track_p;
   std::vector<float>   *Track_phi;
   std::vector<float>   *Track_phiError;
   std::vector<float>   *Track_pt;
   std::vector<float>   *Track_ptError;
   std::vector<float>   *Track_px;
   std::vector<float>   *Track_py;
   std::vector<float>   *Track_pz;
   std::vector<int>     *Track_qualityMask;
   std::vector<float>   *Track_relIso;
   std::vector<float>   *Track_relisolep;
   std::vector<std::vector<int> > *Track_scIndexs;
   std::vector<float>   *Track_sip2D;
   std::vector<float>   *Track_sip3D;
   std::vector<float>   *Track_theta;
   std::vector<float>   *Track_vx;
   std::vector<float>   *Track_vy;
   std::vector<float>   *Track_vz;

   // List of branches
   TBranch        *b_GenVertex_dxy;   //!
   TBranch        *b_GenVertex_eta;   //!
   TBranch        *b_GenVertex_isElectron;   //!
   TBranch        *b_GenVertex_isHadronic;   //!
   TBranch        *b_GenVertex_isMuon;   //!
   TBranch        *b_GenVertex_mass;   //!
   TBranch        *b_GenVertex_nElectron;   //!
   TBranch        *b_GenVertex_nHadronic;   //!
   TBranch        *b_GenVertex_nMuon;   //!
   TBranch        *b_GenVertex_nTotal;   //!
   TBranch        *b_GenVertex_nTracks;   //!
   TBranch        *b_GenVertex_p;   //!
   TBranch        *b_GenVertex_passSelection;   //!
   TBranch        *b_GenVertex_passSelectionAndCuts;   //!
   TBranch        *b_GenVertex_phi;   //!
   TBranch        *b_GenVertex_pt;   //!
   TBranch        *b_GenVertex_px;   //!
   TBranch        *b_GenVertex_py;   //!
   TBranch        *b_GenVertex_pz;   //!
   TBranch        *b_GenVertex_x;   //!
   TBranch        *b_GenVertex_y;   //!
   TBranch        *b_GenVertex_z;   //!
   TBranch        *b_VertexTrack_SCDR;   //!
   TBranch        *b_VertexTrack_energySC;   //!
   TBranch        *b_VertexTrack_isSignalElectron;   //!
   TBranch        *b_VertexTrack_isSignalMuon;   //!
   TBranch        *b_VertexTrack_isSignalTrack;   //!
   TBranch        *b_VertexTrack_ratioPToEnergySC;   //!
   TBranch        *b_VertexTrack_shift3DAfterTrackRemoval;   //!
   TBranch        *b_VertexTrack_shiftDzAfterTrackRemoval;   //!
   TBranch        *b_VertexTrack_trackCompatibility;   //!
   TBranch        *b_VertexTrack_trackCosTheta;   //!
   TBranch        *b_VertexTrack_trackCosThetaAtCM;   //!
   TBranch        *b_VertexTrack_trackIndex;   //!
   TBranch        *b_VertexTrack_vertexIndex;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_cosTheta;   //!
   TBranch        *b_Vertex_cxx;   //!
   TBranch        *b_Vertex_cxy;   //!
   TBranch        *b_Vertex_cxz;   //!
   TBranch        *b_Vertex_cyy;   //!
   TBranch        *b_Vertex_cyz;   //!
   TBranch        *b_Vertex_czz;   //!
   TBranch        *b_Vertex_decayAngle;   //!
   TBranch        *b_Vertex_dxy;   //!
   TBranch        *b_Vertex_dxyError;   //!
   TBranch        *b_Vertex_eta;   //!
   TBranch        *b_Vertex_genVertexIndex;   //!
   TBranch        *b_Vertex_isBronze;   //!
   TBranch        *b_Vertex_isGold;   //!
   TBranch        *b_Vertex_isSilver;   //!
   TBranch        *b_Vertex_isUnique;   //!
   TBranch        *b_Vertex_mass;   //!
   TBranch        *b_Vertex_matchRatio;   //!
   TBranch        *b_Vertex_min3D;   //!
   TBranch        *b_Vertex_nTotal;   //!
   TBranch        *b_Vertex_nTracks;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_nearestGenVertexIndex;   //!
   TBranch        *b_Vertex_normalizedChi2;   //!
   TBranch        *b_Vertex_p;   //!
   TBranch        *b_Vertex_passLooseElectronID;   //!
   TBranch        *b_Vertex_passLooseID;   //!
   TBranch        *b_Vertex_passLooseMuonID;   //!
   TBranch        *b_Vertex_phi;   //!
   TBranch        *b_Vertex_pt;   //!
   TBranch        *b_Vertex_px;   //!
   TBranch        *b_Vertex_py;   //!
   TBranch        *b_Vertex_pz;   //!
   TBranch        *b_Vertex_scMatchRatio;   //!
   TBranch        *b_Vertex_sumCharge;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_ECALRecHit_energy;   //!
   TBranch        *b_ECALRecHit_ID;   //!
   TBranch        *b_ECALRecHit_swCross;   //!
   TBranch        *b_ECALRecHit_0TOF;   //!
   TBranch        *b_ECALRecHit_pvTOF;   //!
   TBranch        *b_ECALRecHit_time;   //!
   TBranch        *b_ECALRecHit_timeError;   //!
   TBranch        *b_ECALRecHit_amplitude;   //!
   TBranch        *b_ECALRecHit_ampres;   //!
   TBranch        *b_ECALRecHit_eta;   //!
   TBranch        *b_ECALRecHit_hasGS1;   //!
   TBranch        *b_ECALRecHit_hasGS6;   //!
   TBranch        *b_ECALRecHit_isDead;   //!
   TBranch        *b_ECALRecHit_isDiWrd;   //!
   TBranch        *b_ECALRecHit_isOOT;   //!
   TBranch        *b_ECALRecHit_isOther;   //!
   TBranch        *b_ECALRecHit_isPoor;   //!
   TBranch        *b_ECALRecHit_isRecov;   //!
   TBranch        *b_ECALRecHit_isTimeValid;   //!
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
   TBranch        *b_SuperCluster_SeedTime;   //!
   TBranch        *b_SuperCluster_TrackIndx;   //!
   TBranch        *b_SuperCluster_x_calo;   //!
   TBranch        *b_SuperCluster_y_calo;   //!
   TBranch        *b_SuperCluster_z_calo;   //!
   TBranch        *b_Electron_DetaSCTV;   //!
   TBranch        *b_Electron_DphiSCTV;   //!
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_genIdx;   //!
   TBranch        *b_Electron_genSigWZId;   //!
   TBranch        *b_Electron_genSigXMomId;   //!
   TBranch        *b_Electron_HOE;   //!
   TBranch        *b_Electron_isLVeto;   //!
   TBranch        *b_Electron_isLoose;   //!
   TBranch        *b_Electron_isMedium;   //!
   TBranch        *b_Electron_isPrompt;   //!
   TBranch        *b_Electron_isVeto;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_px;   //!
   TBranch        *b_Electron_py;   //!
   TBranch        *b_Electron_pz;   //!
   TBranch        *b_Electron_seedTOFTime;   //!
   TBranch        *b_Electron_Sieie;   //!
   TBranch        *b_Electron_trackz;   //!
   TBranch        *b_Electron_Vdetaseed;   //!
   TBranch        *b_Electron_Vdpscvc;   //!
   TBranch        *b_Electron_Vhoec;   //!
   TBranch        *b_Electron_Vmistrack;   //!
   TBranch        *b_Electron_Vooemoop;   //!
   TBranch        *b_Electron_VrelIsoWithEA;   //!
   TBranch        *b_Electron_Vsieiec;   //!
   TBranch        *b_Electron_hasGenMatch;   //!
   TBranch        *b_Electron_isLWZP;   //!
   TBranch        *b_Electron_nElectrons;   //!
   TBranch        *b_Electron_nGenMatch;   //!
   TBranch        *b_Electron_nPrompt;   //!
   TBranch        *b_Electron_nSVMatched;   //!
   TBranch        *b_Electron_nSelElectrons;   //!
   TBranch        *b_Electron_passConversionVeto;   //!
   TBranch        *b_Electron_passPrompt;   //!
   TBranch        *b_Electron_scIndex;   //!
   TBranch        *b_Electron_hasSVMatch;   //!
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
   TBranch        *b_Gen_GMomPdgId;   //!
   TBranch        *b_Gen_isLWZQ;   //!
   TBranch        *b_Gen_isLWZX;   //!
   TBranch        *b_Gen_mass;   //!
   TBranch        *b_Gen_matchDe;   //!
   TBranch        *b_Gen_matchDr;   //!
   TBranch        *b_Gen_momDisplacment;   //!
   TBranch        *b_Gen_motherIdx;   //!
   TBranch        *b_Gen_MomPdgId;   //!
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
   TBranch        *b_Evt_hasLWZQ;   //!
   TBranch        *b_Evt_hasLWZX;   //!
   TBranch        *b_Evt_nLWZX;   //!
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
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_genIdx;   //!
   TBranch        *b_Muon_genSigWZId;   //!
   TBranch        *b_Muon_genSigXMomId;   //!
   TBranch        *b_Muon_isPromptMuon;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_px;   //!
   TBranch        *b_Muon_py;   //!
   TBranch        *b_Muon_pz;   //!
   TBranch        *b_Muon_hasGenMatch;   //!
   TBranch        *b_Muon_isLWZP;   //!
   TBranch        *b_Muon_isLoose;   //!
   TBranch        *b_Muon_isMedium;   //!
   TBranch        *b_Muon_nGenMatch;   //!
   TBranch        *b_Muon_nMuons;   //!
   TBranch        *b_Muon_nPrompt;   //!
   TBranch        *b_Muon_nSelMuons;   //!
   TBranch        *b_Muon_nSVMatched;   //!
   TBranch        *b_Muon_passPrompt;   //!
   TBranch        *b_Muon_hasSVMatch;   //!
   TBranch        *b_Muon_isTight;   //!
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
   TBranch        *b_hasConversionTracks;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_hcalPFClusterIso;   //!
   TBranch        *b_Photon_Hoe_PUcorr;   //!
   TBranch        *b_Photon_nPhoton;   //!
   TBranch        *b_Photon_isSelect;   //!
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
   TBranch        *b_Track_absIso;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_Track_dxy;   //!
   TBranch        *b_Track_ecalDetIds;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_etaError;   //!
   TBranch        *b_Track_hcalDetIds;   //!
   TBranch        *b_Track_ip2D;   //!
   TBranch        *b_Track_ip3D;   //!
   TBranch        *b_Track_isGenTrack;   //!
   TBranch        *b_Track_lambda;   //!
   TBranch        *b_Track_miniRelIso;   //!
   TBranch        *b_Track_pdgId;   //!
   TBranch        *b_Track_mrelisolep;   //!
   TBranch        *b_Track_nIsoTracks;   //!
   TBranch        *b_Track_nLostHits;   //!
   TBranch        *b_Track_nValidHits;   //!
   TBranch        *b_Track_ndof;   //!
   TBranch        *b_Track_normalizedChi2;   //!
   TBranch        *b_Track_p;   //!
   TBranch        *b_Track_phi;   //!
   TBranch        *b_Track_phiError;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_ptError;   //!
   TBranch        *b_Track_px;   //!
   TBranch        *b_Track_py;   //!
   TBranch        *b_Track_pz;   //!
   TBranch        *b_Track_qualityMask;   //!
   TBranch        *b_Track_relIso;   //!
   TBranch        *b_Track_relisolep;   //!
   TBranch        *b_Track_scIndexs;   //!
   TBranch        *b_Track_sip2D;   //!
   TBranch        *b_Track_sip3D;   //!
   TBranch        *b_Track_theta;   //!
   TBranch        *b_Track_vx;   //!
   TBranch        *b_Track_vy;   //!
   TBranch        *b_Track_vz;   //!

   //kuntuple_v30(TTree *tree=0);
   //virtual ~kuntuple_v30();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   //virtual void     Init(TTree *tree);
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);

   void Init( TChain *tree, bool doGenInfo );
   void getBranches( Long64_t entry, bool doGenInfo );

};

/*

kuntuple_v30::kuntuple_v30(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("kucmsntuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("kucmsntuple.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("kucmsntuple.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

kuntuple_v30::~kuntuple_v30()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kuntuple_v30::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kuntuple_v30::LoadTree(Long64_t entry)
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
   GenVertex_dxy = 0;
   GenVertex_eta = 0;
   GenVertex_isElectron = 0;
   GenVertex_isHadronic = 0;
   GenVertex_isMuon = 0;
   GenVertex_mass = 0;
   GenVertex_nTracks = 0;
   GenVertex_p = 0;
   GenVertex_passSelection = 0;
   GenVertex_passSelectionAndCuts = 0;
   GenVertex_phi = 0;
   GenVertex_pt = 0;
   GenVertex_px = 0;
   GenVertex_py = 0;
   GenVertex_pz = 0;
   GenVertex_x = 0;
   GenVertex_y = 0;
   GenVertex_z = 0;
   VertexTrack_SCDR = 0;
   VertexTrack_energySC = 0;
   VertexTrack_isSignalElectron = 0;
   VertexTrack_isSignalMuon = 0;
   VertexTrack_isSignalTrack = 0;
   VertexTrack_ratioPToEnergySC = 0;
   VertexTrack_shift3DAfterTrackRemoval = 0;
   VertexTrack_shiftDzAfterTrackRemoval = 0;
   VertexTrack_trackCompatibility = 0;
   VertexTrack_trackCosTheta = 0;
   VertexTrack_trackCosThetaAtCM = 0;
   VertexTrack_trackIndex = 0;
   VertexTrack_vertexIndex = 0;
   Vertex_chi2 = 0;
   Vertex_cosTheta = 0;
   Vertex_cxx = 0;
   Vertex_cxy = 0;
   Vertex_cxz = 0;
   Vertex_cyy = 0;
   Vertex_cyz = 0;
   Vertex_czz = 0;
   Vertex_decayAngle = 0;
   Vertex_dxy = 0;
   Vertex_dxyError = 0;
   Vertex_eta = 0;
   Vertex_genVertexIndex = 0;
   Vertex_isBronze = 0;
   Vertex_isGold = 0;
   Vertex_isSilver = 0;
   Vertex_isUnique = 0;
   Vertex_mass = 0;
   Vertex_matchRatio = 0;
   Vertex_min3D = 0;
   Vertex_nTracks = 0;
   Vertex_ndof = 0;
   Vertex_nearestGenVertexIndex = 0;
   Vertex_normalizedChi2 = 0;
   Vertex_p = 0;
   Vertex_passLooseElectronID = 0;
   Vertex_passLooseID = 0;
   Vertex_passLooseMuonID = 0;
   Vertex_phi = 0;
   Vertex_pt = 0;
   Vertex_px = 0;
   Vertex_py = 0;
   Vertex_pz = 0;
   Vertex_scMatchRatio = 0;
   Vertex_sumCharge = 0;
   Vertex_x = 0;
   Vertex_y = 0;
   Vertex_z = 0;
   ECALRecHit_energy = 0;
   ECALRecHit_ID = 0;
   ECALRecHit_swCross = 0;
   ECALRecHit_0TOF = 0;
   ECALRecHit_pvTOF = 0;
   ECALRecHit_time = 0;
   ECALRecHit_timeError = 0;
   ECALRecHit_amplitude = 0;
   ECALRecHit_ampres = 0;
   ECALRecHit_eta = 0;
   ECALRecHit_hasGS1 = 0;
   ECALRecHit_hasGS6 = 0;
   ECALRecHit_isDead = 0;
   ECALRecHit_isDiWrd = 0;
   ECALRecHit_isOOT = 0;
   ECALRecHit_isOther = 0;
   ECALRecHit_isPoor = 0;
   ECALRecHit_isRecov = 0;
   ECALRecHit_isTimeValid = 0;
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
   SuperCluster_SeedTime = 0;
   SuperCluster_TrackIndx = 0;
   SuperCluster_x_calo = 0;
   SuperCluster_y_calo = 0;
   SuperCluster_z_calo = 0;
   Electron_DetaSCTV = 0;
   Electron_DphiSCTV = 0;
   Electron_energy = 0;
   Electron_eta = 0;
   Electron_genIdx = 0;
   Electron_genSigWZId = 0;
   Electron_genSigXMomId = 0;
   Electron_HOE = 0;
   Electron_isLVeto = 0;
   Electron_isLoose = 0;
   Electron_isMedium = 0;
   Electron_isPrompt = 0;
   Electron_isVeto = 0;
   Electron_phi = 0;
   Electron_pt = 0;
   Electron_px = 0;
   Electron_py = 0;
   Electron_pz = 0;
   Electron_seedTOFTime = 0;
   Electron_Sieie = 0;
   Electron_trackz = 0;
   Electron_Vdetaseed = 0;
   Electron_Vdpscvc = 0;
   Electron_Vhoec = 0;
   Electron_Vmistrack = 0;
   Electron_Vooemoop = 0;
   Electron_VrelIsoWithEA = 0;
   Electron_Vsieiec = 0;
   Electron_hasGenMatch = 0;
   Electron_isLWZP = 0;
   Electron_passConversionVeto = 0;
   Electron_passPrompt = 0;
   Electron_scIndex = 0;
   Electron_hasSVMatch = 0;
   Gen_charge = 0;
   Gen_energy = 0;
   Gen_eta = 0;
   Gen_GMomPdgId = 0;
   Gen_isLWZQ = 0;
   Gen_isLWZX = 0;
   Gen_mass = 0;
   Gen_matchDe = 0;
   Gen_matchDr = 0;
   Gen_momDisplacment = 0;
   Gen_motherIdx = 0;
   Gen_MomPdgId = 0;
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
   Muon_energy = 0;
   Muon_eta = 0;
   Muon_genIdx = 0;
   Muon_genSigWZId = 0;
   Muon_genSigXMomId = 0;
   Muon_isPromptMuon = 0;
   Muon_phi = 0;
   Muon_pt = 0;
   Muon_px = 0;
   Muon_py = 0;
   Muon_pz = 0;
   Muon_hasGenMatch = 0;
   Muon_isLWZP = 0;
   Muon_isLoose = 0;
   Muon_isMedium = 0;
   Muon_passPrompt = 0;
   Muon_hasSVMatch = 0;
   Muon_isTight = 0;
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
   hasConversionTracks = 0;
   Photon_pixelSeed = 0;
   Photon_hcalPFClusterIso = 0;
   Photon_Hoe_PUcorr = 0;
   Photon_isSelect = 0;
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
   Track_absIso = 0;
   Track_charge = 0;
   Track_chi2 = 0;
   Track_dxy = 0;
   Track_ecalDetIds = 0;
   Track_eta = 0;
   Track_etaError = 0;
   Track_hcalDetIds = 0;
   Track_ip2D = 0;
   Track_ip3D = 0;
   Track_isGenTrack = 0;
   Track_lambda = 0;
   Track_miniRelIso = 0;
   Track_pdgId = 0;
   Track_mrelisolep = 0;
   Track_nLostHits = 0;
   Track_nValidHits = 0;
   Track_ndof = 0;
   Track_normalizedChi2 = 0;
   Track_p = 0;
   Track_phi = 0;
   Track_phiError = 0;
   Track_pt = 0;
   Track_ptError = 0;
   Track_px = 0;
   Track_py = 0;
   Track_pz = 0;
   Track_normalizedChi2 = 0;
   Track_qualityMask = 0;
   Track_relIso = 0;
   Track_relisolep = 0;
   Track_scIndexs = 0;
   Track_sip2D = 0;
   Track_sip3D = 0;
   Track_theta = 0;
   Track_vx = 0;
   Track_vy = 0;
   Track_vz = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);

   if( doGenInfo ){
   fChain->SetBranchAddress("GenVertex_dxy", &GenVertex_dxy, &b_GenVertex_dxy);
   fChain->SetBranchAddress("GenVertex_eta", &GenVertex_eta, &b_GenVertex_eta);
   fChain->SetBranchAddress("GenVertex_isElectron", &GenVertex_isElectron, &b_GenVertex_isElectron);
   fChain->SetBranchAddress("GenVertex_isHadronic", &GenVertex_isHadronic, &b_GenVertex_isHadronic);
   fChain->SetBranchAddress("GenVertex_isMuon", &GenVertex_isMuon, &b_GenVertex_isMuon);
   fChain->SetBranchAddress("GenVertex_mass", &GenVertex_mass, &b_GenVertex_mass);
   fChain->SetBranchAddress("GenVertex_nElectron", &GenVertex_nElectron, &b_GenVertex_nElectron);
   fChain->SetBranchAddress("GenVertex_nHadronic", &GenVertex_nHadronic, &b_GenVertex_nHadronic);
   fChain->SetBranchAddress("GenVertex_nMuon", &GenVertex_nMuon, &b_GenVertex_nMuon);
   fChain->SetBranchAddress("GenVertex_nTotal", &GenVertex_nTotal, &b_GenVertex_nTotal);
   fChain->SetBranchAddress("GenVertex_nTracks", &GenVertex_nTracks, &b_GenVertex_nTracks);
   fChain->SetBranchAddress("GenVertex_p", &GenVertex_p, &b_GenVertex_p);
   fChain->SetBranchAddress("GenVertex_passSelection", &GenVertex_passSelection, &b_GenVertex_passSelection);
   fChain->SetBranchAddress("GenVertex_passSelectionAndCuts", &GenVertex_passSelectionAndCuts, &b_GenVertex_passSelectionAndCuts);
   fChain->SetBranchAddress("GenVertex_phi", &GenVertex_phi, &b_GenVertex_phi);
   fChain->SetBranchAddress("GenVertex_pt", &GenVertex_pt, &b_GenVertex_pt);
   fChain->SetBranchAddress("GenVertex_px", &GenVertex_px, &b_GenVertex_px);
   fChain->SetBranchAddress("GenVertex_py", &GenVertex_py, &b_GenVertex_py);
   fChain->SetBranchAddress("GenVertex_pz", &GenVertex_pz, &b_GenVertex_pz);
   fChain->SetBranchAddress("GenVertex_x", &GenVertex_x, &b_GenVertex_x);
   fChain->SetBranchAddress("GenVertex_y", &GenVertex_y, &b_GenVertex_y);
   fChain->SetBranchAddress("GenVertex_z", &GenVertex_z, &b_GenVertex_z);
   }//<<>>if( doGenInfo )

   fChain->SetBranchAddress("VertexTrack_SCDR", &VertexTrack_SCDR, &b_VertexTrack_SCDR);
   fChain->SetBranchAddress("VertexTrack_energySC", &VertexTrack_energySC, &b_VertexTrack_energySC);
   fChain->SetBranchAddress("VertexTrack_isSignalElectron", &VertexTrack_isSignalElectron, &b_VertexTrack_isSignalElectron);
   fChain->SetBranchAddress("VertexTrack_isSignalMuon", &VertexTrack_isSignalMuon, &b_VertexTrack_isSignalMuon);
   fChain->SetBranchAddress("VertexTrack_isSignalTrack", &VertexTrack_isSignalTrack, &b_VertexTrack_isSignalTrack);
   fChain->SetBranchAddress("VertexTrack_ratioPToEnergySC", &VertexTrack_ratioPToEnergySC, &b_VertexTrack_ratioPToEnergySC);
   fChain->SetBranchAddress("VertexTrack_shift3DAfterTrackRemoval", &VertexTrack_shift3DAfterTrackRemoval, &b_VertexTrack_shift3DAfterTrackRemoval);
   fChain->SetBranchAddress("VertexTrack_shiftDzAfterTrackRemoval", &VertexTrack_shiftDzAfterTrackRemoval, &b_VertexTrack_shiftDzAfterTrackRemoval);
   fChain->SetBranchAddress("VertexTrack_trackCompatibility", &VertexTrack_trackCompatibility, &b_VertexTrack_trackCompatibility);
   fChain->SetBranchAddress("VertexTrack_trackCosTheta", &VertexTrack_trackCosTheta, &b_VertexTrack_trackCosTheta);
   fChain->SetBranchAddress("VertexTrack_trackCosThetaAtCM", &VertexTrack_trackCosThetaAtCM, &b_VertexTrack_trackCosThetaAtCM);
   fChain->SetBranchAddress("VertexTrack_trackIndex", &VertexTrack_trackIndex, &b_VertexTrack_trackIndex);
   fChain->SetBranchAddress("VertexTrack_vertexIndex", &VertexTrack_vertexIndex, &b_VertexTrack_vertexIndex);

   fChain->SetBranchAddress("Vertex_chi2", &Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_cosTheta", &Vertex_cosTheta, &b_Vertex_cosTheta);
   fChain->SetBranchAddress("Vertex_cxx", &Vertex_cxx, &b_Vertex_cxx);
   fChain->SetBranchAddress("Vertex_cxy", &Vertex_cxy, &b_Vertex_cxy);
   fChain->SetBranchAddress("Vertex_cxz", &Vertex_cxz, &b_Vertex_cxz);
   fChain->SetBranchAddress("Vertex_cyy", &Vertex_cyy, &b_Vertex_cyy);
   fChain->SetBranchAddress("Vertex_cyz", &Vertex_cyz, &b_Vertex_cyz);
   fChain->SetBranchAddress("Vertex_czz", &Vertex_czz, &b_Vertex_czz);
   fChain->SetBranchAddress("Vertex_decayAngle", &Vertex_decayAngle, &b_Vertex_decayAngle);
   fChain->SetBranchAddress("Vertex_dxy", &Vertex_dxy, &b_Vertex_dxy);
   fChain->SetBranchAddress("Vertex_dxyError", &Vertex_dxyError, &b_Vertex_dxyError);
   fChain->SetBranchAddress("Vertex_eta", &Vertex_eta, &b_Vertex_eta);
   if( doGenInfo ) fChain->SetBranchAddress("Vertex_genVertexIndex", &Vertex_genVertexIndex, &b_Vertex_genVertexIndex);
   fChain->SetBranchAddress("Vertex_isBronze", &Vertex_isBronze, &b_Vertex_isBronze);
   fChain->SetBranchAddress("Vertex_isGold", &Vertex_isGold, &b_Vertex_isGold);
   fChain->SetBranchAddress("Vertex_isSilver", &Vertex_isSilver, &b_Vertex_isSilver);
   fChain->SetBranchAddress("Vertex_isUnique", &Vertex_isUnique, &b_Vertex_isUnique);
   fChain->SetBranchAddress("Vertex_mass", &Vertex_mass, &b_Vertex_mass);
   fChain->SetBranchAddress("Vertex_matchRatio", &Vertex_matchRatio, &b_Vertex_matchRatio);
   fChain->SetBranchAddress("Vertex_min3D", &Vertex_min3D, &b_Vertex_min3D);
   fChain->SetBranchAddress("Vertex_nTotal", &Vertex_nTotal, &b_Vertex_nTotal);
   fChain->SetBranchAddress("Vertex_nTracks", &Vertex_nTracks, &b_Vertex_nTracks);
   fChain->SetBranchAddress("Vertex_ndof", &Vertex_ndof, &b_Vertex_ndof);
   if( doGenInfo ) fChain->SetBranchAddress("Vertex_nearestGenVertexIndex", &Vertex_nearestGenVertexIndex, &b_Vertex_nearestGenVertexIndex);
   fChain->SetBranchAddress("Vertex_normalizedChi2", &Vertex_normalizedChi2, &b_Vertex_normalizedChi2);
   fChain->SetBranchAddress("Vertex_p", &Vertex_p, &b_Vertex_p);
   fChain->SetBranchAddress("Vertex_passLooseElectronID", &Vertex_passLooseElectronID, &b_Vertex_passLooseElectronID);
   fChain->SetBranchAddress("Vertex_passLooseID", &Vertex_passLooseID, &b_Vertex_passLooseID);
   fChain->SetBranchAddress("Vertex_passLooseMuonID", &Vertex_passLooseMuonID, &b_Vertex_passLooseMuonID);
   fChain->SetBranchAddress("Vertex_phi", &Vertex_phi, &b_Vertex_phi);
   fChain->SetBranchAddress("Vertex_pt", &Vertex_pt, &b_Vertex_pt);
   fChain->SetBranchAddress("Vertex_px", &Vertex_px, &b_Vertex_px);
   fChain->SetBranchAddress("Vertex_py", &Vertex_py, &b_Vertex_py);
   fChain->SetBranchAddress("Vertex_pz", &Vertex_pz, &b_Vertex_pz);
   fChain->SetBranchAddress("Vertex_scMatchRatio", &Vertex_scMatchRatio, &b_Vertex_scMatchRatio);
   fChain->SetBranchAddress("Vertex_sumCharge", &Vertex_sumCharge, &b_Vertex_sumCharge);
   fChain->SetBranchAddress("Vertex_x", &Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", &Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", &Vertex_z, &b_Vertex_z);

   fChain->SetBranchAddress("ECALRecHit_energy", &ECALRecHit_energy, &b_ECALRecHit_energy);
   fChain->SetBranchAddress("ECALRecHit_ID", &ECALRecHit_ID, &b_ECALRecHit_ID);
   fChain->SetBranchAddress("ECALRecHit_swCross", &ECALRecHit_swCross, &b_ECALRecHit_swCross);
   fChain->SetBranchAddress("ECALRecHit_0TOF", &ECALRecHit_0TOF, &b_ECALRecHit_0TOF);
   fChain->SetBranchAddress("ECALRecHit_pvTOF", &ECALRecHit_pvTOF, &b_ECALRecHit_pvTOF);
   fChain->SetBranchAddress("ECALRecHit_time", &ECALRecHit_time, &b_ECALRecHit_time);
   fChain->SetBranchAddress("ECALRecHit_timeError", &ECALRecHit_timeError, &b_ECALRecHit_timeError);
   fChain->SetBranchAddress("ECALRecHit_amplitude", &ECALRecHit_amplitude, &b_ECALRecHit_amplitude);
   fChain->SetBranchAddress("ECALRecHit_ampres", &ECALRecHit_ampres, &b_ECALRecHit_ampres);
   fChain->SetBranchAddress("ECALRecHit_eta", &ECALRecHit_eta, &b_ECALRecHit_eta);
   fChain->SetBranchAddress("ECALRecHit_hasGS1", &ECALRecHit_hasGS1, &b_ECALRecHit_hasGS1);
   fChain->SetBranchAddress("ECALRecHit_hasGS6", &ECALRecHit_hasGS6, &b_ECALRecHit_hasGS6);
   fChain->SetBranchAddress("ECALRecHit_isDead", &ECALRecHit_isDead, &b_ECALRecHit_isDead);
   fChain->SetBranchAddress("ECALRecHit_isDiWrd", &ECALRecHit_isDiWrd, &b_ECALRecHit_isDiWrd);
   fChain->SetBranchAddress("ECALRecHit_isOOT", &ECALRecHit_isOOT, &b_ECALRecHit_isOOT);
   fChain->SetBranchAddress("ECALRecHit_isOther", &ECALRecHit_isOther, &b_ECALRecHit_isOther);
   fChain->SetBranchAddress("ECALRecHit_isPoor", &ECALRecHit_isPoor, &b_ECALRecHit_isPoor);
   fChain->SetBranchAddress("ECALRecHit_isRecov", &ECALRecHit_isRecov, &b_ECALRecHit_isRecov);
   fChain->SetBranchAddress("ECALRecHit_isTimeValid", &ECALRecHit_isTimeValid, &b_ECALRecHit_isTimeValid);
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
   fChain->SetBranchAddress("SuperCluster_ObjectPdgId", &SuperCluster_ObjectPdgId, &b_SuperCluster_ObjectPdgId);
   fChain->SetBranchAddress("SuperCluster_phi", &SuperCluster_phi, &b_SuperCluster_phi);
   fChain->SetBranchAddress("SuperCluster_PhotonIndx", &SuperCluster_PhotonIndx, &b_SuperCluster_PhotonIndx);
   fChain->SetBranchAddress("SuperCluster_rhFracs", &SuperCluster_rhFracs, &b_SuperCluster_rhFracs);
   fChain->SetBranchAddress("SuperCluster_rhIds", &SuperCluster_rhIds, &b_SuperCluster_rhIds);
   fChain->SetBranchAddress("SuperCluster_XtalSeedID", &SuperCluster_XtalSeedID, &b_SuperCluster_XtalSeedID);
   fChain->SetBranchAddress("SuperCluster_nHFXtals", &SuperCluster_nHFXtals, &b_SuperCluster_nHFXtals);
   fChain->SetBranchAddress("SuperCluster_SeedTime", &SuperCluster_SeedTime, &b_SuperCluster_SeedTime);
   fChain->SetBranchAddress("SuperCluster_TrackIndx", &SuperCluster_TrackIndx, &b_SuperCluster_TrackIndx);
   fChain->SetBranchAddress("SuperCluster_x_calo", &SuperCluster_x_calo, &b_SuperCluster_x_calo);
   fChain->SetBranchAddress("SuperCluster_y_calo", &SuperCluster_y_calo, &b_SuperCluster_y_calo);
   fChain->SetBranchAddress("SuperCluster_z_calo", &SuperCluster_z_calo, &b_SuperCluster_z_calo);

   fChain->SetBranchAddress("Electron_DetaSCTV", &Electron_DetaSCTV, &b_Electron_DetaSCTV);
   fChain->SetBranchAddress("Electron_DphiSCTV", &Electron_DphiSCTV, &b_Electron_DphiSCTV);
   fChain->SetBranchAddress("Electron_energy", &Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   if( doGenInfo ){
   fChain->SetBranchAddress("Electron_genIdx", &Electron_genIdx, &b_Electron_genIdx);
   fChain->SetBranchAddress("Electron_genSigWZId", &Electron_genSigWZId, &b_Electron_genSigWZId);
   fChain->SetBranchAddress("Electron_genSigXMomId", &Electron_genSigXMomId, &b_Electron_genSigXMomId);
   fChain->SetBranchAddress("Electron_hasGenMatch", &Electron_hasGenMatch, &b_Electron_hasGenMatch);
   fChain->SetBranchAddress("Electron_isLWZP", &Electron_isLWZP, &b_Electron_isLWZP);
   fChain->SetBranchAddress("Electron_nGenMatch", &Electron_nGenMatch, &b_Electron_nGenMatch);
   }//<<>>if( doGenInfo )
   fChain->SetBranchAddress("Electron_HOE", &Electron_HOE, &b_Electron_HOE);
   fChain->SetBranchAddress("Electron_isLVeto", &Electron_isLVeto, &b_Electron_isLVeto);
   fChain->SetBranchAddress("Electron_isLoose", &Electron_isLoose, &b_Electron_isLoose);
   fChain->SetBranchAddress("Electron_isMedium", &Electron_isMedium, &b_Electron_isMedium);
   fChain->SetBranchAddress("Electron_isPrompt", &Electron_isPrompt, &b_Electron_isPrompt);
   fChain->SetBranchAddress("Electron_isVeto", &Electron_isVeto, &b_Electron_isVeto);
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_px", &Electron_px, &b_Electron_px);
   fChain->SetBranchAddress("Electron_py", &Electron_py, &b_Electron_py);
   fChain->SetBranchAddress("Electron_pz", &Electron_pz, &b_Electron_pz);
   fChain->SetBranchAddress("Electron_seedTOFTime", &Electron_seedTOFTime, &b_Electron_seedTOFTime);
   fChain->SetBranchAddress("Electron_Sieie", &Electron_Sieie, &b_Electron_Sieie);
   fChain->SetBranchAddress("Electron_trackz", &Electron_trackz, &b_Electron_trackz);
   fChain->SetBranchAddress("Electron_Vdetaseed", &Electron_Vdetaseed, &b_Electron_Vdetaseed);
   fChain->SetBranchAddress("Electron_Vdpscvc", &Electron_Vdpscvc, &b_Electron_Vdpscvc);
   fChain->SetBranchAddress("Electron_Vhoec", &Electron_Vhoec, &b_Electron_Vhoec);
   fChain->SetBranchAddress("Electron_Vmistrack", &Electron_Vmistrack, &b_Electron_Vmistrack);
   fChain->SetBranchAddress("Electron_Vooemoop", &Electron_Vooemoop, &b_Electron_Vooemoop);
   fChain->SetBranchAddress("Electron_VrelIsoWithEA", &Electron_VrelIsoWithEA, &b_Electron_VrelIsoWithEA);
   fChain->SetBranchAddress("Electron_Vsieiec", &Electron_Vsieiec, &b_Electron_Vsieiec);
   fChain->SetBranchAddress("Electron_nElectrons", &Electron_nElectrons, &b_Electron_nElectrons);
   fChain->SetBranchAddress("Electron_nPrompt", &Electron_nPrompt, &b_Electron_nPrompt);
   fChain->SetBranchAddress("Electron_nSVMatched", &Electron_nSVMatched, &b_Electron_nSVMatched);
   fChain->SetBranchAddress("Electron_nSelElectrons", &Electron_nSelElectrons, &b_Electron_nSelElectrons);
   fChain->SetBranchAddress("Electron_passConversionVeto", &Electron_passConversionVeto, &b_Electron_passConversionVeto);
   fChain->SetBranchAddress("Electron_passPrompt", &Electron_passPrompt, &b_Electron_passPrompt);
   fChain->SetBranchAddress("Electron_scIndex", &Electron_scIndex, &b_Electron_scIndex);
   fChain->SetBranchAddress("Electron_hasSVMatch", &Electron_hasSVMatch, &b_Electron_hasSVMatch);

   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   //fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
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

   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);

   if( doGenInfo ){
   fChain->SetBranchAddress("Gen_charge", &Gen_charge, &b_Gen_charge);
   fChain->SetBranchAddress("Gen_energy", &Gen_energy, &b_Gen_energy);
   fChain->SetBranchAddress("Gen_eta", &Gen_eta, &b_Gen_eta);
   fChain->SetBranchAddress("Gen_GMomPdgId", &Gen_GMomPdgId, &b_Gen_GMomPdgId);
   fChain->SetBranchAddress("Gen_isLWZQ", &Gen_isLWZQ, &b_Gen_isLWZQ);
   fChain->SetBranchAddress("Gen_isLWZX", &Gen_isLWZX, &b_Gen_isLWZX);
   fChain->SetBranchAddress("Gen_mass", &Gen_mass, &b_Gen_mass);
   fChain->SetBranchAddress("Gen_matchDe", &Gen_matchDe, &b_Gen_matchDe);
   fChain->SetBranchAddress("Gen_matchDr", &Gen_matchDr, &b_Gen_matchDr);
   fChain->SetBranchAddress("Gen_momDisplacment", &Gen_momDisplacment, &b_Gen_momDisplacment);
   fChain->SetBranchAddress("Gen_motherIdx", &Gen_motherIdx, &b_Gen_motherIdx);
   fChain->SetBranchAddress("Gen_MomPdgId", &Gen_MomPdgId, &b_Gen_MomPdgId);
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
   }//<<>>if( doGenInfo )

   if( doGenInfo ){
   fChain->SetBranchAddress("Evt_genWgt", &Evt_genWgt, &b_Evt_genWgt);
   fChain->SetBranchAddress("Evt_hasLWZQ", &Evt_hasLWZQ, &b_Evt_hasLWZQ);
   fChain->SetBranchAddress("Evt_hasLWZX", &Evt_hasLWZX, &b_Evt_hasLWZX);
   fChain->SetBranchAddress("Evt_nLWZX", &Evt_nLWZX, &b_Evt_nLWZX);
   }//<<>>if( doGenInfo )

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

   fChain->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   if( doGenInfo ){
   fChain->SetBranchAddress("Muon_genIdx", &Muon_genIdx, &b_Muon_genIdx);
   fChain->SetBranchAddress("Muon_genSigWZId", &Muon_genSigWZId, &b_Muon_genSigWZId);
   fChain->SetBranchAddress("Muon_genSigXMomId", &Muon_genSigXMomId, &b_Muon_genSigXMomId);
   fChain->SetBranchAddress("Muon_hasGenMatch", &Muon_hasGenMatch, &b_Muon_hasGenMatch);
   }//<<>>if( doGenInfo )
   fChain->SetBranchAddress("Muon_isPromptMuon", &Muon_isPromptMuon, &b_Muon_isPromptMuon);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_px", &Muon_px, &b_Muon_px);
   fChain->SetBranchAddress("Muon_py", &Muon_py, &b_Muon_py);
   fChain->SetBranchAddress("Muon_pz", &Muon_pz, &b_Muon_pz);
   if( doGenInfo ) fChain->SetBranchAddress("Muon_isLWZP", &Muon_isLWZP, &b_Muon_isLWZP);
   fChain->SetBranchAddress("Muon_isLoose", &Muon_isLoose, &b_Muon_isLoose);
   fChain->SetBranchAddress("Muon_isMedium", &Muon_isMedium, &b_Muon_isMedium);
   fChain->SetBranchAddress("Muon_nGenMatch", &Muon_nGenMatch, &b_Muon_nGenMatch);
   fChain->SetBranchAddress("Muon_nMuons", &Muon_nMuons, &b_Muon_nMuons);
   fChain->SetBranchAddress("Muon_nPrompt", &Muon_nPrompt, &b_Muon_nPrompt);
   fChain->SetBranchAddress("Muon_nSelMuons", &Muon_nSelMuons, &b_Muon_nSelMuons);
   fChain->SetBranchAddress("Muon_nSVMatched", &Muon_nSVMatched, &b_Muon_nSVMatched);
   fChain->SetBranchAddress("Muon_passPrompt", &Muon_passPrompt, &b_Muon_passPrompt);
   fChain->SetBranchAddress("Muon_hasSVMatch", &Muon_hasSVMatch, &b_Muon_hasSVMatch);
   fChain->SetBranchAddress("Muon_isTight", &Muon_isTight, &b_Muon_isTight);

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
   fChain->SetBranchAddress("hasConversionTracks", &hasConversionTracks, &b_hasConversionTracks);
   fChain->SetBranchAddress("Photon_pixelSeed", &Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", &Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_Hoe_PUcorr", &Photon_Hoe_PUcorr, &b_Photon_Hoe_PUcorr);
   fChain->SetBranchAddress("Photon_nPhoton", &Photon_nPhoton, &b_Photon_nPhoton);
   fChain->SetBranchAddress("Photon_isSelect", &Photon_isSelect, &b_Photon_isSelect);
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

   fChain->SetBranchAddress("Track_absIso", &Track_absIso, &b_Track_absIso);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("Track_dxy", &Track_dxy, &b_Track_dxy);
   fChain->SetBranchAddress("Track_ecalDetIds", &Track_ecalDetIds, &b_Track_ecalDetIds);
   fChain->SetBranchAddress("Track_eta", &Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_etaError", &Track_etaError, &b_Track_etaError);
   fChain->SetBranchAddress("Track_hcalDetIds", &Track_hcalDetIds, &b_Track_hcalDetIds);
   fChain->SetBranchAddress("Track_ip2D", &Track_ip2D, &b_Track_ip2D);
   fChain->SetBranchAddress("Track_ip3D", &Track_ip3D, &b_Track_ip3D);
   fChain->SetBranchAddress("Track_isGenTrack", &Track_isGenTrack, &b_Track_isGenTrack);
   fChain->SetBranchAddress("Track_lambda", &Track_lambda, &b_Track_lambda);
   fChain->SetBranchAddress("Track_miniRelIso", &Track_miniRelIso, &b_Track_miniRelIso);
   fChain->SetBranchAddress("Track_pdgId", &Track_pdgId, &b_Track_pdgId);
   fChain->SetBranchAddress("Track_mrelisolep", &Track_mrelisolep, &b_Track_mrelisolep);
   fChain->SetBranchAddress("Track_nIsoTracks", &Track_nIsoTracks, &b_Track_nIsoTracks);
   fChain->SetBranchAddress("Track_nLostHits", &Track_nLostHits, &b_Track_nLostHits);
   fChain->SetBranchAddress("Track_nValidHits", &Track_nValidHits, &b_Track_nValidHits);
   fChain->SetBranchAddress("Track_ndof", &Track_ndof, &b_Track_ndof);
   fChain->SetBranchAddress("Track_normalizedChi2", &Track_normalizedChi2, &b_Track_normalizedChi2);
   fChain->SetBranchAddress("Track_p", &Track_p, &b_Track_p);
   fChain->SetBranchAddress("Track_phi", &Track_phi, &b_Track_phi);
   fChain->SetBranchAddress("Track_phiError", &Track_phiError, &b_Track_phiError);
   fChain->SetBranchAddress("Track_pt", &Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_ptError", &Track_ptError, &b_Track_ptError);
   fChain->SetBranchAddress("Track_px", &Track_px, &b_Track_px);
   fChain->SetBranchAddress("Track_py", &Track_py, &b_Track_py);
   fChain->SetBranchAddress("Track_pz", &Track_pz, &b_Track_pz);
//    fChain->SetBranchAddress("Track_normalizedChi2", &Track_normalizedChi2, &b_Track_normalizedChi2);
   fChain->SetBranchAddress("Track_qualityMask", &Track_qualityMask, &b_Track_qualityMask);
   fChain->SetBranchAddress("Track_relIso", &Track_relIso, &b_Track_relIso);
   fChain->SetBranchAddress("Track_relisolep", &Track_relisolep, &b_Track_relisolep);
   fChain->SetBranchAddress("Track_scIndexs", &Track_scIndexs, &b_Track_scIndexs);
   fChain->SetBranchAddress("Track_sip2D", &Track_sip2D, &b_Track_sip2D);
   fChain->SetBranchAddress("Track_sip3D", &Track_sip3D, &b_Track_sip3D);
   fChain->SetBranchAddress("Track_theta", &Track_theta, &b_Track_theta);
   fChain->SetBranchAddress("Track_vx", &Track_vx, &b_Track_vx);
   fChain->SetBranchAddress("Track_vy", &Track_vy, &b_Track_vy);
   fChain->SetBranchAddress("Track_vz", &Track_vz, &b_Track_vz);
   //Notify();
}

void llpgtree::getBranches( Long64_t entry, bool doGenInfo ){

   // List of branches

   if( doGenInfo ){
   b_GenVertex_dxy->GetEntry(entry);   //!
   b_GenVertex_eta->GetEntry(entry);   //!
   b_GenVertex_isElectron->GetEntry(entry);   //!
   b_GenVertex_isHadronic->GetEntry(entry);   //!
   b_GenVertex_isMuon->GetEntry(entry);   //!
   b_GenVertex_mass->GetEntry(entry);   //!
   b_GenVertex_nElectron->GetEntry(entry);   //!
   b_GenVertex_nHadronic->GetEntry(entry);   //!
   b_GenVertex_nMuon->GetEntry(entry);   //!
   b_GenVertex_nTotal->GetEntry(entry);   //!
   b_GenVertex_nTracks->GetEntry(entry);   //!
   b_GenVertex_p->GetEntry(entry);   //!
   b_GenVertex_passSelection->GetEntry(entry);   //!
   b_GenVertex_passSelectionAndCuts->GetEntry(entry);   //!
   b_GenVertex_phi->GetEntry(entry);   //!
   b_GenVertex_pt->GetEntry(entry);   //!
   b_GenVertex_px->GetEntry(entry);   //!
   b_GenVertex_py->GetEntry(entry);   //!
   b_GenVertex_pz->GetEntry(entry);   //!
   b_GenVertex_x->GetEntry(entry);   //!
   b_GenVertex_y->GetEntry(entry);   //!
   b_GenVertex_z->GetEntry(entry);   //!
   }//<<>>if( doGenInfo )

   b_VertexTrack_SCDR->GetEntry(entry);   //!
   b_VertexTrack_energySC->GetEntry(entry);   //!
   b_VertexTrack_isSignalElectron->GetEntry(entry);   //!
   b_VertexTrack_isSignalMuon->GetEntry(entry);   //!
   b_VertexTrack_isSignalTrack->GetEntry(entry);   //!
   b_VertexTrack_ratioPToEnergySC->GetEntry(entry);   //!
   b_VertexTrack_shift3DAfterTrackRemoval->GetEntry(entry);   //!
   b_VertexTrack_shiftDzAfterTrackRemoval->GetEntry(entry);   //!
   b_VertexTrack_trackCompatibility->GetEntry(entry);   //!
   b_VertexTrack_trackCosTheta->GetEntry(entry);   //!
   b_VertexTrack_trackCosThetaAtCM->GetEntry(entry);   //!
   b_VertexTrack_trackIndex->GetEntry(entry);   //!
   b_VertexTrack_vertexIndex->GetEntry(entry);   //!

   b_Vertex_chi2->GetEntry(entry);   //!
   b_Vertex_cosTheta->GetEntry(entry);   //!
   b_Vertex_cxx->GetEntry(entry);   //!
   b_Vertex_cxy->GetEntry(entry);   //!
   b_Vertex_cxz->GetEntry(entry);   //!
   b_Vertex_cyy->GetEntry(entry);   //!
   b_Vertex_cyz->GetEntry(entry);   //!
   b_Vertex_czz->GetEntry(entry);   //!
   b_Vertex_decayAngle->GetEntry(entry);   //!
   b_Vertex_dxy->GetEntry(entry);   //!
   b_Vertex_dxyError->GetEntry(entry);   //!
   b_Vertex_eta->GetEntry(entry);   //!
   if( doGenInfo ) b_Vertex_genVertexIndex->GetEntry(entry);   //!
   b_Vertex_isBronze->GetEntry(entry);   //!
   b_Vertex_isGold->GetEntry(entry);   //!
   b_Vertex_isSilver->GetEntry(entry);   //!
   b_Vertex_isUnique->GetEntry(entry);   //!
   b_Vertex_mass->GetEntry(entry);   //!
   b_Vertex_matchRatio->GetEntry(entry);   //!
   b_Vertex_min3D->GetEntry(entry);   //!
   b_Vertex_nTotal->GetEntry(entry);   //!
   b_Vertex_nTracks->GetEntry(entry);   //!
   b_Vertex_ndof->GetEntry(entry);   //!
   if( doGenInfo ) b_Vertex_nearestGenVertexIndex->GetEntry(entry);   //!
   b_Vertex_normalizedChi2->GetEntry(entry);   //!
   b_Vertex_p->GetEntry(entry);   //!
   b_Vertex_passLooseElectronID->GetEntry(entry);   //!
   b_Vertex_passLooseID->GetEntry(entry);   //!
   b_Vertex_passLooseMuonID->GetEntry(entry);   //!
   b_Vertex_phi->GetEntry(entry);   //!
   b_Vertex_pt->GetEntry(entry);   //!
   b_Vertex_px->GetEntry(entry);   //!
   b_Vertex_py->GetEntry(entry);   //!
   b_Vertex_pz->GetEntry(entry);   //!
   b_Vertex_scMatchRatio->GetEntry(entry);   //!
   b_Vertex_sumCharge->GetEntry(entry);   //!
   b_Vertex_x->GetEntry(entry);   //!
   b_Vertex_y->GetEntry(entry);   //!
   b_Vertex_z->GetEntry(entry);   //!

   b_ECALRecHit_energy->GetEntry(entry);   //!
   b_ECALRecHit_ID->GetEntry(entry);   //!
   b_ECALRecHit_swCross->GetEntry(entry);   //!
   b_ECALRecHit_0TOF->GetEntry(entry);   //!
   b_ECALRecHit_pvTOF->GetEntry(entry);   //!
   b_ECALRecHit_time->GetEntry(entry);   //!
   b_ECALRecHit_timeError->GetEntry(entry);   //!
   b_ECALRecHit_amplitude->GetEntry(entry);   //!
   b_ECALRecHit_ampres->GetEntry(entry);   //!
   b_ECALRecHit_eta->GetEntry(entry);   //!
   b_ECALRecHit_hasGS1->GetEntry(entry);   //!
   b_ECALRecHit_hasGS6->GetEntry(entry);   //!
   b_ECALRecHit_isDead->GetEntry(entry);   //!
   b_ECALRecHit_isDiWrd->GetEntry(entry);   //!
   b_ECALRecHit_isOOT->GetEntry(entry);   //!
   b_ECALRecHit_isOther->GetEntry(entry);   //!
   b_ECALRecHit_isPoor->GetEntry(entry);   //!
   b_ECALRecHit_isRecov->GetEntry(entry);   //!
   b_ECALRecHit_isTimeValid->GetEntry(entry);   //!
   b_ECALRecHit_isWrd->GetEntry(entry);   //!
   b_ECALRecHit_phi->GetEntry(entry);   //!
   b_ECALRecHit_precentUsed->GetEntry(entry);   //!
   b_ECALRecHit_rhx->GetEntry(entry);   //!
   b_ECALRecHit_rhy->GetEntry(entry);   //!
   b_ECALRecHit_rhz->GetEntry(entry);   //!

   b_SuperCluster_covEtaEta->GetEntry(entry);   //!
   b_SuperCluster_covEtaPhi->GetEntry(entry);   //!
   b_SuperCluster_covPhiPhi->GetEntry(entry);   //!
   b_SuperCluster_energyRaw->GetEntry(entry);   //!
   b_SuperCluster_etaWidth->GetEntry(entry);   //!
   b_SuperCluster_excluded->GetEntry(entry);   //!
   b_SuperCluster_seedIsEB->GetEntry(entry);   //!
   b_SuperCluster_isScEtaEB->GetEntry(entry);   //!
   b_SuperCluster_isScEtaEE->GetEntry(entry);   //!
   b_SuperCluster_isOot->GetEntry(entry);   //!
   b_SuperCluster_original->GetEntry(entry);   //!
   b_SuperCluster_phiWidth->GetEntry(entry);   //!
   b_SuperCluster_salp->GetEntry(entry);   //!
   b_SuperCluster_smaj->GetEntry(entry);   //!
   b_SuperCluster_smin->GetEntry(entry);   //!
   b_SuperCluster_seediEtaOriX->GetEntry(entry);   //!
   b_SuperCluster_seediPhiOriY->GetEntry(entry);   //!
   b_SuperCluster_nBasicClusters->GetEntry(entry);   //!
   b_SuperCluster_ElectronIndx->GetEntry(entry);   //!
   b_SuperCluster_energy->GetEntry(entry);   //!
   b_SuperCluster_eta->GetEntry(entry);   //!
   b_SuperCluster_clcx->GetEntry(entry);   //!
   b_SuperCluster_clcy->GetEntry(entry);   //!
   b_SuperCluster_clcz->GetEntry(entry);   //!
   b_SuperCluster_nSuperCluster->GetEntry(entry);   //!
   b_SuperCluster_nRHXtals->GetEntry(entry);   //!
   b_SuperCluster_ObjectPdgId->GetEntry(entry);   //!
   b_SuperCluster_phi->GetEntry(entry);   //!
   b_SuperCluster_PhotonIndx->GetEntry(entry);   //!
   b_SuperCluster_rhFracs->GetEntry(entry);   //!
   b_SuperCluster_rhIds->GetEntry(entry);   //!
   b_SuperCluster_XtalSeedID->GetEntry(entry);   //!
   b_SuperCluster_nHFXtals->GetEntry(entry);   //!
   b_SuperCluster_SeedTime->GetEntry(entry);   //!
   b_SuperCluster_TrackIndx->GetEntry(entry);   //!
   b_SuperCluster_x_calo->GetEntry(entry);   //!
   b_SuperCluster_y_calo->GetEntry(entry);   //!
   b_SuperCluster_z_calo->GetEntry(entry);   //!

   b_Electron_DetaSCTV->GetEntry(entry);   //!
   b_Electron_DphiSCTV->GetEntry(entry);   //!
   b_Electron_energy->GetEntry(entry);   //!
   b_Electron_eta->GetEntry(entry);   //!
   if( doGenInfo ){
   b_Electron_genIdx->GetEntry(entry);   //!
   b_Electron_genSigWZId->GetEntry(entry);   //!
   b_Electron_genSigXMomId->GetEntry(entry);   //!
   }//<<>>if( doGenInfo )
   b_Electron_HOE->GetEntry(entry);   //!
   b_Electron_isLVeto->GetEntry(entry);   //!
   b_Electron_isLoose->GetEntry(entry);   //!
   b_Electron_isMedium->GetEntry(entry);   //!
   b_Electron_isPrompt->GetEntry(entry);   //!
   b_Electron_isVeto->GetEntry(entry);   //!
   b_Electron_phi->GetEntry(entry);   //!
   b_Electron_pt->GetEntry(entry);   //!
   b_Electron_px->GetEntry(entry);   //!
   b_Electron_py->GetEntry(entry);   //!
   b_Electron_pz->GetEntry(entry);   //!
   b_Electron_seedTOFTime->GetEntry(entry);   //!
   b_Electron_Sieie->GetEntry(entry);   //!
   b_Electron_trackz->GetEntry(entry);   //!
   b_Electron_Vdetaseed->GetEntry(entry);   //!
   b_Electron_Vdpscvc->GetEntry(entry);   //!
   b_Electron_Vhoec->GetEntry(entry);   //!
   b_Electron_Vmistrack->GetEntry(entry);   //!
   b_Electron_Vooemoop->GetEntry(entry);   //!
   b_Electron_VrelIsoWithEA->GetEntry(entry);   //!
   b_Electron_Vsieiec->GetEntry(entry);   //!
   if( doGenInfo ) b_Electron_hasGenMatch->GetEntry(entry);   //!
   if( doGenInfo ) b_Electron_isLWZP->GetEntry(entry);   //!
   b_Electron_nElectrons->GetEntry(entry);   //!
   if( doGenInfo ) b_Electron_nGenMatch->GetEntry(entry);   //!
   b_Electron_nPrompt->GetEntry(entry);   //!
   b_Electron_nSVMatched->GetEntry(entry);   //!
   b_Electron_nSelElectrons->GetEntry(entry);   //!
   b_Electron_passConversionVeto->GetEntry(entry);   //!
   b_Electron_passPrompt->GetEntry(entry);   //!
   b_Electron_scIndex->GetEntry(entry);   //!
   b_Electron_hasSVMatch->GetEntry(entry);   //!

   b_Flag_BadChargedCandidateFilter->GetEntry(entry);   //!
   b_Flag_BadPFMuonDzFilter->GetEntry(entry);   //!
   b_Flag_BadPFMuonFilter->GetEntry(entry);   //!
   //b_Flag_EcalDeadCellTriggerPrimitiveFilter->GetEntry(entry);   //!
   b_Flag_HBHENoiseFilter->GetEntry(entry);   //!
   b_Flag_HBHENoiseIsoFilter->GetEntry(entry);   //!
   b_Flag_ecalBadCalibFilter->GetEntry(entry);   //!
   b_Flag_eeBadScFilter->GetEntry(entry);   //!
   b_Flag_globalSuperTightHalo2016Filter->GetEntry(entry);   //!
   b_Flag_goodVertices->GetEntry(entry);   //!
   b_Flag_hfNoisyHitsFilter->GetEntry(entry);   //!

   b_Evt_luminosityBlock->GetEntry(entry);   //!
   b_Evt_run->GetEntry(entry);   //!
   b_Evt_event->GetEntry(entry);   //!

   b_Trigger_hltEG20EtFilterLooseHoverE->GetEntry(entry);   //!
   b_Trigger_hltEG20EtL1TripleEGFilter->GetEntry(entry);   //!
   b_Trigger_hltEG20HEFilterLooseHoverE->GetEntry(entry);   //!
   b_Trigger_hltEG20HEL1TripleEGFilter->GetEntry(entry);   //!
   b_Trigger_hltEG20L1SingleEGLowETOrEtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG20L1SingleEGLowETOrEtaREtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30EBHE10R9Id50b80eHEFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30EBL1SingleAndDoubleEGOrEtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30EBR9Id50b80eR9IdFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30EIso15HE30EcalIsoLastFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30EtFilterLooseHoverE->GetEntry(entry);   //!
   b_Trigger_hltEG30EtL1TripleEGFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30HE30HEFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30HEFilterLooseHoverE->GetEntry(entry);   //!
   b_Trigger_hltEG30HEL1TripleEGFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30L1IsoEGerJetC34drMin0p3EtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30L1SingleAndDoubleEGOrEtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30LHE12R9Id50b80eHEFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30LR9Id50b80eR9IdFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30PVHE10R9Id50b80eHEFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30PVR9Id50b80eR9IdFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30PVrealANDHE10R9Id50b80eHEFilter->GetEntry(entry);   //!
   b_Trigger_hltEG30PVrealANDR9Id50b80eR9IdFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60EtFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60HEFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLClusterShapeFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter->GetEntry(entry);   //!
   b_Trigger_hltEG60R9Id90CaloIdLIsoLR9IdFilter->GetEntry(entry);   //!
   b_Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter->GetEntry(entry);   //!
   b_Trigger_hltHT130Jet30->GetEntry(entry);   //!
   b_Trigger_hltHT175Jet10->GetEntry(entry);   //!
   b_Trigger_hltHT200Jet30->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleAndDoubleEG->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleAndDoubleEGNonIsoOr->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleAndDoubleEGor->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG15->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG18->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG24->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG26->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG34to45->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG34to50->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEG40to50->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS->GetEntry(entry);   //!
   b_Trigger_hltL1sSingleEGor->GetEntry(entry);   //!
   b_Trigger_hltL1sTripleEGOrDoubleEGOrSingleEG->GetEntry(entry);   //!
   b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter->GetEntry(entry);   //!
   b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt->GetEntry(entry);   //!
   b_Trigger_hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt->GetEntry(entry);   //!
   b_Trigger_hltPFHT180Jet30->GetEntry(entry);   //!
   b_Trigger_hltPFHT250Jet30->GetEntry(entry);   //!
   b_Trigger_hltPFHT350Jet15->GetEntry(entry);   //!
   b_Trigger_hltPFMET100->GetEntry(entry);   //!
   b_Trigger_hltPFMET110->GetEntry(entry);   //!
   b_Trigger_hltPFMET120->GetEntry(entry);   //!
   b_Trigger_hltPFMET130->GetEntry(entry);   //!
   b_Trigger_hltPFMET140->GetEntry(entry);   //!
   b_Trigger_hltPFMET200->GetEntry(entry);   //!
   b_Trigger_hltPFMET250->GetEntry(entry);   //!
   b_Trigger_hltPFMET300->GetEntry(entry);   //!
   b_Trigger_hltPFMET50->GetEntry(entry);   //!
   b_Trigger_hltPFMET70->GetEntry(entry);   //!
   b_Trigger_hltPFMET90->GetEntry(entry);   //!
   b_Trigger_hltPFMETNoMu100->GetEntry(entry);   //!
   b_Trigger_hltPFMETNoMu110->GetEntry(entry);   //!
   b_Trigger_hltPFMETNoMu120->GetEntry(entry);   //!
   b_Trigger_hltPFMETNoMu130->GetEntry(entry);   //!
   b_Trigger_hltPFMETNoMu140->GetEntry(entry);   //!
   b_Trigger_hltPFMETNoMu60->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID100->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID110->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID120->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID130->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID140->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID70->GetEntry(entry);   //!
   b_Trigger_hltPFMHTNoMuTightID90->GetEntry(entry);   //!
   b_Trigger_hltPFMHTTightID120->GetEntry(entry);   //!
   b_Trigger_hltPFMHTTightID130->GetEntry(entry);   //!
   b_Trigger_hltPFMHTTightID140->GetEntry(entry);   //!

   b_PV_npvs->GetEntry(entry);   //!
   b_PV_x->GetEntry(entry);   //!
   b_PV_y->GetEntry(entry);   //!
   b_PV_z->GetEntry(entry);   //!

   if( doGenInfo ){
   b_Gen_charge->GetEntry(entry);   //!
   b_Gen_energy->GetEntry(entry);   //!
   b_Gen_eta->GetEntry(entry);   //!
   b_Gen_GMomPdgId->GetEntry(entry);   //!
   b_Gen_isLWZQ->GetEntry(entry);   //!
   b_Gen_isLWZX->GetEntry(entry);   //!
   b_Gen_mass->GetEntry(entry);   //!
   b_Gen_matchDe->GetEntry(entry);   //!
   b_Gen_matchDr->GetEntry(entry);   //!
   b_Gen_momDisplacment->GetEntry(entry);   //!
   b_Gen_motherIdx->GetEntry(entry);   //!
   b_Gen_MomPdgId->GetEntry(entry);   //!
   b_Gen_nTotal->GetEntry(entry);   //!
   b_Gen_pdgId->GetEntry(entry);   //!
   b_Gen_phi->GetEntry(entry);   //!
   b_Gen_pt->GetEntry(entry);   //!
   b_Gen_px->GetEntry(entry);   //!
   b_Gen_py->GetEntry(entry);   //!
   b_Gen_pz->GetEntry(entry);   //!
   b_Gen_status->GetEntry(entry);   //!
   b_Gen_susEvtType->GetEntry(entry);   //!
   b_Gen_susId->GetEntry(entry);   //!
   b_Gen_vx->GetEntry(entry);   //!
   b_Gen_vy->GetEntry(entry);   //!
   b_Gen_vz->GetEntry(entry);   //!
   }//<<>>if( doGenInfo )

   b_Jet_area->GetEntry(entry);   //!
   b_Jet_chEmEF->GetEntry(entry);   //!
   b_Jet_chHEF->GetEntry(entry);   //!
   b_Jet_chHM->GetEntry(entry);   //!
   b_Jet_drRhIds->GetEntry(entry);   //!
   b_Jet_energy->GetEntry(entry);   //!
   b_Jet_eta->GetEntry(entry);   //!
   if( doGenInfo ){
   b_Jet_genJetLlpDr->GetEntry(entry);   //!
   b_Jet_genDptMatch->GetEntry(entry);   //!
   b_Jet_genDrMatch->GetEntry(entry);   //!
   b_Jet_genEnergy->GetEntry(entry);   //!
   b_Jet_genEta->GetEntry(entry);   //!
   b_Jet_genImpactAngle->GetEntry(entry);   //!
   b_Jet_genJetLlpId->GetEntry(entry);   //!
   b_Jet_genJetLlpRp->GetEntry(entry);   //!
   b_Jet_genPhi->GetEntry(entry);   //!
   b_Jet_genPt->GetEntry(entry);   //!
   b_Jet_genQrkLlpDr->GetEntry(entry);   //!
   b_Jet_genQrkLlpId->GetEntry(entry);   //!
   b_Jet_genQrkLlpRp->GetEntry(entry);   //!
   b_Jet_genTOF->GetEntry(entry);   //!
   b_Jet_genTime->GetEntry(entry);   //!
   b_Jet_genTimeLLP->GetEntry(entry);   //!
   }//<<>>if( doGenInfo )
   b_Jet_mass->GetEntry(entry);   //!
   b_Jet_muEF->GetEntry(entry);   //!
   b_Jet_neEmEF->GetEntry(entry);   //!
   b_Jet_neHEF->GetEntry(entry);   //!
   b_Jet_neHM->GetEntry(entry);   //!
   b_Jet_egIndxs->GetEntry(entry);   //!
   b_Jet_phi->GetEntry(entry);   //!
   b_Jet_pt->GetEntry(entry);   //!
   b_Jet_nConstituents->GetEntry(entry);   //!

   b_Muon_energy->GetEntry(entry);   //!
   b_Muon_eta->GetEntry(entry);   //!
   if( doGenInfo ){
   b_Muon_genIdx->GetEntry(entry);   //!
   b_Muon_genSigWZId->GetEntry(entry);   //!
   b_Muon_genSigXMomId->GetEntry(entry);   //!
   }//<<>>if( doGenInfo )
   b_Muon_isPromptMuon->GetEntry(entry);   //!
   b_Muon_phi->GetEntry(entry);   //!
   b_Muon_pt->GetEntry(entry);   //!
   b_Muon_px->GetEntry(entry);   //!
   b_Muon_py->GetEntry(entry);   //!
   b_Muon_pz->GetEntry(entry);   //!
   if( doGenInfo ) b_Muon_hasGenMatch->GetEntry(entry);   //!
   if( doGenInfo ) b_Muon_isLWZP->GetEntry(entry);   //!
   b_Muon_isLoose->GetEntry(entry);   //!
   b_Muon_isMedium->GetEntry(entry);   //!
   b_Muon_nGenMatch->GetEntry(entry);   //!
   b_Muon_nMuons->GetEntry(entry);   //!
   b_Muon_nPrompt->GetEntry(entry);   //!
   b_Muon_nSelMuons->GetEntry(entry);   //!
   b_Muon_nSVMatched->GetEntry(entry);   //!
   b_Muon_passPrompt->GetEntry(entry);   //!
   b_Muon_hasSVMatch->GetEntry(entry);   //!
   b_Muon_isTight->GetEntry(entry);   //!

   b_Met_CPt->GetEntry(entry);   //!
   b_Met_Cpx->GetEntry(entry);   //!
   b_Met_Cpy->GetEntry(entry);   //!
   b_Met_CsumEt->GetEntry(entry);   //!
   b_Met_eta->GetEntry(entry);   //!
   b_Met_phi->GetEntry(entry);   //!
   b_Met_pt->GetEntry(entry);   //!
   b_Met_px->GetEntry(entry);   //!
   b_Met_py->GetEntry(entry);   //!
   b_Met_sumEt->GetEntry(entry);   //!

   b_Photon_ecalRHSumEtConeDR04->GetEntry(entry);   //!
   b_Photon_energy->GetEntry(entry);   //!
   b_Photon_energyErr->GetEntry(entry);   //!
   b_Photon_eta->GetEntry(entry);   //!
   b_Photon_excluded->GetEntry(entry);   //!
   if( doGenInfo ){
   b_Photon_genIdx->GetEntry(entry);   //!
   b_Photon_genSigXMomId->GetEntry(entry);   //!
   }//<<>>if( doGenInfo )
   b_Photon_hadOverEM->GetEntry(entry);   //!
   b_Photon_hadTowOverEM->GetEntry(entry);   //!
   b_Photon_hcalTowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_Photon_hcalTowerSumEtConeDR04->GetEntry(entry);   //!
   b_Photon_isOot->GetEntry(entry);   //!
   b_Photon_nTrkHollowConeDR04->GetEntry(entry);   //!
   b_Photon_nTrkSolidConeDR04->GetEntry(entry);   //!
   b_Photon_phi->GetEntry(entry);   //!
   b_Photon_pt->GetEntry(entry);   //!
   b_Photon_px->GetEntry(entry);   //!
   b_Photon_py->GetEntry(entry);   //!
   b_Photon_pz->GetEntry(entry);   //!
   b_Photon_r9->GetEntry(entry);   //!
   b_Photon_s4->GetEntry(entry);   //!
   b_Photon_seedTOFTime->GetEntry(entry);   //!
   b_Photon_SigmaIEtaIEta->GetEntry(entry);   //!
   b_Photon_trkSumPtHollowConeDR03->GetEntry(entry);   //!
   b_Photon_trkSumPtHollowConeDR04->GetEntry(entry);   //!
   b_Photon_trkSumPtSolidConeDR04->GetEntry(entry);   //!
   b_Photon_ecalPFClusterIso->GetEntry(entry);   //!
   b_Photon_electronVeto->GetEntry(entry);   //!
   b_hasConversionTracks->GetEntry(entry);   //!
   b_Photon_pixelSeed->GetEntry(entry);   //!
   b_Photon_hcalPFClusterIso->GetEntry(entry);   //!
   b_Photon_Hoe_PUcorr->GetEntry(entry);   //!
   b_Photon_nPhoton->GetEntry(entry);   //!
   b_Photon_isSelect->GetEntry(entry);   //!
   b_Photon_pfChargedIso->GetEntry(entry);   //!
   b_Photon_pfChargedIsoPFPV->GetEntry(entry);   //!
   b_Photon_pfChargedIsoWorstVtx->GetEntry(entry);   //!
   b_Photon_pfPhoIso03->GetEntry(entry);   //!
   b_pfRelIso03_all_quadratic->GetEntry(entry);   //!
   b_pfRelIso03_chg_quadratic->GetEntry(entry);   //!
   b_Photon_scIndex->GetEntry(entry);   //!
   b_Photon_sieie->GetEntry(entry);   //!
   b_Photon_sieip->GetEntry(entry);   //!
   b_Photon_sipip->GetEntry(entry);   //!

   b_Track_absIso->GetEntry(entry);   //!
   b_Track_charge->GetEntry(entry);   //!
   b_Track_chi2->GetEntry(entry);   //!
   b_Track_dxy->GetEntry(entry);   //!
   b_Track_ecalDetIds->GetEntry(entry);   //!
   b_Track_eta->GetEntry(entry);   //!
   b_Track_etaError->GetEntry(entry);   //!
   b_Track_hcalDetIds->GetEntry(entry);   //!
   b_Track_ip2D->GetEntry(entry);   //!
   b_Track_ip3D->GetEntry(entry);   //!
   b_Track_isGenTrack->GetEntry(entry);   //!
   b_Track_lambda->GetEntry(entry);   //!
   b_Track_miniRelIso->GetEntry(entry);   //!
   b_Track_pdgId->GetEntry(entry);   //!
   b_Track_mrelisolep->GetEntry(entry);   //!
   b_Track_nIsoTracks->GetEntry(entry);   //!
   b_Track_nLostHits->GetEntry(entry);   //!
   b_Track_nValidHits->GetEntry(entry);   //!
   b_Track_ndof->GetEntry(entry);   //!
   b_Track_normalizedChi2->GetEntry(entry);   //!
   b_Track_p->GetEntry(entry);   //!
   b_Track_phi->GetEntry(entry);   //!
   b_Track_phiError->GetEntry(entry);   //!
   b_Track_pt->GetEntry(entry);   //!
   b_Track_ptError->GetEntry(entry);   //!
   b_Track_px->GetEntry(entry);   //!
   b_Track_py->GetEntry(entry);   //!
   b_Track_pz->GetEntry(entry);   //!
   b_Track_normalizedChi2->GetEntry(entry);   //!
   b_Track_qualityMask->GetEntry(entry);   //!
   b_Track_relIso->GetEntry(entry);   //!
   b_Track_relisolep->GetEntry(entry);   //!
   b_Track_scIndexs->GetEntry(entry);   //!
   b_Track_sip2D->GetEntry(entry);   //!
   b_Track_sip3D->GetEntry(entry);   //!
   b_Track_theta->GetEntry(entry);   //!
   b_Track_vx->GetEntry(entry);   //!
   b_Track_vy->GetEntry(entry);   //!
   b_Track_vz->GetEntry(entry);   //!

}//<<>>void root_base::getBranches(Long64_t entry)

/*
Bool_t kuntuple_v30::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kuntuple_v30::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kuntuple_v30::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/

#endif // #ifdef kuntuple_v30_cxx
