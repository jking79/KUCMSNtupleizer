//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 14 13:09:36 2024 by ROOT version 6.26/11
// from TTree llpgtree/KUCMSNtuple
// found on file: gmsb_AODSIM_KUCMSNtuplizer_Objectified_v20_oottrue.root
//////////////////////////////////////////////////////////

#ifndef ntuplev20_h
#define ntuplev20_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class ntuplev20 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *DisplacedElectron_charge;
   Float_t         DisplacedElectron_costGen;
   Float_t         DisplacedElectron_costSC;
   vector<float>   *DisplacedElectron_dRGen;
   vector<float>   *DisplacedElectron_dRSC;
   vector<float>   *DisplacedElectron_dxy;
   vector<float>   *DisplacedElectron_energy;
   vector<float>   *DisplacedElectron_eta;
   vector<float>   *DisplacedElectron_etaECAL;
   vector<int>     *DisplacedElectron_genSigIndex;
   vector<int>     *DisplacedElectron_indexGen;
   vector<int>     *DisplacedElectron_indexSC;
   vector<int>     *DisplacedElectron_indexTrack;
   vector<float>   *DisplacedElectron_ip2D;
   vector<float>   *DisplacedElectron_ip2Dbs;
   vector<float>   *DisplacedElectron_ip3D;
   vector<float>   *DisplacedElectron_ip3Dbs;
   vector<bool>    *DisplacedElectron_isBarrel;
   vector<bool>    *DisplacedElectron_isConversion;
   vector<bool>    *DisplacedElectron_isEndCap;
   vector<bool>    *DisplacedElectron_isGeneralTrack;
   vector<bool>    *DisplacedElectron_isGsfTrack;
   vector<bool>    *DisplacedElectron_isHeavy;
   vector<bool>    *DisplacedElectron_isLight;
   vector<bool>    *DisplacedElectron_isOOTSC;
   vector<bool>    *DisplacedElectron_isSC;
   vector<bool>    *DisplacedElectron_isSignal;
   vector<bool>    *DisplacedElectron_isTau;
   vector<bool>    *DisplacedElectron_isUnmatched;
   Int_t           DisplacedElectron_nGenSig;
   vector<int>     *DisplacedElectron_nHits;
   Int_t           DisplacedElectron_nOOTSCs;
   Int_t           DisplacedElectron_nSCs;
   Int_t           DisplacedElectron_nTotal;
   vector<float>   *DisplacedElectron_p;
   vector<int>     *DisplacedElectron_pdgId;
   vector<float>   *DisplacedElectron_pfIso0p4;
   vector<float>   *DisplacedElectron_phi;
   vector<float>   *DisplacedElectron_phiECAL;
   vector<float>   *DisplacedElectron_pt;
   vector<float>   *DisplacedElectron_scIso0p4;
   vector<float>   *DisplacedElectron_sip2D;
   vector<float>   *DisplacedElectron_sip2Dbs;
   vector<float>   *DisplacedElectron_sip3D;
   vector<float>   *DisplacedElectron_sip3Dbs;
   vector<float>   *DisplacedElectron_trackIso0p4;
   vector<float>   *DisplacedElectron_x;
   vector<float>   *DisplacedElectron_xECAL;
   vector<float>   *DisplacedElectron_y;
   vector<float>   *DisplacedElectron_yECAL;
   vector<float>   *DisplacedElectron_z;
   vector<float>   *DisplacedElectron_zECAL;
   vector<float>   *ECALRecHit_energy;
   vector<unsigned int> *ECALRecHit_ID;
   vector<float>   *ECALRecHit_swCross;
   vector<float>   *ECALRecHit_0TOF;
   vector<float>   *ECALRecHit_pvTOF;
   vector<float>   *ECALRecHit_time;
   vector<float>   *ECALRecHit_amplitude;
   vector<float>   *ECALRecHit_ampres;
   vector<float>   *ECALRecHit_eta;
   vector<bool>    *ECALRecHit_isDead;
   vector<bool>    *ECALRecHit_isDiWrd;
   vector<bool>    *ECALRecHit_isOOT;
   vector<bool>    *ECALRecHit_isOther;
   vector<bool>    *ECALRecHit_isPoor;
   vector<bool>    *ECALRecHit_isRecov;
   vector<bool>    *ECALRecHit_isWrd;
   vector<float>   *ECALRecHit_phi;
   vector<float>   *ECALRecHit_precentUsed;
   vector<float>   *ECALRecHit_rhx;
   vector<float>   *ECALRecHit_rhy;
   vector<float>   *ECALRecHit_rhz;
   vector<float>   *SuperCluster_covEtaEta;
   vector<float>   *SuperCluster_covEtaPhi;
   vector<float>   *SuperCluster_covPhiPhi;
   vector<float>   *SuperCluster_energyRaw;
   vector<float>   *SuperCluster_etaWidth;
   vector<bool>    *SuperCluster_excluded;
   vector<bool>    *SuperCluster_seedIsEB;
   vector<bool>    *SuperCluster_isScEtaEB;
   vector<bool>    *SuperCluster_isScEtaEE;
   vector<bool>    *SuperCluster_isOot;
   vector<float>   *SuperCluster_phiWidth;
   vector<float>   *SuperCluster_salp;
   vector<float>   *SuperCluster_smaj;
   vector<float>   *SuperCluster_smin;
   vector<int>     *SuperCluster_seediEtaOriX;
   vector<int>     *SuperCluster_seediPhiOriY;
   vector<unsigned int> *SuperCluster_nBasicClusters;
   vector<int>     *SuperCluster_ElectronIndx;
   vector<float>   *SuperCluster_energy;
   vector<float>   *SuperCluster_eta;
   vector<float>   *SuperCluster_clcx;
   vector<float>   *SuperCluster_clcy;
   vector<float>   *SuperCluster_clcz;
   vector<float>   *OSuperCluster_diffEnergy;
   vector<float>   *OSuperCluster_XdiffEnergy;
   vector<float>   *OSuperCluster_dR;
   vector<float>   *OSuperCluster_XdR;
   vector<float>   *OSuperCluster_perDiffEnergy;
   vector<float>   *OSuperCluster_XPerDiffEnergy;
   vector<unsigned int> *OSuperCluster_otherMatchSeedID;
   vector<unsigned int> *OSuperCluster_XotherMatchSeedID;
   Int_t           OSuperCluster_nOther;
   Int_t           OSuperCluster_nOtherEx;
   Int_t           OSuperCluster_nOtherIn;
   vector<unsigned int> *OSuperCluster_otherSeedID;
   vector<int>     *OSuperCluster_nXtalOverlap;
   vector<int>     *OSuperCluster_XnXtalOverlap;
   Int_t           SuperCluster_nSuperCluster;
   vector<unsigned int> *SuperCluster_nRHXtals;
   vector<int>     *SuperCluster_ObjectPdgId;
   vector<float>   *OSuperCluster_precentXtalOverlap;
   vector<float>   *OSuperCluster_XprecentXtalOverlap;
   vector<float>   *SuperCluster_phi;
   vector<int>     *SuperCluster_PhotonIndx;
   vector<vector<float> > *SuperCluster_rhFracs;
   vector<vector<unsigned int> > *SuperCluster_rhIds;
   vector<unsigned int> *SuperCluster_XtalSeedID;
   vector<unsigned int> *SuperCluster_nHFXtals;
   vector<float>   *SuperCluster_x_calo;
   vector<float>   *SuperCluster_y_calo;
   vector<float>   *SuperCluster_z_calo;
   vector<unsigned int> *ECALTrackDetID_detId;
   vector<bool>    *ECALTrackDetID_isECAL;
   vector<bool>    *ECALTrackDetID_isHCAL;
   vector<unsigned int> *ECALTrackDetID_trackIndex;
   vector<int>     *ECALTrack_charge;
   vector<float>   *ECALTrack_chi2;
   vector<float>   *ECALTrack_dxy;
   vector<float>   *ECALTrack_eta;
   vector<float>   *ECALTrack_etaError;
   vector<float>   *ECALTrack_genDeltaR;
   vector<int>     *ECALTrack_genIndex;
   vector<bool>    *ECALTrack_isGeneral;
   vector<bool>    *ECALTrack_isGsf;
   vector<float>   *ECALTrack_lambda;
   Int_t           ECALTrack_nGeneral;
   Int_t           ECALTrack_nGsf;
   vector<int>     *ECALTrack_nLostHits;
   Int_t           ECALTrack_nTracks;
   vector<int>     *ECALTrack_nValidHits;
   vector<float>   *ECALTrack_ndof;
   vector<float>   *ECALTrack_normalizedChi2;
   vector<float>   *ECALTrack_p;
   vector<int>     *ECALTrack_pdgId;
   vector<float>   *ECALTrack_phi;
   vector<float>   *ECALTrack_phiError;
   vector<float>   *ECALTrack_pt;
   vector<float>   *ECALTrack_ptError;
   vector<float>   *ECALTrack_px;
   vector<float>   *ECALTrack_py;
   vector<float>   *ECALTrack_pz;
   vector<float>   *ECALTrack_qoverp;
   vector<int>     *ECALTrack_qualityMask;
   vector<float>   *ECALTrack_theta;
   vector<float>   *ECALTrack_vx;
   vector<float>   *ECALTrack_vy;
   vector<float>   *ECALTrack_vz;
   vector<float>   *Electron_energy;
   vector<float>   *Electron_eta;
   vector<int>     *Electron_genIdx;
   vector<int>     *Electron_genSigMomId;
   vector<int>     *Electron_genSigWZId;
   vector<float>   *Electron_phi;
   vector<float>   *Electron_pt;
   vector<float>   *Electron_px;
   vector<float>   *Electron_py;
   vector<float>   *Electron_pz;
   vector<float>   *Electron_seedTOFTime;
   vector<float>   *Electron_trackz;
   Int_t           Electron_nElectrons;
   vector<int>     *Electron_scIndex;
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
   UInt_t          PV_npvs;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   vector<int>     *Gen_charge;
   vector<float>   *Gen_energy;
   vector<float>   *Gen_eta;
   vector<float>   *Gen_mass;
   UInt_t          Gen_nTotal;
   vector<unsigned int> *Gen_pdgId;
   vector<float>   *Gen_phi;
   vector<float>   *Gen_pt;
   vector<float>   *Gen_px;
   vector<float>   *Gen_py;
   vector<float>   *Gen_pz;
   vector<bool>    *Gen_status;
   vector<int>     *Gen_susEvtType;
   vector<int>     *Gen_susId;
   vector<float>   *Gen_vx;
   vector<float>   *Gen_vy;
   vector<float>   *Gen_vz;
   Float_t         Evt_genWgt;
   vector<float>   *Jet_area;
   vector<float>   *Jet_chEmEF;
   vector<float>   *Jet_chHEF;
   vector<float>   *Jet_chHM;
   vector<vector<unsigned int> > *Jet_drRhIds;
   vector<float>   *Jet_energy;
   vector<float>   *Jet_eta;
   vector<float>   *Jet_genDptMatch;
   vector<float>   *Jet_genDrMatch;
   vector<float>   *Jet_genEnergy;
   vector<float>   *Jet_genEta;
   vector<float>   *Jet_genImpactAngle;
   vector<float>   *Jet_genLlpDp;
   vector<float>   *Jet_genLlpDr;
   vector<float>   *Jet_genLlpId;
   vector<float>   *Jet_genPhi;
   vector<float>   *Jet_genPt;
   vector<float>   *Jet_genTOF;
   vector<float>   *Jet_genTime;
   vector<float>   *Jet_genTimeLLP;
   vector<float>   *Jet_mass;
   vector<float>   *Jet_muEF;
   vector<float>   *Jet_neEmEF;
   vector<float>   *Jet_neHEF;
   vector<float>   *Jet_neHM;
   vector<vector<unsigned int> > *Jet_egIndxs;
   vector<float>   *Jet_phi;
   vector<float>   *Jet_pt;
   vector<int>     *Jet_nConstituents;
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
   vector<float>   *Photon_ecalRHSumEtConeDR04;
   vector<float>   *Photon_energy;
   vector<float>   *Photon_energyErr;
   vector<float>   *Photon_eta;
   vector<bool>    *Photon_excluded;
   vector<int>     *Photon_genIdx;
   vector<int>     *Photon_genSigMomId;
   vector<float>   *Photon_hadOverEM;
   vector<float>   *Photon_hadTowOverEM;
   vector<float>   *Photon_hcalTowerSumEtBcConeDR04;
   vector<float>   *Photon_hcalTowerSumEtConeDR04;
   vector<bool>    *Photon_isOot;
   vector<float>   *Photon_nTrkHollowConeDR04;
   vector<float>   *Photon_nTrkSolidConeDR04;
   vector<float>   *Photon_phi;
   vector<float>   *Photon_pt;
   vector<float>   *Photon_px;
   vector<float>   *Photon_py;
   vector<float>   *Photon_pz;
   vector<float>   *Photon_r9;
   vector<float>   *Photon_s4;
   vector<float>   *Photon_seedTOFTime;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_trkSumPtHollowConeDR03;
   vector<float>   *Photon_trkSumPtHollowConeDR04;
   vector<float>   *Photon_trkSumPtSolidConeDR04;
   vector<float>   *Photon_ecalPFClusterIso;
   vector<bool>    *Photon_electronVeto;
   vector<vector<unsigned int> > *Photon_gloResRhId;
   vector<bool>    *hasConversionTracks;
   vector<bool>    *Photon_pixelSeed;
   vector<float>   *Photon_hcalPFClusterIso;
   vector<float>   *Photon_Hoe_PUcorr;
   vector<vector<unsigned int> > *Photon_locResRhId;
   Int_t           Photon_nPhoton;
   vector<float>   *Photon_pfChargedIso;
   vector<float>   *Photon_pfChargedIsoPFPV;
   vector<float>   *Photon_pfChargedIsoWorstVtx;
   vector<float>   *Photon_pfPhoIso03;
   vector<float>   *pfRelIso03_all_quadratic;
   vector<float>   *pfRelIso03_chg_quadratic;
   vector<int>     *Photon_scIndex;
   vector<float>   *Photon_sieie;
   vector<float>   *Photon_sieip;
   vector<float>   *Photon_sipip;

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
   TBranch        *b_OSuperCluster_diffEnergy;   //!
   TBranch        *b_OSuperCluster_XdiffEnergy;   //!
   TBranch        *b_OSuperCluster_dR;   //!
   TBranch        *b_OSuperCluster_XdR;   //!
   TBranch        *b_OSuperCluster_perDiffEnergy;   //!
   TBranch        *b_OSuperCluster_XPerDiffEnergy;   //!
   TBranch        *b_OSuperCluster_otherMatchSeedID;   //!
   TBranch        *b_OSuperCluster_XotherMatchSeedID;   //!
   TBranch        *b_OSuperCluster_nOther;   //!
   TBranch        *b_OSuperCluster_nOtherEx;   //!
   TBranch        *b_OSuperCluster_nOtherIn;   //!
   TBranch        *b_OSuperCluster_otherSeedID;   //!
   TBranch        *b_OSuperCluster_nXtalOverlap;   //!
   TBranch        *b_OSuperCluster_XnXtalOverlap;   //!
   TBranch        *b_SuperCluster_nSuperCluster;   //!
   TBranch        *b_SuperCluster_nRHXtals;   //!
   TBranch        *b_SuperCluster_ObjectPdgId;   //!
   TBranch        *b_OSuperCluster_precentXtalOverlap;   //!
   TBranch        *b_OSuperCluster_XprecentXtalOverlap;   //!
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
   TBranch        *b_ECALTrack_isGeneral;   //!
   TBranch        *b_ECALTrack_isGsf;   //!
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
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_genIdx;   //!
   TBranch        *b_Electron_genSigMomId;   //!
   TBranch        *b_Electron_genSigWZId;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_px;   //!
   TBranch        *b_Electron_py;   //!
   TBranch        *b_Electron_pz;   //!
   TBranch        *b_Electron_seedTOFTime;   //!
   TBranch        *b_Electron_trackz;   //!
   TBranch        *b_Electron_nElectrons;   //!
   TBranch        *b_Electron_scIndex;   //!
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
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_Gen_charge;   //!
   TBranch        *b_Gen_energy;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_mass;   //!
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
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_chHM;   //!
   TBranch        *b_Jet_drRhIds;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_genDptMatch;   //!
   TBranch        *b_Jet_genDrMatch;   //!
   TBranch        *b_Jet_genEnergy;   //!
   TBranch        *b_Jet_genEta;   //!
   TBranch        *b_Jet_genImpactAngle;   //!
   TBranch        *b_Jet_genLlpDp;   //!
   TBranch        *b_Jet_genLlpDr;   //!
   TBranch        *b_Jet_genLlpId;   //!
   TBranch        *b_Jet_genPhi;   //!
   TBranch        *b_Jet_genPt;   //!
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
   TBranch        *b_Photon_genSigMomId;   //!
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

   ntuplev20(TTree *tree=0);
   virtual ~ntuplev20();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntuplev20_cxx
ntuplev20::ntuplev20(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v20_oottrue.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v20_oottrue.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v20_oottrue.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

ntuplev20::~ntuplev20()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntuplev20::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntuplev20::LoadTree(Long64_t entry)
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

void ntuplev20::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

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
   OSuperCluster_diffEnergy = 0;
   OSuperCluster_XdiffEnergy = 0;
   OSuperCluster_dR = 0;
   OSuperCluster_XdR = 0;
   OSuperCluster_perDiffEnergy = 0;
   OSuperCluster_XPerDiffEnergy = 0;
   OSuperCluster_otherMatchSeedID = 0;
   OSuperCluster_XotherMatchSeedID = 0;
   OSuperCluster_otherSeedID = 0;
   OSuperCluster_nXtalOverlap = 0;
   OSuperCluster_XnXtalOverlap = 0;
   SuperCluster_nRHXtals = 0;
   SuperCluster_ObjectPdgId = 0;
   OSuperCluster_precentXtalOverlap = 0;
   OSuperCluster_XprecentXtalOverlap = 0;
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
   ECALTrack_isGeneral = 0;
   ECALTrack_isGsf = 0;
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
   Electron_energy = 0;
   Electron_eta = 0;
   Electron_genIdx = 0;
   Electron_genSigMomId = 0;
   Electron_genSigWZId = 0;
   Electron_phi = 0;
   Electron_pt = 0;
   Electron_px = 0;
   Electron_py = 0;
   Electron_pz = 0;
   Electron_seedTOFTime = 0;
   Electron_trackz = 0;
   Electron_scIndex = 0;
   Gen_charge = 0;
   Gen_energy = 0;
   Gen_eta = 0;
   Gen_mass = 0;
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
   Jet_genDptMatch = 0;
   Jet_genDrMatch = 0;
   Jet_genEnergy = 0;
   Jet_genEta = 0;
   Jet_genImpactAngle = 0;
   Jet_genLlpDp = 0;
   Jet_genLlpDr = 0;
   Jet_genLlpId = 0;
   Jet_genPhi = 0;
   Jet_genPt = 0;
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
   Photon_genSigMomId = 0;
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

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
   fChain->SetBranchAddress("DisplacedElectron_nHits", &DisplacedElectron_nHits, &b_DisplacedElectron_nHits);
   fChain->SetBranchAddress("DisplacedElectron_nOOTSCs", &DisplacedElectron_nOOTSCs, &b_DisplacedElectron_nOOTSCs);
   fChain->SetBranchAddress("DisplacedElectron_nSCs", &DisplacedElectron_nSCs, &b_DisplacedElectron_nSCs);
   fChain->SetBranchAddress("DisplacedElectron_nTotal", &DisplacedElectron_nTotal, &b_DisplacedElectron_nTotal);
   fChain->SetBranchAddress("DisplacedElectron_p", &DisplacedElectron_p, &b_DisplacedElectron_p);
   fChain->SetBranchAddress("DisplacedElectron_pdgId", &DisplacedElectron_pdgId, &b_DisplacedElectron_pdgId);
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
   fChain->SetBranchAddress("ECALRecHit_energy", &ECALRecHit_energy, &b_ECALRecHit_energy);
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
   fChain->SetBranchAddress("OSuperCluster_diffEnergy", &OSuperCluster_diffEnergy, &b_OSuperCluster_diffEnergy);
   fChain->SetBranchAddress("OSuperCluster_XdiffEnergy", &OSuperCluster_XdiffEnergy, &b_OSuperCluster_XdiffEnergy);
   fChain->SetBranchAddress("OSuperCluster_dR", &OSuperCluster_dR, &b_OSuperCluster_dR);
   fChain->SetBranchAddress("OSuperCluster_XdR", &OSuperCluster_XdR, &b_OSuperCluster_XdR);
   fChain->SetBranchAddress("OSuperCluster_perDiffEnergy", &OSuperCluster_perDiffEnergy, &b_OSuperCluster_perDiffEnergy);
   fChain->SetBranchAddress("OSuperCluster_XPerDiffEnergy", &OSuperCluster_XPerDiffEnergy, &b_OSuperCluster_XPerDiffEnergy);
   fChain->SetBranchAddress("OSuperCluster_otherMatchSeedID", &OSuperCluster_otherMatchSeedID, &b_OSuperCluster_otherMatchSeedID);
   fChain->SetBranchAddress("OSuperCluster_XotherMatchSeedID", &OSuperCluster_XotherMatchSeedID, &b_OSuperCluster_XotherMatchSeedID);
   fChain->SetBranchAddress("OSuperCluster_nOther", &OSuperCluster_nOther, &b_OSuperCluster_nOther);
   fChain->SetBranchAddress("OSuperCluster_nOtherEx", &OSuperCluster_nOtherEx, &b_OSuperCluster_nOtherEx);
   fChain->SetBranchAddress("OSuperCluster_nOtherIn", &OSuperCluster_nOtherIn, &b_OSuperCluster_nOtherIn);
   fChain->SetBranchAddress("OSuperCluster_otherSeedID", &OSuperCluster_otherSeedID, &b_OSuperCluster_otherSeedID);
   fChain->SetBranchAddress("OSuperCluster_nXtalOverlap", &OSuperCluster_nXtalOverlap, &b_OSuperCluster_nXtalOverlap);
   fChain->SetBranchAddress("OSuperCluster_XnXtalOverlap", &OSuperCluster_XnXtalOverlap, &b_OSuperCluster_XnXtalOverlap);
   fChain->SetBranchAddress("SuperCluster_nSuperCluster", &SuperCluster_nSuperCluster, &b_SuperCluster_nSuperCluster);
   fChain->SetBranchAddress("SuperCluster_nRHXtals", &SuperCluster_nRHXtals, &b_SuperCluster_nRHXtals);
   fChain->SetBranchAddress("SuperCluster_ObjectPdgId", &SuperCluster_ObjectPdgId, &b_SuperCluster_ObjectPdgId);
   fChain->SetBranchAddress("OSuperCluster_precentXtalOverlap", &OSuperCluster_precentXtalOverlap, &b_OSuperCluster_precentXtalOverlap);
   fChain->SetBranchAddress("OSuperCluster_XprecentXtalOverlap", &OSuperCluster_XprecentXtalOverlap, &b_OSuperCluster_XprecentXtalOverlap);
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
   fChain->SetBranchAddress("ECALTrack_genDeltaR", &ECALTrack_genDeltaR, &b_ECALTrack_genDeltaR);
   fChain->SetBranchAddress("ECALTrack_genIndex", &ECALTrack_genIndex, &b_ECALTrack_genIndex);
   fChain->SetBranchAddress("ECALTrack_isGeneral", &ECALTrack_isGeneral, &b_ECALTrack_isGeneral);
   fChain->SetBranchAddress("ECALTrack_isGsf", &ECALTrack_isGsf, &b_ECALTrack_isGsf);
   fChain->SetBranchAddress("ECALTrack_lambda", &ECALTrack_lambda, &b_ECALTrack_lambda);
   fChain->SetBranchAddress("ECALTrack_nGeneral", &ECALTrack_nGeneral, &b_ECALTrack_nGeneral);
   fChain->SetBranchAddress("ECALTrack_nGsf", &ECALTrack_nGsf, &b_ECALTrack_nGsf);
   fChain->SetBranchAddress("ECALTrack_nLostHits", &ECALTrack_nLostHits, &b_ECALTrack_nLostHits);
   fChain->SetBranchAddress("ECALTrack_nTracks", &ECALTrack_nTracks, &b_ECALTrack_nTracks);
   fChain->SetBranchAddress("ECALTrack_nValidHits", &ECALTrack_nValidHits, &b_ECALTrack_nValidHits);
   fChain->SetBranchAddress("ECALTrack_ndof", &ECALTrack_ndof, &b_ECALTrack_ndof);
   fChain->SetBranchAddress("ECALTrack_normalizedChi2", &ECALTrack_normalizedChi2, &b_ECALTrack_normalizedChi2);
   fChain->SetBranchAddress("ECALTrack_p", &ECALTrack_p, &b_ECALTrack_p);
   fChain->SetBranchAddress("ECALTrack_pdgId", &ECALTrack_pdgId, &b_ECALTrack_pdgId);
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
   fChain->SetBranchAddress("Electron_energy", &Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_genIdx", &Electron_genIdx, &b_Electron_genIdx);
   fChain->SetBranchAddress("Electron_genSigMomId", &Electron_genSigMomId, &b_Electron_genSigMomId);
   fChain->SetBranchAddress("Electron_genSigWZId", &Electron_genSigWZId, &b_Electron_genSigWZId);
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_px", &Electron_px, &b_Electron_px);
   fChain->SetBranchAddress("Electron_py", &Electron_py, &b_Electron_py);
   fChain->SetBranchAddress("Electron_pz", &Electron_pz, &b_Electron_pz);
   fChain->SetBranchAddress("Electron_seedTOFTime", &Electron_seedTOFTime, &b_Electron_seedTOFTime);
   fChain->SetBranchAddress("Electron_trackz", &Electron_trackz, &b_Electron_trackz);
   fChain->SetBranchAddress("Electron_nElectrons", &Electron_nElectrons, &b_Electron_nElectrons);
   fChain->SetBranchAddress("Electron_scIndex", &Electron_scIndex, &b_Electron_scIndex);
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
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("Gen_charge", &Gen_charge, &b_Gen_charge);
   fChain->SetBranchAddress("Gen_energy", &Gen_energy, &b_Gen_energy);
   fChain->SetBranchAddress("Gen_eta", &Gen_eta, &b_Gen_eta);
   fChain->SetBranchAddress("Gen_mass", &Gen_mass, &b_Gen_mass);
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
   fChain->SetBranchAddress("Jet_area", &Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", &Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_chHM", &Jet_chHM, &b_Jet_chHM);
   fChain->SetBranchAddress("Jet_drRhIds", &Jet_drRhIds, &b_Jet_drRhIds);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_genDptMatch", &Jet_genDptMatch, &b_Jet_genDptMatch);
   fChain->SetBranchAddress("Jet_genDrMatch", &Jet_genDrMatch, &b_Jet_genDrMatch);
   fChain->SetBranchAddress("Jet_genEnergy", &Jet_genEnergy, &b_Jet_genEnergy);
   fChain->SetBranchAddress("Jet_genEta", &Jet_genEta, &b_Jet_genEta);
   fChain->SetBranchAddress("Jet_genImpactAngle", &Jet_genImpactAngle, &b_Jet_genImpactAngle);
   fChain->SetBranchAddress("Jet_genLlpDp", &Jet_genLlpDp, &b_Jet_genLlpDp);
   fChain->SetBranchAddress("Jet_genLlpDr", &Jet_genLlpDr, &b_Jet_genLlpDr);
   fChain->SetBranchAddress("Jet_genLlpId", &Jet_genLlpId, &b_Jet_genLlpId);
   fChain->SetBranchAddress("Jet_genPhi", &Jet_genPhi, &b_Jet_genPhi);
   fChain->SetBranchAddress("Jet_genPt", &Jet_genPt, &b_Jet_genPt);
   fChain->SetBranchAddress("Jet_genTOF", &Jet_genTOF, &b_Jet_genTOF);
   fChain->SetBranchAddress("Jet_genTime", &Jet_genTime, &b_Jet_genTime);
   fChain->SetBranchAddress("Jet_genTimeLLP", &Jet_genTimeLLP, &b_Jet_genTimeLLP);
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
   fChain->SetBranchAddress("Photon_genIdx", &Photon_genIdx, &b_Photon_genIdx);
   fChain->SetBranchAddress("Photon_genSigMomId", &Photon_genSigMomId, &b_Photon_genSigMomId);
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
   Notify();
}

Bool_t ntuplev20::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntuplev20::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntuplev20::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntuplev20_cxx
