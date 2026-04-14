//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  7 17:10:33 2026 by ROOT version 6.26/11
// from TTree kuSkimTree/output root file for kUCMSSkimmer
// found on file: SMS_v34p1_FULLMINI_SVHPM100_mN2-1300_mGl-2300_mN2-1300_mN1-1000_ct0p1_rjrskim.root
//////////////////////////////////////////////////////////

#ifndef _h
#define _h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class  {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_MetFilters;
   Bool_t          Flag_eLHemVeto;
   Bool_t          Flag_eMHemVeto;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_hemRun;
   Bool_t          Flag_hemVeto;
   Bool_t          Flag_hfNoisyHitsFilter;
   Bool_t          Flag_jLHemVeto;
   Bool_t          Flag_jMHemVeto;
   Bool_t          Flag_mLHemVeto;
   Bool_t          Flag_mMHemVeto;
   Bool_t          Flag_pLHemVeto;
   Bool_t          Flag_pMHemVeto;
   Float_t         PVx;
   Float_t         PVy;
   Float_t         PVz;
   Bool_t          Trigger_PFMET120_PFMHT120_IDTight;
   Bool_t          Trigger_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t          Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   string          *DataSetKey;
   Float_t         evtFillWgt;
   Float_t         evtGenWgt;
   Float_t         evtXSection;
   Float_t         selCMet;
   Float_t         selCMetPx;
   Float_t         selCMetPy;
   Float_t         selMet;
   Float_t         selMetPx;
   Float_t         selMetPy;
   Int_t           Trk_nIsoTrack;
   Bool_t          Evt_isGG;
   Bool_t          Evt_isGZ;
   Bool_t          Evt_isZZ;
   Int_t           Evt_nXs;
   Float_t         Xa_Displacment;
   Float_t         Xa_beta;
   Float_t         Xa_ctau;
   Float_t         Xa_energy;
   Float_t         Xa_mass;
   Float_t         Xa_p;
   UInt_t          Xa_pdgId;
   Float_t         Xa_phi;
   Float_t         Xa_pt;
   Float_t         Xa_vx;
   Float_t         Xa_vy;
   Float_t         Xa_vz;
   Float_t         Xb_Displacment;
   Float_t         Xb_beta;
   Float_t         Xb_ctau;
   Float_t         Xb_energy;
   Float_t         Xb_mass;
   Float_t         Xb_p;
   UInt_t          Xb_pdgId;
   Float_t         Xb_phi;
   Float_t         Xb_pt;
   Float_t         Xb_vx;
   Float_t         Xb_vy;
   Float_t         Xb_vz;
   vector<int>     *genCharge;
   vector<float>   *genGrvtinoMass;
   vector<float>   *genLLPMass;
   vector<float>   *genLSPMass;
   vector<float>   *genMass;
   vector<float>   *genMomCTau;
   vector<float>   *genPartEnergy;
   vector<float>   *genPartEta;
   vector<int>     *genPartMomIdx;
   vector<unsigned int> *genPartPdgId;
   vector<float>   *genPartPhi;
   vector<float>   *genPartPt;
   vector<int>     *genPartSusId;
   vector<float>   *genPx;
   vector<float>   *genPy;
   vector<float>   *genPz;
   vector<float>   *genSGMass;
   vector<float>   *genSQMass;
   Bool_t          genSTFlagQP;
   Bool_t          genSTFlagQQP;
   Bool_t          genSTFlagQQZ;
   Bool_t          genSTFlagQQZL;
   Bool_t          genSTFlagQQZLL;
   Int_t           genSigType;
   vector<int>     *genStatus;
   vector<float>   *genVx;
   vector<float>   *genVy;
   vector<float>   *genVz;
   vector<float>   *genXMomCTau;
   vector<bool>    *barrel_photon_baseline;
   vector<float>   *barrel_photon_energy;
   vector<float>   *barrel_photon_erhsecdr4;
   vector<float>   *barrel_photon_eta;
   vector<float>   *barrel_photon_htoem;
   vector<float>   *barrel_photon_isoANNScore;
   vector<float>   *barrel_photon_phi;
   vector<float>   *barrel_photon_pt;
   vector<float>   *barrel_photon_tspscdr4;
   vector<float>   *baseLinePhoton_ClstrRn;
   vector<float>   *baseLinePhoton_CorEnergy;
   vector<float>   *baseLinePhoton_CorPt;
   vector<float>   *baseLinePhoton_CovEtaEta;
   vector<float>   *baseLinePhoton_CovEtaPhi;
   vector<float>   *baseLinePhoton_CovPhiPhi;
   vector<bool>    *baseLinePhoton_DiJetsCR;
   vector<float>   *baseLinePhoton_EcalPFClusterIso;
   vector<float>   *baseLinePhoton_EcalRHSumEtConeDR04;
   UInt_t          baseLinePhoton_EleVeto;
   vector<float>   *baseLinePhoton_Energy;
   vector<float>   *baseLinePhoton_Eta;
   vector<float>   *baseLinePhoton_EtaWidth;
   vector<bool>    *baseLinePhoton_GJetsCR;
   vector<int>     *baseLinePhoton_GenIdx;
   vector<float>   *baseLinePhoton_GenLabTime;
   vector<float>   *baseLinePhoton_GenLabTimeSig;
   vector<float>   *baseLinePhoton_GenPt;
   vector<float>   *baseLinePhoton_GenSigMomEnergy;
   vector<float>   *baseLinePhoton_GenSigMomEta;
   vector<float>   *baseLinePhoton_GenSigMomMass;
   vector<float>   *baseLinePhoton_GenSigMomPhi;
   vector<float>   *baseLinePhoton_GenSigMomPt;
   vector<float>   *baseLinePhoton_GenSigMomPx;
   vector<float>   *baseLinePhoton_GenSigMomPy;
   vector<float>   *baseLinePhoton_GenSigMomPz;
   vector<float>   *baseLinePhoton_GenSigMomVx;
   vector<float>   *baseLinePhoton_GenSigMomVy;
   vector<float>   *baseLinePhoton_GenSigMomVz;
   vector<float>   *baseLinePhoton_GenTime;
   vector<float>   *baseLinePhoton_GenTimeSig;
   vector<float>   *baseLinePhoton_HadOverEM;
   vector<float>   *baseLinePhoton_HadTowOverEM;
   vector<bool>    *baseLinePhoton_HasConversionTracks;
   vector<float>   *baseLinePhoton_HcalPFClusterIso;
   vector<float>   *baseLinePhoton_HcalTowerSumEtBcConeDR04;
   vector<float>   *baseLinePhoton_HcalTowerSumEtConeDR04;
   vector<float>   *baseLinePhoton_Hoe_PUcorr;
   vector<float>   *baseLinePhoton_LSCross;
   vector<float>   *baseLinePhoton_LTRes;
   vector<float>   *baseLinePhoton_LTime;
   vector<float>   *baseLinePhoton_LTimeSig;
   vector<bool>    *baseLinePhoton_LWeird;
   vector<float>   *baseLinePhoton_MBetaEql;
   vector<float>   *baseLinePhoton_MBetaPmt;
   vector<float>   *baseLinePhoton_NTrkHollowConeDR04;
   vector<float>   *baseLinePhoton_NTrkSolidConeDR04;
   vector<unsigned int> *baseLinePhoton_Nrh;
   vector<bool>    *baseLinePhoton_OOT;
   vector<float>   *baseLinePhoton_PVtoRH;
   vector<float>   *baseLinePhoton_PfChargedIso;
   vector<float>   *baseLinePhoton_PfChargedIsoPFPV;
   vector<float>   *baseLinePhoton_PfChargedIsoWorstVtx;
   vector<float>   *baseLinePhoton_PfPhoIso03;
   vector<float>   *baseLinePhoton_PfRelIso03_all_quadratic;
   vector<float>   *baseLinePhoton_PfRelIso03_chg_quadratic;
   vector<float>   *baseLinePhoton_Phi;
   vector<float>   *baseLinePhoton_PhiWidth;
   vector<float>   *baseLinePhoton_PhoIsoDr;
   vector<bool>    *baseLinePhoton_PixelSeed;
   vector<float>   *baseLinePhoton_Pt;
   vector<int>     *baseLinePhoton_Quality;
   vector<float>   *baseLinePhoton_R9;
   vector<float>   *baseLinePhoton_RawTime;
   vector<float>   *baseLinePhoton_S4;
   vector<float>   *baseLinePhoton_SAlp;
   vector<float>   *baseLinePhoton_SCx;
   vector<float>   *baseLinePhoton_SCy;
   vector<float>   *baseLinePhoton_SCz;
   vector<float>   *baseLinePhoton_SMaj;
   vector<float>   *baseLinePhoton_SMin;
   vector<float>   *baseLinePhoton_SSCross;
   vector<float>   *baseLinePhoton_STRes;
   vector<float>   *baseLinePhoton_STime;
   vector<float>   *baseLinePhoton_STimeSig;
   vector<bool>    *baseLinePhoton_SWeird;
   vector<bool>    *baseLinePhoton_ShasGS;
   vector<float>   *baseLinePhoton_Sieie;
   vector<float>   *baseLinePhoton_Sieip;
   vector<float>   *baseLinePhoton_SigmaIEtaIEta;
   vector<float>   *baseLinePhoton_Sipip;
   vector<float>   *baseLinePhoton_SusyId;
   vector<float>   *baseLinePhoton_Time;
   vector<float>   *baseLinePhoton_TrkSumPtHollowConeDR03;
   vector<float>   *baseLinePhoton_TrkSumPtHollowConeDR04;
   vector<float>   *baseLinePhoton_TrkSumPtSolidConeDR04;
   vector<float>   *baseLinePhoton_WTRes;
   vector<float>   *baseLinePhoton_WTResOld;
   vector<float>   *baseLinePhoton_WTime;
   vector<float>   *baseLinePhoton_WTimeOld;
   vector<float>   *baseLinePhoton_WTimeSig;
   vector<float>   *baseLinePhoton_WTimeSigOld;
   vector<float>   *baseLinePhoton_beamHaloCNNScore;
   vector<bool>    *baseLinePhoton_beamHaloCR;
   vector<bool>    *baseLinePhoton_beamHaloID;
   vector<float>   *baseLinePhoton_isoANNScore;
   vector<bool>    *baseLinePhoton_looseNotTightIso1ID;
   vector<bool>    *baseLinePhoton_looseNotTightIso2ID;
   vector<bool>    *baseLinePhoton_looseNotTightIsoID;
   vector<float>   *baseLinePhoton_minJetDr;
   vector<float>   *baseLinePhoton_nonIsoANNScore;
   vector<bool>    *baseLinePhoton_notBeamHaloID;
   vector<float>   *baseLinePhoton_physBkgCNNScore;
   vector<bool>    *baseLinePhoton_spikeCR;
   vector<bool>    *baseLinePhoton_tightIsoID;
   vector<bool>    *endcap_photon_baseline;
   vector<float>   *endcap_photon_energy;
   vector<float>   *endcap_photon_erhsecdr4;
   vector<float>   *endcap_photon_eta;
   vector<float>   *endcap_photon_htoem;
   vector<float>   *endcap_photon_isoANNScore;
   vector<float>   *endcap_photon_phi;
   vector<float>   *endcap_photon_pt;
   vector<float>   *endcap_photon_tspscdr4;
   UInt_t          nBaseLinePhotons;
   UInt_t          nPhotons;
   UInt_t          nblBHPhotons;
   UInt_t          nblLooseNotTightIsoPhotons;
   UInt_t          nblNonpromptPhotons;
   UInt_t          nblNotBHPhotons;
   UInt_t          nblPromptPhotons;
   UInt_t          nblTightIsoPhotons;
   Bool_t          passNPhoEq1SelectionBeamHaloCR;
   Bool_t          passNPhoEq1SelectionEarlyBeamHaloCR;
   Bool_t          passNPhoEq1SelectionEarlyNotBHCR;
   Bool_t          passNPhoEq1SelectionLateBeamHaloCR;
   Bool_t          passNPhoEq1SelectionLateNotBHSR;
   Bool_t          passNPhoEq1SelectionPromptLooseNotTightIso1CR;
   Bool_t          passNPhoEq1SelectionPromptLooseNotTightIso2CR;
   Bool_t          passNPhoEq1SelectionPromptLooseNotTightIsoCR;
   Bool_t          passNPhoEq1SelectionPromptTightIsoSR;
   Bool_t          passNPhoEq2SelectionBeamHaloCR;
   Bool_t          passNPhoEq2SelectionEarlyBeamHaloCR;
   Bool_t          passNPhoEq2SelectionEarlyNotBHCR;
   Bool_t          passNPhoEq2SelectionLateBeamHaloCR;
   Bool_t          passNPhoEq2SelectionLateNotBHSR;
   Bool_t          passNPhoEq2SelectionPromptLooseNotTightIso1CR;
   Bool_t          passNPhoEq2SelectionPromptLooseNotTightIso2CR;
   Bool_t          passNPhoEq2SelectionPromptLooseNotTightIsoCR;
   Bool_t          passNPhoEq2SelectionPromptTightIsoSR;
   Bool_t          passNPhoGe1SelectionBeamHaloCR;
   Bool_t          passNPhoGe1SelectionEarlyBeamHaloCR;
   Bool_t          passNPhoGe1SelectionEarlyNotBHCR;
   Bool_t          passNPhoGe1SelectionLateBeamHaloCR;
   Bool_t          passNPhoGe1SelectionLateNotBHSR;
   Bool_t          passNPhoGe1SelectionNotBHCR;
   Bool_t          passNPhoGe1SelectionNotBHPrompt;
   Bool_t          passNPhoGe1SelectionPromptLooseNotTightIso1CR;
   Bool_t          passNPhoGe1SelectionPromptLooseNotTightIso2CR;
   Bool_t          passNPhoGe1SelectionPromptLooseNotTightIsoCR;
   Bool_t          passNPhoGe1SelectionPromptTightIsoSR;
   vector<float>   *photon_E;
   vector<float>   *photon_Eta;
   vector<float>   *photon_LeadTime;
   vector<float>   *photon_LeadTimeSig;
   vector<float>   *photon_Phi;
   vector<bool>    *photon_PixSeed;
   vector<float>   *photon_Pt;
   vector<float>   *photon_WTime;
   vector<float>   *photon_WTimeSig;
   vector<bool>    *photon_baseline;
   vector<float>   *photon_beamHaloCNNScore;
   vector<float>   *photon_isoANNScore;
   vector<float>   *photon_minJetDr;
   Float_t         tagged_lead_timesig;
   vector<bool>    *rjrABSide;
   vector<float>   *rjr_Rdiff;
   vector<bool>    *rjrCleaningVeto0;
   vector<bool>    *rjrCleaningVeto1;
   vector<bool>    *rjrCleaningVeto2;
   vector<float>   *rjrDPhiSI;
   vector<float>   *rjrDphiMETV;
   vector<float>   *rjrEVa;
   vector<float>   *rjrEVb;
   vector<int>     *rjrLeadPhoLocation;
   vector<float>   *rjrMET;
   vector<float>   *rjrN2Px;
   vector<float>   *rjrN2Py;
   vector<float>   *rjrN2Pz;
   vector<int>     *rjrNJets;
   vector<int>     *rjrNJetsJa;
   vector<int>     *rjrNJetsJb;
   vector<int>     *rjrNPhotons;
   vector<int>     *rjrNVisObjects;
   vector<float>   *rjrPS;
   vector<float>   *rjrPTS;
   vector<float>   *rjrPVa;
   vector<float>   *rjrPVb;
   vector<float>   *rjrPVlab;
   vector<float>   *rjrPZS;
   vector<float>   *rjrSCosA;
   vector<float>   *rjrSMass;
   vector<float>   *rjrSdphiBV;
   vector<float>   *rjrSdphiDA;
   vector<int>     *rjrSubPhoLocation;
   vector<int>     *rjrType;
   vector<bool>    *rjrVisJetSide;
   vector<bool>    *rjrVisPhoSide;
   vector<float>   *rjrX1aCosA;
   vector<float>   *rjrX1aMass;
   vector<float>   *rjrX1a_Eta;
   vector<float>   *rjrX1a_Phi;
   vector<float>   *rjrX1a_Pt;
   vector<float>   *rjrX1bCosA;
   vector<float>   *rjrX1bMass;
   vector<float>   *rjrX1b_Eta;
   vector<float>   *rjrX1b_Phi;
   vector<float>   *rjrX1b_Pt;
   vector<float>   *rjrX2aCosA;
   vector<float>   *rjrX2aMass;
   vector<float>   *rjrX2aPs;
   vector<float>   *rjrX2bCosA;
   vector<float>   *rjrX2bMass;
   vector<float>   *rjrX2bPs;
   vector<float>   *rjr_Mr;
   vector<float>   *rjr_Ms;
   vector<float>   *rjr_Mva;
   vector<float>   *rjr_Mvb;
   vector<float>   *rjr_R;
   vector<float>   *rjr_Rm;
   vector<float>   *rjr_Rs;
   vector<float>   *rjr_Rv;
   vector<float>   *rjr_Rx;
   vector<float>   *rjr_Rxa;
   vector<float>   *rjr_Rxb;
   vector<float>   *rjr_p1Pts;
   vector<float>   *rjr_p1Ptxa;
   vector<float>   *rjr_p1Ptxb;
   vector<float>   *rjr_p1Pxa;
   vector<float>   *rjr_p1Pxb;
   vector<float>   *rjr_p2Pts;
   vector<float>   *rjr_p2Ptxa;
   vector<float>   *rjr_p2Ptxb;
   vector<float>   *rjr_p2Pxa;
   vector<float>   *rjr_p2Pxb;
   vector<float>   *rjr_pHs11;
   vector<float>   *rjr_pHs11a;
   vector<float>   *rjr_pHs11b;
   vector<float>   *rjr_pHs20;
   vector<float>   *rjr_pHs21;
   vector<float>   *rjr_pHs21a;
   vector<float>   *rjr_pHs21b;
   vector<float>   *rjr_pHs22;
   vector<float>   *rjr_pHs40;
   vector<float>   *rjr_pHs41;
   vector<float>   *rjr_pHts11;
   vector<float>   *rjr_pHts11a;
   vector<float>   *rjr_pHts11b;
   vector<float>   *rjr_pHts20;
   vector<float>   *rjr_pHts21;
   vector<float>   *rjr_pHts21a;
   vector<float>   *rjr_pHts21b;
   vector<float>   *rjr_pHts22;
   vector<float>   *rjr_pHts40;
   vector<float>   *rjr_pHts41;
   vector<float>   *rjr_pHtxa10;
   vector<float>   *rjr_pHtxa11;
   vector<float>   *rjr_pHtxa20;
   vector<float>   *rjr_pHtxa21;
   vector<float>   *rjr_pHtxb10;
   vector<float>   *rjr_pHtxb11;
   vector<float>   *rjr_pHtxb20;
   vector<float>   *rjr_pHtxb21;
   vector<float>   *rjr_pHxa10;
   vector<float>   *rjr_pHxa11;
   vector<float>   *rjr_pHxa20;
   vector<float>   *rjr_pHxa21;
   vector<float>   *rjr_pHxb10;
   vector<float>   *rjr_pHxb11;
   vector<float>   *rjr_pHxb20;
   vector<float>   *rjr_pHxb21;
   vector<float>   *rjr_sJ1a_e;
   vector<float>   *rjr_sJ1a_px;
   vector<float>   *rjr_sJ1a_py;
   vector<float>   *rjr_sJ1a_pz;
   vector<float>   *rjr_sJ1b_e;
   vector<float>   *rjr_sJ1b_px;
   vector<float>   *rjr_sJ1b_py;
   vector<float>   *rjr_sJ1b_pz;
   vector<float>   *rjr_sJ2a_e;
   vector<float>   *rjr_sJ2a_px;
   vector<float>   *rjr_sJ2a_py;
   vector<float>   *rjr_sJ2a_pz;
   vector<float>   *rjr_sJ2b_e;
   vector<float>   *rjr_sJ2b_px;
   vector<float>   *rjr_sJ2b_py;
   vector<float>   *rjr_sJ2b_pz;
   Bool_t          rjrIsrABSide;
   Float_t         rjr_Rdiff;
   Bool_t          rjrIsrCleaningVeto0;
   Bool_t          rjrIsrCleaningVeto1;
   Bool_t          rjrIsrCleaningVeto2;
   Float_t         rjrIsrDPhiSI;
   Float_t         rjrIsrDphiMETV;
   Float_t         rjrIsrEVa;
   Float_t         rjrIsrEVb;
   Int_t           rjrIsrLeadPhoLocation;
   Float_t         rjrIsrMET;
   Int_t           rjrIsrNISRPhos;
   Int_t           rjrIsrNJets;
   Int_t           rjrIsrNJetsISR;
   Int_t           rjrIsrNJetsJa;
   Int_t           rjrIsrNJetsJb;
   Int_t           rjrIsrNPhotons;
   Int_t           rjrIsrNVisObjects;
   Float_t         rjrIsrPS;
   Float_t         rjrIsrPTS;
   Float_t         rjrIsrPVa;
   Float_t         rjrIsrPVb;
   Float_t         rjrIsrPVlab;
   Float_t         rjrIsrPZS;
   Float_t         rjrIsrSCosA;
   Float_t         rjrIsrSMass;
   Float_t         rjrIsrSdphiBV;
   Float_t         rjrIsrSdphiDA;
   Int_t           rjrIsrSubPhoLocation;
   Int_t           rjrIsrVisJetSide;
   Int_t           rjrIsrVisPhoSide;
   Float_t         rjrIsrX1aCosA;
   Float_t         rjrIsrX1aMass;
   Float_t         rjrIsrX1a_Eta;
   Float_t         rjrIsrX1a_Phi;
   Float_t         rjrIsrX1a_Pt;
   Float_t         rjrIsrX1bCosA;
   Float_t         rjrIsrX1bMass;
   Float_t         rjrIsrX1b_Eta;
   Float_t         rjrIsrX1b_Phi;
   Float_t         rjrIsrX1b_Pt;
   Float_t         rjrIsrX2aCosA;
   Float_t         rjrIsrX2aMass;
   Float_t         rjrIsrX2aPs;
   Float_t         rjrIsrX2bCosA;
   Float_t         rjrIsrX2bMass;
   Float_t         rjrIsrX2bPs;
   Float_t         rjrIsr_Mr;
   Float_t         rjrIsr_Ms;
   Float_t         rjrIsr_MsPerp;
   Float_t         rjrIsr_Mva;
   Float_t         rjrIsr_Mvb;
   Float_t         rjrIsr_PIsr;
   Float_t         rjrIsr_PtIsr;
   Float_t         rjrIsr_R;
   Float_t         rjrIsr_RIsr;
   Float_t         rjrIsr_Rm;
   Float_t         rjrIsr_Rs;
   Float_t         rjrIsr_Rv;
   Float_t         rjrIsr_Rx;
   Float_t         rjrIsr_Rxa;
   Float_t         rjrIsr_Rxb;
   Int_t           rjrIsr_nIsrVisObjects;
   Int_t           rjrIsr_nSVisObjects;
   Float_t         rjrIsr_p1Pts;
   Float_t         rjrIsr_p1Ptxa;
   Float_t         rjrIsr_p1Ptxb;
   Float_t         rjrIsr_p1Pxa;
   Float_t         rjrIsr_p1Pxb;
   Float_t         rjrIsr_p2Pts;
   Float_t         rjrIsr_p2Ptxa;
   Float_t         rjrIsr_p2Ptxb;
   Float_t         rjrIsr_p2Pxa;
   Float_t         rjrIsr_p2Pxb;
   Float_t         rjrIsr_pHs11;
   Float_t         rjrIsr_pHs11a;
   Float_t         rjrIsr_pHs11b;
   Float_t         rjrIsr_pHs20;
   Float_t         rjrIsr_pHs21;
   Float_t         rjrIsr_pHs21a;
   Float_t         rjrIsr_pHs21b;
   Float_t         rjrIsr_pHs22;
   Float_t         rjrIsr_pHs40;
   Float_t         rjrIsr_pHs41;
   Float_t         rjrIsr_pHts11;
   Float_t         rjrIsr_pHts11a;
   Float_t         rjrIsr_pHts11b;
   Float_t         rjrIsr_pHts20;
   Float_t         rjrIsr_pHts21;
   Float_t         rjrIsr_pHts21a;
   Float_t         rjrIsr_pHts21b;
   Float_t         rjrIsr_pHts22;
   Float_t         rjrIsr_pHts40;
   Float_t         rjrIsr_pHts41;
   Float_t         rjrIsr_pHtxa10;
   Float_t         rjrIsr_pHtxa11;
   Float_t         rjrIsr_pHtxa20;
   Float_t         rjrIsr_pHtxa21;
   Float_t         rjrIsr_pHtxb10;
   Float_t         rjrIsr_pHtxb11;
   Float_t         rjrIsr_pHtxb20;
   Float_t         rjrIsr_pHtxb21;
   Float_t         rjrIsr_pHxa10;
   Float_t         rjrIsr_pHxa11;
   Float_t         rjrIsr_pHxa20;
   Float_t         rjrIsr_pHxa21;
   Float_t         rjrIsr_pHxb10;
   Float_t         rjrIsr_pHxb11;
   Float_t         rjrIsr_pHxb20;
   Float_t         rjrIsr_pHxb21;
   Float_t         rjrIsr_sJ1a_e;
   Float_t         rjrIsr_sJ1a_px;
   Float_t         rjrIsr_sJ1a_py;
   Float_t         rjrIsr_sJ1a_pz;
   Float_t         rjrIsr_sJ1b_e;
   Float_t         rjrIsr_sJ1b_px;
   Float_t         rjrIsr_sJ1b_py;
   Float_t         rjrIsr_sJ1b_pz;
   Float_t         rjrIsr_sJ2a_e;
   Float_t         rjrIsr_sJ2a_px;
   Float_t         rjrIsr_sJ2a_py;
   Float_t         rjrIsr_sJ2a_pz;
   Float_t         rjrIsr_sJ2b_e;
   Float_t         rjrIsr_sJ2b_px;
   Float_t         rjrIsr_sJ2b_py;
   Float_t         rjrIsr_sJ2b_pz;
   vector<bool>    *eleHasSVMatch;
   vector<bool>    *eleIsLoose;
   vector<float>   *elePhoMinDr;
   vector<float>   *epDrMatch;
   vector<float>   *epSeedIdMatch;
   vector<float>   *epSeedIdUnique;
   Int_t           nEleSVMatched;
   UInt_t          nElectrons;
   UInt_t          nSelElectrons;
   UInt_t          nSelIsoElectrons;
   vector<float>   *allJetEnergy;
   vector<float>   *allJetEta;
   vector<float>   *allJetMass;
   vector<int>     *allJetNPho;
   vector<float>   *allJetPhi;
   vector<int>     *allJetPixNPho;
   vector<float>   *allJetPixWTime;
   vector<float>   *allJetPixWTimeRes;
   vector<float>   *allJetPixWTimeSig;
   vector<float>   *allJetPt;
   vector<vector<unsigned int> > *allJetSCIndexs;
   vector<float>   *allJetWTime;
   vector<float>   *allJetWTimeRes;
   vector<float>   *allJetWTimeSig;
   vector<int>     *allJetsWithWPixTime;
   vector<int>     *allJetsWithWTime;
   UInt_t          nJets;
   UInt_t          nSelJets;
   vector<int>     *nSelJetsWithWPixTime;
   vector<int>     *nSelJetsWithWTime;
   vector<float>   *pv_dDiJetPEnergy;
   vector<float>   *pv_dDiJetTime;
   vector<float>   *pv_dGJMinDr;
   vector<float>   *pv_dGJPEnergy;
   vector<float>   *pv_dGJTime;
   vector<float>   *pv_diJet1Time;
   vector<float>   *pv_diJet2Time;
   vector<float>   *pv_gjGammaTime;
   vector<float>   *pv_gjJetTime;
   vector<int>     *pv_gjPhoIndex;
   vector<int>     *pv_haspvtime;
   vector<float>   *pv_wpenergy;
   vector<float>   *pv_wphartres;
   vector<float>   *pv_wtime;
   vector<float>   *pv_wtimeres;
   vector<float>   *selGenJetDpt;
   vector<float>   *selGenJetEnergy;
   vector<float>   *selGenJetImpAng;
   vector<float>   *selGenJetLlpTime;
   vector<float>   *selGenJetPt;
   vector<float>   *selGenJetTime;
   vector<float>   *selGenJetTof;
   vector<float>   *selGenJetdr;
   vector<float>   *selGenJeteta;
   vector<float>   *selJetArea;
   vector<float>   *selJetChEmEF;
   vector<float>   *selJetChHM;
   vector<float>   *selJetEnergy;
   vector<float>   *selJetEta;
   vector<int>     *selJetGenLlpId;
   vector<float>   *selJetLlpDp;
   vector<float>   *selJetLlpDr;
   vector<float>   *selJetMass;
   vector<float>   *selJetMuEF;
   vector<float>   *selJetNeEmEF;
   vector<float>   *selJetNeHEF;
   vector<float>   *selJetNeHM;
   vector<float>   *selJetPhi;
   vector<float>   *selJetPixWTime;
   vector<float>   *selJetPixWTimeRes;
   vector<float>   *selJetPixWTimeSig;
   vector<float>   *selJetPt;
   vector<int>     *selJetQrkLlpId;
   vector<int>     *selJetQuality;
   vector<float>   *selJetTime;
   vector<float>   *selJetWTime;
   vector<float>   *selJetWTimeRes;
   vector<float>   *selJetWTimeSig;
   vector<float>   *selJetchHEF;
   vector<bool>    *muonIsLoose;
   vector<bool>    *muonIsMedium;
   UInt_t          nLooseMuons;
   Int_t           nMuonSVMatched;
   UInt_t          nMuons;
   vector<float>   *HadronicSV_cosTheta;
   vector<float>   *HadronicSV_decayAngle;
   vector<float>   *HadronicSV_dxy;
   vector<float>   *HadronicSV_dxySig;
   vector<float>   *HadronicSV_mass;
   vector<float>   *HadronicSV_massOverNtracks;
   vector<float>   *HadronicSV_matchRatio;
   UInt_t          HadronicSV_nTotal;
   vector<int>     *HadronicSV_nTracks;
   vector<float>   *HadronicSV_pOverE;
   vector<bool>    *HadronicSV_passMaterialVeto;
   vector<bool>    *HadronicSV_tightZWindowSelection;
   vector<float>   *HadronicSV_x;
   vector<float>   *HadronicSV_y;
   vector<float>   *HadronicSV_z;
   vector<float>   *LeptonicSV_cosTheta;
   vector<float>   *LeptonicSV_decayAngle;
   vector<float>   *LeptonicSV_dxy;
   vector<float>   *LeptonicSV_dxySig;
   vector<int>     *LeptonicSV_electronIndex;
   vector<bool>    *LeptonicSV_isGold;
   vector<float>   *LeptonicSV_mass;
   vector<int>     *LeptonicSV_muonIndex;
   UInt_t          LeptonicSV_nTotal;
   vector<int>     *LeptonicSV_nTracks;
   vector<float>   *LeptonicSV_pOverE;
   vector<bool>    *LeptonicSV_passElectronID;
   vector<bool>    *LeptonicSV_passMaterialVeto;
   vector<bool>    *LeptonicSV_passMuonID;
   vector<bool>    *LeptonicSV_tightZWindowSelection;
   vector<float>   *LeptonicSV_x;
   vector<float>   *LeptonicSV_y;
   vector<float>   *LeptonicSV_z;
   Int_t           SV_nElectron;
   UInt_t          SV_nElectronPassZWindow;
   Int_t           SV_nHadCmpDxyGeq2;
   Int_t           SV_nHadCmpDxySigGeq1000;
   Int_t           SV_nHadStdDxyGeq2;
   Int_t           SV_nHadStdDxySigGeq1000;
   UInt_t          SV_nHadronPassZWindow;
   Int_t           SV_nHadronic;
   Int_t           SV_nLepCmpDxyGeq2;
   Int_t           SV_nLepCmpDxySigGeq1000;
   Int_t           SV_nLepStdDxyGeq2;
   Int_t           SV_nLepStdDxySigGeq1000;
   UInt_t          SV_nLeptonPassZWindow;
   Int_t           SV_nLeptonic;
   Int_t           SV_nMuon;
   UInt_t          SV_nMuonPassZWindow;
   Bool_t          passGe1HadAndLepSelectionCRCmp;
   Bool_t          passGe1HadAndLepSelectionCRLoose;
   Bool_t          passGe1HadAndLepSelectionCRTight;
   Bool_t          passGe1HadAndLepSelectionSRCmp;
   Bool_t          passGe1HadAndLepSelectionSRLoose;
   Bool_t          passGe1HadAndLepSelectionSRTight;
   Bool_t          passNHad1SelectionCRCmp;
   Bool_t          passNHad1SelectionCRLoose;
   Bool_t          passNHad1SelectionCRTight;
   Bool_t          passNHad1SelectionSRCmp;
   Bool_t          passNHad1SelectionSRLoose;
   Bool_t          passNHad1SelectionSRTight;
   Bool_t          passNHadGe2SelectionCRCmp;
   Bool_t          passNHadGe2SelectionCRLoose;
   Bool_t          passNHadGe2SelectionCRTight;
   Bool_t          passNHadGe2SelectionSRCmp;
   Bool_t          passNHadGe2SelectionSRLoose;
   Bool_t          passNHadGe2SelectionSRTight;
   Bool_t          passNLep1SelectionCRCmp;
   Bool_t          passNLep1SelectionCRLoose;
   Bool_t          passNLep1SelectionCRTight;
   Bool_t          passNLep1SelectionSRCmp;
   Bool_t          passNLep1SelectionSRLoose;
   Bool_t          passNLep1SelectionSRTight;
   Bool_t          passNLepGe2SelectionCRCmp;
   Bool_t          passNLepGe2SelectionCRLoose;
   Bool_t          passNLepGe2SelectionCRTight;
   Bool_t          passNLepGe2SelectionSRCmp;
   Bool_t          passNLepGe2SelectionSRLoose;
   Bool_t          passNLepGe2SelectionSRTight;

   // List of branches
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_MetFilters;   //!
   TBranch        *b_Flag_eLHemVeto;   //!
   TBranch        *b_Flag_eMHemVeto;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_hemRun;   //!
   TBranch        *b_Flag_hemVeto;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_Flag_jLHemVeto;   //!
   TBranch        *b_Flag_jMHemVeto;   //!
   TBranch        *b_Flag_mLHemVeto;   //!
   TBranch        *b_Flag_mMHemVeto;   //!
   TBranch        *b_Flag_pLHemVeto;   //!
   TBranch        *b_Flag_pMHemVeto;   //!
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_Trigger_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_Trigger_PFMET120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_DataSetKey;   //!
   TBranch        *b_evtFillWgt;   //!
   TBranch        *b_evtGenWgt;   //!
   TBranch        *b_evtXSection;   //!
   TBranch        *b_selCMet;   //!
   TBranch        *b_selCMetPx;   //!
   TBranch        *b_selCMetPy;   //!
   TBranch        *b_selMet;   //!
   TBranch        *b_selMetPx;   //!
   TBranch        *b_selMetPy;   //!
   TBranch        *b_Trk_nIsoTrack;   //!
   TBranch        *b_Evt_isGG;   //!
   TBranch        *b_Evt_isGZ;   //!
   TBranch        *b_Evt_isZZ;   //!
   TBranch        *b_Evt_nXs;   //!
   TBranch        *b_Xa_Displacment;   //!
   TBranch        *b_Xa_beta;   //!
   TBranch        *b_Xa_ctau;   //!
   TBranch        *b_Xa_energy;   //!
   TBranch        *b_Xa_mass;   //!
   TBranch        *b_Xa_p;   //!
   TBranch        *b_Xa_pdgId;   //!
   TBranch        *b_Xa_phi;   //!
   TBranch        *b_Xa_pt;   //!
   TBranch        *b_Xa_vx;   //!
   TBranch        *b_Xa_vy;   //!
   TBranch        *b_Xa_vz;   //!
   TBranch        *b_Xb_Displacment;   //!
   TBranch        *b_Xb_beta;   //!
   TBranch        *b_Xb_ctau;   //!
   TBranch        *b_Xb_energy;   //!
   TBranch        *b_Xb_mass;   //!
   TBranch        *b_Xb_p;   //!
   TBranch        *b_Xb_pdgId;   //!
   TBranch        *b_Xb_phi;   //!
   TBranch        *b_Xb_pt;   //!
   TBranch        *b_Xb_vx;   //!
   TBranch        *b_Xb_vy;   //!
   TBranch        *b_Xb_vz;   //!
   TBranch        *b_genCharge;   //!
   TBranch        *b_genGrvtinoMass;   //!
   TBranch        *b_genLLPMass;   //!
   TBranch        *b_genLSPMass;   //!
   TBranch        *b_genMass;   //!
   TBranch        *b_genMomCTau;   //!
   TBranch        *b_genPartEnergy;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartMomIdx;   //!
   TBranch        *b_genPartPdgId;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartSusId;   //!
   TBranch        *b_genPx;   //!
   TBranch        *b_genPy;   //!
   TBranch        *b_genPz;   //!
   TBranch        *b_genSGMass;   //!
   TBranch        *b_genSQMass;   //!
   TBranch        *b_genSTFlagQP;   //!
   TBranch        *b_genSTFlagQQP;   //!
   TBranch        *b_genSTFlagQQZ;   //!
   TBranch        *b_genSTFlagQQZL;   //!
   TBranch        *b_genSTFlagQQZLL;   //!
   TBranch        *b_genSigType;   //!
   TBranch        *b_genStatus;   //!
   TBranch        *b_genVx;   //!
   TBranch        *b_genVy;   //!
   TBranch        *b_genVz;   //!
   TBranch        *b_genXMomCTau;   //!
   TBranch        *b_barrel_photon_baseline;   //!
   TBranch        *b_barrel_photon_energy;   //!
   TBranch        *b_barrel_photon_erhsecdr4;   //!
   TBranch        *b_barrel_photon_eta;   //!
   TBranch        *b_barrel_photon_htoem;   //!
   TBranch        *b_barrel_photon_isoANNScore;   //!
   TBranch        *b_barrel_photon_phi;   //!
   TBranch        *b_barrel_photon_pt;   //!
   TBranch        *b_barrel_photon_tspscdr4;   //!
   TBranch        *b_baseLinePhoton_ClstrRn;   //!
   TBranch        *b_baseLinePhoton_CorEnergy;   //!
   TBranch        *b_baseLinePhoton_CorPt;   //!
   TBranch        *b_baseLinePhoton_CovEtaEta;   //!
   TBranch        *b_baseLinePhoton_CovEtaPhi;   //!
   TBranch        *b_baseLinePhoton_CovPhiPhi;   //!
   TBranch        *b_baseLinePhoton_DiJetsCR;   //!
   TBranch        *b_baseLinePhoton_EcalPFClusterIso;   //!
   TBranch        *b_baseLinePhoton_EcalRHSumEtConeDR04;   //!
   TBranch        *b_baseLinePhoton_EleVeto;   //!
   TBranch        *b_baseLinePhoton_Energy;   //!
   TBranch        *b_baseLinePhoton_Eta;   //!
   TBranch        *b_baseLinePhoton_EtaWidth;   //!
   TBranch        *b_baseLinePhoton_GJetsCR;   //!
   TBranch        *b_baseLinePhoton_GenIdx;   //!
   TBranch        *b_baseLinePhoton_GenLabTime;   //!
   TBranch        *b_baseLinePhoton_GenLabTimeSig;   //!
   TBranch        *b_baseLinePhoton_GenPt;   //!
   TBranch        *b_baseLinePhoton_GenSigMomEnergy;   //!
   TBranch        *b_baseLinePhoton_GenSigMomEta;   //!
   TBranch        *b_baseLinePhoton_GenSigMomMass;   //!
   TBranch        *b_baseLinePhoton_GenSigMomPhi;   //!
   TBranch        *b_baseLinePhoton_GenSigMomPt;   //!
   TBranch        *b_baseLinePhoton_GenSigMomPx;   //!
   TBranch        *b_baseLinePhoton_GenSigMomPy;   //!
   TBranch        *b_baseLinePhoton_GenSigMomPz;   //!
   TBranch        *b_baseLinePhoton_GenSigMomVx;   //!
   TBranch        *b_baseLinePhoton_GenSigMomVy;   //!
   TBranch        *b_baseLinePhoton_GenSigMomVz;   //!
   TBranch        *b_baseLinePhoton_GenTime;   //!
   TBranch        *b_baseLinePhoton_GenTimeSig;   //!
   TBranch        *b_baseLinePhoton_HadOverEM;   //!
   TBranch        *b_baseLinePhoton_HadTowOverEM;   //!
   TBranch        *b_baseLinePhoton_HasConversionTracks;   //!
   TBranch        *b_baseLinePhoton_HcalPFClusterIso;   //!
   TBranch        *b_baseLinePhoton_HcalTowerSumEtBcConeDR04;   //!
   TBranch        *b_baseLinePhoton_HcalTowerSumEtConeDR04;   //!
   TBranch        *b_baseLinePhoton_Hoe_PUcorr;   //!
   TBranch        *b_baseLinePhoton_LSCross;   //!
   TBranch        *b_baseLinePhoton_LTRes;   //!
   TBranch        *b_baseLinePhoton_LTime;   //!
   TBranch        *b_baseLinePhoton_LTimeSig;   //!
   TBranch        *b_baseLinePhoton_LWeird;   //!
   TBranch        *b_baseLinePhoton_MBetaEql;   //!
   TBranch        *b_baseLinePhoton_MBetaPmt;   //!
   TBranch        *b_baseLinePhoton_NTrkHollowConeDR04;   //!
   TBranch        *b_baseLinePhoton_NTrkSolidConeDR04;   //!
   TBranch        *b_baseLinePhoton_Nrh;   //!
   TBranch        *b_baseLinePhoton_OOT;   //!
   TBranch        *b_baseLinePhoton_PVtoRH;   //!
   TBranch        *b_baseLinePhoton_PfChargedIso;   //!
   TBranch        *b_baseLinePhoton_PfChargedIsoPFPV;   //!
   TBranch        *b_baseLinePhoton_PfChargedIsoWorstVtx;   //!
   TBranch        *b_baseLinePhoton_PfPhoIso03;   //!
   TBranch        *b_baseLinePhoton_PfRelIso03_all_quadratic;   //!
   TBranch        *b_baseLinePhoton_PfRelIso03_chg_quadratic;   //!
   TBranch        *b_baseLinePhoton_Phi;   //!
   TBranch        *b_baseLinePhoton_PhiWidth;   //!
   TBranch        *b_baseLinePhoton_PhoIsoDr;   //!
   TBranch        *b_baseLinePhoton_PixelSeed;   //!
   TBranch        *b_baseLinePhoton_Pt;   //!
   TBranch        *b_baseLinePhoton_Quality;   //!
   TBranch        *b_baseLinePhoton_R9;   //!
   TBranch        *b_baseLinePhoton_RawTime;   //!
   TBranch        *b_baseLinePhoton_S4;   //!
   TBranch        *b_baseLinePhoton_SAlp;   //!
   TBranch        *b_baseLinePhoton_SCx;   //!
   TBranch        *b_baseLinePhoton_SCy;   //!
   TBranch        *b_baseLinePhoton_SCz;   //!
   TBranch        *b_baseLinePhoton_SMaj;   //!
   TBranch        *b_baseLinePhoton_SMin;   //!
   TBranch        *b_baseLinePhoton_SSCross;   //!
   TBranch        *b_baseLinePhoton_STRes;   //!
   TBranch        *b_baseLinePhoton_STime;   //!
   TBranch        *b_baseLinePhoton_STimeSig;   //!
   TBranch        *b_baseLinePhoton_SWeird;   //!
   TBranch        *b_baseLinePhoton_ShasGS;   //!
   TBranch        *b_baseLinePhoton_Sieie;   //!
   TBranch        *b_baseLinePhoton_Sieip;   //!
   TBranch        *b_baseLinePhoton_SigmaIEtaIEta;   //!
   TBranch        *b_baseLinePhoton_Sipip;   //!
   TBranch        *b_baseLinePhoton_SusyId;   //!
   TBranch        *b_baseLinePhoton_Time;   //!
   TBranch        *b_baseLinePhoton_TrkSumPtHollowConeDR03;   //!
   TBranch        *b_baseLinePhoton_TrkSumPtHollowConeDR04;   //!
   TBranch        *b_baseLinePhoton_TrkSumPtSolidConeDR04;   //!
   TBranch        *b_baseLinePhoton_WTRes;   //!
   TBranch        *b_baseLinePhoton_WTResOld;   //!
   TBranch        *b_baseLinePhoton_WTime;   //!
   TBranch        *b_baseLinePhoton_WTimeOld;   //!
   TBranch        *b_baseLinePhoton_WTimeSig;   //!
   TBranch        *b_baseLinePhoton_WTimeSigOld;   //!
   TBranch        *b_baseLinePhoton_beamHaloCNNScore;   //!
   TBranch        *b_baseLinePhoton_beamHaloCR;   //!
   TBranch        *b_baseLinePhoton_beamHaloID;   //!
   TBranch        *b_baseLinePhoton_isoANNScore;   //!
   TBranch        *b_baseLinePhoton_looseNotTightIso1ID;   //!
   TBranch        *b_baseLinePhoton_looseNotTightIso2ID;   //!
   TBranch        *b_baseLinePhoton_looseNotTightIsoID;   //!
   TBranch        *b_baseLinePhoton_minJetDr;   //!
   TBranch        *b_baseLinePhoton_nonIsoANNScore;   //!
   TBranch        *b_baseLinePhoton_notBeamHaloID;   //!
   TBranch        *b_baseLinePhoton_physBkgCNNScore;   //!
   TBranch        *b_baseLinePhoton_spikeCR;   //!
   TBranch        *b_baseLinePhoton_tightIsoID;   //!
   TBranch        *b_endcap_photon_baseline;   //!
   TBranch        *b_endcap_photon_energy;   //!
   TBranch        *b_endcap_photon_erhsecdr4;   //!
   TBranch        *b_endcap_photon_eta;   //!
   TBranch        *b_endcap_photon_htoem;   //!
   TBranch        *b_endcap_photon_isoANNScore;   //!
   TBranch        *b_endcap_photon_phi;   //!
   TBranch        *b_endcap_photon_pt;   //!
   TBranch        *b_endcap_photon_tspscdr4;   //!
   TBranch        *b_nBaseLinePhotons;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_nblBHPhotons;   //!
   TBranch        *b_nblLooseNotTightIsoPhotons;   //!
   TBranch        *b_nblNonpromptPhotons;   //!
   TBranch        *b_nblNotBHPhotons;   //!
   TBranch        *b_nblPromptPhotons;   //!
   TBranch        *b_nblTightIsoPhotons;   //!
   TBranch        *b_passNPhoEq1SelectionBeamHaloCR;   //!
   TBranch        *b_passNPhoEq1SelectionEarlyBeamHaloCR;   //!
   TBranch        *b_passNPhoEq1SelectionEarlyNotBHCR;   //!
   TBranch        *b_passNPhoEq1SelectionLateBeamHaloCR;   //!
   TBranch        *b_passNPhoEq1SelectionLateNotBHSR;   //!
   TBranch        *b_passNPhoEq1SelectionPromptLooseNotTightIso1CR;   //!
   TBranch        *b_passNPhoEq1SelectionPromptLooseNotTightIso2CR;   //!
   TBranch        *b_passNPhoEq1SelectionPromptLooseNotTightIsoCR;   //!
   TBranch        *b_passNPhoEq1SelectionPromptTightIsoSR;   //!
   TBranch        *b_passNPhoEq2SelectionBeamHaloCR;   //!
   TBranch        *b_passNPhoEq2SelectionEarlyBeamHaloCR;   //!
   TBranch        *b_passNPhoEq2SelectionEarlyNotBHCR;   //!
   TBranch        *b_passNPhoEq2SelectionLateBeamHaloCR;   //!
   TBranch        *b_passNPhoEq2SelectionLateNotBHSR;   //!
   TBranch        *b_passNPhoEq2SelectionPromptLooseNotTightIso1CR;   //!
   TBranch        *b_passNPhoEq2SelectionPromptLooseNotTightIso2CR;   //!
   TBranch        *b_passNPhoEq2SelectionPromptLooseNotTightIsoCR;   //!
   TBranch        *b_passNPhoEq2SelectionPromptTightIsoSR;   //!
   TBranch        *b_passNPhoGe1SelectionBeamHaloCR;   //!
   TBranch        *b_passNPhoGe1SelectionEarlyBeamHaloCR;   //!
   TBranch        *b_passNPhoGe1SelectionEarlyNotBHCR;   //!
   TBranch        *b_passNPhoGe1SelectionLateBeamHaloCR;   //!
   TBranch        *b_passNPhoGe1SelectionLateNotBHSR;   //!
   TBranch        *b_passNPhoGe1SelectionNotBHCR;   //!
   TBranch        *b_passNPhoGe1SelectionNotBHPrompt;   //!
   TBranch        *b_passNPhoGe1SelectionPromptLooseNotTightIso1CR;   //!
   TBranch        *b_passNPhoGe1SelectionPromptLooseNotTightIso2CR;   //!
   TBranch        *b_passNPhoGe1SelectionPromptLooseNotTightIsoCR;   //!
   TBranch        *b_passNPhoGe1SelectionPromptTightIsoSR;   //!
   TBranch        *b_photon_E;   //!
   TBranch        *b_photon_Eta;   //!
   TBranch        *b_photon_LeadTime;   //!
   TBranch        *b_photon_LeadTimeSig;   //!
   TBranch        *b_photon_Phi;   //!
   TBranch        *b_photon_PixSeed;   //!
   TBranch        *b_photon_Pt;   //!
   TBranch        *b_photon_WTime;   //!
   TBranch        *b_photon_WTimeSig;   //!
   TBranch        *b_photon_baseline;   //!
   TBranch        *b_photon_beamHaloCNNScore;   //!
   TBranch        *b_photon_isoANNScore;   //!
   TBranch        *b_photon_minJetDr;   //!
   TBranch        *b_tagged_lead_timesig;   //!
   TBranch        *b_rjrABSide;   //!
   TBranch        *b_rjr_Rdiff;   //!
   TBranch        *b_rjrCleaningVeto0;   //!
   TBranch        *b_rjrCleaningVeto1;   //!
   TBranch        *b_rjrCleaningVeto2;   //!
   TBranch        *b_rjrDPhiSI;   //!
   TBranch        *b_rjrDphiMETV;   //!
   TBranch        *b_rjrEVa;   //!
   TBranch        *b_rjrEVb;   //!
   TBranch        *b_rjrLeadPhoLocation;   //!
   TBranch        *b_rjrMET;   //!
   TBranch        *b_rjrN2Px;   //!
   TBranch        *b_rjrN2Py;   //!
   TBranch        *b_rjrN2Pz;   //!
   TBranch        *b_rjrNJets;   //!
   TBranch        *b_rjrNJetsJa;   //!
   TBranch        *b_rjrNJetsJb;   //!
   TBranch        *b_rjrNPhotons;   //!
   TBranch        *b_rjrNVisObjects;   //!
   TBranch        *b_rjrPS;   //!
   TBranch        *b_rjrPTS;   //!
   TBranch        *b_rjrPVa;   //!
   TBranch        *b_rjrPVb;   //!
   TBranch        *b_rjrPVlab;   //!
   TBranch        *b_rjrPZS;   //!
   TBranch        *b_rjrSCosA;   //!
   TBranch        *b_rjrSMass;   //!
   TBranch        *b_rjrSdphiBV;   //!
   TBranch        *b_rjrSdphiDA;   //!
   TBranch        *b_rjrSubPhoLocation;   //!
   TBranch        *b_rjrType;   //!
   TBranch        *b_rjrVisJetSide;   //!
   TBranch        *b_rjrVisPhoSide;   //!
   TBranch        *b_rjrX1aCosA;   //!
   TBranch        *b_rjrX1aMass;   //!
   TBranch        *b_rjrX1a_Eta;   //!
   TBranch        *b_rjrX1a_Phi;   //!
   TBranch        *b_rjrX1a_Pt;   //!
   TBranch        *b_rjrX1bCosA;   //!
   TBranch        *b_rjrX1bMass;   //!
   TBranch        *b_rjrX1b_Eta;   //!
   TBranch        *b_rjrX1b_Phi;   //!
   TBranch        *b_rjrX1b_Pt;   //!
   TBranch        *b_rjrX2aCosA;   //!
   TBranch        *b_rjrX2aMass;   //!
   TBranch        *b_rjrX2aPs;   //!
   TBranch        *b_rjrX2bCosA;   //!
   TBranch        *b_rjrX2bMass;   //!
   TBranch        *b_rjrX2bPs;   //!
   TBranch        *b_rjr_Mr;   //!
   TBranch        *b_rjr_Ms;   //!
   TBranch        *b_rjr_Mva;   //!
   TBranch        *b_rjr_Mvb;   //!
   TBranch        *b_rjr_R;   //!
   TBranch        *b_rjr_Rm;   //!
   TBranch        *b_rjr_Rs;   //!
   TBranch        *b_rjr_Rv;   //!
   TBranch        *b_rjr_Rx;   //!
   TBranch        *b_rjr_Rxa;   //!
   TBranch        *b_rjr_Rxb;   //!
   TBranch        *b_rjr_p1Pts;   //!
   TBranch        *b_rjr_p1Ptxa;   //!
   TBranch        *b_rjr_p1Ptxb;   //!
   TBranch        *b_rjr_p1Pxa;   //!
   TBranch        *b_rjr_p1Pxb;   //!
   TBranch        *b_rjr_p2Pts;   //!
   TBranch        *b_rjr_p2Ptxa;   //!
   TBranch        *b_rjr_p2Ptxb;   //!
   TBranch        *b_rjr_p2Pxa;   //!
   TBranch        *b_rjr_p2Pxb;   //!
   TBranch        *b_rjr_pHs11;   //!
   TBranch        *b_rjr_pHs11a;   //!
   TBranch        *b_rjr_pHs11b;   //!
   TBranch        *b_rjr_pHs20;   //!
   TBranch        *b_rjr_pHs21;   //!
   TBranch        *b_rjr_pHs21a;   //!
   TBranch        *b_rjr_pHs21b;   //!
   TBranch        *b_rjr_pHs22;   //!
   TBranch        *b_rjr_pHs40;   //!
   TBranch        *b_rjr_pHs41;   //!
   TBranch        *b_rjr_pHts11;   //!
   TBranch        *b_rjr_pHts11a;   //!
   TBranch        *b_rjr_pHts11b;   //!
   TBranch        *b_rjr_pHts20;   //!
   TBranch        *b_rjr_pHts21;   //!
   TBranch        *b_rjr_pHts21a;   //!
   TBranch        *b_rjr_pHts21b;   //!
   TBranch        *b_rjr_pHts22;   //!
   TBranch        *b_rjr_pHts40;   //!
   TBranch        *b_rjr_pHts41;   //!
   TBranch        *b_rjr_pHtxa10;   //!
   TBranch        *b_rjr_pHtxa11;   //!
   TBranch        *b_rjr_pHtxa20;   //!
   TBranch        *b_rjr_pHtxa21;   //!
   TBranch        *b_rjr_pHtxb10;   //!
   TBranch        *b_rjr_pHtxb11;   //!
   TBranch        *b_rjr_pHtxb20;   //!
   TBranch        *b_rjr_pHtxb21;   //!
   TBranch        *b_rjr_pHxa10;   //!
   TBranch        *b_rjr_pHxa11;   //!
   TBranch        *b_rjr_pHxa20;   //!
   TBranch        *b_rjr_pHxa21;   //!
   TBranch        *b_rjr_pHxb10;   //!
   TBranch        *b_rjr_pHxb11;   //!
   TBranch        *b_rjr_pHxb20;   //!
   TBranch        *b_rjr_pHxb21;   //!
   TBranch        *b_rjr_sJ1a_e;   //!
   TBranch        *b_rjr_sJ1a_px;   //!
   TBranch        *b_rjr_sJ1a_py;   //!
   TBranch        *b_rjr_sJ1a_pz;   //!
   TBranch        *b_rjr_sJ1b_e;   //!
   TBranch        *b_rjr_sJ1b_px;   //!
   TBranch        *b_rjr_sJ1b_py;   //!
   TBranch        *b_rjr_sJ1b_pz;   //!
   TBranch        *b_rjr_sJ2a_e;   //!
   TBranch        *b_rjr_sJ2a_px;   //!
   TBranch        *b_rjr_sJ2a_py;   //!
   TBranch        *b_rjr_sJ2a_pz;   //!
   TBranch        *b_rjr_sJ2b_e;   //!
   TBranch        *b_rjr_sJ2b_px;   //!
   TBranch        *b_rjr_sJ2b_py;   //!
   TBranch        *b_rjr_sJ2b_pz;   //!
   TBranch        *b_rjrIsrABSide;   //!
   TBranch        *b_rjr_Rdiff;   //!
   TBranch        *b_rjrIsrCleaningVeto0;   //!
   TBranch        *b_rjrIsrCleaningVeto1;   //!
   TBranch        *b_rjrIsrCleaningVeto2;   //!
   TBranch        *b_rjrIsrDPhiSI;   //!
   TBranch        *b_rjrIsrDphiMETV;   //!
   TBranch        *b_rjrIsrEVa;   //!
   TBranch        *b_rjrIsrEVb;   //!
   TBranch        *b_rjrIsrLeadPhoLocation;   //!
   TBranch        *b_rjrIsrMET;   //!
   TBranch        *b_rjrIsrNISRPhos;   //!
   TBranch        *b_rjrIsrNJets;   //!
   TBranch        *b_rjrIsrNJetsISR;   //!
   TBranch        *b_rjrIsrNJetsJa;   //!
   TBranch        *b_rjrIsrNJetsJb;   //!
   TBranch        *b_rjrIsrNPhotons;   //!
   TBranch        *b_rjrIsrNVisObjects;   //!
   TBranch        *b_rjrIsrPS;   //!
   TBranch        *b_rjrIsrPTS;   //!
   TBranch        *b_rjrIsrPVa;   //!
   TBranch        *b_rjrIsrPVb;   //!
   TBranch        *b_rjrIsrPVlab;   //!
   TBranch        *b_rjrIsrPZS;   //!
   TBranch        *b_rjrIsrSCosA;   //!
   TBranch        *b_rjrIsrSMass;   //!
   TBranch        *b_rjrIsrSdphiBV;   //!
   TBranch        *b_rjrIsrSdphiDA;   //!
   TBranch        *b_rjrIsrSubPhoLocation;   //!
   TBranch        *b_rjrIsrVisJetSide;   //!
   TBranch        *b_rjrIsrVisPhoSide;   //!
   TBranch        *b_rjrIsrX1aCosA;   //!
   TBranch        *b_rjrIsrX1aMass;   //!
   TBranch        *b_rjrIsrX1a_Eta;   //!
   TBranch        *b_rjrIsrX1a_Phi;   //!
   TBranch        *b_rjrIsrX1a_Pt;   //!
   TBranch        *b_rjrIsrX1bCosA;   //!
   TBranch        *b_rjrIsrX1bMass;   //!
   TBranch        *b_rjrIsrX1b_Eta;   //!
   TBranch        *b_rjrIsrX1b_Phi;   //!
   TBranch        *b_rjrIsrX1b_Pt;   //!
   TBranch        *b_rjrIsrX2aCosA;   //!
   TBranch        *b_rjrIsrX2aMass;   //!
   TBranch        *b_rjrIsrX2aPs;   //!
   TBranch        *b_rjrIsrX2bCosA;   //!
   TBranch        *b_rjrIsrX2bMass;   //!
   TBranch        *b_rjrIsrX2bPs;   //!
   TBranch        *b_rjrIsr_Mr;   //!
   TBranch        *b_rjrIsr_Ms;   //!
   TBranch        *b_rjrIsr_MsPerp;   //!
   TBranch        *b_rjrIsr_Mva;   //!
   TBranch        *b_rjrIsr_Mvb;   //!
   TBranch        *b_rjrIsr_PIsr;   //!
   TBranch        *b_rjrIsr_PtIsr;   //!
   TBranch        *b_rjrIsr_R;   //!
   TBranch        *b_rjrIsr_RIsr;   //!
   TBranch        *b_rjrIsr_Rm;   //!
   TBranch        *b_rjrIsr_Rs;   //!
   TBranch        *b_rjrIsr_Rv;   //!
   TBranch        *b_rjrIsr_Rx;   //!
   TBranch        *b_rjrIsr_Rxa;   //!
   TBranch        *b_rjrIsr_Rxb;   //!
   TBranch        *b_rjrIsr_nIsrVisObjects;   //!
   TBranch        *b_rjrIsr_nSVisObjects;   //!
   TBranch        *b_rjrIsr_p1Pts;   //!
   TBranch        *b_rjrIsr_p1Ptxa;   //!
   TBranch        *b_rjrIsr_p1Ptxb;   //!
   TBranch        *b_rjrIsr_p1Pxa;   //!
   TBranch        *b_rjrIsr_p1Pxb;   //!
   TBranch        *b_rjrIsr_p2Pts;   //!
   TBranch        *b_rjrIsr_p2Ptxa;   //!
   TBranch        *b_rjrIsr_p2Ptxb;   //!
   TBranch        *b_rjrIsr_p2Pxa;   //!
   TBranch        *b_rjrIsr_p2Pxb;   //!
   TBranch        *b_rjrIsr_pHs11;   //!
   TBranch        *b_rjrIsr_pHs11a;   //!
   TBranch        *b_rjrIsr_pHs11b;   //!
   TBranch        *b_rjrIsr_pHs20;   //!
   TBranch        *b_rjrIsr_pHs21;   //!
   TBranch        *b_rjrIsr_pHs21a;   //!
   TBranch        *b_rjrIsr_pHs21b;   //!
   TBranch        *b_rjrIsr_pHs22;   //!
   TBranch        *b_rjrIsr_pHs40;   //!
   TBranch        *b_rjrIsr_pHs41;   //!
   TBranch        *b_rjrIsr_pHts11;   //!
   TBranch        *b_rjrIsr_pHts11a;   //!
   TBranch        *b_rjrIsr_pHts11b;   //!
   TBranch        *b_rjrIsr_pHts20;   //!
   TBranch        *b_rjrIsr_pHts21;   //!
   TBranch        *b_rjrIsr_pHts21a;   //!
   TBranch        *b_rjrIsr_pHts21b;   //!
   TBranch        *b_rjrIsr_pHts22;   //!
   TBranch        *b_rjrIsr_pHts40;   //!
   TBranch        *b_rjrIsr_pHts41;   //!
   TBranch        *b_rjrIsr_pHtxa10;   //!
   TBranch        *b_rjrIsr_pHtxa11;   //!
   TBranch        *b_rjrIsr_pHtxa20;   //!
   TBranch        *b_rjrIsr_pHtxa21;   //!
   TBranch        *b_rjrIsr_pHtxb10;   //!
   TBranch        *b_rjrIsr_pHtxb11;   //!
   TBranch        *b_rjrIsr_pHtxb20;   //!
   TBranch        *b_rjrIsr_pHtxb21;   //!
   TBranch        *b_rjrIsr_pHxa10;   //!
   TBranch        *b_rjrIsr_pHxa11;   //!
   TBranch        *b_rjrIsr_pHxa20;   //!
   TBranch        *b_rjrIsr_pHxa21;   //!
   TBranch        *b_rjrIsr_pHxb10;   //!
   TBranch        *b_rjrIsr_pHxb11;   //!
   TBranch        *b_rjrIsr_pHxb20;   //!
   TBranch        *b_rjrIsr_pHxb21;   //!
   TBranch        *b_rjrIsr_sJ1a_e;   //!
   TBranch        *b_rjrIsr_sJ1a_px;   //!
   TBranch        *b_rjrIsr_sJ1a_py;   //!
   TBranch        *b_rjrIsr_sJ1a_pz;   //!
   TBranch        *b_rjrIsr_sJ1b_e;   //!
   TBranch        *b_rjrIsr_sJ1b_px;   //!
   TBranch        *b_rjrIsr_sJ1b_py;   //!
   TBranch        *b_rjrIsr_sJ1b_pz;   //!
   TBranch        *b_rjrIsr_sJ2a_e;   //!
   TBranch        *b_rjrIsr_sJ2a_px;   //!
   TBranch        *b_rjrIsr_sJ2a_py;   //!
   TBranch        *b_rjrIsr_sJ2a_pz;   //!
   TBranch        *b_rjrIsr_sJ2b_e;   //!
   TBranch        *b_rjrIsr_sJ2b_px;   //!
   TBranch        *b_rjrIsr_sJ2b_py;   //!
   TBranch        *b_rjrIsr_sJ2b_pz;   //!
   TBranch        *b_eleHasSVMatch;   //!
   TBranch        *b_eleIsLoose;   //!
   TBranch        *b_elePhoMinDr;   //!
   TBranch        *b_epDrMatch;   //!
   TBranch        *b_epSeedIdMatch;   //!
   TBranch        *b_epSeedIdUnique;   //!
   TBranch        *b_nEleSVMatched;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nSelElectrons;   //!
   TBranch        *b_nSelIsoElectrons;   //!
   TBranch        *b_allJetEnergy;   //!
   TBranch        *b_allJetEta;   //!
   TBranch        *b_allJetMass;   //!
   TBranch        *b_allJetNPho;   //!
   TBranch        *b_allJetPhi;   //!
   TBranch        *b_allJetPixNPho;   //!
   TBranch        *b_allJetPixWTime;   //!
   TBranch        *b_allJetPixWTimeRes;   //!
   TBranch        *b_allJetPixWTimeSig;   //!
   TBranch        *b_allJetPt;   //!
   TBranch        *b_allJetSCIndexs;   //!
   TBranch        *b_allJetWTime;   //!
   TBranch        *b_allJetWTimeRes;   //!
   TBranch        *b_allJetWTimeSig;   //!
   TBranch        *b_allJetsWithWPixTime;   //!
   TBranch        *b_allJetsWithWTime;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_nSelJets;   //!
   TBranch        *b_nSelJetsWithWPixTime;   //!
   TBranch        *b_nSelJetsWithWTime;   //!
   TBranch        *b_pv_dDiJetPEnergy;   //!
   TBranch        *b_pv_dDiJetTime;   //!
   TBranch        *b_pv_dGJMinDr;   //!
   TBranch        *b_pv_dGJPEnergy;   //!
   TBranch        *b_pv_dGJTime;   //!
   TBranch        *b_pv_diJet1Time;   //!
   TBranch        *b_pv_diJet2Time;   //!
   TBranch        *b_pv_gjGammaTime;   //!
   TBranch        *b_pv_gjJetTime;   //!
   TBranch        *b_pv_gjPhoIndex;   //!
   TBranch        *b_pv_haspvtime;   //!
   TBranch        *b_pv_wpenergy;   //!
   TBranch        *b_pv_wphartres;   //!
   TBranch        *b_pv_wtime;   //!
   TBranch        *b_pv_wtimeres;   //!
   TBranch        *b_selGenJetDpt;   //!
   TBranch        *b_selGenJetEnergy;   //!
   TBranch        *b_selGenJetImpAng;   //!
   TBranch        *b_selGenJetLlpTime;   //!
   TBranch        *b_selGenJetPt;   //!
   TBranch        *b_selGenJetTime;   //!
   TBranch        *b_selGenJetTof;   //!
   TBranch        *b_selGenJetdr;   //!
   TBranch        *b_selGenJeteta;   //!
   TBranch        *b_selJetArea;   //!
   TBranch        *b_selJetChEmEF;   //!
   TBranch        *b_selJetChHM;   //!
   TBranch        *b_selJetEnergy;   //!
   TBranch        *b_selJetEta;   //!
   TBranch        *b_selJetGenLlpId;   //!
   TBranch        *b_selJetLlpDp;   //!
   TBranch        *b_selJetLlpDr;   //!
   TBranch        *b_selJetMass;   //!
   TBranch        *b_selJetMuEF;   //!
   TBranch        *b_selJetNeEmEF;   //!
   TBranch        *b_selJetNeHEF;   //!
   TBranch        *b_selJetNeHM;   //!
   TBranch        *b_selJetPhi;   //!
   TBranch        *b_selJetPixWTime;   //!
   TBranch        *b_selJetPixWTimeRes;   //!
   TBranch        *b_selJetPixWTimeSig;   //!
   TBranch        *b_selJetPt;   //!
   TBranch        *b_selJetQrkLlpId;   //!
   TBranch        *b_selJetQuality;   //!
   TBranch        *b_selJetTime;   //!
   TBranch        *b_selJetWTime;   //!
   TBranch        *b_selJetWTimeRes;   //!
   TBranch        *b_selJetWTimeSig;   //!
   TBranch        *b_selJetchHEF;   //!
   TBranch        *b_muonIsLoose;   //!
   TBranch        *b_muonIsMedium;   //!
   TBranch        *b_nLooseMuons;   //!
   TBranch        *b_nMuonSVMatched;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_HadronicSV_cosTheta;   //!
   TBranch        *b_HadronicSV_decayAngle;   //!
   TBranch        *b_HadronicSV_dxy;   //!
   TBranch        *b_HadronicSV_dxySig;   //!
   TBranch        *b_HadronicSV_mass;   //!
   TBranch        *b_HadronicSV_massOverNtracks;   //!
   TBranch        *b_HadronicSV_matchRatio;   //!
   TBranch        *b_HadronicSV_nTotal;   //!
   TBranch        *b_HadronicSV_nTracks;   //!
   TBranch        *b_HadronicSV_pOverE;   //!
   TBranch        *b_HadronicSV_passMaterialVeto;   //!
   TBranch        *b_HadronicSV_tightZWindowSelection;   //!
   TBranch        *b_HadronicSV_x;   //!
   TBranch        *b_HadronicSV_y;   //!
   TBranch        *b_HadronicSV_z;   //!
   TBranch        *b_LeptonicSV_cosTheta;   //!
   TBranch        *b_LeptonicSV_decayAngle;   //!
   TBranch        *b_LeptonicSV_dxy;   //!
   TBranch        *b_LeptonicSV_dxySig;   //!
   TBranch        *b_LeptonicSV_electronIndex;   //!
   TBranch        *b_LeptonicSV_isGold;   //!
   TBranch        *b_LeptonicSV_mass;   //!
   TBranch        *b_LeptonicSV_muonIndex;   //!
   TBranch        *b_LeptonicSV_nTotal;   //!
   TBranch        *b_LeptonicSV_nTracks;   //!
   TBranch        *b_LeptonicSV_pOverE;   //!
   TBranch        *b_LeptonicSV_passElectronID;   //!
   TBranch        *b_LeptonicSV_passMaterialVeto;   //!
   TBranch        *b_LeptonicSV_passMuonID;   //!
   TBranch        *b_LeptonicSV_tightZWindowSelection;   //!
   TBranch        *b_LeptonicSV_x;   //!
   TBranch        *b_LeptonicSV_y;   //!
   TBranch        *b_LeptonicSV_z;   //!
   TBranch        *b_SV_nElectron;   //!
   TBranch        *b_SV_nElectronPassZWindow;   //!
   TBranch        *b_SV_nHadCmpDxyGeq2;   //!
   TBranch        *b_SV_nHadCmpDxySigGeq1000;   //!
   TBranch        *b_SV_nHadStdDxyGeq2;   //!
   TBranch        *b_SV_nHadStdDxySigGeq1000;   //!
   TBranch        *b_SV_nHadronPassZWindow;   //!
   TBranch        *b_SV_nHadronic;   //!
   TBranch        *b_SV_nLepCmpDxyGeq2;   //!
   TBranch        *b_SV_nLepCmpDxySigGeq1000;   //!
   TBranch        *b_SV_nLepStdDxyGeq2;   //!
   TBranch        *b_SV_nLepStdDxySigGeq1000;   //!
   TBranch        *b_SV_nLeptonPassZWindow;   //!
   TBranch        *b_SV_nLeptonic;   //!
   TBranch        *b_SV_nMuon;   //!
   TBranch        *b_SV_nMuonPassZWindow;   //!
   TBranch        *b_passGe1HadAndLepSelectionCRCmp;   //!
   TBranch        *b_passGe1HadAndLepSelectionCRLoose;   //!
   TBranch        *b_passGe1HadAndLepSelectionCRTight;   //!
   TBranch        *b_passGe1HadAndLepSelectionSRCmp;   //!
   TBranch        *b_passGe1HadAndLepSelectionSRLoose;   //!
   TBranch        *b_passGe1HadAndLepSelectionSRTight;   //!
   TBranch        *b_passNHad1SelectionCRCmp;   //!
   TBranch        *b_passNHad1SelectionCRLoose;   //!
   TBranch        *b_passNHad1SelectionCRTight;   //!
   TBranch        *b_passNHad1SelectionSRCmp;   //!
   TBranch        *b_passNHad1SelectionSRLoose;   //!
   TBranch        *b_passNHad1SelectionSRTight;   //!
   TBranch        *b_passNHadGe2SelectionCRCmp;   //!
   TBranch        *b_passNHadGe2SelectionCRLoose;   //!
   TBranch        *b_passNHadGe2SelectionCRTight;   //!
   TBranch        *b_passNHadGe2SelectionSRCmp;   //!
   TBranch        *b_passNHadGe2SelectionSRLoose;   //!
   TBranch        *b_passNHadGe2SelectionSRTight;   //!
   TBranch        *b_passNLep1SelectionCRCmp;   //!
   TBranch        *b_passNLep1SelectionCRLoose;   //!
   TBranch        *b_passNLep1SelectionCRTight;   //!
   TBranch        *b_passNLep1SelectionSRCmp;   //!
   TBranch        *b_passNLep1SelectionSRLoose;   //!
   TBranch        *b_passNLep1SelectionSRTight;   //!
   TBranch        *b_passNLepGe2SelectionCRCmp;   //!
   TBranch        *b_passNLepGe2SelectionCRLoose;   //!
   TBranch        *b_passNLepGe2SelectionCRTight;   //!
   TBranch        *b_passNLepGe2SelectionSRCmp;   //!
   TBranch        *b_passNLepGe2SelectionSRLoose;   //!
   TBranch        *b_passNLepGe2SelectionSRTight;   //!

   (TTree *tree=0);
   virtual ~();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef _cxx
::(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SMS_v34p1_FULLMINI_SVHPM100_mN2-1300_mGl-2300_mN2-1300_mN1-1000_ct0p1_rjrskim.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SMS_v34p1_FULLMINI_SVHPM100_mN2-1300_mGl-2300_mN2-1300_mN1-1000_ct0p1_rjrskim.root");
      }
      f->GetObject("kuSkimTree",tree);

   }
   Init(tree);
}

::~()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ::LoadTree(Long64_t entry)
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

void ::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   DataSetKey = 0;
   genCharge = 0;
   genGrvtinoMass = 0;
   genLLPMass = 0;
   genLSPMass = 0;
   genMass = 0;
   genMomCTau = 0;
   genPartEnergy = 0;
   genPartEta = 0;
   genPartMomIdx = 0;
   genPartPdgId = 0;
   genPartPhi = 0;
   genPartPt = 0;
   genPartSusId = 0;
   genPx = 0;
   genPy = 0;
   genPz = 0;
   genSGMass = 0;
   genSQMass = 0;
   genStatus = 0;
   genVx = 0;
   genVy = 0;
   genVz = 0;
   genXMomCTau = 0;
   barrel_photon_baseline = 0;
   barrel_photon_energy = 0;
   barrel_photon_erhsecdr4 = 0;
   barrel_photon_eta = 0;
   barrel_photon_htoem = 0;
   barrel_photon_isoANNScore = 0;
   barrel_photon_phi = 0;
   barrel_photon_pt = 0;
   barrel_photon_tspscdr4 = 0;
   baseLinePhoton_ClstrRn = 0;
   baseLinePhoton_CorEnergy = 0;
   baseLinePhoton_CorPt = 0;
   baseLinePhoton_CovEtaEta = 0;
   baseLinePhoton_CovEtaPhi = 0;
   baseLinePhoton_CovPhiPhi = 0;
   baseLinePhoton_DiJetsCR = 0;
   baseLinePhoton_EcalPFClusterIso = 0;
   baseLinePhoton_EcalRHSumEtConeDR04 = 0;
   baseLinePhoton_Energy = 0;
   baseLinePhoton_Eta = 0;
   baseLinePhoton_EtaWidth = 0;
   baseLinePhoton_GJetsCR = 0;
   baseLinePhoton_GenIdx = 0;
   baseLinePhoton_GenLabTime = 0;
   baseLinePhoton_GenLabTimeSig = 0;
   baseLinePhoton_GenPt = 0;
   baseLinePhoton_GenSigMomEnergy = 0;
   baseLinePhoton_GenSigMomEta = 0;
   baseLinePhoton_GenSigMomMass = 0;
   baseLinePhoton_GenSigMomPhi = 0;
   baseLinePhoton_GenSigMomPt = 0;
   baseLinePhoton_GenSigMomPx = 0;
   baseLinePhoton_GenSigMomPy = 0;
   baseLinePhoton_GenSigMomPz = 0;
   baseLinePhoton_GenSigMomVx = 0;
   baseLinePhoton_GenSigMomVy = 0;
   baseLinePhoton_GenSigMomVz = 0;
   baseLinePhoton_GenTime = 0;
   baseLinePhoton_GenTimeSig = 0;
   baseLinePhoton_HadOverEM = 0;
   baseLinePhoton_HadTowOverEM = 0;
   baseLinePhoton_HasConversionTracks = 0;
   baseLinePhoton_HcalPFClusterIso = 0;
   baseLinePhoton_HcalTowerSumEtBcConeDR04 = 0;
   baseLinePhoton_HcalTowerSumEtConeDR04 = 0;
   baseLinePhoton_Hoe_PUcorr = 0;
   baseLinePhoton_LSCross = 0;
   baseLinePhoton_LTRes = 0;
   baseLinePhoton_LTime = 0;
   baseLinePhoton_LTimeSig = 0;
   baseLinePhoton_LWeird = 0;
   baseLinePhoton_MBetaEql = 0;
   baseLinePhoton_MBetaPmt = 0;
   baseLinePhoton_NTrkHollowConeDR04 = 0;
   baseLinePhoton_NTrkSolidConeDR04 = 0;
   baseLinePhoton_Nrh = 0;
   baseLinePhoton_OOT = 0;
   baseLinePhoton_PVtoRH = 0;
   baseLinePhoton_PfChargedIso = 0;
   baseLinePhoton_PfChargedIsoPFPV = 0;
   baseLinePhoton_PfChargedIsoWorstVtx = 0;
   baseLinePhoton_PfPhoIso03 = 0;
   baseLinePhoton_PfRelIso03_all_quadratic = 0;
   baseLinePhoton_PfRelIso03_chg_quadratic = 0;
   baseLinePhoton_Phi = 0;
   baseLinePhoton_PhiWidth = 0;
   baseLinePhoton_PhoIsoDr = 0;
   baseLinePhoton_PixelSeed = 0;
   baseLinePhoton_Pt = 0;
   baseLinePhoton_Quality = 0;
   baseLinePhoton_R9 = 0;
   baseLinePhoton_RawTime = 0;
   baseLinePhoton_S4 = 0;
   baseLinePhoton_SAlp = 0;
   baseLinePhoton_SCx = 0;
   baseLinePhoton_SCy = 0;
   baseLinePhoton_SCz = 0;
   baseLinePhoton_SMaj = 0;
   baseLinePhoton_SMin = 0;
   baseLinePhoton_SSCross = 0;
   baseLinePhoton_STRes = 0;
   baseLinePhoton_STime = 0;
   baseLinePhoton_STimeSig = 0;
   baseLinePhoton_SWeird = 0;
   baseLinePhoton_ShasGS = 0;
   baseLinePhoton_Sieie = 0;
   baseLinePhoton_Sieip = 0;
   baseLinePhoton_SigmaIEtaIEta = 0;
   baseLinePhoton_Sipip = 0;
   baseLinePhoton_SusyId = 0;
   baseLinePhoton_Time = 0;
   baseLinePhoton_TrkSumPtHollowConeDR03 = 0;
   baseLinePhoton_TrkSumPtHollowConeDR04 = 0;
   baseLinePhoton_TrkSumPtSolidConeDR04 = 0;
   baseLinePhoton_WTRes = 0;
   baseLinePhoton_WTResOld = 0;
   baseLinePhoton_WTime = 0;
   baseLinePhoton_WTimeOld = 0;
   baseLinePhoton_WTimeSig = 0;
   baseLinePhoton_WTimeSigOld = 0;
   baseLinePhoton_beamHaloCNNScore = 0;
   baseLinePhoton_beamHaloCR = 0;
   baseLinePhoton_beamHaloID = 0;
   baseLinePhoton_isoANNScore = 0;
   baseLinePhoton_looseNotTightIso1ID = 0;
   baseLinePhoton_looseNotTightIso2ID = 0;
   baseLinePhoton_looseNotTightIsoID = 0;
   baseLinePhoton_minJetDr = 0;
   baseLinePhoton_nonIsoANNScore = 0;
   baseLinePhoton_notBeamHaloID = 0;
   baseLinePhoton_physBkgCNNScore = 0;
   baseLinePhoton_spikeCR = 0;
   baseLinePhoton_tightIsoID = 0;
   endcap_photon_baseline = 0;
   endcap_photon_energy = 0;
   endcap_photon_erhsecdr4 = 0;
   endcap_photon_eta = 0;
   endcap_photon_htoem = 0;
   endcap_photon_isoANNScore = 0;
   endcap_photon_phi = 0;
   endcap_photon_pt = 0;
   endcap_photon_tspscdr4 = 0;
   photon_E = 0;
   photon_Eta = 0;
   photon_LeadTime = 0;
   photon_LeadTimeSig = 0;
   photon_Phi = 0;
   photon_PixSeed = 0;
   photon_Pt = 0;
   photon_WTime = 0;
   photon_WTimeSig = 0;
   photon_baseline = 0;
   photon_beamHaloCNNScore = 0;
   photon_isoANNScore = 0;
   photon_minJetDr = 0;
   rjrABSide = 0;
   rjr_Rdiff = 0;
   rjrCleaningVeto0 = 0;
   rjrCleaningVeto1 = 0;
   rjrCleaningVeto2 = 0;
   rjrDPhiSI = 0;
   rjrDphiMETV = 0;
   rjrEVa = 0;
   rjrEVb = 0;
   rjrLeadPhoLocation = 0;
   rjrMET = 0;
   rjrN2Px = 0;
   rjrN2Py = 0;
   rjrN2Pz = 0;
   rjrNJets = 0;
   rjrNJetsJa = 0;
   rjrNJetsJb = 0;
   rjrNPhotons = 0;
   rjrNVisObjects = 0;
   rjrPS = 0;
   rjrPTS = 0;
   rjrPVa = 0;
   rjrPVb = 0;
   rjrPVlab = 0;
   rjrPZS = 0;
   rjrSCosA = 0;
   rjrSMass = 0;
   rjrSdphiBV = 0;
   rjrSdphiDA = 0;
   rjrSubPhoLocation = 0;
   rjrType = 0;
   rjrVisJetSide = 0;
   rjrVisPhoSide = 0;
   rjrX1aCosA = 0;
   rjrX1aMass = 0;
   rjrX1a_Eta = 0;
   rjrX1a_Phi = 0;
   rjrX1a_Pt = 0;
   rjrX1bCosA = 0;
   rjrX1bMass = 0;
   rjrX1b_Eta = 0;
   rjrX1b_Phi = 0;
   rjrX1b_Pt = 0;
   rjrX2aCosA = 0;
   rjrX2aMass = 0;
   rjrX2aPs = 0;
   rjrX2bCosA = 0;
   rjrX2bMass = 0;
   rjrX2bPs = 0;
   rjr_Mr = 0;
   rjr_Ms = 0;
   rjr_Mva = 0;
   rjr_Mvb = 0;
   rjr_R = 0;
   rjr_Rm = 0;
   rjr_Rs = 0;
   rjr_Rv = 0;
   rjr_Rx = 0;
   rjr_Rxa = 0;
   rjr_Rxb = 0;
   rjr_p1Pts = 0;
   rjr_p1Ptxa = 0;
   rjr_p1Ptxb = 0;
   rjr_p1Pxa = 0;
   rjr_p1Pxb = 0;
   rjr_p2Pts = 0;
   rjr_p2Ptxa = 0;
   rjr_p2Ptxb = 0;
   rjr_p2Pxa = 0;
   rjr_p2Pxb = 0;
   rjr_pHs11 = 0;
   rjr_pHs11a = 0;
   rjr_pHs11b = 0;
   rjr_pHs20 = 0;
   rjr_pHs21 = 0;
   rjr_pHs21a = 0;
   rjr_pHs21b = 0;
   rjr_pHs22 = 0;
   rjr_pHs40 = 0;
   rjr_pHs41 = 0;
   rjr_pHts11 = 0;
   rjr_pHts11a = 0;
   rjr_pHts11b = 0;
   rjr_pHts20 = 0;
   rjr_pHts21 = 0;
   rjr_pHts21a = 0;
   rjr_pHts21b = 0;
   rjr_pHts22 = 0;
   rjr_pHts40 = 0;
   rjr_pHts41 = 0;
   rjr_pHtxa10 = 0;
   rjr_pHtxa11 = 0;
   rjr_pHtxa20 = 0;
   rjr_pHtxa21 = 0;
   rjr_pHtxb10 = 0;
   rjr_pHtxb11 = 0;
   rjr_pHtxb20 = 0;
   rjr_pHtxb21 = 0;
   rjr_pHxa10 = 0;
   rjr_pHxa11 = 0;
   rjr_pHxa20 = 0;
   rjr_pHxa21 = 0;
   rjr_pHxb10 = 0;
   rjr_pHxb11 = 0;
   rjr_pHxb20 = 0;
   rjr_pHxb21 = 0;
   rjr_sJ1a_e = 0;
   rjr_sJ1a_px = 0;
   rjr_sJ1a_py = 0;
   rjr_sJ1a_pz = 0;
   rjr_sJ1b_e = 0;
   rjr_sJ1b_px = 0;
   rjr_sJ1b_py = 0;
   rjr_sJ1b_pz = 0;
   rjr_sJ2a_e = 0;
   rjr_sJ2a_px = 0;
   rjr_sJ2a_py = 0;
   rjr_sJ2a_pz = 0;
   rjr_sJ2b_e = 0;
   rjr_sJ2b_px = 0;
   rjr_sJ2b_py = 0;
   rjr_sJ2b_pz = 0;
   eleHasSVMatch = 0;
   eleIsLoose = 0;
   elePhoMinDr = 0;
   epDrMatch = 0;
   epSeedIdMatch = 0;
   epSeedIdUnique = 0;
   allJetEnergy = 0;
   allJetEta = 0;
   allJetMass = 0;
   allJetNPho = 0;
   allJetPhi = 0;
   allJetPixNPho = 0;
   allJetPixWTime = 0;
   allJetPixWTimeRes = 0;
   allJetPixWTimeSig = 0;
   allJetPt = 0;
   allJetSCIndexs = 0;
   allJetWTime = 0;
   allJetWTimeRes = 0;
   allJetWTimeSig = 0;
   allJetsWithWPixTime = 0;
   allJetsWithWTime = 0;
   nSelJetsWithWPixTime = 0;
   nSelJetsWithWTime = 0;
   pv_dDiJetPEnergy = 0;
   pv_dDiJetTime = 0;
   pv_dGJMinDr = 0;
   pv_dGJPEnergy = 0;
   pv_dGJTime = 0;
   pv_diJet1Time = 0;
   pv_diJet2Time = 0;
   pv_gjGammaTime = 0;
   pv_gjJetTime = 0;
   pv_gjPhoIndex = 0;
   pv_haspvtime = 0;
   pv_wpenergy = 0;
   pv_wphartres = 0;
   pv_wtime = 0;
   pv_wtimeres = 0;
   selGenJetDpt = 0;
   selGenJetEnergy = 0;
   selGenJetImpAng = 0;
   selGenJetLlpTime = 0;
   selGenJetPt = 0;
   selGenJetTime = 0;
   selGenJetTof = 0;
   selGenJetdr = 0;
   selGenJeteta = 0;
   selJetArea = 0;
   selJetChEmEF = 0;
   selJetChHM = 0;
   selJetEnergy = 0;
   selJetEta = 0;
   selJetGenLlpId = 0;
   selJetLlpDp = 0;
   selJetLlpDr = 0;
   selJetMass = 0;
   selJetMuEF = 0;
   selJetNeEmEF = 0;
   selJetNeHEF = 0;
   selJetNeHM = 0;
   selJetPhi = 0;
   selJetPixWTime = 0;
   selJetPixWTimeRes = 0;
   selJetPixWTimeSig = 0;
   selJetPt = 0;
   selJetQrkLlpId = 0;
   selJetQuality = 0;
   selJetTime = 0;
   selJetWTime = 0;
   selJetWTimeRes = 0;
   selJetWTimeSig = 0;
   selJetchHEF = 0;
   muonIsLoose = 0;
   muonIsMedium = 0;
   HadronicSV_cosTheta = 0;
   HadronicSV_decayAngle = 0;
   HadronicSV_dxy = 0;
   HadronicSV_dxySig = 0;
   HadronicSV_mass = 0;
   HadronicSV_massOverNtracks = 0;
   HadronicSV_matchRatio = 0;
   HadronicSV_nTracks = 0;
   HadronicSV_pOverE = 0;
   HadronicSV_passMaterialVeto = 0;
   HadronicSV_tightZWindowSelection = 0;
   HadronicSV_x = 0;
   HadronicSV_y = 0;
   HadronicSV_z = 0;
   LeptonicSV_cosTheta = 0;
   LeptonicSV_decayAngle = 0;
   LeptonicSV_dxy = 0;
   LeptonicSV_dxySig = 0;
   LeptonicSV_electronIndex = 0;
   LeptonicSV_isGold = 0;
   LeptonicSV_mass = 0;
   LeptonicSV_muonIndex = 0;
   LeptonicSV_nTracks = 0;
   LeptonicSV_pOverE = 0;
   LeptonicSV_passElectronID = 0;
   LeptonicSV_passMaterialVeto = 0;
   LeptonicSV_passMuonID = 0;
   LeptonicSV_tightZWindowSelection = 0;
   LeptonicSV_x = 0;
   LeptonicSV_y = 0;
   LeptonicSV_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_MetFilters", &Flag_MetFilters, &b_Flag_MetFilters);
   fChain->SetBranchAddress("Flag_eLHemVeto", &Flag_eLHemVeto, &b_Flag_eLHemVeto);
   fChain->SetBranchAddress("Flag_eMHemVeto", &Flag_eMHemVeto, &b_Flag_eMHemVeto);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_hemRun", &Flag_hemRun, &b_Flag_hemRun);
   fChain->SetBranchAddress("Flag_hemVeto", &Flag_hemVeto, &b_Flag_hemVeto);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("Flag_jLHemVeto", &Flag_jLHemVeto, &b_Flag_jLHemVeto);
   fChain->SetBranchAddress("Flag_jMHemVeto", &Flag_jMHemVeto, &b_Flag_jMHemVeto);
   fChain->SetBranchAddress("Flag_mLHemVeto", &Flag_mLHemVeto, &b_Flag_mLHemVeto);
   fChain->SetBranchAddress("Flag_mMHemVeto", &Flag_mMHemVeto, &b_Flag_mMHemVeto);
   fChain->SetBranchAddress("Flag_pLHemVeto", &Flag_pLHemVeto, &b_Flag_pLHemVeto);
   fChain->SetBranchAddress("Flag_pMHemVeto", &Flag_pMHemVeto, &b_Flag_pMHemVeto);
   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("Trigger_PFMET120_PFMHT120_IDTight", &Trigger_PFMET120_PFMHT120_IDTight, &b_Trigger_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("Trigger_PFMET120_PFMHT120_IDTight_PFHT60", &Trigger_PFMET120_PFMHT120_IDTight_PFHT60, &b_Trigger_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight", &Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("DataSetKey", &DataSetKey, &b_DataSetKey);
   fChain->SetBranchAddress("evtFillWgt", &evtFillWgt, &b_evtFillWgt);
   fChain->SetBranchAddress("evtGenWgt", &evtGenWgt, &b_evtGenWgt);
   fChain->SetBranchAddress("evtXSection", &evtXSection, &b_evtXSection);
   fChain->SetBranchAddress("selCMet", &selCMet, &b_selCMet);
   fChain->SetBranchAddress("selCMetPx", &selCMetPx, &b_selCMetPx);
   fChain->SetBranchAddress("selCMetPy", &selCMetPy, &b_selCMetPy);
   fChain->SetBranchAddress("selMet", &selMet, &b_selMet);
   fChain->SetBranchAddress("selMetPx", &selMetPx, &b_selMetPx);
   fChain->SetBranchAddress("selMetPy", &selMetPy, &b_selMetPy);
   fChain->SetBranchAddress("Trk_nIsoTrack", &Trk_nIsoTrack, &b_Trk_nIsoTrack);
   fChain->SetBranchAddress("Evt_isGG", &Evt_isGG, &b_Evt_isGG);
   fChain->SetBranchAddress("Evt_isGZ", &Evt_isGZ, &b_Evt_isGZ);
   fChain->SetBranchAddress("Evt_isZZ", &Evt_isZZ, &b_Evt_isZZ);
   fChain->SetBranchAddress("Evt_nXs", &Evt_nXs, &b_Evt_nXs);
   fChain->SetBranchAddress("Xa_Displacment", &Xa_Displacment, &b_Xa_Displacment);
   fChain->SetBranchAddress("Xa_beta", &Xa_beta, &b_Xa_beta);
   fChain->SetBranchAddress("Xa_ctau", &Xa_ctau, &b_Xa_ctau);
   fChain->SetBranchAddress("Xa_energy", &Xa_energy, &b_Xa_energy);
   fChain->SetBranchAddress("Xa_mass", &Xa_mass, &b_Xa_mass);
   fChain->SetBranchAddress("Xa_p", &Xa_p, &b_Xa_p);
   fChain->SetBranchAddress("Xa_pdgId", &Xa_pdgId, &b_Xa_pdgId);
   fChain->SetBranchAddress("Xa_phi", &Xa_phi, &b_Xa_phi);
   fChain->SetBranchAddress("Xa_pt", &Xa_pt, &b_Xa_pt);
   fChain->SetBranchAddress("Xa_vx", &Xa_vx, &b_Xa_vx);
   fChain->SetBranchAddress("Xa_vy", &Xa_vy, &b_Xa_vy);
   fChain->SetBranchAddress("Xa_vz", &Xa_vz, &b_Xa_vz);
   fChain->SetBranchAddress("Xb_Displacment", &Xb_Displacment, &b_Xb_Displacment);
   fChain->SetBranchAddress("Xb_beta", &Xb_beta, &b_Xb_beta);
   fChain->SetBranchAddress("Xb_ctau", &Xb_ctau, &b_Xb_ctau);
   fChain->SetBranchAddress("Xb_energy", &Xb_energy, &b_Xb_energy);
   fChain->SetBranchAddress("Xb_mass", &Xb_mass, &b_Xb_mass);
   fChain->SetBranchAddress("Xb_p", &Xb_p, &b_Xb_p);
   fChain->SetBranchAddress("Xb_pdgId", &Xb_pdgId, &b_Xb_pdgId);
   fChain->SetBranchAddress("Xb_phi", &Xb_phi, &b_Xb_phi);
   fChain->SetBranchAddress("Xb_pt", &Xb_pt, &b_Xb_pt);
   fChain->SetBranchAddress("Xb_vx", &Xb_vx, &b_Xb_vx);
   fChain->SetBranchAddress("Xb_vy", &Xb_vy, &b_Xb_vy);
   fChain->SetBranchAddress("Xb_vz", &Xb_vz, &b_Xb_vz);
   fChain->SetBranchAddress("genCharge", &genCharge, &b_genCharge);
   fChain->SetBranchAddress("genGrvtinoMass", &genGrvtinoMass, &b_genGrvtinoMass);
   fChain->SetBranchAddress("genLLPMass", &genLLPMass, &b_genLLPMass);
   fChain->SetBranchAddress("genLSPMass", &genLSPMass, &b_genLSPMass);
   fChain->SetBranchAddress("genMass", &genMass, &b_genMass);
   fChain->SetBranchAddress("genMomCTau", &genMomCTau, &b_genMomCTau);
   fChain->SetBranchAddress("genPartEnergy", &genPartEnergy, &b_genPartEnergy);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartMomIdx", &genPartMomIdx, &b_genPartMomIdx);
   fChain->SetBranchAddress("genPartPdgId", &genPartPdgId, &b_genPartPdgId);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartSusId", &genPartSusId, &b_genPartSusId);
   fChain->SetBranchAddress("genPx", &genPx, &b_genPx);
   fChain->SetBranchAddress("genPy", &genPy, &b_genPy);
   fChain->SetBranchAddress("genPz", &genPz, &b_genPz);
   fChain->SetBranchAddress("genSGMass", &genSGMass, &b_genSGMass);
   fChain->SetBranchAddress("genSQMass", &genSQMass, &b_genSQMass);
   fChain->SetBranchAddress("genSTFlagQP", &genSTFlagQP, &b_genSTFlagQP);
   fChain->SetBranchAddress("genSTFlagQQP", &genSTFlagQQP, &b_genSTFlagQQP);
   fChain->SetBranchAddress("genSTFlagQQZ", &genSTFlagQQZ, &b_genSTFlagQQZ);
   fChain->SetBranchAddress("genSTFlagQQZL", &genSTFlagQQZL, &b_genSTFlagQQZL);
   fChain->SetBranchAddress("genSTFlagQQZLL", &genSTFlagQQZLL, &b_genSTFlagQQZLL);
   fChain->SetBranchAddress("genSigType", &genSigType, &b_genSigType);
   fChain->SetBranchAddress("genStatus", &genStatus, &b_genStatus);
   fChain->SetBranchAddress("genVx", &genVx, &b_genVx);
   fChain->SetBranchAddress("genVy", &genVy, &b_genVy);
   fChain->SetBranchAddress("genVz", &genVz, &b_genVz);
   fChain->SetBranchAddress("genXMomCTau", &genXMomCTau, &b_genXMomCTau);
   fChain->SetBranchAddress("barrel_photon_baseline", &barrel_photon_baseline, &b_barrel_photon_baseline);
   fChain->SetBranchAddress("barrel_photon_energy", &barrel_photon_energy, &b_barrel_photon_energy);
   fChain->SetBranchAddress("barrel_photon_erhsecdr4", &barrel_photon_erhsecdr4, &b_barrel_photon_erhsecdr4);
   fChain->SetBranchAddress("barrel_photon_eta", &barrel_photon_eta, &b_barrel_photon_eta);
   fChain->SetBranchAddress("barrel_photon_htoem", &barrel_photon_htoem, &b_barrel_photon_htoem);
   fChain->SetBranchAddress("barrel_photon_isoANNScore", &barrel_photon_isoANNScore, &b_barrel_photon_isoANNScore);
   fChain->SetBranchAddress("barrel_photon_phi", &barrel_photon_phi, &b_barrel_photon_phi);
   fChain->SetBranchAddress("barrel_photon_pt", &barrel_photon_pt, &b_barrel_photon_pt);
   fChain->SetBranchAddress("barrel_photon_tspscdr4", &barrel_photon_tspscdr4, &b_barrel_photon_tspscdr4);
   fChain->SetBranchAddress("baseLinePhoton_ClstrRn", &baseLinePhoton_ClstrRn, &b_baseLinePhoton_ClstrRn);
   fChain->SetBranchAddress("baseLinePhoton_CorEnergy", &baseLinePhoton_CorEnergy, &b_baseLinePhoton_CorEnergy);
   fChain->SetBranchAddress("baseLinePhoton_CorPt", &baseLinePhoton_CorPt, &b_baseLinePhoton_CorPt);
   fChain->SetBranchAddress("baseLinePhoton_CovEtaEta", &baseLinePhoton_CovEtaEta, &b_baseLinePhoton_CovEtaEta);
   fChain->SetBranchAddress("baseLinePhoton_CovEtaPhi", &baseLinePhoton_CovEtaPhi, &b_baseLinePhoton_CovEtaPhi);
   fChain->SetBranchAddress("baseLinePhoton_CovPhiPhi", &baseLinePhoton_CovPhiPhi, &b_baseLinePhoton_CovPhiPhi);
   fChain->SetBranchAddress("baseLinePhoton_DiJetsCR", &baseLinePhoton_DiJetsCR, &b_baseLinePhoton_DiJetsCR);
   fChain->SetBranchAddress("baseLinePhoton_EcalPFClusterIso", &baseLinePhoton_EcalPFClusterIso, &b_baseLinePhoton_EcalPFClusterIso);
   fChain->SetBranchAddress("baseLinePhoton_EcalRHSumEtConeDR04", &baseLinePhoton_EcalRHSumEtConeDR04, &b_baseLinePhoton_EcalRHSumEtConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_EleVeto", &baseLinePhoton_EleVeto, &b_baseLinePhoton_EleVeto);
   fChain->SetBranchAddress("baseLinePhoton_Energy", &baseLinePhoton_Energy, &b_baseLinePhoton_Energy);
   fChain->SetBranchAddress("baseLinePhoton_Eta", &baseLinePhoton_Eta, &b_baseLinePhoton_Eta);
   fChain->SetBranchAddress("baseLinePhoton_EtaWidth", &baseLinePhoton_EtaWidth, &b_baseLinePhoton_EtaWidth);
   fChain->SetBranchAddress("baseLinePhoton_GJetsCR", &baseLinePhoton_GJetsCR, &b_baseLinePhoton_GJetsCR);
   fChain->SetBranchAddress("baseLinePhoton_GenIdx", &baseLinePhoton_GenIdx, &b_baseLinePhoton_GenIdx);
   fChain->SetBranchAddress("baseLinePhoton_GenLabTime", &baseLinePhoton_GenLabTime, &b_baseLinePhoton_GenLabTime);
   fChain->SetBranchAddress("baseLinePhoton_GenLabTimeSig", &baseLinePhoton_GenLabTimeSig, &b_baseLinePhoton_GenLabTimeSig);
   fChain->SetBranchAddress("baseLinePhoton_GenPt", &baseLinePhoton_GenPt, &b_baseLinePhoton_GenPt);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomEnergy", &baseLinePhoton_GenSigMomEnergy, &b_baseLinePhoton_GenSigMomEnergy);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomEta", &baseLinePhoton_GenSigMomEta, &b_baseLinePhoton_GenSigMomEta);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomMass", &baseLinePhoton_GenSigMomMass, &b_baseLinePhoton_GenSigMomMass);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomPhi", &baseLinePhoton_GenSigMomPhi, &b_baseLinePhoton_GenSigMomPhi);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomPt", &baseLinePhoton_GenSigMomPt, &b_baseLinePhoton_GenSigMomPt);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomPx", &baseLinePhoton_GenSigMomPx, &b_baseLinePhoton_GenSigMomPx);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomPy", &baseLinePhoton_GenSigMomPy, &b_baseLinePhoton_GenSigMomPy);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomPz", &baseLinePhoton_GenSigMomPz, &b_baseLinePhoton_GenSigMomPz);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomVx", &baseLinePhoton_GenSigMomVx, &b_baseLinePhoton_GenSigMomVx);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomVy", &baseLinePhoton_GenSigMomVy, &b_baseLinePhoton_GenSigMomVy);
   fChain->SetBranchAddress("baseLinePhoton_GenSigMomVz", &baseLinePhoton_GenSigMomVz, &b_baseLinePhoton_GenSigMomVz);
   fChain->SetBranchAddress("baseLinePhoton_GenTime", &baseLinePhoton_GenTime, &b_baseLinePhoton_GenTime);
   fChain->SetBranchAddress("baseLinePhoton_GenTimeSig", &baseLinePhoton_GenTimeSig, &b_baseLinePhoton_GenTimeSig);
   fChain->SetBranchAddress("baseLinePhoton_HadOverEM", &baseLinePhoton_HadOverEM, &b_baseLinePhoton_HadOverEM);
   fChain->SetBranchAddress("baseLinePhoton_HadTowOverEM", &baseLinePhoton_HadTowOverEM, &b_baseLinePhoton_HadTowOverEM);
   fChain->SetBranchAddress("baseLinePhoton_HasConversionTracks", &baseLinePhoton_HasConversionTracks, &b_baseLinePhoton_HasConversionTracks);
   fChain->SetBranchAddress("baseLinePhoton_HcalPFClusterIso", &baseLinePhoton_HcalPFClusterIso, &b_baseLinePhoton_HcalPFClusterIso);
   fChain->SetBranchAddress("baseLinePhoton_HcalTowerSumEtBcConeDR04", &baseLinePhoton_HcalTowerSumEtBcConeDR04, &b_baseLinePhoton_HcalTowerSumEtBcConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_HcalTowerSumEtConeDR04", &baseLinePhoton_HcalTowerSumEtConeDR04, &b_baseLinePhoton_HcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_Hoe_PUcorr", &baseLinePhoton_Hoe_PUcorr, &b_baseLinePhoton_Hoe_PUcorr);
   fChain->SetBranchAddress("baseLinePhoton_LSCross", &baseLinePhoton_LSCross, &b_baseLinePhoton_LSCross);
   fChain->SetBranchAddress("baseLinePhoton_LTRes", &baseLinePhoton_LTRes, &b_baseLinePhoton_LTRes);
   fChain->SetBranchAddress("baseLinePhoton_LTime", &baseLinePhoton_LTime, &b_baseLinePhoton_LTime);
   fChain->SetBranchAddress("baseLinePhoton_LTimeSig", &baseLinePhoton_LTimeSig, &b_baseLinePhoton_LTimeSig);
   fChain->SetBranchAddress("baseLinePhoton_LWeird", &baseLinePhoton_LWeird, &b_baseLinePhoton_LWeird);
   fChain->SetBranchAddress("baseLinePhoton_MBetaEql", &baseLinePhoton_MBetaEql, &b_baseLinePhoton_MBetaEql);
   fChain->SetBranchAddress("baseLinePhoton_MBetaPmt", &baseLinePhoton_MBetaPmt, &b_baseLinePhoton_MBetaPmt);
   fChain->SetBranchAddress("baseLinePhoton_NTrkHollowConeDR04", &baseLinePhoton_NTrkHollowConeDR04, &b_baseLinePhoton_NTrkHollowConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_NTrkSolidConeDR04", &baseLinePhoton_NTrkSolidConeDR04, &b_baseLinePhoton_NTrkSolidConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_Nrh", &baseLinePhoton_Nrh, &b_baseLinePhoton_Nrh);
   fChain->SetBranchAddress("baseLinePhoton_OOT", &baseLinePhoton_OOT, &b_baseLinePhoton_OOT);
   fChain->SetBranchAddress("baseLinePhoton_PVtoRH", &baseLinePhoton_PVtoRH, &b_baseLinePhoton_PVtoRH);
   fChain->SetBranchAddress("baseLinePhoton_PfChargedIso", &baseLinePhoton_PfChargedIso, &b_baseLinePhoton_PfChargedIso);
   fChain->SetBranchAddress("baseLinePhoton_PfChargedIsoPFPV", &baseLinePhoton_PfChargedIsoPFPV, &b_baseLinePhoton_PfChargedIsoPFPV);
   fChain->SetBranchAddress("baseLinePhoton_PfChargedIsoWorstVtx", &baseLinePhoton_PfChargedIsoWorstVtx, &b_baseLinePhoton_PfChargedIsoWorstVtx);
   fChain->SetBranchAddress("baseLinePhoton_PfPhoIso03", &baseLinePhoton_PfPhoIso03, &b_baseLinePhoton_PfPhoIso03);
   fChain->SetBranchAddress("baseLinePhoton_PfRelIso03_all_quadratic", &baseLinePhoton_PfRelIso03_all_quadratic, &b_baseLinePhoton_PfRelIso03_all_quadratic);
   fChain->SetBranchAddress("baseLinePhoton_PfRelIso03_chg_quadratic", &baseLinePhoton_PfRelIso03_chg_quadratic, &b_baseLinePhoton_PfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("baseLinePhoton_Phi", &baseLinePhoton_Phi, &b_baseLinePhoton_Phi);
   fChain->SetBranchAddress("baseLinePhoton_PhiWidth", &baseLinePhoton_PhiWidth, &b_baseLinePhoton_PhiWidth);
   fChain->SetBranchAddress("baseLinePhoton_PhoIsoDr", &baseLinePhoton_PhoIsoDr, &b_baseLinePhoton_PhoIsoDr);
   fChain->SetBranchAddress("baseLinePhoton_PixelSeed", &baseLinePhoton_PixelSeed, &b_baseLinePhoton_PixelSeed);
   fChain->SetBranchAddress("baseLinePhoton_Pt", &baseLinePhoton_Pt, &b_baseLinePhoton_Pt);
   fChain->SetBranchAddress("baseLinePhoton_Quality", &baseLinePhoton_Quality, &b_baseLinePhoton_Quality);
   fChain->SetBranchAddress("baseLinePhoton_R9", &baseLinePhoton_R9, &b_baseLinePhoton_R9);
   fChain->SetBranchAddress("baseLinePhoton_RawTime", &baseLinePhoton_RawTime, &b_baseLinePhoton_RawTime);
   fChain->SetBranchAddress("baseLinePhoton_S4", &baseLinePhoton_S4, &b_baseLinePhoton_S4);
   fChain->SetBranchAddress("baseLinePhoton_SAlp", &baseLinePhoton_SAlp, &b_baseLinePhoton_SAlp);
   fChain->SetBranchAddress("baseLinePhoton_SCx", &baseLinePhoton_SCx, &b_baseLinePhoton_SCx);
   fChain->SetBranchAddress("baseLinePhoton_SCy", &baseLinePhoton_SCy, &b_baseLinePhoton_SCy);
   fChain->SetBranchAddress("baseLinePhoton_SCz", &baseLinePhoton_SCz, &b_baseLinePhoton_SCz);
   fChain->SetBranchAddress("baseLinePhoton_SMaj", &baseLinePhoton_SMaj, &b_baseLinePhoton_SMaj);
   fChain->SetBranchAddress("baseLinePhoton_SMin", &baseLinePhoton_SMin, &b_baseLinePhoton_SMin);
   fChain->SetBranchAddress("baseLinePhoton_SSCross", &baseLinePhoton_SSCross, &b_baseLinePhoton_SSCross);
   fChain->SetBranchAddress("baseLinePhoton_STRes", &baseLinePhoton_STRes, &b_baseLinePhoton_STRes);
   fChain->SetBranchAddress("baseLinePhoton_STime", &baseLinePhoton_STime, &b_baseLinePhoton_STime);
   fChain->SetBranchAddress("baseLinePhoton_STimeSig", &baseLinePhoton_STimeSig, &b_baseLinePhoton_STimeSig);
   fChain->SetBranchAddress("baseLinePhoton_SWeird", &baseLinePhoton_SWeird, &b_baseLinePhoton_SWeird);
   fChain->SetBranchAddress("baseLinePhoton_ShasGS", &baseLinePhoton_ShasGS, &b_baseLinePhoton_ShasGS);
   fChain->SetBranchAddress("baseLinePhoton_Sieie", &baseLinePhoton_Sieie, &b_baseLinePhoton_Sieie);
   fChain->SetBranchAddress("baseLinePhoton_Sieip", &baseLinePhoton_Sieip, &b_baseLinePhoton_Sieip);
   fChain->SetBranchAddress("baseLinePhoton_SigmaIEtaIEta", &baseLinePhoton_SigmaIEtaIEta, &b_baseLinePhoton_SigmaIEtaIEta);
   fChain->SetBranchAddress("baseLinePhoton_Sipip", &baseLinePhoton_Sipip, &b_baseLinePhoton_Sipip);
   fChain->SetBranchAddress("baseLinePhoton_SusyId", &baseLinePhoton_SusyId, &b_baseLinePhoton_SusyId);
   fChain->SetBranchAddress("baseLinePhoton_Time", &baseLinePhoton_Time, &b_baseLinePhoton_Time);
   fChain->SetBranchAddress("baseLinePhoton_TrkSumPtHollowConeDR03", &baseLinePhoton_TrkSumPtHollowConeDR03, &b_baseLinePhoton_TrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("baseLinePhoton_TrkSumPtHollowConeDR04", &baseLinePhoton_TrkSumPtHollowConeDR04, &b_baseLinePhoton_TrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_TrkSumPtSolidConeDR04", &baseLinePhoton_TrkSumPtSolidConeDR04, &b_baseLinePhoton_TrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("baseLinePhoton_WTRes", &baseLinePhoton_WTRes, &b_baseLinePhoton_WTRes);
   fChain->SetBranchAddress("baseLinePhoton_WTResOld", &baseLinePhoton_WTResOld, &b_baseLinePhoton_WTResOld);
   fChain->SetBranchAddress("baseLinePhoton_WTime", &baseLinePhoton_WTime, &b_baseLinePhoton_WTime);
   fChain->SetBranchAddress("baseLinePhoton_WTimeOld", &baseLinePhoton_WTimeOld, &b_baseLinePhoton_WTimeOld);
   fChain->SetBranchAddress("baseLinePhoton_WTimeSig", &baseLinePhoton_WTimeSig, &b_baseLinePhoton_WTimeSig);
   fChain->SetBranchAddress("baseLinePhoton_WTimeSigOld", &baseLinePhoton_WTimeSigOld, &b_baseLinePhoton_WTimeSigOld);
   fChain->SetBranchAddress("baseLinePhoton_beamHaloCNNScore", &baseLinePhoton_beamHaloCNNScore, &b_baseLinePhoton_beamHaloCNNScore);
   fChain->SetBranchAddress("baseLinePhoton_beamHaloCR", &baseLinePhoton_beamHaloCR, &b_baseLinePhoton_beamHaloCR);
   fChain->SetBranchAddress("baseLinePhoton_beamHaloID", &baseLinePhoton_beamHaloID, &b_baseLinePhoton_beamHaloID);
   fChain->SetBranchAddress("baseLinePhoton_isoANNScore", &baseLinePhoton_isoANNScore, &b_baseLinePhoton_isoANNScore);
   fChain->SetBranchAddress("baseLinePhoton_looseNotTightIso1ID", &baseLinePhoton_looseNotTightIso1ID, &b_baseLinePhoton_looseNotTightIso1ID);
   fChain->SetBranchAddress("baseLinePhoton_looseNotTightIso2ID", &baseLinePhoton_looseNotTightIso2ID, &b_baseLinePhoton_looseNotTightIso2ID);
   fChain->SetBranchAddress("baseLinePhoton_looseNotTightIsoID", &baseLinePhoton_looseNotTightIsoID, &b_baseLinePhoton_looseNotTightIsoID);
   fChain->SetBranchAddress("baseLinePhoton_minJetDr", &baseLinePhoton_minJetDr, &b_baseLinePhoton_minJetDr);
   fChain->SetBranchAddress("baseLinePhoton_nonIsoANNScore", &baseLinePhoton_nonIsoANNScore, &b_baseLinePhoton_nonIsoANNScore);
   fChain->SetBranchAddress("baseLinePhoton_notBeamHaloID", &baseLinePhoton_notBeamHaloID, &b_baseLinePhoton_notBeamHaloID);
   fChain->SetBranchAddress("baseLinePhoton_physBkgCNNScore", &baseLinePhoton_physBkgCNNScore, &b_baseLinePhoton_physBkgCNNScore);
   fChain->SetBranchAddress("baseLinePhoton_spikeCR", &baseLinePhoton_spikeCR, &b_baseLinePhoton_spikeCR);
   fChain->SetBranchAddress("baseLinePhoton_tightIsoID", &baseLinePhoton_tightIsoID, &b_baseLinePhoton_tightIsoID);
   fChain->SetBranchAddress("endcap_photon_baseline", &endcap_photon_baseline, &b_endcap_photon_baseline);
   fChain->SetBranchAddress("endcap_photon_energy", &endcap_photon_energy, &b_endcap_photon_energy);
   fChain->SetBranchAddress("endcap_photon_erhsecdr4", &endcap_photon_erhsecdr4, &b_endcap_photon_erhsecdr4);
   fChain->SetBranchAddress("endcap_photon_eta", &endcap_photon_eta, &b_endcap_photon_eta);
   fChain->SetBranchAddress("endcap_photon_htoem", &endcap_photon_htoem, &b_endcap_photon_htoem);
   fChain->SetBranchAddress("endcap_photon_isoANNScore", &endcap_photon_isoANNScore, &b_endcap_photon_isoANNScore);
   fChain->SetBranchAddress("endcap_photon_phi", &endcap_photon_phi, &b_endcap_photon_phi);
   fChain->SetBranchAddress("endcap_photon_pt", &endcap_photon_pt, &b_endcap_photon_pt);
   fChain->SetBranchAddress("endcap_photon_tspscdr4", &endcap_photon_tspscdr4, &b_endcap_photon_tspscdr4);
   fChain->SetBranchAddress("nBaseLinePhotons", &nBaseLinePhotons, &b_nBaseLinePhotons);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("nblBHPhotons", &nblBHPhotons, &b_nblBHPhotons);
   fChain->SetBranchAddress("nblLooseNotTightIsoPhotons", &nblLooseNotTightIsoPhotons, &b_nblLooseNotTightIsoPhotons);
   fChain->SetBranchAddress("nblNonpromptPhotons", &nblNonpromptPhotons, &b_nblNonpromptPhotons);
   fChain->SetBranchAddress("nblNotBHPhotons", &nblNotBHPhotons, &b_nblNotBHPhotons);
   fChain->SetBranchAddress("nblPromptPhotons", &nblPromptPhotons, &b_nblPromptPhotons);
   fChain->SetBranchAddress("nblTightIsoPhotons", &nblTightIsoPhotons, &b_nblTightIsoPhotons);
   fChain->SetBranchAddress("passNPhoEq1SelectionBeamHaloCR", &passNPhoEq1SelectionBeamHaloCR, &b_passNPhoEq1SelectionBeamHaloCR);
   fChain->SetBranchAddress("passNPhoEq1SelectionEarlyBeamHaloCR", &passNPhoEq1SelectionEarlyBeamHaloCR, &b_passNPhoEq1SelectionEarlyBeamHaloCR);
   fChain->SetBranchAddress("passNPhoEq1SelectionEarlyNotBHCR", &passNPhoEq1SelectionEarlyNotBHCR, &b_passNPhoEq1SelectionEarlyNotBHCR);
   fChain->SetBranchAddress("passNPhoEq1SelectionLateBeamHaloCR", &passNPhoEq1SelectionLateBeamHaloCR, &b_passNPhoEq1SelectionLateBeamHaloCR);
   fChain->SetBranchAddress("passNPhoEq1SelectionLateNotBHSR", &passNPhoEq1SelectionLateNotBHSR, &b_passNPhoEq1SelectionLateNotBHSR);
   fChain->SetBranchAddress("passNPhoEq1SelectionPromptLooseNotTightIso1CR", &passNPhoEq1SelectionPromptLooseNotTightIso1CR, &b_passNPhoEq1SelectionPromptLooseNotTightIso1CR);
   fChain->SetBranchAddress("passNPhoEq1SelectionPromptLooseNotTightIso2CR", &passNPhoEq1SelectionPromptLooseNotTightIso2CR, &b_passNPhoEq1SelectionPromptLooseNotTightIso2CR);
   fChain->SetBranchAddress("passNPhoEq1SelectionPromptLooseNotTightIsoCR", &passNPhoEq1SelectionPromptLooseNotTightIsoCR, &b_passNPhoEq1SelectionPromptLooseNotTightIsoCR);
   fChain->SetBranchAddress("passNPhoEq1SelectionPromptTightIsoSR", &passNPhoEq1SelectionPromptTightIsoSR, &b_passNPhoEq1SelectionPromptTightIsoSR);
   fChain->SetBranchAddress("passNPhoEq2SelectionBeamHaloCR", &passNPhoEq2SelectionBeamHaloCR, &b_passNPhoEq2SelectionBeamHaloCR);
   fChain->SetBranchAddress("passNPhoEq2SelectionEarlyBeamHaloCR", &passNPhoEq2SelectionEarlyBeamHaloCR, &b_passNPhoEq2SelectionEarlyBeamHaloCR);
   fChain->SetBranchAddress("passNPhoEq2SelectionEarlyNotBHCR", &passNPhoEq2SelectionEarlyNotBHCR, &b_passNPhoEq2SelectionEarlyNotBHCR);
   fChain->SetBranchAddress("passNPhoEq2SelectionLateBeamHaloCR", &passNPhoEq2SelectionLateBeamHaloCR, &b_passNPhoEq2SelectionLateBeamHaloCR);
   fChain->SetBranchAddress("passNPhoEq2SelectionLateNotBHSR", &passNPhoEq2SelectionLateNotBHSR, &b_passNPhoEq2SelectionLateNotBHSR);
   fChain->SetBranchAddress("passNPhoEq2SelectionPromptLooseNotTightIso1CR", &passNPhoEq2SelectionPromptLooseNotTightIso1CR, &b_passNPhoEq2SelectionPromptLooseNotTightIso1CR);
   fChain->SetBranchAddress("passNPhoEq2SelectionPromptLooseNotTightIso2CR", &passNPhoEq2SelectionPromptLooseNotTightIso2CR, &b_passNPhoEq2SelectionPromptLooseNotTightIso2CR);
   fChain->SetBranchAddress("passNPhoEq2SelectionPromptLooseNotTightIsoCR", &passNPhoEq2SelectionPromptLooseNotTightIsoCR, &b_passNPhoEq2SelectionPromptLooseNotTightIsoCR);
   fChain->SetBranchAddress("passNPhoEq2SelectionPromptTightIsoSR", &passNPhoEq2SelectionPromptTightIsoSR, &b_passNPhoEq2SelectionPromptTightIsoSR);
   fChain->SetBranchAddress("passNPhoGe1SelectionBeamHaloCR", &passNPhoGe1SelectionBeamHaloCR, &b_passNPhoGe1SelectionBeamHaloCR);
   fChain->SetBranchAddress("passNPhoGe1SelectionEarlyBeamHaloCR", &passNPhoGe1SelectionEarlyBeamHaloCR, &b_passNPhoGe1SelectionEarlyBeamHaloCR);
   fChain->SetBranchAddress("passNPhoGe1SelectionEarlyNotBHCR", &passNPhoGe1SelectionEarlyNotBHCR, &b_passNPhoGe1SelectionEarlyNotBHCR);
   fChain->SetBranchAddress("passNPhoGe1SelectionLateBeamHaloCR", &passNPhoGe1SelectionLateBeamHaloCR, &b_passNPhoGe1SelectionLateBeamHaloCR);
   fChain->SetBranchAddress("passNPhoGe1SelectionLateNotBHSR", &passNPhoGe1SelectionLateNotBHSR, &b_passNPhoGe1SelectionLateNotBHSR);
   fChain->SetBranchAddress("passNPhoGe1SelectionNotBHCR", &passNPhoGe1SelectionNotBHCR, &b_passNPhoGe1SelectionNotBHCR);
   fChain->SetBranchAddress("passNPhoGe1SelectionNotBHPrompt", &passNPhoGe1SelectionNotBHPrompt, &b_passNPhoGe1SelectionNotBHPrompt);
   fChain->SetBranchAddress("passNPhoGe1SelectionPromptLooseNotTightIso1CR", &passNPhoGe1SelectionPromptLooseNotTightIso1CR, &b_passNPhoGe1SelectionPromptLooseNotTightIso1CR);
   fChain->SetBranchAddress("passNPhoGe1SelectionPromptLooseNotTightIso2CR", &passNPhoGe1SelectionPromptLooseNotTightIso2CR, &b_passNPhoGe1SelectionPromptLooseNotTightIso2CR);
   fChain->SetBranchAddress("passNPhoGe1SelectionPromptLooseNotTightIsoCR", &passNPhoGe1SelectionPromptLooseNotTightIsoCR, &b_passNPhoGe1SelectionPromptLooseNotTightIsoCR);
   fChain->SetBranchAddress("passNPhoGe1SelectionPromptTightIsoSR", &passNPhoGe1SelectionPromptTightIsoSR, &b_passNPhoGe1SelectionPromptTightIsoSR);
   fChain->SetBranchAddress("photon_E", &photon_E, &b_photon_E);
   fChain->SetBranchAddress("photon_Eta", &photon_Eta, &b_photon_Eta);
   fChain->SetBranchAddress("photon_LeadTime", &photon_LeadTime, &b_photon_LeadTime);
   fChain->SetBranchAddress("photon_LeadTimeSig", &photon_LeadTimeSig, &b_photon_LeadTimeSig);
   fChain->SetBranchAddress("photon_Phi", &photon_Phi, &b_photon_Phi);
   fChain->SetBranchAddress("photon_PixSeed", &photon_PixSeed, &b_photon_PixSeed);
   fChain->SetBranchAddress("photon_Pt", &photon_Pt, &b_photon_Pt);
   fChain->SetBranchAddress("photon_WTime", &photon_WTime, &b_photon_WTime);
   fChain->SetBranchAddress("photon_WTimeSig", &photon_WTimeSig, &b_photon_WTimeSig);
   fChain->SetBranchAddress("photon_baseline", &photon_baseline, &b_photon_baseline);
   fChain->SetBranchAddress("photon_beamHaloCNNScore", &photon_beamHaloCNNScore, &b_photon_beamHaloCNNScore);
   fChain->SetBranchAddress("photon_isoANNScore", &photon_isoANNScore, &b_photon_isoANNScore);
   fChain->SetBranchAddress("photon_minJetDr", &photon_minJetDr, &b_photon_minJetDr);
   fChain->SetBranchAddress("tagged_lead_timesig", &tagged_lead_timesig, &b_tagged_lead_timesig);
   fChain->SetBranchAddress("rjrABSide", &rjrABSide, &b_rjrABSide);
   fChain->SetBranchAddress("rjr_Rdiff", &rjr_Rdiff, &b_rjr_Rdiff);
   fChain->SetBranchAddress("rjrCleaningVeto0", &rjrCleaningVeto0, &b_rjrCleaningVeto0);
   fChain->SetBranchAddress("rjrCleaningVeto1", &rjrCleaningVeto1, &b_rjrCleaningVeto1);
   fChain->SetBranchAddress("rjrCleaningVeto2", &rjrCleaningVeto2, &b_rjrCleaningVeto2);
   fChain->SetBranchAddress("rjrDPhiSI", &rjrDPhiSI, &b_rjrDPhiSI);
   fChain->SetBranchAddress("rjrDphiMETV", &rjrDphiMETV, &b_rjrDphiMETV);
   fChain->SetBranchAddress("rjrEVa", &rjrEVa, &b_rjrEVa);
   fChain->SetBranchAddress("rjrEVb", &rjrEVb, &b_rjrEVb);
   fChain->SetBranchAddress("rjrLeadPhoLocation", &rjrLeadPhoLocation, &b_rjrLeadPhoLocation);
   fChain->SetBranchAddress("rjrMET", &rjrMET, &b_rjrMET);
   fChain->SetBranchAddress("rjrN2Px", &rjrN2Px, &b_rjrN2Px);
   fChain->SetBranchAddress("rjrN2Py", &rjrN2Py, &b_rjrN2Py);
   fChain->SetBranchAddress("rjrN2Pz", &rjrN2Pz, &b_rjrN2Pz);
   fChain->SetBranchAddress("rjrNJets", &rjrNJets, &b_rjrNJets);
   fChain->SetBranchAddress("rjrNJetsJa", &rjrNJetsJa, &b_rjrNJetsJa);
   fChain->SetBranchAddress("rjrNJetsJb", &rjrNJetsJb, &b_rjrNJetsJb);
   fChain->SetBranchAddress("rjrNPhotons", &rjrNPhotons, &b_rjrNPhotons);
   fChain->SetBranchAddress("rjrNVisObjects", &rjrNVisObjects, &b_rjrNVisObjects);
   fChain->SetBranchAddress("rjrPS", &rjrPS, &b_rjrPS);
   fChain->SetBranchAddress("rjrPTS", &rjrPTS, &b_rjrPTS);
   fChain->SetBranchAddress("rjrPVa", &rjrPVa, &b_rjrPVa);
   fChain->SetBranchAddress("rjrPVb", &rjrPVb, &b_rjrPVb);
   fChain->SetBranchAddress("rjrPVlab", &rjrPVlab, &b_rjrPVlab);
   fChain->SetBranchAddress("rjrPZS", &rjrPZS, &b_rjrPZS);
   fChain->SetBranchAddress("rjrSCosA", &rjrSCosA, &b_rjrSCosA);
   fChain->SetBranchAddress("rjrSMass", &rjrSMass, &b_rjrSMass);
   fChain->SetBranchAddress("rjrSdphiBV", &rjrSdphiBV, &b_rjrSdphiBV);
   fChain->SetBranchAddress("rjrSdphiDA", &rjrSdphiDA, &b_rjrSdphiDA);
   fChain->SetBranchAddress("rjrSubPhoLocation", &rjrSubPhoLocation, &b_rjrSubPhoLocation);
   fChain->SetBranchAddress("rjrType", &rjrType, &b_rjrType);
   fChain->SetBranchAddress("rjrVisJetSide", &rjrVisJetSide, &b_rjrVisJetSide);
   fChain->SetBranchAddress("rjrVisPhoSide", &rjrVisPhoSide, &b_rjrVisPhoSide);
   fChain->SetBranchAddress("rjrX1aCosA", &rjrX1aCosA, &b_rjrX1aCosA);
   fChain->SetBranchAddress("rjrX1aMass", &rjrX1aMass, &b_rjrX1aMass);
   fChain->SetBranchAddress("rjrX1a_Eta", &rjrX1a_Eta, &b_rjrX1a_Eta);
   fChain->SetBranchAddress("rjrX1a_Phi", &rjrX1a_Phi, &b_rjrX1a_Phi);
   fChain->SetBranchAddress("rjrX1a_Pt", &rjrX1a_Pt, &b_rjrX1a_Pt);
   fChain->SetBranchAddress("rjrX1bCosA", &rjrX1bCosA, &b_rjrX1bCosA);
   fChain->SetBranchAddress("rjrX1bMass", &rjrX1bMass, &b_rjrX1bMass);
   fChain->SetBranchAddress("rjrX1b_Eta", &rjrX1b_Eta, &b_rjrX1b_Eta);
   fChain->SetBranchAddress("rjrX1b_Phi", &rjrX1b_Phi, &b_rjrX1b_Phi);
   fChain->SetBranchAddress("rjrX1b_Pt", &rjrX1b_Pt, &b_rjrX1b_Pt);
   fChain->SetBranchAddress("rjrX2aCosA", &rjrX2aCosA, &b_rjrX2aCosA);
   fChain->SetBranchAddress("rjrX2aMass", &rjrX2aMass, &b_rjrX2aMass);
   fChain->SetBranchAddress("rjrX2aPs", &rjrX2aPs, &b_rjrX2aPs);
   fChain->SetBranchAddress("rjrX2bCosA", &rjrX2bCosA, &b_rjrX2bCosA);
   fChain->SetBranchAddress("rjrX2bMass", &rjrX2bMass, &b_rjrX2bMass);
   fChain->SetBranchAddress("rjrX2bPs", &rjrX2bPs, &b_rjrX2bPs);
   fChain->SetBranchAddress("rjr_Mr", &rjr_Mr, &b_rjr_Mr);
   fChain->SetBranchAddress("rjr_Ms", &rjr_Ms, &b_rjr_Ms);
   fChain->SetBranchAddress("rjr_Mva", &rjr_Mva, &b_rjr_Mva);
   fChain->SetBranchAddress("rjr_Mvb", &rjr_Mvb, &b_rjr_Mvb);
   fChain->SetBranchAddress("rjr_R", &rjr_R, &b_rjr_R);
   fChain->SetBranchAddress("rjr_Rm", &rjr_Rm, &b_rjr_Rm);
   fChain->SetBranchAddress("rjr_Rs", &rjr_Rs, &b_rjr_Rs);
   fChain->SetBranchAddress("rjr_Rv", &rjr_Rv, &b_rjr_Rv);
   fChain->SetBranchAddress("rjr_Rx", &rjr_Rx, &b_rjr_Rx);
   fChain->SetBranchAddress("rjr_Rxa", &rjr_Rxa, &b_rjr_Rxa);
   fChain->SetBranchAddress("rjr_Rxb", &rjr_Rxb, &b_rjr_Rxb);
   fChain->SetBranchAddress("rjr_p1Pts", &rjr_p1Pts, &b_rjr_p1Pts);
   fChain->SetBranchAddress("rjr_p1Ptxa", &rjr_p1Ptxa, &b_rjr_p1Ptxa);
   fChain->SetBranchAddress("rjr_p1Ptxb", &rjr_p1Ptxb, &b_rjr_p1Ptxb);
   fChain->SetBranchAddress("rjr_p1Pxa", &rjr_p1Pxa, &b_rjr_p1Pxa);
   fChain->SetBranchAddress("rjr_p1Pxb", &rjr_p1Pxb, &b_rjr_p1Pxb);
   fChain->SetBranchAddress("rjr_p2Pts", &rjr_p2Pts, &b_rjr_p2Pts);
   fChain->SetBranchAddress("rjr_p2Ptxa", &rjr_p2Ptxa, &b_rjr_p2Ptxa);
   fChain->SetBranchAddress("rjr_p2Ptxb", &rjr_p2Ptxb, &b_rjr_p2Ptxb);
   fChain->SetBranchAddress("rjr_p2Pxa", &rjr_p2Pxa, &b_rjr_p2Pxa);
   fChain->SetBranchAddress("rjr_p2Pxb", &rjr_p2Pxb, &b_rjr_p2Pxb);
   fChain->SetBranchAddress("rjr_pHs11", &rjr_pHs11, &b_rjr_pHs11);
   fChain->SetBranchAddress("rjr_pHs11a", &rjr_pHs11a, &b_rjr_pHs11a);
   fChain->SetBranchAddress("rjr_pHs11b", &rjr_pHs11b, &b_rjr_pHs11b);
   fChain->SetBranchAddress("rjr_pHs20", &rjr_pHs20, &b_rjr_pHs20);
   fChain->SetBranchAddress("rjr_pHs21", &rjr_pHs21, &b_rjr_pHs21);
   fChain->SetBranchAddress("rjr_pHs21a", &rjr_pHs21a, &b_rjr_pHs21a);
   fChain->SetBranchAddress("rjr_pHs21b", &rjr_pHs21b, &b_rjr_pHs21b);
   fChain->SetBranchAddress("rjr_pHs22", &rjr_pHs22, &b_rjr_pHs22);
   fChain->SetBranchAddress("rjr_pHs40", &rjr_pHs40, &b_rjr_pHs40);
   fChain->SetBranchAddress("rjr_pHs41", &rjr_pHs41, &b_rjr_pHs41);
   fChain->SetBranchAddress("rjr_pHts11", &rjr_pHts11, &b_rjr_pHts11);
   fChain->SetBranchAddress("rjr_pHts11a", &rjr_pHts11a, &b_rjr_pHts11a);
   fChain->SetBranchAddress("rjr_pHts11b", &rjr_pHts11b, &b_rjr_pHts11b);
   fChain->SetBranchAddress("rjr_pHts20", &rjr_pHts20, &b_rjr_pHts20);
   fChain->SetBranchAddress("rjr_pHts21", &rjr_pHts21, &b_rjr_pHts21);
   fChain->SetBranchAddress("rjr_pHts21a", &rjr_pHts21a, &b_rjr_pHts21a);
   fChain->SetBranchAddress("rjr_pHts21b", &rjr_pHts21b, &b_rjr_pHts21b);
   fChain->SetBranchAddress("rjr_pHts22", &rjr_pHts22, &b_rjr_pHts22);
   fChain->SetBranchAddress("rjr_pHts40", &rjr_pHts40, &b_rjr_pHts40);
   fChain->SetBranchAddress("rjr_pHts41", &rjr_pHts41, &b_rjr_pHts41);
   fChain->SetBranchAddress("rjr_pHtxa10", &rjr_pHtxa10, &b_rjr_pHtxa10);
   fChain->SetBranchAddress("rjr_pHtxa11", &rjr_pHtxa11, &b_rjr_pHtxa11);
   fChain->SetBranchAddress("rjr_pHtxa20", &rjr_pHtxa20, &b_rjr_pHtxa20);
   fChain->SetBranchAddress("rjr_pHtxa21", &rjr_pHtxa21, &b_rjr_pHtxa21);
   fChain->SetBranchAddress("rjr_pHtxb10", &rjr_pHtxb10, &b_rjr_pHtxb10);
   fChain->SetBranchAddress("rjr_pHtxb11", &rjr_pHtxb11, &b_rjr_pHtxb11);
   fChain->SetBranchAddress("rjr_pHtxb20", &rjr_pHtxb20, &b_rjr_pHtxb20);
   fChain->SetBranchAddress("rjr_pHtxb21", &rjr_pHtxb21, &b_rjr_pHtxb21);
   fChain->SetBranchAddress("rjr_pHxa10", &rjr_pHxa10, &b_rjr_pHxa10);
   fChain->SetBranchAddress("rjr_pHxa11", &rjr_pHxa11, &b_rjr_pHxa11);
   fChain->SetBranchAddress("rjr_pHxa20", &rjr_pHxa20, &b_rjr_pHxa20);
   fChain->SetBranchAddress("rjr_pHxa21", &rjr_pHxa21, &b_rjr_pHxa21);
   fChain->SetBranchAddress("rjr_pHxb10", &rjr_pHxb10, &b_rjr_pHxb10);
   fChain->SetBranchAddress("rjr_pHxb11", &rjr_pHxb11, &b_rjr_pHxb11);
   fChain->SetBranchAddress("rjr_pHxb20", &rjr_pHxb20, &b_rjr_pHxb20);
   fChain->SetBranchAddress("rjr_pHxb21", &rjr_pHxb21, &b_rjr_pHxb21);
   fChain->SetBranchAddress("rjr_sJ1a_e", &rjr_sJ1a_e, &b_rjr_sJ1a_e);
   fChain->SetBranchAddress("rjr_sJ1a_px", &rjr_sJ1a_px, &b_rjr_sJ1a_px);
   fChain->SetBranchAddress("rjr_sJ1a_py", &rjr_sJ1a_py, &b_rjr_sJ1a_py);
   fChain->SetBranchAddress("rjr_sJ1a_pz", &rjr_sJ1a_pz, &b_rjr_sJ1a_pz);
   fChain->SetBranchAddress("rjr_sJ1b_e", &rjr_sJ1b_e, &b_rjr_sJ1b_e);
   fChain->SetBranchAddress("rjr_sJ1b_px", &rjr_sJ1b_px, &b_rjr_sJ1b_px);
   fChain->SetBranchAddress("rjr_sJ1b_py", &rjr_sJ1b_py, &b_rjr_sJ1b_py);
   fChain->SetBranchAddress("rjr_sJ1b_pz", &rjr_sJ1b_pz, &b_rjr_sJ1b_pz);
   fChain->SetBranchAddress("rjr_sJ2a_e", &rjr_sJ2a_e, &b_rjr_sJ2a_e);
   fChain->SetBranchAddress("rjr_sJ2a_px", &rjr_sJ2a_px, &b_rjr_sJ2a_px);
   fChain->SetBranchAddress("rjr_sJ2a_py", &rjr_sJ2a_py, &b_rjr_sJ2a_py);
   fChain->SetBranchAddress("rjr_sJ2a_pz", &rjr_sJ2a_pz, &b_rjr_sJ2a_pz);
   fChain->SetBranchAddress("rjr_sJ2b_e", &rjr_sJ2b_e, &b_rjr_sJ2b_e);
   fChain->SetBranchAddress("rjr_sJ2b_px", &rjr_sJ2b_px, &b_rjr_sJ2b_px);
   fChain->SetBranchAddress("rjr_sJ2b_py", &rjr_sJ2b_py, &b_rjr_sJ2b_py);
   fChain->SetBranchAddress("rjr_sJ2b_pz", &rjr_sJ2b_pz, &b_rjr_sJ2b_pz);
   fChain->SetBranchAddress("rjrIsrABSide", &rjrIsrABSide, &b_rjrIsrABSide);
//    fChain->SetBranchAddress("rjr_Rdiff", &rjr_Rdiff, &b_rjr_Rdiff);
   fChain->SetBranchAddress("rjrIsrCleaningVeto0", &rjrIsrCleaningVeto0, &b_rjrIsrCleaningVeto0);
   fChain->SetBranchAddress("rjrIsrCleaningVeto1", &rjrIsrCleaningVeto1, &b_rjrIsrCleaningVeto1);
   fChain->SetBranchAddress("rjrIsrCleaningVeto2", &rjrIsrCleaningVeto2, &b_rjrIsrCleaningVeto2);
   fChain->SetBranchAddress("rjrIsrDPhiSI", &rjrIsrDPhiSI, &b_rjrIsrDPhiSI);
   fChain->SetBranchAddress("rjrIsrDphiMETV", &rjrIsrDphiMETV, &b_rjrIsrDphiMETV);
   fChain->SetBranchAddress("rjrIsrEVa", &rjrIsrEVa, &b_rjrIsrEVa);
   fChain->SetBranchAddress("rjrIsrEVb", &rjrIsrEVb, &b_rjrIsrEVb);
   fChain->SetBranchAddress("rjrIsrLeadPhoLocation", &rjrIsrLeadPhoLocation, &b_rjrIsrLeadPhoLocation);
   fChain->SetBranchAddress("rjrIsrMET", &rjrIsrMET, &b_rjrIsrMET);
   fChain->SetBranchAddress("rjrIsrNISRPhos", &rjrIsrNISRPhos, &b_rjrIsrNISRPhos);
   fChain->SetBranchAddress("rjrIsrNJets", &rjrIsrNJets, &b_rjrIsrNJets);
   fChain->SetBranchAddress("rjrIsrNJetsISR", &rjrIsrNJetsISR, &b_rjrIsrNJetsISR);
   fChain->SetBranchAddress("rjrIsrNJetsJa", &rjrIsrNJetsJa, &b_rjrIsrNJetsJa);
   fChain->SetBranchAddress("rjrIsrNJetsJb", &rjrIsrNJetsJb, &b_rjrIsrNJetsJb);
   fChain->SetBranchAddress("rjrIsrNPhotons", &rjrIsrNPhotons, &b_rjrIsrNPhotons);
   fChain->SetBranchAddress("rjrIsrNVisObjects", &rjrIsrNVisObjects, &b_rjrIsrNVisObjects);
   fChain->SetBranchAddress("rjrIsrPS", &rjrIsrPS, &b_rjrIsrPS);
   fChain->SetBranchAddress("rjrIsrPTS", &rjrIsrPTS, &b_rjrIsrPTS);
   fChain->SetBranchAddress("rjrIsrPVa", &rjrIsrPVa, &b_rjrIsrPVa);
   fChain->SetBranchAddress("rjrIsrPVb", &rjrIsrPVb, &b_rjrIsrPVb);
   fChain->SetBranchAddress("rjrIsrPVlab", &rjrIsrPVlab, &b_rjrIsrPVlab);
   fChain->SetBranchAddress("rjrIsrPZS", &rjrIsrPZS, &b_rjrIsrPZS);
   fChain->SetBranchAddress("rjrIsrSCosA", &rjrIsrSCosA, &b_rjrIsrSCosA);
   fChain->SetBranchAddress("rjrIsrSMass", &rjrIsrSMass, &b_rjrIsrSMass);
   fChain->SetBranchAddress("rjrIsrSdphiBV", &rjrIsrSdphiBV, &b_rjrIsrSdphiBV);
   fChain->SetBranchAddress("rjrIsrSdphiDA", &rjrIsrSdphiDA, &b_rjrIsrSdphiDA);
   fChain->SetBranchAddress("rjrIsrSubPhoLocation", &rjrIsrSubPhoLocation, &b_rjrIsrSubPhoLocation);
   fChain->SetBranchAddress("rjrIsrVisJetSide", &rjrIsrVisJetSide, &b_rjrIsrVisJetSide);
   fChain->SetBranchAddress("rjrIsrVisPhoSide", &rjrIsrVisPhoSide, &b_rjrIsrVisPhoSide);
   fChain->SetBranchAddress("rjrIsrX1aCosA", &rjrIsrX1aCosA, &b_rjrIsrX1aCosA);
   fChain->SetBranchAddress("rjrIsrX1aMass", &rjrIsrX1aMass, &b_rjrIsrX1aMass);
   fChain->SetBranchAddress("rjrIsrX1a_Eta", &rjrIsrX1a_Eta, &b_rjrIsrX1a_Eta);
   fChain->SetBranchAddress("rjrIsrX1a_Phi", &rjrIsrX1a_Phi, &b_rjrIsrX1a_Phi);
   fChain->SetBranchAddress("rjrIsrX1a_Pt", &rjrIsrX1a_Pt, &b_rjrIsrX1a_Pt);
   fChain->SetBranchAddress("rjrIsrX1bCosA", &rjrIsrX1bCosA, &b_rjrIsrX1bCosA);
   fChain->SetBranchAddress("rjrIsrX1bMass", &rjrIsrX1bMass, &b_rjrIsrX1bMass);
   fChain->SetBranchAddress("rjrIsrX1b_Eta", &rjrIsrX1b_Eta, &b_rjrIsrX1b_Eta);
   fChain->SetBranchAddress("rjrIsrX1b_Phi", &rjrIsrX1b_Phi, &b_rjrIsrX1b_Phi);
   fChain->SetBranchAddress("rjrIsrX1b_Pt", &rjrIsrX1b_Pt, &b_rjrIsrX1b_Pt);
   fChain->SetBranchAddress("rjrIsrX2aCosA", &rjrIsrX2aCosA, &b_rjrIsrX2aCosA);
   fChain->SetBranchAddress("rjrIsrX2aMass", &rjrIsrX2aMass, &b_rjrIsrX2aMass);
   fChain->SetBranchAddress("rjrIsrX2aPs", &rjrIsrX2aPs, &b_rjrIsrX2aPs);
   fChain->SetBranchAddress("rjrIsrX2bCosA", &rjrIsrX2bCosA, &b_rjrIsrX2bCosA);
   fChain->SetBranchAddress("rjrIsrX2bMass", &rjrIsrX2bMass, &b_rjrIsrX2bMass);
   fChain->SetBranchAddress("rjrIsrX2bPs", &rjrIsrX2bPs, &b_rjrIsrX2bPs);
   fChain->SetBranchAddress("rjrIsr_Mr", &rjrIsr_Mr, &b_rjrIsr_Mr);
   fChain->SetBranchAddress("rjrIsr_Ms", &rjrIsr_Ms, &b_rjrIsr_Ms);
   fChain->SetBranchAddress("rjrIsr_MsPerp", &rjrIsr_MsPerp, &b_rjrIsr_MsPerp);
   fChain->SetBranchAddress("rjrIsr_Mva", &rjrIsr_Mva, &b_rjrIsr_Mva);
   fChain->SetBranchAddress("rjrIsr_Mvb", &rjrIsr_Mvb, &b_rjrIsr_Mvb);
   fChain->SetBranchAddress("rjrIsr_PIsr", &rjrIsr_PIsr, &b_rjrIsr_PIsr);
   fChain->SetBranchAddress("rjrIsr_PtIsr", &rjrIsr_PtIsr, &b_rjrIsr_PtIsr);
   fChain->SetBranchAddress("rjrIsr_R", &rjrIsr_R, &b_rjrIsr_R);
   fChain->SetBranchAddress("rjrIsr_RIsr", &rjrIsr_RIsr, &b_rjrIsr_RIsr);
   fChain->SetBranchAddress("rjrIsr_Rm", &rjrIsr_Rm, &b_rjrIsr_Rm);
   fChain->SetBranchAddress("rjrIsr_Rs", &rjrIsr_Rs, &b_rjrIsr_Rs);
   fChain->SetBranchAddress("rjrIsr_Rv", &rjrIsr_Rv, &b_rjrIsr_Rv);
   fChain->SetBranchAddress("rjrIsr_Rx", &rjrIsr_Rx, &b_rjrIsr_Rx);
   fChain->SetBranchAddress("rjrIsr_Rxa", &rjrIsr_Rxa, &b_rjrIsr_Rxa);
   fChain->SetBranchAddress("rjrIsr_Rxb", &rjrIsr_Rxb, &b_rjrIsr_Rxb);
   fChain->SetBranchAddress("rjrIsr_nIsrVisObjects", &rjrIsr_nIsrVisObjects, &b_rjrIsr_nIsrVisObjects);
   fChain->SetBranchAddress("rjrIsr_nSVisObjects", &rjrIsr_nSVisObjects, &b_rjrIsr_nSVisObjects);
   fChain->SetBranchAddress("rjrIsr_p1Pts", &rjrIsr_p1Pts, &b_rjrIsr_p1Pts);
   fChain->SetBranchAddress("rjrIsr_p1Ptxa", &rjrIsr_p1Ptxa, &b_rjrIsr_p1Ptxa);
   fChain->SetBranchAddress("rjrIsr_p1Ptxb", &rjrIsr_p1Ptxb, &b_rjrIsr_p1Ptxb);
   fChain->SetBranchAddress("rjrIsr_p1Pxa", &rjrIsr_p1Pxa, &b_rjrIsr_p1Pxa);
   fChain->SetBranchAddress("rjrIsr_p1Pxb", &rjrIsr_p1Pxb, &b_rjrIsr_p1Pxb);
   fChain->SetBranchAddress("rjrIsr_p2Pts", &rjrIsr_p2Pts, &b_rjrIsr_p2Pts);
   fChain->SetBranchAddress("rjrIsr_p2Ptxa", &rjrIsr_p2Ptxa, &b_rjrIsr_p2Ptxa);
   fChain->SetBranchAddress("rjrIsr_p2Ptxb", &rjrIsr_p2Ptxb, &b_rjrIsr_p2Ptxb);
   fChain->SetBranchAddress("rjrIsr_p2Pxa", &rjrIsr_p2Pxa, &b_rjrIsr_p2Pxa);
   fChain->SetBranchAddress("rjrIsr_p2Pxb", &rjrIsr_p2Pxb, &b_rjrIsr_p2Pxb);
   fChain->SetBranchAddress("rjrIsr_pHs11", &rjrIsr_pHs11, &b_rjrIsr_pHs11);
   fChain->SetBranchAddress("rjrIsr_pHs11a", &rjrIsr_pHs11a, &b_rjrIsr_pHs11a);
   fChain->SetBranchAddress("rjrIsr_pHs11b", &rjrIsr_pHs11b, &b_rjrIsr_pHs11b);
   fChain->SetBranchAddress("rjrIsr_pHs20", &rjrIsr_pHs20, &b_rjrIsr_pHs20);
   fChain->SetBranchAddress("rjrIsr_pHs21", &rjrIsr_pHs21, &b_rjrIsr_pHs21);
   fChain->SetBranchAddress("rjrIsr_pHs21a", &rjrIsr_pHs21a, &b_rjrIsr_pHs21a);
   fChain->SetBranchAddress("rjrIsr_pHs21b", &rjrIsr_pHs21b, &b_rjrIsr_pHs21b);
   fChain->SetBranchAddress("rjrIsr_pHs22", &rjrIsr_pHs22, &b_rjrIsr_pHs22);
   fChain->SetBranchAddress("rjrIsr_pHs40", &rjrIsr_pHs40, &b_rjrIsr_pHs40);
   fChain->SetBranchAddress("rjrIsr_pHs41", &rjrIsr_pHs41, &b_rjrIsr_pHs41);
   fChain->SetBranchAddress("rjrIsr_pHts11", &rjrIsr_pHts11, &b_rjrIsr_pHts11);
   fChain->SetBranchAddress("rjrIsr_pHts11a", &rjrIsr_pHts11a, &b_rjrIsr_pHts11a);
   fChain->SetBranchAddress("rjrIsr_pHts11b", &rjrIsr_pHts11b, &b_rjrIsr_pHts11b);
   fChain->SetBranchAddress("rjrIsr_pHts20", &rjrIsr_pHts20, &b_rjrIsr_pHts20);
   fChain->SetBranchAddress("rjrIsr_pHts21", &rjrIsr_pHts21, &b_rjrIsr_pHts21);
   fChain->SetBranchAddress("rjrIsr_pHts21a", &rjrIsr_pHts21a, &b_rjrIsr_pHts21a);
   fChain->SetBranchAddress("rjrIsr_pHts21b", &rjrIsr_pHts21b, &b_rjrIsr_pHts21b);
   fChain->SetBranchAddress("rjrIsr_pHts22", &rjrIsr_pHts22, &b_rjrIsr_pHts22);
   fChain->SetBranchAddress("rjrIsr_pHts40", &rjrIsr_pHts40, &b_rjrIsr_pHts40);
   fChain->SetBranchAddress("rjrIsr_pHts41", &rjrIsr_pHts41, &b_rjrIsr_pHts41);
   fChain->SetBranchAddress("rjrIsr_pHtxa10", &rjrIsr_pHtxa10, &b_rjrIsr_pHtxa10);
   fChain->SetBranchAddress("rjrIsr_pHtxa11", &rjrIsr_pHtxa11, &b_rjrIsr_pHtxa11);
   fChain->SetBranchAddress("rjrIsr_pHtxa20", &rjrIsr_pHtxa20, &b_rjrIsr_pHtxa20);
   fChain->SetBranchAddress("rjrIsr_pHtxa21", &rjrIsr_pHtxa21, &b_rjrIsr_pHtxa21);
   fChain->SetBranchAddress("rjrIsr_pHtxb10", &rjrIsr_pHtxb10, &b_rjrIsr_pHtxb10);
   fChain->SetBranchAddress("rjrIsr_pHtxb11", &rjrIsr_pHtxb11, &b_rjrIsr_pHtxb11);
   fChain->SetBranchAddress("rjrIsr_pHtxb20", &rjrIsr_pHtxb20, &b_rjrIsr_pHtxb20);
   fChain->SetBranchAddress("rjrIsr_pHtxb21", &rjrIsr_pHtxb21, &b_rjrIsr_pHtxb21);
   fChain->SetBranchAddress("rjrIsr_pHxa10", &rjrIsr_pHxa10, &b_rjrIsr_pHxa10);
   fChain->SetBranchAddress("rjrIsr_pHxa11", &rjrIsr_pHxa11, &b_rjrIsr_pHxa11);
   fChain->SetBranchAddress("rjrIsr_pHxa20", &rjrIsr_pHxa20, &b_rjrIsr_pHxa20);
   fChain->SetBranchAddress("rjrIsr_pHxa21", &rjrIsr_pHxa21, &b_rjrIsr_pHxa21);
   fChain->SetBranchAddress("rjrIsr_pHxb10", &rjrIsr_pHxb10, &b_rjrIsr_pHxb10);
   fChain->SetBranchAddress("rjrIsr_pHxb11", &rjrIsr_pHxb11, &b_rjrIsr_pHxb11);
   fChain->SetBranchAddress("rjrIsr_pHxb20", &rjrIsr_pHxb20, &b_rjrIsr_pHxb20);
   fChain->SetBranchAddress("rjrIsr_pHxb21", &rjrIsr_pHxb21, &b_rjrIsr_pHxb21);
   fChain->SetBranchAddress("rjrIsr_sJ1a_e", &rjrIsr_sJ1a_e, &b_rjrIsr_sJ1a_e);
   fChain->SetBranchAddress("rjrIsr_sJ1a_px", &rjrIsr_sJ1a_px, &b_rjrIsr_sJ1a_px);
   fChain->SetBranchAddress("rjrIsr_sJ1a_py", &rjrIsr_sJ1a_py, &b_rjrIsr_sJ1a_py);
   fChain->SetBranchAddress("rjrIsr_sJ1a_pz", &rjrIsr_sJ1a_pz, &b_rjrIsr_sJ1a_pz);
   fChain->SetBranchAddress("rjrIsr_sJ1b_e", &rjrIsr_sJ1b_e, &b_rjrIsr_sJ1b_e);
   fChain->SetBranchAddress("rjrIsr_sJ1b_px", &rjrIsr_sJ1b_px, &b_rjrIsr_sJ1b_px);
   fChain->SetBranchAddress("rjrIsr_sJ1b_py", &rjrIsr_sJ1b_py, &b_rjrIsr_sJ1b_py);
   fChain->SetBranchAddress("rjrIsr_sJ1b_pz", &rjrIsr_sJ1b_pz, &b_rjrIsr_sJ1b_pz);
   fChain->SetBranchAddress("rjrIsr_sJ2a_e", &rjrIsr_sJ2a_e, &b_rjrIsr_sJ2a_e);
   fChain->SetBranchAddress("rjrIsr_sJ2a_px", &rjrIsr_sJ2a_px, &b_rjrIsr_sJ2a_px);
   fChain->SetBranchAddress("rjrIsr_sJ2a_py", &rjrIsr_sJ2a_py, &b_rjrIsr_sJ2a_py);
   fChain->SetBranchAddress("rjrIsr_sJ2a_pz", &rjrIsr_sJ2a_pz, &b_rjrIsr_sJ2a_pz);
   fChain->SetBranchAddress("rjrIsr_sJ2b_e", &rjrIsr_sJ2b_e, &b_rjrIsr_sJ2b_e);
   fChain->SetBranchAddress("rjrIsr_sJ2b_px", &rjrIsr_sJ2b_px, &b_rjrIsr_sJ2b_px);
   fChain->SetBranchAddress("rjrIsr_sJ2b_py", &rjrIsr_sJ2b_py, &b_rjrIsr_sJ2b_py);
   fChain->SetBranchAddress("rjrIsr_sJ2b_pz", &rjrIsr_sJ2b_pz, &b_rjrIsr_sJ2b_pz);
   fChain->SetBranchAddress("eleHasSVMatch", &eleHasSVMatch, &b_eleHasSVMatch);
   fChain->SetBranchAddress("eleIsLoose", &eleIsLoose, &b_eleIsLoose);
   fChain->SetBranchAddress("elePhoMinDr", &elePhoMinDr, &b_elePhoMinDr);
   fChain->SetBranchAddress("epDrMatch", &epDrMatch, &b_epDrMatch);
   fChain->SetBranchAddress("epSeedIdMatch", &epSeedIdMatch, &b_epSeedIdMatch);
   fChain->SetBranchAddress("epSeedIdUnique", &epSeedIdUnique, &b_epSeedIdUnique);
   fChain->SetBranchAddress("nEleSVMatched", &nEleSVMatched, &b_nEleSVMatched);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nSelElectrons", &nSelElectrons, &b_nSelElectrons);
   fChain->SetBranchAddress("nSelIsoElectrons", &nSelIsoElectrons, &b_nSelIsoElectrons);
   fChain->SetBranchAddress("allJetEnergy", &allJetEnergy, &b_allJetEnergy);
   fChain->SetBranchAddress("allJetEta", &allJetEta, &b_allJetEta);
   fChain->SetBranchAddress("allJetMass", &allJetMass, &b_allJetMass);
   fChain->SetBranchAddress("allJetNPho", &allJetNPho, &b_allJetNPho);
   fChain->SetBranchAddress("allJetPhi", &allJetPhi, &b_allJetPhi);
   fChain->SetBranchAddress("allJetPixNPho", &allJetPixNPho, &b_allJetPixNPho);
   fChain->SetBranchAddress("allJetPixWTime", &allJetPixWTime, &b_allJetPixWTime);
   fChain->SetBranchAddress("allJetPixWTimeRes", &allJetPixWTimeRes, &b_allJetPixWTimeRes);
   fChain->SetBranchAddress("allJetPixWTimeSig", &allJetPixWTimeSig, &b_allJetPixWTimeSig);
   fChain->SetBranchAddress("allJetPt", &allJetPt, &b_allJetPt);
   fChain->SetBranchAddress("allJetSCIndexs", &allJetSCIndexs, &b_allJetSCIndexs);
   fChain->SetBranchAddress("allJetWTime", &allJetWTime, &b_allJetWTime);
   fChain->SetBranchAddress("allJetWTimeRes", &allJetWTimeRes, &b_allJetWTimeRes);
   fChain->SetBranchAddress("allJetWTimeSig", &allJetWTimeSig, &b_allJetWTimeSig);
   fChain->SetBranchAddress("allJetsWithWPixTime", &allJetsWithWPixTime, &b_allJetsWithWPixTime);
   fChain->SetBranchAddress("allJetsWithWTime", &allJetsWithWTime, &b_allJetsWithWTime);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("nSelJets", &nSelJets, &b_nSelJets);
   fChain->SetBranchAddress("nSelJetsWithWPixTime", &nSelJetsWithWPixTime, &b_nSelJetsWithWPixTime);
   fChain->SetBranchAddress("nSelJetsWithWTime", &nSelJetsWithWTime, &b_nSelJetsWithWTime);
   fChain->SetBranchAddress("pv_dDiJetPEnergy", &pv_dDiJetPEnergy, &b_pv_dDiJetPEnergy);
   fChain->SetBranchAddress("pv_dDiJetTime", &pv_dDiJetTime, &b_pv_dDiJetTime);
   fChain->SetBranchAddress("pv_dGJMinDr", &pv_dGJMinDr, &b_pv_dGJMinDr);
   fChain->SetBranchAddress("pv_dGJPEnergy", &pv_dGJPEnergy, &b_pv_dGJPEnergy);
   fChain->SetBranchAddress("pv_dGJTime", &pv_dGJTime, &b_pv_dGJTime);
   fChain->SetBranchAddress("pv_diJet1Time", &pv_diJet1Time, &b_pv_diJet1Time);
   fChain->SetBranchAddress("pv_diJet2Time", &pv_diJet2Time, &b_pv_diJet2Time);
   fChain->SetBranchAddress("pv_gjGammaTime", &pv_gjGammaTime, &b_pv_gjGammaTime);
   fChain->SetBranchAddress("pv_gjJetTime", &pv_gjJetTime, &b_pv_gjJetTime);
   fChain->SetBranchAddress("pv_gjPhoIndex", &pv_gjPhoIndex, &b_pv_gjPhoIndex);
   fChain->SetBranchAddress("pv_haspvtime", &pv_haspvtime, &b_pv_haspvtime);
   fChain->SetBranchAddress("pv_wpenergy", &pv_wpenergy, &b_pv_wpenergy);
   fChain->SetBranchAddress("pv_wphartres", &pv_wphartres, &b_pv_wphartres);
   fChain->SetBranchAddress("pv_wtime", &pv_wtime, &b_pv_wtime);
   fChain->SetBranchAddress("pv_wtimeres", &pv_wtimeres, &b_pv_wtimeres);
   fChain->SetBranchAddress("selGenJetDpt", &selGenJetDpt, &b_selGenJetDpt);
   fChain->SetBranchAddress("selGenJetEnergy", &selGenJetEnergy, &b_selGenJetEnergy);
   fChain->SetBranchAddress("selGenJetImpAng", &selGenJetImpAng, &b_selGenJetImpAng);
   fChain->SetBranchAddress("selGenJetLlpTime", &selGenJetLlpTime, &b_selGenJetLlpTime);
   fChain->SetBranchAddress("selGenJetPt", &selGenJetPt, &b_selGenJetPt);
   fChain->SetBranchAddress("selGenJetTime", &selGenJetTime, &b_selGenJetTime);
   fChain->SetBranchAddress("selGenJetTof", &selGenJetTof, &b_selGenJetTof);
   fChain->SetBranchAddress("selGenJetdr", &selGenJetdr, &b_selGenJetdr);
   fChain->SetBranchAddress("selGenJeteta", &selGenJeteta, &b_selGenJeteta);
   fChain->SetBranchAddress("selJetArea", &selJetArea, &b_selJetArea);
   fChain->SetBranchAddress("selJetChEmEF", &selJetChEmEF, &b_selJetChEmEF);
   fChain->SetBranchAddress("selJetChHM", &selJetChHM, &b_selJetChHM);
   fChain->SetBranchAddress("selJetEnergy", &selJetEnergy, &b_selJetEnergy);
   fChain->SetBranchAddress("selJetEta", &selJetEta, &b_selJetEta);
   fChain->SetBranchAddress("selJetGenLlpId", &selJetGenLlpId, &b_selJetGenLlpId);
   fChain->SetBranchAddress("selJetLlpDp", &selJetLlpDp, &b_selJetLlpDp);
   fChain->SetBranchAddress("selJetLlpDr", &selJetLlpDr, &b_selJetLlpDr);
   fChain->SetBranchAddress("selJetMass", &selJetMass, &b_selJetMass);
   fChain->SetBranchAddress("selJetMuEF", &selJetMuEF, &b_selJetMuEF);
   fChain->SetBranchAddress("selJetNeEmEF", &selJetNeEmEF, &b_selJetNeEmEF);
   fChain->SetBranchAddress("selJetNeHEF", &selJetNeHEF, &b_selJetNeHEF);
   fChain->SetBranchAddress("selJetNeHM", &selJetNeHM, &b_selJetNeHM);
   fChain->SetBranchAddress("selJetPhi", &selJetPhi, &b_selJetPhi);
   fChain->SetBranchAddress("selJetPixWTime", &selJetPixWTime, &b_selJetPixWTime);
   fChain->SetBranchAddress("selJetPixWTimeRes", &selJetPixWTimeRes, &b_selJetPixWTimeRes);
   fChain->SetBranchAddress("selJetPixWTimeSig", &selJetPixWTimeSig, &b_selJetPixWTimeSig);
   fChain->SetBranchAddress("selJetPt", &selJetPt, &b_selJetPt);
   fChain->SetBranchAddress("selJetQrkLlpId", &selJetQrkLlpId, &b_selJetQrkLlpId);
   fChain->SetBranchAddress("selJetQuality", &selJetQuality, &b_selJetQuality);
   fChain->SetBranchAddress("selJetTime", &selJetTime, &b_selJetTime);
   fChain->SetBranchAddress("selJetWTime", &selJetWTime, &b_selJetWTime);
   fChain->SetBranchAddress("selJetWTimeRes", &selJetWTimeRes, &b_selJetWTimeRes);
   fChain->SetBranchAddress("selJetWTimeSig", &selJetWTimeSig, &b_selJetWTimeSig);
   fChain->SetBranchAddress("selJetchHEF", &selJetchHEF, &b_selJetchHEF);
   fChain->SetBranchAddress("muonIsLoose", &muonIsLoose, &b_muonIsLoose);
   fChain->SetBranchAddress("muonIsMedium", &muonIsMedium, &b_muonIsMedium);
   fChain->SetBranchAddress("nLooseMuons", &nLooseMuons, &b_nLooseMuons);
   fChain->SetBranchAddress("nMuonSVMatched", &nMuonSVMatched, &b_nMuonSVMatched);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("HadronicSV_cosTheta", &HadronicSV_cosTheta, &b_HadronicSV_cosTheta);
   fChain->SetBranchAddress("HadronicSV_decayAngle", &HadronicSV_decayAngle, &b_HadronicSV_decayAngle);
   fChain->SetBranchAddress("HadronicSV_dxy", &HadronicSV_dxy, &b_HadronicSV_dxy);
   fChain->SetBranchAddress("HadronicSV_dxySig", &HadronicSV_dxySig, &b_HadronicSV_dxySig);
   fChain->SetBranchAddress("HadronicSV_mass", &HadronicSV_mass, &b_HadronicSV_mass);
   fChain->SetBranchAddress("HadronicSV_massOverNtracks", &HadronicSV_massOverNtracks, &b_HadronicSV_massOverNtracks);
   fChain->SetBranchAddress("HadronicSV_matchRatio", &HadronicSV_matchRatio, &b_HadronicSV_matchRatio);
   fChain->SetBranchAddress("HadronicSV_nTotal", &HadronicSV_nTotal, &b_HadronicSV_nTotal);
   fChain->SetBranchAddress("HadronicSV_nTracks", &HadronicSV_nTracks, &b_HadronicSV_nTracks);
   fChain->SetBranchAddress("HadronicSV_pOverE", &HadronicSV_pOverE, &b_HadronicSV_pOverE);
   fChain->SetBranchAddress("HadronicSV_passMaterialVeto", &HadronicSV_passMaterialVeto, &b_HadronicSV_passMaterialVeto);
   fChain->SetBranchAddress("HadronicSV_tightZWindowSelection", &HadronicSV_tightZWindowSelection, &b_HadronicSV_tightZWindowSelection);
   fChain->SetBranchAddress("HadronicSV_x", &HadronicSV_x, &b_HadronicSV_x);
   fChain->SetBranchAddress("HadronicSV_y", &HadronicSV_y, &b_HadronicSV_y);
   fChain->SetBranchAddress("HadronicSV_z", &HadronicSV_z, &b_HadronicSV_z);
   fChain->SetBranchAddress("LeptonicSV_cosTheta", &LeptonicSV_cosTheta, &b_LeptonicSV_cosTheta);
   fChain->SetBranchAddress("LeptonicSV_decayAngle", &LeptonicSV_decayAngle, &b_LeptonicSV_decayAngle);
   fChain->SetBranchAddress("LeptonicSV_dxy", &LeptonicSV_dxy, &b_LeptonicSV_dxy);
   fChain->SetBranchAddress("LeptonicSV_dxySig", &LeptonicSV_dxySig, &b_LeptonicSV_dxySig);
   fChain->SetBranchAddress("LeptonicSV_electronIndex", &LeptonicSV_electronIndex, &b_LeptonicSV_electronIndex);
   fChain->SetBranchAddress("LeptonicSV_isGold", &LeptonicSV_isGold, &b_LeptonicSV_isGold);
   fChain->SetBranchAddress("LeptonicSV_mass", &LeptonicSV_mass, &b_LeptonicSV_mass);
   fChain->SetBranchAddress("LeptonicSV_muonIndex", &LeptonicSV_muonIndex, &b_LeptonicSV_muonIndex);
   fChain->SetBranchAddress("LeptonicSV_nTotal", &LeptonicSV_nTotal, &b_LeptonicSV_nTotal);
   fChain->SetBranchAddress("LeptonicSV_nTracks", &LeptonicSV_nTracks, &b_LeptonicSV_nTracks);
   fChain->SetBranchAddress("LeptonicSV_pOverE", &LeptonicSV_pOverE, &b_LeptonicSV_pOverE);
   fChain->SetBranchAddress("LeptonicSV_passElectronID", &LeptonicSV_passElectronID, &b_LeptonicSV_passElectronID);
   fChain->SetBranchAddress("LeptonicSV_passMaterialVeto", &LeptonicSV_passMaterialVeto, &b_LeptonicSV_passMaterialVeto);
   fChain->SetBranchAddress("LeptonicSV_passMuonID", &LeptonicSV_passMuonID, &b_LeptonicSV_passMuonID);
   fChain->SetBranchAddress("LeptonicSV_tightZWindowSelection", &LeptonicSV_tightZWindowSelection, &b_LeptonicSV_tightZWindowSelection);
   fChain->SetBranchAddress("LeptonicSV_x", &LeptonicSV_x, &b_LeptonicSV_x);
   fChain->SetBranchAddress("LeptonicSV_y", &LeptonicSV_y, &b_LeptonicSV_y);
   fChain->SetBranchAddress("LeptonicSV_z", &LeptonicSV_z, &b_LeptonicSV_z);
   fChain->SetBranchAddress("SV_nElectron", &SV_nElectron, &b_SV_nElectron);
   fChain->SetBranchAddress("SV_nElectronPassZWindow", &SV_nElectronPassZWindow, &b_SV_nElectronPassZWindow);
   fChain->SetBranchAddress("SV_nHadCmpDxyGeq2", &SV_nHadCmpDxyGeq2, &b_SV_nHadCmpDxyGeq2);
   fChain->SetBranchAddress("SV_nHadCmpDxySigGeq1000", &SV_nHadCmpDxySigGeq1000, &b_SV_nHadCmpDxySigGeq1000);
   fChain->SetBranchAddress("SV_nHadStdDxyGeq2", &SV_nHadStdDxyGeq2, &b_SV_nHadStdDxyGeq2);
   fChain->SetBranchAddress("SV_nHadStdDxySigGeq1000", &SV_nHadStdDxySigGeq1000, &b_SV_nHadStdDxySigGeq1000);
   fChain->SetBranchAddress("SV_nHadronPassZWindow", &SV_nHadronPassZWindow, &b_SV_nHadronPassZWindow);
   fChain->SetBranchAddress("SV_nHadronic", &SV_nHadronic, &b_SV_nHadronic);
   fChain->SetBranchAddress("SV_nLepCmpDxyGeq2", &SV_nLepCmpDxyGeq2, &b_SV_nLepCmpDxyGeq2);
   fChain->SetBranchAddress("SV_nLepCmpDxySigGeq1000", &SV_nLepCmpDxySigGeq1000, &b_SV_nLepCmpDxySigGeq1000);
   fChain->SetBranchAddress("SV_nLepStdDxyGeq2", &SV_nLepStdDxyGeq2, &b_SV_nLepStdDxyGeq2);
   fChain->SetBranchAddress("SV_nLepStdDxySigGeq1000", &SV_nLepStdDxySigGeq1000, &b_SV_nLepStdDxySigGeq1000);
   fChain->SetBranchAddress("SV_nLeptonPassZWindow", &SV_nLeptonPassZWindow, &b_SV_nLeptonPassZWindow);
   fChain->SetBranchAddress("SV_nLeptonic", &SV_nLeptonic, &b_SV_nLeptonic);
   fChain->SetBranchAddress("SV_nMuon", &SV_nMuon, &b_SV_nMuon);
   fChain->SetBranchAddress("SV_nMuonPassZWindow", &SV_nMuonPassZWindow, &b_SV_nMuonPassZWindow);
   fChain->SetBranchAddress("passGe1HadAndLepSelectionCRCmp", &passGe1HadAndLepSelectionCRCmp, &b_passGe1HadAndLepSelectionCRCmp);
   fChain->SetBranchAddress("passGe1HadAndLepSelectionCRLoose", &passGe1HadAndLepSelectionCRLoose, &b_passGe1HadAndLepSelectionCRLoose);
   fChain->SetBranchAddress("passGe1HadAndLepSelectionCRTight", &passGe1HadAndLepSelectionCRTight, &b_passGe1HadAndLepSelectionCRTight);
   fChain->SetBranchAddress("passGe1HadAndLepSelectionSRCmp", &passGe1HadAndLepSelectionSRCmp, &b_passGe1HadAndLepSelectionSRCmp);
   fChain->SetBranchAddress("passGe1HadAndLepSelectionSRLoose", &passGe1HadAndLepSelectionSRLoose, &b_passGe1HadAndLepSelectionSRLoose);
   fChain->SetBranchAddress("passGe1HadAndLepSelectionSRTight", &passGe1HadAndLepSelectionSRTight, &b_passGe1HadAndLepSelectionSRTight);
   fChain->SetBranchAddress("passNHad1SelectionCRCmp", &passNHad1SelectionCRCmp, &b_passNHad1SelectionCRCmp);
   fChain->SetBranchAddress("passNHad1SelectionCRLoose", &passNHad1SelectionCRLoose, &b_passNHad1SelectionCRLoose);
   fChain->SetBranchAddress("passNHad1SelectionCRTight", &passNHad1SelectionCRTight, &b_passNHad1SelectionCRTight);
   fChain->SetBranchAddress("passNHad1SelectionSRCmp", &passNHad1SelectionSRCmp, &b_passNHad1SelectionSRCmp);
   fChain->SetBranchAddress("passNHad1SelectionSRLoose", &passNHad1SelectionSRLoose, &b_passNHad1SelectionSRLoose);
   fChain->SetBranchAddress("passNHad1SelectionSRTight", &passNHad1SelectionSRTight, &b_passNHad1SelectionSRTight);
   fChain->SetBranchAddress("passNHadGe2SelectionCRCmp", &passNHadGe2SelectionCRCmp, &b_passNHadGe2SelectionCRCmp);
   fChain->SetBranchAddress("passNHadGe2SelectionCRLoose", &passNHadGe2SelectionCRLoose, &b_passNHadGe2SelectionCRLoose);
   fChain->SetBranchAddress("passNHadGe2SelectionCRTight", &passNHadGe2SelectionCRTight, &b_passNHadGe2SelectionCRTight);
   fChain->SetBranchAddress("passNHadGe2SelectionSRCmp", &passNHadGe2SelectionSRCmp, &b_passNHadGe2SelectionSRCmp);
   fChain->SetBranchAddress("passNHadGe2SelectionSRLoose", &passNHadGe2SelectionSRLoose, &b_passNHadGe2SelectionSRLoose);
   fChain->SetBranchAddress("passNHadGe2SelectionSRTight", &passNHadGe2SelectionSRTight, &b_passNHadGe2SelectionSRTight);
   fChain->SetBranchAddress("passNLep1SelectionCRCmp", &passNLep1SelectionCRCmp, &b_passNLep1SelectionCRCmp);
   fChain->SetBranchAddress("passNLep1SelectionCRLoose", &passNLep1SelectionCRLoose, &b_passNLep1SelectionCRLoose);
   fChain->SetBranchAddress("passNLep1SelectionCRTight", &passNLep1SelectionCRTight, &b_passNLep1SelectionCRTight);
   fChain->SetBranchAddress("passNLep1SelectionSRCmp", &passNLep1SelectionSRCmp, &b_passNLep1SelectionSRCmp);
   fChain->SetBranchAddress("passNLep1SelectionSRLoose", &passNLep1SelectionSRLoose, &b_passNLep1SelectionSRLoose);
   fChain->SetBranchAddress("passNLep1SelectionSRTight", &passNLep1SelectionSRTight, &b_passNLep1SelectionSRTight);
   fChain->SetBranchAddress("passNLepGe2SelectionCRCmp", &passNLepGe2SelectionCRCmp, &b_passNLepGe2SelectionCRCmp);
   fChain->SetBranchAddress("passNLepGe2SelectionCRLoose", &passNLepGe2SelectionCRLoose, &b_passNLepGe2SelectionCRLoose);
   fChain->SetBranchAddress("passNLepGe2SelectionCRTight", &passNLepGe2SelectionCRTight, &b_passNLepGe2SelectionCRTight);
   fChain->SetBranchAddress("passNLepGe2SelectionSRCmp", &passNLepGe2SelectionSRCmp, &b_passNLepGe2SelectionSRCmp);
   fChain->SetBranchAddress("passNLepGe2SelectionSRLoose", &passNLepGe2SelectionSRLoose, &b_passNLepGe2SelectionSRLoose);
   fChain->SetBranchAddress("passNLepGe2SelectionSRTight", &passNLepGe2SelectionSRTight, &b_passNLepGe2SelectionSRTight);
   Notify();
}

Bool_t ::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef _cxx
