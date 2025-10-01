//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep  7 12:52:33 2025 by ROOT version 6.26/11
// from TTree kuSkimTree/output root file for kUCMSSkimmer
// found on file: /uscms/home/jaking/nobackup/llpana_skims/GJets_R18_SVIPM100_v30_GJets_HT-200To400_rjrskim_v36.root
//////////////////////////////////////////////////////////

#ifndef skimtuple_v36_h
#define skimtuple_v36_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"

//#define sbDEBUG true
#define sbDEBUG false

class skimtuple_v36 {
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
   Bool_t          Flag_MetFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_hfNoisyHitsFilter;

   Float_t         PVx;
   Float_t         PVy;
   Float_t         PVz;

   std::vector<float>   *SV_cosTheta;
   std::vector<float>   *SV_decayAngle;
   std::vector<float>   *SV_dxy;
   std::vector<float>   *SV_dxySig;
   std::vector<float>   *SV_mass;
   std::vector<float>   *SV_massOverNtracks;
   Int_t           SV_nElectron;
   Int_t           SV_nHadronic;
   Int_t           SV_nLeptonic;
   Int_t           SV_nMuon;
   std::vector<int>     *SV_nTracks;
   std::vector<float>   *SV_pOverE;
   std::string          *DataSetKey;

   Float_t         evtFillWgt;
   Float_t         evtGenWgt;
   Float_t         evtXSection;
   Float_t         Evt_genWgt;

   Float_t         selCMet;
   Float_t         selCMetPx;
   Float_t         selCMetPy;
   Float_t         selMet;
   Float_t         selMetPx;
   Float_t         selMetPy;

   std::vector<int>     *genCharge;
   std::vector<float>   *genGrvtinoMass;
   std::vector<float>   *genLLPMass;
   std::vector<float>   *genLSPMass;
   std::vector<float>   *genMass;
   std::vector<float>   *genMomCTau;
   std::vector<float>   *genPartEnergy;
   std::vector<float>   *genPartEta;
   std::vector<int>     *genPartMomIdx;
   std::vector<unsigned int> *genPartPdgId;
   std::vector<float>   *genPartPhi;
   std::vector<float>   *genPartPt;
   std::vector<int>     *genPartSusId;
   std::vector<float>   *genPx;
   std::vector<float>   *genPy;
   std::vector<float>   *genPz;
   std::vector<float>   *genSGMass;
   std::vector<float>   *genSQMass;
   Bool_t          genSTFlagQP;
   Bool_t          genSTFlagQQP;
   Bool_t          genSTFlagQQZ;
   Bool_t          genSTFlagQQZL;
   Bool_t          genSTFlagQQZLL;
   Int_t           genSigType;
   std::vector<int>     *genStatus;
   std::vector<float>   *genVx;
   std::vector<float>   *genVy;
   std::vector<float>   *genVz;
   std::vector<float>   *genXMomCTau;

   std::vector<bool>    *eleHasSVMatch;
   std::vector<bool>    *eleIsLoose;
   std::vector<float>   *elePhoMinDr;
   std::vector<float>   *epDrMatch;
   std::vector<float>   *epSeedIdMatch;
   std::vector<float>   *epSeedIdUnique;
   Int_t           nEleSVMatched;
   UInt_t          nElectrons;
   UInt_t          nSelElectrons;
   UInt_t          nSelIsoElectrons;

   UInt_t          nPhotons;
   UInt_t          nSelPhotons;
   std::vector<float>   *selPhoClstrRn;
   std::vector<float>   *selPhoCorEnergy;
   std::vector<float>   *selPhoCorPt;
   std::vector<float>   *selPhoCovEtaEta;
   std::vector<float>   *selPhoCovEtaPhi;
   std::vector<float>   *selPhoCovPhiPhi;
   std::vector<float>   *selPhoEcalPFClusterIso;
   std::vector<float>   *selPhoEcalRHSumEtConeDR04;
   UInt_t          selPhoEleVeto;
   std::vector<float>   *selPhoEnergy;
   std::vector<float>   *selPhoEta;
   std::vector<float>   *selPhoEtaWidth;
   std::vector<int>     *selPhoGenIdx;
   std::vector<float>   *selPhoGenPt;
   std::vector<float>   *selPhoGenSigMomEnergy;
   std::vector<float>   *selPhoGenSigMomEta;
   std::vector<float>   *selPhoGenSigMomMass;
   std::vector<float>   *selPhoGenSigMomPhi;
   std::vector<float>   *selPhoGenSigMomPt;
   std::vector<float>   *selPhoGenSigMomPx;
   std::vector<float>   *selPhoGenSigMomPy;
   std::vector<float>   *selPhoGenSigMomPz;
   std::vector<float>   *selPhoGenSigMomVx;
   std::vector<float>   *selPhoGenSigMomVy;
   std::vector<float>   *selPhoGenSigMomVz;
   std::vector<float>   *selPhoHadOverEM;
   std::vector<float>   *selPhoHadTowOverEM;
   std::vector<bool>    *selPhoHasConversionTracks;
   std::vector<float>   *selPhoHcalPFClusterIso;
   std::vector<float>   *selPhoHcalTowerSumEtBcConeDR04;
   std::vector<float>   *selPhoHcalTowerSumEtConeDR04;
   std::vector<float>   *selPhoHoe_PUcorr;
   std::vector<float>   *selPhoMBetaEql;
   std::vector<float>   *selPhoMBetaPmt;
   std::vector<float>   *selPhoNTrkHollowConeDR04;
   std::vector<float>   *selPhoNTrkSolidConeDR04;
   std::vector<unsigned int> *selPhoNrh;
   std::vector<bool>    *selPhoOOT;
   std::vector<float>   *selPhoPVtoRH;
   std::vector<float>   *selPhoPfChargedIso;
   std::vector<float>   *selPhoPfChargedIsoPFPV;
   std::vector<float>   *selPhoPfChargedIsoWorstVtx;
   std::vector<float>   *selPhoPfPhoIso03;
   std::vector<float>   *selPhoPfRelIso03_all_quadratic;
   std::vector<float>   *selPhoPfRelIso03_chg_quadratic;
   std::vector<float>   *selPhoPhi;
   std::vector<float>   *selPhoPhiWidth;
   std::vector<float>   *selPhoPhoIsoDr;
   std::vector<bool>    *selPhoPixelSeed;
   std::vector<float>   *selPhoPt;
   std::vector<int>     *selPhoQuality;
   std::vector<float>   *selPhoR9;
   std::vector<float>   *selPhoS4;
   std::vector<float>   *selPhoSAlp;
   std::vector<float>   *selPhoSCx;
   std::vector<float>   *selPhoSCy;
   std::vector<float>   *selPhoSCz;
   std::vector<float>   *selPhoSMaj;
   std::vector<float>   *selPhoSMin;
   std::vector<float>   *selPhoSieie;
   std::vector<float>   *selPhoSieip;
   std::vector<float>   *selPhoSigmaIEtaIEta;
   std::vector<float>   *selPhoSipip;
   std::vector<float>   *selPhoSusyId;
   std::vector<float>   *selPhoTime;
   std::vector<float>   *selPhoTrkSumPtHollowConeDR03;
   std::vector<float>   *selPhoTrkSumPtHollowConeDR04;
   std::vector<float>   *selPhoTrkSumPtSolidConeDR04;
   std::vector<float>   *selPhoSTimeSig;
   std::vector<float>   *selPhoLTimeSig;
   std::vector<float>   *selPhoWTimeSig;

   std::vector<bool>    *muonIsLoose;
   std::vector<bool>    *muonIsMedium;
   UInt_t          nLooseMuons;
   Int_t           nMuonSVMatched;
   UInt_t          nMuons;

   UInt_t          nJets;
   UInt_t          nSelJets;
   std::vector<float>   *selGenJetDpt;
   std::vector<float>   *selGenJetEnergy;
   std::vector<float>   *selGenJetImpAng;
   std::vector<float>   *selGenJetLlpTime;
   std::vector<float>   *selGenJetPt;
   std::vector<float>   *selGenJetTime;
   std::vector<float>   *selGenJetTof;
   std::vector<float>   *selGenJetdr;
   std::vector<float>   *selGenJeteta;
   std::vector<float>   *selJetArea;
   std::vector<float>   *selJetChEmEF;
   std::vector<float>   *selJetChHM;
   std::vector<float>   *selJetEnergy;
   std::vector<float>   *selJetEta;
   std::vector<int>     *selJetGenLlpId;
   std::vector<float>   *selJetLlpDp;
   std::vector<float>   *selJetLlpDr;
   std::vector<float>   *selJetMass;
   std::vector<float>   *selJetMuEF;
   std::vector<float>   *selJetNeEmEF;
   std::vector<float>   *selJetNeHEF;
   std::vector<float>   *selJetNeHM;
   std::vector<float>   *selJetPhi;
   std::vector<float>   *selJetPt;
   std::vector<int>     *selJetQrkLlpId;
   std::vector<int>     *selJetQuality;
   std::vector<float>   *selJetTime;
   std::vector<float>   *selJetchHEF;

   std::vector<bool>    *rjrABSide;
   std::vector<float>   *rjr_Rdiff;
   std::vector<float>   *rjrDphiMETV;
   std::vector<float>   *rjrEVa;
   std::vector<float>   *rjrEVb;
   std::vector<float>   *rjrMET;
   std::vector<float>   *rjrN2Px;
   std::vector<float>   *rjrN2Py;
   std::vector<float>   *rjrN2Pz;
   std::vector<int>     *rjrNJets;
   std::vector<int>     *rjrNJetsJa;
   std::vector<int>     *rjrNJetsJb;
   std::vector<int>     *rjrNPhotons;
   std::vector<int>     *rjrNVisObjects;
   std::vector<float>   *rjrPTS;
   std::vector<float>   *rjrPVa;
   std::vector<float>   *rjrPVb;
   std::vector<float>   *rjrPVlab;
   std::vector<float>   *rjrPZS;
   std::vector<float>   *rjrSCosA;
   std::vector<float>   *rjrSMass;
   std::vector<float>   *rjrSdphiBV;
   std::vector<float>   *rjrSdphiDA;
   std::vector<int>     *rjrSubPhoLocation;
   std::vector<int>     *rjrType;
   std::vector<bool>    *rjrVisJetSide;
   std::vector<bool>    *rjrVisPhoSide;
   std::vector<float>   *rjrX1aCosA;
   std::vector<float>   *rjrX1aMass;
   std::vector<float>   *rjrX1a_Eta;
   std::vector<float>   *rjrX1a_Phi;
   std::vector<float>   *rjrX1a_Pt;
   std::vector<float>   *rjrX1bCosA;
   std::vector<float>   *rjrX1bMass;
   std::vector<float>   *rjrX1b_Eta;
   std::vector<float>   *rjrX1b_Phi;
   std::vector<float>   *rjrX1b_Pt;
   std::vector<float>   *rjrX2aCosA;
   std::vector<float>   *rjrX2aMass;
   std::vector<float>   *rjrX2aPs;
   std::vector<float>   *rjrX2bCosA;
   std::vector<float>   *rjrX2bMass;
   std::vector<float>   *rjrX2bPs;
   std::vector<float>   *rjr_Mr;
   std::vector<float>   *rjr_Ms;
   std::vector<float>   *rjr_Mva;
   std::vector<float>   *rjr_Mvb;
   std::vector<float>   *rjr_R;
   std::vector<float>   *rjr_Rm;
   std::vector<float>   *rjr_Rs;
   std::vector<float>   *rjr_Rv;
   std::vector<float>   *rjr_Rx;
   std::vector<float>   *rjr_Rxa;
   std::vector<float>   *rjr_Rxb;
   std::vector<float>   *rjr_p1Pts;
   std::vector<float>   *rjr_p1Ptxa;
   std::vector<float>   *rjr_p1Ptxb;
   std::vector<float>   *rjr_p1Pxa;
   std::vector<float>   *rjr_p1Pxb;
   std::vector<float>   *rjr_p2Pts;
   std::vector<float>   *rjr_p2Ptxa;
   std::vector<float>   *rjr_p2Ptxb;
   std::vector<float>   *rjr_p2Pxa;
   std::vector<float>   *rjr_p2Pxb;
   std::vector<float>   *rjr_pHs11;
   std::vector<float>   *rjr_pHs11a;
   std::vector<float>   *rjr_pHs11b;
   std::vector<float>   *rjr_pHs20;
   std::vector<float>   *rjr_pHs21;
   std::vector<float>   *rjr_pHs21a;
   std::vector<float>   *rjr_pHs21b;
   std::vector<float>   *rjr_pHs22;
   std::vector<float>   *rjr_pHs40;
   std::vector<float>   *rjr_pHs41;
   std::vector<float>   *rjr_pHts11;
   std::vector<float>   *rjr_pHts11a;
   std::vector<float>   *rjr_pHts11b;
   std::vector<float>   *rjr_pHts20;
   std::vector<float>   *rjr_pHts21;
   std::vector<float>   *rjr_pHts21a;
   std::vector<float>   *rjr_pHts21b;
   std::vector<float>   *rjr_pHts22;
   std::vector<float>   *rjr_pHts40;
   std::vector<float>   *rjr_pHts41;
   std::vector<float>   *rjr_pHtxa10;
   std::vector<float>   *rjr_pHtxa11;
   std::vector<float>   *rjr_pHtxa20;
   std::vector<float>   *rjr_pHtxa21;
   std::vector<float>   *rjr_pHtxb10;
   std::vector<float>   *rjr_pHtxb11;
   std::vector<float>   *rjr_pHtxb20;
   std::vector<float>   *rjr_pHtxb21;
   std::vector<float>   *rjr_pHxa10;
   std::vector<float>   *rjr_pHxa11;
   std::vector<float>   *rjr_pHxa20;
   std::vector<float>   *rjr_pHxa21;
   std::vector<float>   *rjr_pHxb10;
   std::vector<float>   *rjr_pHxb11;
   std::vector<float>   *rjr_pHxb20;
   std::vector<float>   *rjr_pHxb21;

   // List of branches
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_MetFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!

   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!

   TBranch        *b_SV_cosTheta;   //!
   TBranch        *b_SV_decayAngle;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_massOverNtracks;   //!
   TBranch        *b_SV_nElectron;   //!
   TBranch        *b_SV_nHadronic;   //!
   TBranch        *b_SV_nLeptonic;   //!
   TBranch        *b_SV_nMuon;   //!
   TBranch        *b_SV_nTracks;   //!
   TBranch        *b_SV_pOverE;   //!

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

   TBranch        *b_nPhotons;   //!
   TBranch        *b_nSelPhotons;   //!
   TBranch        *b_selPhoClstrRn;   //!
   TBranch        *b_selPhoCorEnergy;   //!
   TBranch        *b_selPhoCorPt;   //!
   TBranch        *b_selPhoCovEtaEta;   //!
   TBranch        *b_selPhoCovEtaPhi;   //!
   TBranch        *b_selPhoCovPhiPhi;   //!
   TBranch        *b_selPhoEcalPFClusterIso;   //!
   TBranch        *b_selPhoEcalRHSumEtConeDR04;   //!
   TBranch        *b_selPhoEleVeto;   //!
   TBranch        *b_selPhoEnergy;   //!
   TBranch        *b_selPhoEta;   //!
   TBranch        *b_selPhoEtaWidth;   //!
   TBranch        *b_selPhoGenIdx;   //!
   TBranch        *b_selPhoGenPt;   //!
   TBranch        *b_selPhoGenSigMomEnergy;   //!
   TBranch        *b_selPhoGenSigMomEta;   //!
   TBranch        *b_selPhoGenSigMomMass;   //!
   TBranch        *b_selPhoGenSigMomPhi;   //!
   TBranch        *b_selPhoGenSigMomPt;   //!
   TBranch        *b_selPhoGenSigMomPx;   //!
   TBranch        *b_selPhoGenSigMomPy;   //!
   TBranch        *b_selPhoGenSigMomPz;   //!
   TBranch        *b_selPhoGenSigMomVx;   //!
   TBranch        *b_selPhoGenSigMomVy;   //!
   TBranch        *b_selPhoGenSigMomVz;   //!
   TBranch        *b_selPhoHadOverEM;   //!
   TBranch        *b_selPhoHadTowOverEM;   //!
   TBranch        *b_selPhoHasConversionTracks;   //!
   TBranch        *b_selPhoHcalPFClusterIso;   //!
   TBranch        *b_selPhoHcalTowerSumEtBcConeDR04;   //!
   TBranch        *b_selPhoHcalTowerSumEtConeDR04;   //!
   TBranch        *b_selPhoHoe_PUcorr;   //!
   TBranch        *b_selPhoMBetaEql;   //!
   TBranch        *b_selPhoMBetaPmt;   //!
   TBranch        *b_selPhoNTrkHollowConeDR04;   //!
   TBranch        *b_selPhoNTrkSolidConeDR04;   //!
   TBranch        *b_selPhoNrh;   //!
   TBranch        *b_selPhoOOT;   //!
   TBranch        *b_selPhoPVtoRH;   //!
   TBranch        *b_selPhoPfChargedIso;   //!
   TBranch        *b_selPhoPfChargedIsoPFPV;   //!
   TBranch        *b_selPhoPfChargedIsoWorstVtx;   //!
   TBranch        *b_selPhoPfPhoIso03;   //!
   TBranch        *b_selPhoPfRelIso03_all_quadratic;   //!
   TBranch        *b_selPhoPfRelIso03_chg_quadratic;   //!
   TBranch        *b_selPhoPhi;   //!
   TBranch        *b_selPhoPhiWidth;   //!
   TBranch        *b_selPhoPhoIsoDr;   //!
   TBranch        *b_selPhoPixelSeed;   //!
   TBranch        *b_selPhoPt;   //!
   TBranch        *b_selPhoQuality;   //!
   TBranch        *b_selPhoR9;   //!
   TBranch        *b_selPhoS4;   //!
   TBranch        *b_selPhoSAlp;   //!
   TBranch        *b_selPhoSCx;   //!
   TBranch        *b_selPhoSCy;   //!
   TBranch        *b_selPhoSCz;   //!
   TBranch        *b_selPhoSMaj;   //!
   TBranch        *b_selPhoSMin;   //!
   TBranch        *b_selPhoSieie;   //!
   TBranch        *b_selPhoSieip;   //!
   TBranch        *b_selPhoSigmaIEtaIEta;   //!
   TBranch        *b_selPhoSipip;   //!
   TBranch        *b_selPhoSusyId;   //!
   TBranch        *b_selPhoTime;   //!
   TBranch        *b_selPhoTrkSumPtHollowConeDR03;   //!
   TBranch        *b_selPhoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_selPhoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_selPhoSTimeSig;   //!
   TBranch        *b_selPhoLTimeSig;   //!
   TBranch        *b_selPhoWTimeSig;   //!

   TBranch        *b_muonIsLoose;   //!
   TBranch        *b_muonIsMedium;   //!
   TBranch        *b_nLooseMuons;   //!
   TBranch        *b_nMuonSVMatched;   //!
   TBranch        *b_nMuons;   //!

   TBranch        *b_nJets;   //!
   TBranch        *b_nSelJets;   //!
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
   TBranch        *b_selJetPt;   //!
   TBranch        *b_selJetQrkLlpId;   //!
   TBranch        *b_selJetQuality;   //!
   TBranch        *b_selJetTime;   //!
   TBranch        *b_selJetchHEF;   //!

   TBranch        *b_rjrABSide;   //!
   TBranch        *b_rjr_Rdiff;   //!
   TBranch        *b_rjrDphiMETV;   //!
   TBranch        *b_rjrEVa;   //!
   TBranch        *b_rjrEVb;   //!
   TBranch        *b_rjrMET;   //!
   TBranch        *b_rjrN2Px;   //!
   TBranch        *b_rjrN2Py;   //!
   TBranch        *b_rjrN2Pz;   //!
   TBranch        *b_rjrNJets;   //!
   TBranch        *b_rjrNJetsJa;   //!
   TBranch        *b_rjrNJetsJb;   //!
   TBranch        *b_rjrNPhotons;   //!
   TBranch        *b_rjrNVisObjects;   //!
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

   //skimtuple_v36(TTree *tree=0);
   //virtual ~skimtuple_v36();
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
skimtuple_v36::skimtuple_v36(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms/home/jaking/nobackup/llpana_skims/GJets_R18_SVIPM100_v30_GJets_HT-200To400_rjrskim_v36.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uscms/home/jaking/nobackup/llpana_skims/GJets_R18_SVIPM100_v30_GJets_HT-200To400_rjrskim_v36.root");
      }
      f->GetObject("kuSkimTree",tree);

   }
   Init(tree);
}

skimtuple_v36::~skimtuple_v36()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimtuple_v36::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimtuple_v36::LoadTree(Long64_t entry)
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

void skimtuple_v36::Init( TChain *tree, bool doGenInfo ){
//void skimtuple_v36::Init(TTree *tree){

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   Evt_genWgt = 1;

   // Set object pointer
   SV_cosTheta = 0;
   SV_decayAngle = 0;
   SV_dxy = 0;
   SV_dxySig = 0;
   SV_mass = 0;
   SV_massOverNtracks = 0;
   SV_nTracks = 0;
   SV_pOverE = 0;
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
   eleHasSVMatch = 0;
   eleIsLoose = 0;
   elePhoMinDr = 0;
   epDrMatch = 0;
   epSeedIdMatch = 0;
   epSeedIdUnique = 0;
   selPhoClstrRn = 0;
   selPhoCorEnergy = 0;
   selPhoCorPt = 0;
   selPhoCovEtaEta = 0;
   selPhoCovEtaPhi = 0;
   selPhoCovPhiPhi = 0;
   selPhoEcalPFClusterIso = 0;
   selPhoEcalRHSumEtConeDR04 = 0;
   selPhoEnergy = 0;
   selPhoEta = 0;
   selPhoEtaWidth = 0;
   selPhoGenIdx = 0;
   selPhoGenPt = 0;
   selPhoGenSigMomEnergy = 0;
   selPhoGenSigMomEta = 0;
   selPhoGenSigMomMass = 0;
   selPhoGenSigMomPhi = 0;
   selPhoGenSigMomPt = 0;
   selPhoGenSigMomPx = 0;
   selPhoGenSigMomPy = 0;
   selPhoGenSigMomPz = 0;
   selPhoGenSigMomVx = 0;
   selPhoGenSigMomVy = 0;
   selPhoGenSigMomVz = 0;
   selPhoHadOverEM = 0;
   selPhoHadTowOverEM = 0;
   selPhoHasConversionTracks = 0;
   selPhoHcalPFClusterIso = 0;
   selPhoHcalTowerSumEtBcConeDR04 = 0;
   selPhoHcalTowerSumEtConeDR04 = 0;
   selPhoHoe_PUcorr = 0;
   selPhoMBetaEql = 0;
   selPhoMBetaPmt = 0;
   selPhoNTrkHollowConeDR04 = 0;
   selPhoNTrkSolidConeDR04 = 0;
   selPhoNrh = 0;
   selPhoOOT = 0;
   selPhoPVtoRH = 0;
   selPhoPfChargedIso = 0;
   selPhoPfChargedIsoPFPV = 0;
   selPhoPfChargedIsoWorstVtx = 0;
   selPhoPfPhoIso03 = 0;
   selPhoPfRelIso03_all_quadratic = 0;
   selPhoPfRelIso03_chg_quadratic = 0;
   selPhoPhi = 0;
   selPhoPhiWidth = 0;
   selPhoPhoIsoDr = 0;
   selPhoPixelSeed = 0;
   selPhoPt = 0;
   selPhoQuality = 0;
   selPhoR9 = 0;
   selPhoS4 = 0;
   selPhoSAlp = 0;
   selPhoSCx = 0;
   selPhoSCy = 0;
   selPhoSCz = 0;
   selPhoSMaj = 0;
   selPhoSMin = 0;
   selPhoSieie = 0;
   selPhoSieip = 0;
   selPhoSigmaIEtaIEta = 0;
   selPhoSipip = 0;
   selPhoSusyId = 0;
   selPhoTime = 0;
   selPhoSTimeSig = 0;
   selPhoWTimeSig = 0;
   selPhoLTimeSig = 0;
   selPhoTrkSumPtHollowConeDR03 = 0;
   selPhoTrkSumPtHollowConeDR04 = 0;
   selPhoTrkSumPtSolidConeDR04 = 0;
   muonIsLoose = 0;
   muonIsMedium = 0;
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
   selJetPt = 0;
   selJetQrkLlpId = 0;
   selJetQuality = 0;
   selJetTime = 0;
   selJetchHEF = 0;
   rjrABSide = 0;
   rjr_Rdiff = 0;
   rjrDphiMETV = 0;
   rjrEVa = 0;
   rjrEVb = 0;
   rjrMET = 0;
   rjrN2Px = 0;
   rjrN2Py = 0;
   rjrN2Pz = 0;
   rjrNJets = 0;
   rjrNJetsJa = 0;
   rjrNJetsJb = 0;
   rjrNPhotons = 0;
   rjrNVisObjects = 0;
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

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   //fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_MetFilter", &Flag_MetFilter, &b_Flag_MetFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);

   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);

/*
   fChain->SetBranchAddress("SV_cosTheta", &SV_cosTheta, &b_SV_cosTheta);
   fChain->SetBranchAddress("SV_decayAngle", &SV_decayAngle, &b_SV_decayAngle);
   fChain->SetBranchAddress("SV_dxy", &SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", &SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_mass", &SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_massOverNtracks", &SV_massOverNtracks, &b_SV_massOverNtracks);
   fChain->SetBranchAddress("SV_nElectron", &SV_nElectron, &b_SV_nElectron);
   fChain->SetBranchAddress("SV_nHadronic", &SV_nHadronic, &b_SV_nHadronic);
   fChain->SetBranchAddress("SV_nLeptonic", &SV_nLeptonic, &b_SV_nLeptonic);
   fChain->SetBranchAddress("SV_nMuon", &SV_nMuon, &b_SV_nMuon);
   fChain->SetBranchAddress("SV_nTracks", &SV_nTracks, &b_SV_nTracks);
   fChain->SetBranchAddress("SV_pOverE", &SV_pOverE, &b_SV_pOverE);
*/

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

   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("nSelPhotons", &nSelPhotons, &b_nSelPhotons);
   fChain->SetBranchAddress("selPhoClstrRn", &selPhoClstrRn, &b_selPhoClstrRn);
   fChain->SetBranchAddress("selPhoCorEnergy", &selPhoCorEnergy, &b_selPhoCorEnergy);
   fChain->SetBranchAddress("selPhoCorPt", &selPhoCorPt, &b_selPhoCorPt);
   fChain->SetBranchAddress("selPhoCovEtaEta", &selPhoCovEtaEta, &b_selPhoCovEtaEta);
   fChain->SetBranchAddress("selPhoCovEtaPhi", &selPhoCovEtaPhi, &b_selPhoCovEtaPhi);
   fChain->SetBranchAddress("selPhoCovPhiPhi", &selPhoCovPhiPhi, &b_selPhoCovPhiPhi);
   fChain->SetBranchAddress("selPhoEcalPFClusterIso", &selPhoEcalPFClusterIso, &b_selPhoEcalPFClusterIso);
   fChain->SetBranchAddress("selPhoEcalRHSumEtConeDR04", &selPhoEcalRHSumEtConeDR04, &b_selPhoEcalRHSumEtConeDR04);
   fChain->SetBranchAddress("selPhoEleVeto", &selPhoEleVeto, &b_selPhoEleVeto);
   fChain->SetBranchAddress("selPhoEnergy", &selPhoEnergy, &b_selPhoEnergy);
   fChain->SetBranchAddress("selPhoEta", &selPhoEta, &b_selPhoEta);
   fChain->SetBranchAddress("selPhoEtaWidth", &selPhoEtaWidth, &b_selPhoEtaWidth);
   fChain->SetBranchAddress("selPhoGenIdx", &selPhoGenIdx, &b_selPhoGenIdx);
   fChain->SetBranchAddress("selPhoGenPt", &selPhoGenPt, &b_selPhoGenPt);
   fChain->SetBranchAddress("selPhoGenSigMomEnergy", &selPhoGenSigMomEnergy, &b_selPhoGenSigMomEnergy);
   fChain->SetBranchAddress("selPhoGenSigMomEta", &selPhoGenSigMomEta, &b_selPhoGenSigMomEta);
   fChain->SetBranchAddress("selPhoGenSigMomMass", &selPhoGenSigMomMass, &b_selPhoGenSigMomMass);
   fChain->SetBranchAddress("selPhoGenSigMomPhi", &selPhoGenSigMomPhi, &b_selPhoGenSigMomPhi);
   fChain->SetBranchAddress("selPhoGenSigMomPt", &selPhoGenSigMomPt, &b_selPhoGenSigMomPt);
   fChain->SetBranchAddress("selPhoGenSigMomPx", &selPhoGenSigMomPx, &b_selPhoGenSigMomPx);
   fChain->SetBranchAddress("selPhoGenSigMomPy", &selPhoGenSigMomPy, &b_selPhoGenSigMomPy);
   fChain->SetBranchAddress("selPhoGenSigMomPz", &selPhoGenSigMomPz, &b_selPhoGenSigMomPz);
   fChain->SetBranchAddress("selPhoGenSigMomVx", &selPhoGenSigMomVx, &b_selPhoGenSigMomVx);
   fChain->SetBranchAddress("selPhoGenSigMomVy", &selPhoGenSigMomVy, &b_selPhoGenSigMomVy);
   fChain->SetBranchAddress("selPhoGenSigMomVz", &selPhoGenSigMomVz, &b_selPhoGenSigMomVz);
   fChain->SetBranchAddress("selPhoHadOverEM", &selPhoHadOverEM, &b_selPhoHadOverEM);
   fChain->SetBranchAddress("selPhoHadTowOverEM", &selPhoHadTowOverEM, &b_selPhoHadTowOverEM);
   fChain->SetBranchAddress("selPhoHasConversionTracks", &selPhoHasConversionTracks, &b_selPhoHasConversionTracks);
   fChain->SetBranchAddress("selPhoHcalPFClusterIso", &selPhoHcalPFClusterIso, &b_selPhoHcalPFClusterIso);
   fChain->SetBranchAddress("selPhoHcalTowerSumEtBcConeDR04", &selPhoHcalTowerSumEtBcConeDR04, &b_selPhoHcalTowerSumEtBcConeDR04);
   fChain->SetBranchAddress("selPhoHcalTowerSumEtConeDR04", &selPhoHcalTowerSumEtConeDR04, &b_selPhoHcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("selPhoHoe_PUcorr", &selPhoHoe_PUcorr, &b_selPhoHoe_PUcorr);
   fChain->SetBranchAddress("selPhoMBetaEql", &selPhoMBetaEql, &b_selPhoMBetaEql);
   fChain->SetBranchAddress("selPhoMBetaPmt", &selPhoMBetaPmt, &b_selPhoMBetaPmt);
   fChain->SetBranchAddress("selPhoNTrkHollowConeDR04", &selPhoNTrkHollowConeDR04, &b_selPhoNTrkHollowConeDR04);
   fChain->SetBranchAddress("selPhoNTrkSolidConeDR04", &selPhoNTrkSolidConeDR04, &b_selPhoNTrkSolidConeDR04);
   fChain->SetBranchAddress("selPhoNrh", &selPhoNrh, &b_selPhoNrh);
   fChain->SetBranchAddress("selPhoOOT", &selPhoOOT, &b_selPhoOOT);
   fChain->SetBranchAddress("selPhoPVtoRH", &selPhoPVtoRH, &b_selPhoPVtoRH);
   fChain->SetBranchAddress("selPhoPfChargedIso", &selPhoPfChargedIso, &b_selPhoPfChargedIso);
   fChain->SetBranchAddress("selPhoPfChargedIsoPFPV", &selPhoPfChargedIsoPFPV, &b_selPhoPfChargedIsoPFPV);
   fChain->SetBranchAddress("selPhoPfChargedIsoWorstVtx", &selPhoPfChargedIsoWorstVtx, &b_selPhoPfChargedIsoWorstVtx);
   fChain->SetBranchAddress("selPhoPfPhoIso03", &selPhoPfPhoIso03, &b_selPhoPfPhoIso03);
   fChain->SetBranchAddress("selPhoPfRelIso03_all_quadratic", &selPhoPfRelIso03_all_quadratic, &b_selPhoPfRelIso03_all_quadratic);
   fChain->SetBranchAddress("selPhoPfRelIso03_chg_quadratic", &selPhoPfRelIso03_chg_quadratic, &b_selPhoPfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("selPhoPhi", &selPhoPhi, &b_selPhoPhi);
   fChain->SetBranchAddress("selPhoPhiWidth", &selPhoPhiWidth, &b_selPhoPhiWidth);
   fChain->SetBranchAddress("selPhoPhoIsoDr", &selPhoPhoIsoDr, &b_selPhoPhoIsoDr);
   fChain->SetBranchAddress("selPhoPixelSeed", &selPhoPixelSeed, &b_selPhoPixelSeed);
   fChain->SetBranchAddress("selPhoPt", &selPhoPt, &b_selPhoPt);
   fChain->SetBranchAddress("selPhoQuality", &selPhoQuality, &b_selPhoQuality);
   fChain->SetBranchAddress("selPhoR9", &selPhoR9, &b_selPhoR9);
   fChain->SetBranchAddress("selPhoS4", &selPhoS4, &b_selPhoS4);
   fChain->SetBranchAddress("selPhoSAlp", &selPhoSAlp, &b_selPhoSAlp);
   fChain->SetBranchAddress("selPhoSCx", &selPhoSCx, &b_selPhoSCx);
   fChain->SetBranchAddress("selPhoSCy", &selPhoSCy, &b_selPhoSCy);
   fChain->SetBranchAddress("selPhoSCz", &selPhoSCz, &b_selPhoSCz);
   fChain->SetBranchAddress("selPhoSMaj", &selPhoSMaj, &b_selPhoSMaj);
   fChain->SetBranchAddress("selPhoSMin", &selPhoSMin, &b_selPhoSMin);
   fChain->SetBranchAddress("selPhoSieie", &selPhoSieie, &b_selPhoSieie);
   fChain->SetBranchAddress("selPhoSieip", &selPhoSieip, &b_selPhoSieip);
   fChain->SetBranchAddress("selPhoSigmaIEtaIEta", &selPhoSigmaIEtaIEta, &b_selPhoSigmaIEtaIEta);
   fChain->SetBranchAddress("selPhoSipip", &selPhoSipip, &b_selPhoSipip);
   fChain->SetBranchAddress("selPhoSusyId", &selPhoSusyId, &b_selPhoSusyId);
   fChain->SetBranchAddress("selPhoTime", &selPhoTime, &b_selPhoTime);
   fChain->SetBranchAddress("selPhoTrkSumPtHollowConeDR03", &selPhoTrkSumPtHollowConeDR03, &b_selPhoTrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("selPhoTrkSumPtHollowConeDR04", &selPhoTrkSumPtHollowConeDR04, &b_selPhoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("selPhoTrkSumPtSolidConeDR04", &selPhoTrkSumPtSolidConeDR04, &b_selPhoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("selPhoSTimeSig", &selPhoSTimeSig, &b_selPhoSTimeSig);
   fChain->SetBranchAddress("selPhoLTimeSig", &selPhoLTimeSig, &b_selPhoLTimeSig);
   fChain->SetBranchAddress("selPhoWTimeSig", &selPhoWTimeSig, &b_selPhoWTimeSig);


   fChain->SetBranchAddress("muonIsLoose", &muonIsLoose, &b_muonIsLoose);
   fChain->SetBranchAddress("muonIsMedium", &muonIsMedium, &b_muonIsMedium);
   fChain->SetBranchAddress("nLooseMuons", &nLooseMuons, &b_nLooseMuons);
   fChain->SetBranchAddress("nMuonSVMatched", &nMuonSVMatched, &b_nMuonSVMatched);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);

   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("nSelJets", &nSelJets, &b_nSelJets);
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
   fChain->SetBranchAddress("selJetPt", &selJetPt, &b_selJetPt);
   fChain->SetBranchAddress("selJetQrkLlpId", &selJetQrkLlpId, &b_selJetQrkLlpId);
   fChain->SetBranchAddress("selJetQuality", &selJetQuality, &b_selJetQuality);
   fChain->SetBranchAddress("selJetTime", &selJetTime, &b_selJetTime);
   fChain->SetBranchAddress("selJetchHEF", &selJetchHEF, &b_selJetchHEF);

   fChain->SetBranchAddress("rjrABSide", &rjrABSide, &b_rjrABSide);
   fChain->SetBranchAddress("rjr_Rdiff", &rjr_Rdiff, &b_rjr_Rdiff);
   fChain->SetBranchAddress("rjrDphiMETV", &rjrDphiMETV, &b_rjrDphiMETV);
   fChain->SetBranchAddress("rjrEVa", &rjrEVa, &b_rjrEVa);
   fChain->SetBranchAddress("rjrEVb", &rjrEVb, &b_rjrEVb);
   fChain->SetBranchAddress("rjrMET", &rjrMET, &b_rjrMET);
   fChain->SetBranchAddress("rjrN2Px", &rjrN2Px, &b_rjrN2Px);
   fChain->SetBranchAddress("rjrN2Py", &rjrN2Py, &b_rjrN2Py);
   fChain->SetBranchAddress("rjrN2Pz", &rjrN2Pz, &b_rjrN2Pz);
   fChain->SetBranchAddress("rjrNJets", &rjrNJets, &b_rjrNJets);
   fChain->SetBranchAddress("rjrNJetsJa", &rjrNJetsJa, &b_rjrNJetsJa);
   fChain->SetBranchAddress("rjrNJetsJb", &rjrNJetsJb, &b_rjrNJetsJb);
   fChain->SetBranchAddress("rjrNPhotons", &rjrNPhotons, &b_rjrNPhotons);
   fChain->SetBranchAddress("rjrNVisObjects", &rjrNVisObjects, &b_rjrNVisObjects);
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
   //fChain->SetBranchAddress("rjr_p1Pxa", &rjr_p1Pxa, &b_rjr_p1Pxa);
   //fChain->SetBranchAddress("rjr_p1Pxb", &rjr_p1Pxb, &b_rjr_p1Pxb);
   fChain->SetBranchAddress("rjr_p2Pts", &rjr_p2Pts, &b_rjr_p2Pts);
   fChain->SetBranchAddress("rjr_p2Ptxa", &rjr_p2Ptxa, &b_rjr_p2Ptxa);
   fChain->SetBranchAddress("rjr_p2Ptxb", &rjr_p2Ptxb, &b_rjr_p2Ptxb);
   //fChain->SetBranchAddress("rjr_p2Pxa", &rjr_p2Pxa, &b_rjr_p2Pxa);
   //fChain->SetBranchAddress("rjr_p2Pxb", &rjr_p2Pxb, &b_rjr_p2Pxb);
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
   //Notify();

}//<<>>void llpgtree::Init( TTree *tree, bool doGenInfo )

void skimtuple_v36::getBranches( Long64_t entry, bool doGenInfo ){

   b_Flag_BadChargedCandidateFilter->GetEntry(entry);  //   Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   b_Flag_BadPFMuonDzFilter->GetEntry(entry);  //   Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   b_Flag_BadPFMuonFilter->GetEntry(entry);  //   Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   //b_Flag_EcalDeadCellTriggerPrimitiveFilter->GetEntry(entry);  //   Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   b_Flag_HBHENoiseFilter->GetEntry(entry);  //   Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   b_Flag_HBHENoiseIsoFilter->GetEntry(entry);  //   Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   b_Flag_MetFilter->GetEntry(entry);  //   Flag_MetFilter, &b_Flag_MetFilter);
   b_Flag_ecalBadCalibFilter->GetEntry(entry);  //   Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   b_Flag_eeBadScFilter->GetEntry(entry);  //   Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   b_Flag_globalSuperTightHalo2016Filter->GetEntry(entry);  //   Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   b_Flag_goodVertices->GetEntry(entry);  //   Flag_goodVertices, &b_Flag_goodVertices);
   b_Flag_hfNoisyHitsFilter->GetEntry(entry);  //   Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);

   b_PVx->GetEntry(entry);  //   PVx, &b_PVx);
   b_PVy->GetEntry(entry);  //   PVy, &b_PVy);
   b_PVz->GetEntry(entry);  //   PVz, &b_PVz);

/*
   b_SV_cosTheta->GetEntry(entry);  //   SV_cosTheta, &b_SV_cosTheta);
   b_SV_decayAngle->GetEntry(entry);  //   SV_decayAngle, &b_SV_decayAngle);
   b_SV_dxy->GetEntry(entry);  //   SV_dxy, &b_SV_dxy);
   b_SV_dxySig->GetEntry(entry);  //   SV_dxySig, &b_SV_dxySig);
   b_SV_mass->GetEntry(entry);  //   SV_mass, &b_SV_mass);
   b_SV_massOverNtracks->GetEntry(entry);  //   SV_massOverNtracks, &b_SV_massOverNtracks);
   b_SV_nElectron->GetEntry(entry);  //   SV_nElectron, &b_SV_nElectron);
   b_SV_nHadronic->GetEntry(entry);  //   SV_nHadronic, &b_SV_nHadronic);
   b_SV_nLeptonic->GetEntry(entry);  //   SV_nLeptonic, &b_SV_nLeptonic);
   b_SV_nMuon->GetEntry(entry);  //   SV_nMuon, &b_SV_nMuon);
   b_SV_nTracks->GetEntry(entry);  //   SV_nTracks, &b_SV_nTracks);
   b_SV_pOverE->GetEntry(entry);  //   SV_pOverE, &b_SV_pOverE);
*/

   b_DataSetKey->GetEntry(entry);  //   DataSetKey, &b_DataSetKey);
   b_evtFillWgt->GetEntry(entry);  //   evtFillWgt, &b_evtFillWgt);
   b_evtGenWgt->GetEntry(entry);  //   evtGenWgt, &b_evtGenWgt);
   b_evtXSection->GetEntry(entry);  //   evtXSection, &b_evtXSection);

   b_selCMet->GetEntry(entry);  //   selCMet, &b_selCMet);
   b_selCMetPx->GetEntry(entry);  //   selCMetPx, &b_selCMetPx);
   b_selCMetPy->GetEntry(entry);  //   selCMetPy, &b_selCMetPy);
   b_selMet->GetEntry(entry);  //   selMet, &b_selMet);
   b_selMetPx->GetEntry(entry);  //   selMetPx, &b_selMetPx);
   b_selMetPy->GetEntry(entry);  //   selMetPy, &b_selMetPy);

   b_genCharge->GetEntry(entry);  //   genCharge, &b_genCharge);
   b_genGrvtinoMass->GetEntry(entry);  //   genGrvtinoMass, &b_genGrvtinoMass);
   b_genLLPMass->GetEntry(entry);  //   genLLPMass, &b_genLLPMass);
   b_genLSPMass->GetEntry(entry);  //   genLSPMass, &b_genLSPMass);
   b_genMass->GetEntry(entry);  //   genMass, &b_genMass);
   b_genMomCTau->GetEntry(entry);  //   genMomCTau, &b_genMomCTau);
   b_genPartEnergy->GetEntry(entry);  //   genPartEnergy, &b_genPartEnergy);
   b_genPartEta->GetEntry(entry);  //   genPartEta, &b_genPartEta);
   b_genPartMomIdx->GetEntry(entry);  //   genPartMomIdx, &b_genPartMomIdx);
   b_genPartPdgId->GetEntry(entry);  //   genPartPdgId, &b_genPartPdgId);
   b_genPartPhi->GetEntry(entry);  //   genPartPhi, &b_genPartPhi);
   b_genPartPt->GetEntry(entry);  //   genPartPt, &b_genPartPt);
   b_genPartSusId->GetEntry(entry);  //   genPartSusId, &b_genPartSusId);
   b_genPx->GetEntry(entry);  //   genPx, &b_genPx);
   b_genPy->GetEntry(entry);  //   genPy, &b_genPy);
   b_genPz->GetEntry(entry);  //   genPz, &b_genPz);
   b_genSGMass->GetEntry(entry);  //   genSGMass, &b_genSGMass);
   b_genSQMass->GetEntry(entry);  //   genSQMass, &b_genSQMass);
   b_genSTFlagQP->GetEntry(entry);  //   genSTFlagQP, &b_genSTFlagQP);
   b_genSTFlagQQP->GetEntry(entry);  //   genSTFlagQQP, &b_genSTFlagQQP);
   b_genSTFlagQQZ->GetEntry(entry);  //   genSTFlagQQZ, &b_genSTFlagQQZ);
   b_genSTFlagQQZL->GetEntry(entry);  //   genSTFlagQQZL, &b_genSTFlagQQZL);
   b_genSTFlagQQZLL->GetEntry(entry);  //   genSTFlagQQZLL, &b_genSTFlagQQZLL);
   b_genSigType->GetEntry(entry);  //   genSigType, &b_genSigType);
   b_genStatus->GetEntry(entry);  //   genStatus, &b_genStatus);
   b_genVx->GetEntry(entry);  //   genVx, &b_genVx);
   b_genVy->GetEntry(entry);  //   genVy, &b_genVy);
   b_genVz->GetEntry(entry);  //   genVz, &b_genVz);
   b_genXMomCTau->GetEntry(entry);  //   genXMomCTau, &b_genXMomCTau);

   b_eleHasSVMatch->GetEntry(entry);  //   eleHasSVMatch, &b_eleHasSVMatch);
   b_eleIsLoose->GetEntry(entry);  //   eleIsLoose, &b_eleIsLoose);
   b_elePhoMinDr->GetEntry(entry);  //   elePhoMinDr, &b_elePhoMinDr);
   b_epDrMatch->GetEntry(entry);  //   epDrMatch, &b_epDrMatch);
   b_epSeedIdMatch->GetEntry(entry);  //   epSeedIdMatch, &b_epSeedIdMatch);
   b_epSeedIdUnique->GetEntry(entry);  //   epSeedIdUnique, &b_epSeedIdUnique);
   b_nEleSVMatched->GetEntry(entry);  //   nEleSVMatched, &b_nEleSVMatched);
   b_nElectrons->GetEntry(entry);  //   nElectrons, &b_nElectrons);
   b_nSelElectrons->GetEntry(entry);  //   nSelElectrons, &b_nSelElectrons);
   b_nSelIsoElectrons->GetEntry(entry);  //   nSelIsoElectrons, &b_nSelIsoElectrons);

   b_nPhotons->GetEntry(entry);  //   nPhotons, &b_nPhotons);
   b_nSelPhotons->GetEntry(entry);  //   nSelPhotons, &b_nSelPhotons);
   b_selPhoClstrRn->GetEntry(entry);  //   selPhoClstrRn, &b_selPhoClstrRn);
   b_selPhoCorEnergy->GetEntry(entry);  //   selPhoCorEnergy, &b_selPhoCorEnergy);
   b_selPhoCorPt->GetEntry(entry);  //   selPhoCorPt, &b_selPhoCorPt);
   b_selPhoCovEtaEta->GetEntry(entry);  //   selPhoCovEtaEta, &b_selPhoCovEtaEta);
   b_selPhoCovEtaPhi->GetEntry(entry);  //   selPhoCovEtaPhi, &b_selPhoCovEtaPhi);
   b_selPhoCovPhiPhi->GetEntry(entry);  //   selPhoCovPhiPhi, &b_selPhoCovPhiPhi);
   b_selPhoEcalPFClusterIso->GetEntry(entry);  //   selPhoEcalPFClusterIso, &b_selPhoEcalPFClusterIso);
   b_selPhoEcalRHSumEtConeDR04->GetEntry(entry);  //   selPhoEcalRHSumEtConeDR04, &b_selPhoEcalRHSumEtConeDR04);
   b_selPhoEleVeto->GetEntry(entry);  //   selPhoEleVeto, &b_selPhoEleVeto);
   b_selPhoEnergy->GetEntry(entry);  //   selPhoEnergy, &b_selPhoEnergy);
   b_selPhoEta->GetEntry(entry);  //   selPhoEta, &b_selPhoEta);
   b_selPhoEtaWidth->GetEntry(entry);  //   selPhoEtaWidth, &b_selPhoEtaWidth);
   b_selPhoGenIdx->GetEntry(entry);  //   selPhoGenIdx, &b_selPhoGenIdx);
   b_selPhoGenPt->GetEntry(entry);  //   selPhoGenPt, &b_selPhoGenPt);
   b_selPhoGenSigMomEnergy->GetEntry(entry);  //   selPhoGenSigMomEnergy, &b_selPhoGenSigMomEnergy);
   b_selPhoGenSigMomEta->GetEntry(entry);  //   selPhoGenSigMomEta, &b_selPhoGenSigMomEta);
   b_selPhoGenSigMomMass->GetEntry(entry);  //   selPhoGenSigMomMass, &b_selPhoGenSigMomMass);
   b_selPhoGenSigMomPhi->GetEntry(entry);  //   selPhoGenSigMomPhi, &b_selPhoGenSigMomPhi);
   b_selPhoGenSigMomPt->GetEntry(entry);  //   selPhoGenSigMomPt, &b_selPhoGenSigMomPt);
   b_selPhoGenSigMomPx->GetEntry(entry);  //   selPhoGenSigMomPx, &b_selPhoGenSigMomPx);
   b_selPhoGenSigMomPy->GetEntry(entry);  //   selPhoGenSigMomPy, &b_selPhoGenSigMomPy);
   b_selPhoGenSigMomPz->GetEntry(entry);  //   selPhoGenSigMomPz, &b_selPhoGenSigMomPz);
   b_selPhoGenSigMomVx->GetEntry(entry);  //   selPhoGenSigMomVx, &b_selPhoGenSigMomVx);
   b_selPhoGenSigMomVy->GetEntry(entry);  //   selPhoGenSigMomVy, &b_selPhoGenSigMomVy);
   b_selPhoGenSigMomVz->GetEntry(entry);  //   selPhoGenSigMomVz, &b_selPhoGenSigMomVz);
   b_selPhoHadOverEM->GetEntry(entry);  //   selPhoHadOverEM, &b_selPhoHadOverEM);
   b_selPhoHadTowOverEM->GetEntry(entry);  //   selPhoHadTowOverEM, &b_selPhoHadTowOverEM);
   b_selPhoHasConversionTracks->GetEntry(entry);  //   selPhoHasConversionTracks, &b_selPhoHasConversionTracks);
   b_selPhoHcalPFClusterIso->GetEntry(entry);  //   selPhoHcalPFClusterIso, &b_selPhoHcalPFClusterIso);
   b_selPhoHcalTowerSumEtBcConeDR04->GetEntry(entry);  //   selPhoHcalTowerSumEtBcConeDR04, &b_selPhoHcalTowerSumEtBcConeDR04);
   b_selPhoHcalTowerSumEtConeDR04->GetEntry(entry);  //   selPhoHcalTowerSumEtConeDR04, &b_selPhoHcalTowerSumEtConeDR04);
   b_selPhoHoe_PUcorr->GetEntry(entry);  //   selPhoHoe_PUcorr, &b_selPhoHoe_PUcorr);
   b_selPhoMBetaEql->GetEntry(entry);  //   selPhoMBetaEql, &b_selPhoMBetaEql);
   b_selPhoMBetaPmt->GetEntry(entry);  //   selPhoMBetaPmt, &b_selPhoMBetaPmt);
   b_selPhoNTrkHollowConeDR04->GetEntry(entry);  //   selPhoNTrkHollowConeDR04, &b_selPhoNTrkHollowConeDR04);
   b_selPhoNTrkSolidConeDR04->GetEntry(entry);  //   selPhoNTrkSolidConeDR04, &b_selPhoNTrkSolidConeDR04);
   b_selPhoNrh->GetEntry(entry);  //   selPhoNrh, &b_selPhoNrh);
   b_selPhoOOT->GetEntry(entry);  //   selPhoOOT, &b_selPhoOOT);
   b_selPhoPVtoRH->GetEntry(entry);  //   selPhoPVtoRH, &b_selPhoPVtoRH);
   b_selPhoPfChargedIso->GetEntry(entry);  //   selPhoPfChargedIso, &b_selPhoPfChargedIso);
   b_selPhoPfChargedIsoPFPV->GetEntry(entry);  //   selPhoPfChargedIsoPFPV, &b_selPhoPfChargedIsoPFPV);
   b_selPhoPfChargedIsoWorstVtx->GetEntry(entry);  //   selPhoPfChargedIsoWorstVtx, &b_selPhoPfChargedIsoWorstVtx);
   b_selPhoPfPhoIso03->GetEntry(entry);  //   selPhoPfPhoIso03, &b_selPhoPfPhoIso03);
   b_selPhoPfRelIso03_all_quadratic->GetEntry(entry);  //   selPhoPfRelIso03_all_quadratic, &b_selPhoPfRelIso03_all_quadratic);
   b_selPhoPfRelIso03_chg_quadratic->GetEntry(entry);  //   selPhoPfRelIso03_chg_quadratic, &b_selPhoPfRelIso03_chg_quadratic);
   b_selPhoPhi->GetEntry(entry);  //   selPhoPhi, &b_selPhoPhi);
   b_selPhoPhiWidth->GetEntry(entry);  //   selPhoPhiWidth, &b_selPhoPhiWidth);
   b_selPhoPhoIsoDr->GetEntry(entry);  //   selPhoPhoIsoDr, &b_selPhoPhoIsoDr);
   b_selPhoPixelSeed->GetEntry(entry);  //   selPhoPixelSeed, &b_selPhoPixelSeed);
   b_selPhoPt->GetEntry(entry);  //   selPhoPt, &b_selPhoPt);
   b_selPhoQuality->GetEntry(entry);  //   selPhoQuality, &b_selPhoQuality);
   b_selPhoR9->GetEntry(entry);  //   selPhoR9, &b_selPhoR9);
   b_selPhoS4->GetEntry(entry);  //   selPhoS4, &b_selPhoS4);
   b_selPhoSAlp->GetEntry(entry);  //   selPhoSAlp, &b_selPhoSAlp);
   b_selPhoSCx->GetEntry(entry);  //   selPhoSCx, &b_selPhoSCx);
   b_selPhoSCy->GetEntry(entry);  //   selPhoSCy, &b_selPhoSCy);
   b_selPhoSCz->GetEntry(entry);  //   selPhoSCz, &b_selPhoSCz);
   b_selPhoSMaj->GetEntry(entry);  //   selPhoSMaj, &b_selPhoSMaj);
   b_selPhoSMin->GetEntry(entry);  //   selPhoSMin, &b_selPhoSMin);
   b_selPhoSieie->GetEntry(entry);  //   selPhoSieie, &b_selPhoSieie);
   b_selPhoSieip->GetEntry(entry);  //   selPhoSieip, &b_selPhoSieip);
   b_selPhoSigmaIEtaIEta->GetEntry(entry);  //   selPhoSigmaIEtaIEta, &b_selPhoSigmaIEtaIEta);
   b_selPhoSipip->GetEntry(entry);  //   selPhoSipip, &b_selPhoSipip);
   b_selPhoSusyId->GetEntry(entry);  //   selPhoSusyId, &b_selPhoSusyId);
   b_selPhoTime->GetEntry(entry);  //   selPhoTime, &b_selPhoTime);
   b_selPhoTrkSumPtHollowConeDR03->GetEntry(entry);  //   selPhoTrkSumPtHollowConeDR03, &b_selPhoTrkSumPtHollowConeDR03);
   b_selPhoTrkSumPtHollowConeDR04->GetEntry(entry);  //   selPhoTrkSumPtHollowConeDR04, &b_selPhoTrkSumPtHollowConeDR04);
   b_selPhoTrkSumPtSolidConeDR04->GetEntry(entry);  //   selPhoTrkSumPtSolidConeDR04, &b_selPhoTrkSumPtSolidConeDR04);
   b_selPhoSTimeSig->GetEntry(entry);  //   selPhoSTimeSig, &b_selPhoSTimeSig);
   b_selPhoLTimeSig->GetEntry(entry);  //   selPhoLTimeSig, &b_selPhoLTimeSig);
   b_selPhoWTimeSig->GetEntry(entry);  //   selPhoWTimeSig, &b_selPhoWTimeSig);

   b_muonIsLoose->GetEntry(entry);  //   muonIsLoose, &b_muonIsLoose);
   b_muonIsMedium->GetEntry(entry);  //   muonIsMedium, &b_muonIsMedium);
   b_nLooseMuons->GetEntry(entry);  //   nLooseMuons, &b_nLooseMuons);
   b_nMuonSVMatched->GetEntry(entry);  //   nMuonSVMatched, &b_nMuonSVMatched);
   b_nMuons->GetEntry(entry);  //   nMuons, &b_nMuons);

   b_nJets->GetEntry(entry);  //   nJets, &b_nJets);
   b_nSelJets->GetEntry(entry);  //   nSelJets, &b_nSelJets);
   b_selGenJetDpt->GetEntry(entry);  //   selGenJetDpt, &b_selGenJetDpt);
   b_selGenJetEnergy->GetEntry(entry);  //   selGenJetEnergy, &b_selGenJetEnergy);
   b_selGenJetImpAng->GetEntry(entry);  //   selGenJetImpAng, &b_selGenJetImpAng);
   b_selGenJetLlpTime->GetEntry(entry);  //   selGenJetLlpTime, &b_selGenJetLlpTime);
   b_selGenJetPt->GetEntry(entry);  //   selGenJetPt, &b_selGenJetPt);
   b_selGenJetTime->GetEntry(entry);  //   selGenJetTime, &b_selGenJetTime);
   b_selGenJetTof->GetEntry(entry);  //   selGenJetTof, &b_selGenJetTof);
   b_selGenJetdr->GetEntry(entry);  //   selGenJetdr, &b_selGenJetdr);
   b_selGenJeteta->GetEntry(entry);  //   selGenJeteta, &b_selGenJeteta);
   b_selJetArea->GetEntry(entry);  //   selJetArea, &b_selJetArea);
   b_selJetChEmEF->GetEntry(entry);  //   selJetChEmEF, &b_selJetChEmEF);
   b_selJetChHM->GetEntry(entry);  //   selJetChHM, &b_selJetChHM);
   b_selJetEnergy->GetEntry(entry);  //   selJetEnergy, &b_selJetEnergy);
   b_selJetEta->GetEntry(entry);  //   selJetEta, &b_selJetEta);
   b_selJetGenLlpId->GetEntry(entry);  //   selJetGenLlpId, &b_selJetGenLlpId);
   b_selJetLlpDp->GetEntry(entry);  //   selJetLlpDp, &b_selJetLlpDp);
   b_selJetLlpDr->GetEntry(entry);  //   selJetLlpDr, &b_selJetLlpDr);
   b_selJetMass->GetEntry(entry);  //   selJetMass, &b_selJetMass);
   b_selJetMuEF->GetEntry(entry);  //   selJetMuEF, &b_selJetMuEF);
   b_selJetNeEmEF->GetEntry(entry);  //   selJetNeEmEF, &b_selJetNeEmEF);
   b_selJetNeHEF->GetEntry(entry);  //   selJetNeHEF, &b_selJetNeHEF);
   b_selJetNeHM->GetEntry(entry);  //   selJetNeHM, &b_selJetNeHM);
   b_selJetPhi->GetEntry(entry);  //   selJetPhi, &b_selJetPhi);
   b_selJetPt->GetEntry(entry);  //   selJetPt, &b_selJetPt);
   b_selJetQrkLlpId->GetEntry(entry);  //   selJetQrkLlpId, &b_selJetQrkLlpId);
   b_selJetQuality->GetEntry(entry);  //   selJetQuality, &b_selJetQuality);
   b_selJetTime->GetEntry(entry);  //   selJetTime, &b_selJetTime);
   b_selJetchHEF->GetEntry(entry);  //   selJetchHEF, &b_selJetchHEF);

   b_rjrABSide->GetEntry(entry);  //   rjrABSide, &b_rjrABSide);
   b_rjr_Rdiff->GetEntry(entry);  //   rjr_Rdiff, &b_rjr_Rdiff);
   b_rjrDphiMETV->GetEntry(entry);  //   rjrDphiMETV, &b_rjrDphiMETV);
   b_rjrEVa->GetEntry(entry);  //   rjrEVa, &b_rjrEVa);
   b_rjrEVb->GetEntry(entry);  //   rjrEVb, &b_rjrEVb);
   b_rjrMET->GetEntry(entry);  //   rjrMET, &b_rjrMET);
   b_rjrN2Px->GetEntry(entry);  //   rjrN2Px, &b_rjrN2Px);
   b_rjrN2Py->GetEntry(entry);  //   rjrN2Py, &b_rjrN2Py);
   b_rjrN2Pz->GetEntry(entry);  //   rjrN2Pz, &b_rjrN2Pz);
   b_rjrNJets->GetEntry(entry);  //   rjrNJets, &b_rjrNJets);
   b_rjrNJetsJa->GetEntry(entry);  //   rjrNJetsJa, &b_rjrNJetsJa);
   b_rjrNJetsJb->GetEntry(entry);  //   rjrNJetsJb, &b_rjrNJetsJb);
   b_rjrNPhotons->GetEntry(entry);  //   rjrNPhotons, &b_rjrNPhotons);
   b_rjrNVisObjects->GetEntry(entry);  //   rjrNVisObjects, &b_rjrNVisObjects);
   b_rjrPTS->GetEntry(entry);  //   rjrPTS, &b_rjrPTS);
   b_rjrPVa->GetEntry(entry);  //   rjrPVa, &b_rjrPVa);
   b_rjrPVb->GetEntry(entry);  //   rjrPVb, &b_rjrPVb);
   b_rjrPVlab->GetEntry(entry);  //   rjrPVlab, &b_rjrPVlab);
   b_rjrPZS->GetEntry(entry);  //   rjrPZS, &b_rjrPZS);
   b_rjrSCosA->GetEntry(entry);  //   rjrSCosA, &b_rjrSCosA);
   b_rjrSMass->GetEntry(entry);  //   rjrSMass, &b_rjrSMass);
   b_rjrSdphiBV->GetEntry(entry);  //   rjrSdphiBV, &b_rjrSdphiBV);
   b_rjrSdphiDA->GetEntry(entry);  //   rjrSdphiDA, &b_rjrSdphiDA);
   b_rjrSubPhoLocation->GetEntry(entry);  //   rjrSubPhoLocation, &b_rjrSubPhoLocation);
   b_rjrType->GetEntry(entry);  //   rjrType, &b_rjrType);
   b_rjrVisJetSide->GetEntry(entry);  //   rjrVisJetSide, &b_rjrVisJetSide);
   b_rjrVisPhoSide->GetEntry(entry);  //   rjrVisPhoSide, &b_rjrVisPhoSide);
   b_rjrX1aCosA->GetEntry(entry);  //   rjrX1aCosA, &b_rjrX1aCosA);
   b_rjrX1aMass->GetEntry(entry);  //   rjrX1aMass, &b_rjrX1aMass);
   b_rjrX1a_Eta->GetEntry(entry);  //   rjrX1a_Eta, &b_rjrX1a_Eta);
   b_rjrX1a_Phi->GetEntry(entry);  //   rjrX1a_Phi, &b_rjrX1a_Phi);
   b_rjrX1a_Pt->GetEntry(entry);  //   rjrX1a_Pt, &b_rjrX1a_Pt);
   b_rjrX1bCosA->GetEntry(entry);  //   rjrX1bCosA, &b_rjrX1bCosA);
   b_rjrX1bMass->GetEntry(entry);  //   rjrX1bMass, &b_rjrX1bMass);
   b_rjrX1b_Eta->GetEntry(entry);  //   rjrX1b_Eta, &b_rjrX1b_Eta);
   b_rjrX1b_Phi->GetEntry(entry);  //   rjrX1b_Phi, &b_rjrX1b_Phi);
   b_rjrX1b_Pt->GetEntry(entry);  //   rjrX1b_Pt, &b_rjrX1b_Pt);
   b_rjrX2aCosA->GetEntry(entry);  //   rjrX2aCosA, &b_rjrX2aCosA);
   b_rjrX2aMass->GetEntry(entry);  //   rjrX2aMass, &b_rjrX2aMass);
   b_rjrX2aPs->GetEntry(entry);  //   rjrX2aPs, &b_rjrX2aPs);
   b_rjrX2bCosA->GetEntry(entry);  //   rjrX2bCosA, &b_rjrX2bCosA);
   b_rjrX2bMass->GetEntry(entry);  //   rjrX2bMass, &b_rjrX2bMass);
   b_rjrX2bPs->GetEntry(entry);  //   rjrX2bPs, &b_rjrX2bPs);
   b_rjr_Mr->GetEntry(entry);  //   rjr_Mr, &b_rjr_Mr);
   b_rjr_Ms->GetEntry(entry);  //   rjr_Ms, &b_rjr_Ms);
   b_rjr_Mva->GetEntry(entry);  //   rjr_Mva, &b_rjr_Mva);
   b_rjr_Mvb->GetEntry(entry);  //   rjr_Mvb, &b_rjr_Mvb);
   b_rjr_R->GetEntry(entry);  //   rjr_R, &b_rjr_R);
   b_rjr_Rm->GetEntry(entry);  //   rjr_Rm, &b_rjr_Rm);
   b_rjr_Rs->GetEntry(entry);  //   rjr_Rs, &b_rjr_Rs);
   b_rjr_Rv->GetEntry(entry);  //   rjr_Rv, &b_rjr_Rv);
   b_rjr_Rx->GetEntry(entry);  //   rjr_Rx, &b_rjr_Rx);
   b_rjr_Rxa->GetEntry(entry);  //   rjr_Rxa, &b_rjr_Rxa);
   b_rjr_Rxb->GetEntry(entry);  //   rjr_Rxb, &b_rjr_Rxb);
   b_rjr_p1Pts->GetEntry(entry);  //   rjr_p1Pts, &b_rjr_p1Pts);
   b_rjr_p1Ptxa->GetEntry(entry);  //   rjr_p1Ptxa, &b_rjr_p1Ptxa);
   b_rjr_p1Ptxb->GetEntry(entry);  //   rjr_p1Ptxb, &b_rjr_p1Ptxb);
   //b_rjr_p1Pxa->GetEntry(entry);  //   rjr_p1Ptxa, &b_rjr_p1Ptxa);
   //b_rjr_p1Pxb->GetEntry(entry);  //   rjr_p1Ptxb, &b_rjr_p1Ptxb);
   b_rjr_p2Pts->GetEntry(entry);  //   rjr_p2Pts, &b_rjr_p2Pts);
   b_rjr_p2Ptxa->GetEntry(entry);  //   rjr_p2Ptxa, &b_rjr_p2Ptxa);
   b_rjr_p2Ptxb->GetEntry(entry);  //   rjr_p2Ptxb, &b_rjr_p2Ptxb);
   //b_rjr_p2Pxa->GetEntry(entry);  //   rjr_p2Ptxa, &b_rjr_p2Ptxa);
   //b_rjr_p2Pxb->GetEntry(entry);  //   rjr_p2Ptxb, &b_rjr_p2Ptxb);
   b_rjr_pHs11->GetEntry(entry);  //   rjr_pHs11, &b_rjr_pHs11);
   b_rjr_pHs11a->GetEntry(entry);  //   rjr_pHs11a, &b_rjr_pHs11a);
   b_rjr_pHs11b->GetEntry(entry);  //   rjr_pHs11b, &b_rjr_pHs11b);
   b_rjr_pHs20->GetEntry(entry);  //   rjr_pHs20, &b_rjr_pHs20);
   b_rjr_pHs21->GetEntry(entry);  //   rjr_pHs21, &b_rjr_pHs21);
   b_rjr_pHs21a->GetEntry(entry);  //   rjr_pHs21a, &b_rjr_pHs21a);
   b_rjr_pHs21b->GetEntry(entry);  //   rjr_pHs21b, &b_rjr_pHs21b);
   b_rjr_pHs22->GetEntry(entry);  //   rjr_pHs22, &b_rjr_pHs22);
   b_rjr_pHs40->GetEntry(entry);  //   rjr_pHs40, &b_rjr_pHs40);
   b_rjr_pHs41->GetEntry(entry);  //   rjr_pHs41, &b_rjr_pHs41);
   b_rjr_pHts11->GetEntry(entry);  //   rjr_pHts11, &b_rjr_pHts11);
   b_rjr_pHts11a->GetEntry(entry);  //   rjr_pHts11a, &b_rjr_pHts11a);
   b_rjr_pHts11b->GetEntry(entry);  //   rjr_pHts11b, &b_rjr_pHts11b);
   b_rjr_pHts20->GetEntry(entry);  //   rjr_pHts20, &b_rjr_pHts20);
   b_rjr_pHts21->GetEntry(entry);  //   rjr_pHts21, &b_rjr_pHts21);
   b_rjr_pHts21a->GetEntry(entry);  //   rjr_pHts21a, &b_rjr_pHts21a);
   b_rjr_pHts21b->GetEntry(entry);  //   rjr_pHts21b, &b_rjr_pHts21b);
   b_rjr_pHts22->GetEntry(entry);  //   rjr_pHts22, &b_rjr_pHts22);
   b_rjr_pHts40->GetEntry(entry);  //   rjr_pHts40, &b_rjr_pHts40);
   b_rjr_pHts41->GetEntry(entry);  //   rjr_pHts41, &b_rjr_pHts41);
   b_rjr_pHtxa10->GetEntry(entry);  //   rjr_pHtxa10, &b_rjr_pHtxa10);
   b_rjr_pHtxa11->GetEntry(entry);  //   rjr_pHtxa11, &b_rjr_pHtxa11);
   b_rjr_pHtxa20->GetEntry(entry);  //   rjr_pHtxa20, &b_rjr_pHtxa20);
   b_rjr_pHtxa21->GetEntry(entry);  //   rjr_pHtxa21, &b_rjr_pHtxa21);
   b_rjr_pHtxb10->GetEntry(entry);  //   rjr_pHtxb10, &b_rjr_pHtxb10);
   b_rjr_pHtxb11->GetEntry(entry);  //   rjr_pHtxb11, &b_rjr_pHtxb11);
   b_rjr_pHtxb20->GetEntry(entry);  //   rjr_pHtxb20, &b_rjr_pHtxb20);
   b_rjr_pHtxb21->GetEntry(entry);  //   rjr_pHtxb21, &b_rjr_pHtxb21);
   b_rjr_pHxa10->GetEntry(entry);  //   rjr_pHxa10, &b_rjr_pHxa10);
   b_rjr_pHxa11->GetEntry(entry);  //   rjr_pHxa11, &b_rjr_pHxa11);
   b_rjr_pHxa20->GetEntry(entry);  //   rjr_pHxa20, &b_rjr_pHxa20);
   b_rjr_pHxa21->GetEntry(entry);  //   rjr_pHxa21, &b_rjr_pHxa21);
   b_rjr_pHxb10->GetEntry(entry);  //   rjr_pHxb10, &b_rjr_pHxb10);
   b_rjr_pHxb11->GetEntry(entry);  //   rjr_pHxb11, &b_rjr_pHxb11);
   b_rjr_pHxb20->GetEntry(entry);  //   rjr_pHxb20, &b_rjr_pHxb20);
   b_rjr_pHxb21->GetEntry(entry);  //   rjr_pHxb21, &b_rjr_pHxb21);

}//<<>>void root_base::getBranches(Long64_t entry)

/*
Bool_t skimtuple_v36::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimtuple_v36::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimtuple_v36::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/

#endif // #ifdef skimtuple_v36_cxx
