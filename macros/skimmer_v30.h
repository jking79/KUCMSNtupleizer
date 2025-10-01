//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 21 13:10:14 2025 by ROOT version 6.26/11
// from TTree kuSkimTree/output root file for kUCMSSkimmer
// found on file: /uscms/home/jaking/nobackup/llpana_skims/SMS_SVIPM100_v30_SqSq_AODSIM_mGl-1850_mN2-1650_mN1-100_ct0p1_rjrskim_v30.root
//////////////////////////////////////////////////////////

#ifndef skimmer_v30_h
#define skimmer_v30_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "string"
#include "vector"
#include "vector"

class skimmer_v30 {
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
   vector<float>   *SV_cosTheta;
   vector<float>   *SV_decayAngle;
   vector<float>   *SV_dxy;
   vector<float>   *SV_dxySig;
   vector<float>   *SV_mass;
   vector<float>   *SV_massOverNtracks;
   Int_t           SV_nElectron;
   Int_t           SV_nHadronic;
   Int_t           SV_nLeptonic;
   Int_t           SV_nMuon;
   vector<int>     *SV_nTracks;
   vector<float>   *SV_pOverE;
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
   UInt_t          nPhotons;
   UInt_t          nSelPhotons;
   vector<float>   *selPhoClstrRn;
   vector<float>   *selPhoCorEnergy;
   vector<float>   *selPhoCorPt;
   vector<float>   *selPhoCovEtaEta;
   vector<float>   *selPhoCovEtaPhi;
   vector<float>   *selPhoCovPhiPhi;
   vector<float>   *selPhoEcalPFClusterIso;
   vector<float>   *selPhoEcalRHSumEtConeDR04;
   UInt_t          selPhoEleVeto;
   vector<float>   *selPhoEnergy;
   vector<float>   *selPhoEta;
   vector<float>   *selPhoEtaWidth;
   vector<int>     *selPhoGenIdx;
   vector<float>   *selPhoGenPt;
   vector<float>   *selPhoGenSigMomEnergy;
   vector<float>   *selPhoGenSigMomEta;
   vector<float>   *selPhoGenSigMomMass;
   vector<float>   *selPhoGenSigMomPhi;
   vector<float>   *selPhoGenSigMomPt;
   vector<float>   *selPhoGenSigMomPx;
   vector<float>   *selPhoGenSigMomPy;
   vector<float>   *selPhoGenSigMomPz;
   vector<float>   *selPhoGenSigMomVx;
   vector<float>   *selPhoGenSigMomVy;
   vector<float>   *selPhoGenSigMomVz;
   vector<float>   *selPhoHadOverEM;
   vector<float>   *selPhoHadTowOverEM;
   vector<bool>    *selPhoHasConversionTracks;
   vector<float>   *selPhoHcalPFClusterIso;
   vector<float>   *selPhoHcalTowerSumEtBcConeDR04;
   vector<float>   *selPhoHcalTowerSumEtConeDR04;
   vector<float>   *selPhoHoe_PUcorr;
   vector<float>   *selPhoMBetaEql;
   vector<float>   *selPhoMBetaPmt;
   vector<float>   *selPhoNTrkHollowConeDR04;
   vector<float>   *selPhoNTrkSolidConeDR04;
   vector<unsigned int> *selPhoNrh;
   vector<bool>    *selPhoOOT;
   vector<float>   *selPhoPVtoRH;
   vector<float>   *selPhoPfChargedIso;
   vector<float>   *selPhoPfChargedIsoPFPV;
   vector<float>   *selPhoPfChargedIsoWorstVtx;
   vector<float>   *selPhoPfPhoIso03;
   vector<float>   *selPhoPfRelIso03_all_quadratic;
   vector<float>   *selPhoPfRelIso03_chg_quadratic;
   vector<float>   *selPhoPhi;
   vector<float>   *selPhoPhiWidth;
   vector<float>   *selPhoPhoIsoDr;
   vector<bool>    *selPhoPixelSeed;
   vector<float>   *selPhoPt;
   vector<int>     *selPhoQuality;
   vector<float>   *selPhoR9;
   vector<float>   *selPhoS4;
   vector<float>   *selPhoSAlp;
   vector<float>   *selPhoSCx;
   vector<float>   *selPhoSCy;
   vector<float>   *selPhoSCz;
   vector<float>   *selPhoSMaj;
   vector<float>   *selPhoSMin;
   vector<float>   *selPhoSieie;
   vector<float>   *selPhoSieip;
   vector<float>   *selPhoSigmaIEtaIEta;
   vector<float>   *selPhoSipip;
   vector<float>   *selPhoSusyId;
   vector<float>   *selPhoTime;
   vector<float>   *selPhoTrkSumPtHollowConeDR03;
   vector<float>   *selPhoTrkSumPtHollowConeDR04;
   vector<float>   *selPhoTrkSumPtSolidConeDR04;
   vector<bool>    *muonIsLoose;
   vector<bool>    *muonIsMedium;
   UInt_t          nLooseMuons;
   Int_t           nMuonSVMatched;
   UInt_t          nMuons;
   UInt_t          nJets;
   UInt_t          nSelJets;
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
   vector<float>   *selJetPt;
   vector<int>     *selJetQrkLlpId;
   vector<int>     *selJetQuality;
   vector<float>   *selJetTime;
   vector<float>   *selJetchHEF;
   vector<bool>    *rjrABSide;
   vector<float>   *rjr_Rdiff;
   vector<float>   *rjrDphiMETV;
   vector<float>   *rjrEVa;
   vector<float>   *rjrEVb;
   vector<float>   *rjrMET;
   vector<float>   *rjrN2Px;
   vector<float>   *rjrN2Py;
   vector<float>   *rjrN2Pz;
   vector<int>     *rjrNJets;
   vector<int>     *rjrNJetsJa;
   vector<int>     *rjrNJetsJb;
   vector<int>     *rjrNPhotons;
   vector<int>     *rjrNVisObjects;
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
   vector<float>   *rjrX2aPtS;
   vector<float>   *rjrX2bCosA;
   vector<float>   *rjrX2bMass;
   vector<float>   *rjrX2bPtS;
   vector<float>   *rjr_Mr;
   vector<float>   *rjr_R;
   vector<float>   *rjr_Rv;

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
   TBranch        *b_rjrX2aPtS;   //!
   TBranch        *b_rjrX2bCosA;   //!
   TBranch        *b_rjrX2bMass;   //!
   TBranch        *b_rjrX2bPtS;   //!
   TBranch        *b_rjr_Mr;   //!
   TBranch        *b_rjr_R;   //!
   TBranch        *b_rjr_Rv;   //!

   skimmer_v30(TTree *tree=0);
   virtual ~skimmer_v30();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef skimmer_v30_cxx
skimmer_v30::skimmer_v30(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms/home/jaking/nobackup/llpana_skims/SMS_SVIPM100_v30_SqSq_AODSIM_mGl-1850_mN2-1650_mN1-100_ct0p1_rjrskim_v30.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uscms/home/jaking/nobackup/llpana_skims/SMS_SVIPM100_v30_SqSq_AODSIM_mGl-1850_mN2-1650_mN1-100_ct0p1_rjrskim_v30.root");
      }
      f->GetObject("kuSkimTree",tree);

   }
   Init(tree);
}

skimmer_v30::~skimmer_v30()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimmer_v30::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimmer_v30::LoadTree(Long64_t entry)
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

void skimmer_v30::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

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
   rjrX2aPtS = 0;
   rjrX2bCosA = 0;
   rjrX2bMass = 0;
   rjrX2bPtS = 0;
   rjr_Mr = 0;
   rjr_R = 0;
   rjr_Rv = 0;
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
   fChain->SetBranchAddress("Flag_MetFilter", &Flag_MetFilter, &b_Flag_MetFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
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
   fChain->SetBranchAddress("rjrX2aPtS", &rjrX2aPtS, &b_rjrX2aPtS);
   fChain->SetBranchAddress("rjrX2bCosA", &rjrX2bCosA, &b_rjrX2bCosA);
   fChain->SetBranchAddress("rjrX2bMass", &rjrX2bMass, &b_rjrX2bMass);
   fChain->SetBranchAddress("rjrX2bPtS", &rjrX2bPtS, &b_rjrX2bPtS);
   fChain->SetBranchAddress("rjr_Mr", &rjr_Mr, &b_rjr_Mr);
   fChain->SetBranchAddress("rjr_R", &rjr_R, &b_rjr_R);
   fChain->SetBranchAddress("rjr_Rv", &rjr_Rv, &b_rjr_Rv);
   Notify();
}

Bool_t skimmer_v30::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimmer_v30::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimmer_v30::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skimmer_v30_cxx
