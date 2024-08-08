//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 17 05:42:16 2024 by ROOT version 6.26/11
// from TTree kuSkimTree/output root file for kUCMSSkimmer
// found on file: rjr_skim_files/kuntuple_GJets_HT-100To200_Met0_R17_v19_Ntuple_v19_LLPgama_Skim_v19_Met0_rjrcombo.root
//////////////////////////////////////////////////////////

#ifndef skim_mod2_h
#define skim_mod2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class skim_mod2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         PVx;
   Float_t         PVy;
   Float_t         PVz;
   string          *DataSetKey;
   Float_t         evtGenWgt;
   Float_t         evtXSection;
   Float_t         selCMet;
   Float_t         selCMetPx;
   Float_t         selCMetPy;
   Float_t         selMet;
   Float_t         selMetPx;
   Float_t         selMetPy;
   vector<int>     *genCharge;
   vector<float>   *genMass;
   vector<float>   *genPartEnergy;
   vector<float>   *genPartEta;
   vector<unsigned int> *genPartPdgId;
   vector<float>   *genPartPhi;
   vector<float>   *genPartPt;
   vector<int>     *genPartSusId;
   vector<float>   *genPx;
   vector<float>   *genPy;
   vector<float>   *genPz;
   vector<float>   *genSGMass;
   vector<float>   *genSQMass;
   Int_t           genSigType;
   Bool_t          genStatus;
   vector<float>   *genVx;
   vector<float>   *genVy;
   vector<float>   *genVz;
   vector<float>   *elePhoMinDr;
   vector<float>   *epDrMatch;
   vector<float>   *epSeedIdMatch;
   vector<float>   *epSeedIdUnique;
   UInt_t          nElectrons;
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
   vector<float>   *selJetLlpDp;
   vector<float>   *selJetLlpDr;
   vector<float>   *selJetMass;
   vector<float>   *selJetMuEF;
   vector<float>   *selJetNeEmEF;
   vector<float>   *selJetNeHEF;
   vector<float>   *selJetNeHM;
   vector<float>   *selJetPhi;
   vector<float>   *selJetPt;
   vector<int>     *selJetQuality;
   vector<float>   *selJetSusyId;
   vector<float>   *selJetTime;
   vector<float>   *selJetchHEF;
   vector<float>   *rjrASMass;
   vector<float>   *rjrAX2Ave;
   vector<float>   *rjrAX2Diff;
   vector<float>   *rjrAX2GMean;
   vector<float>   *rjrAX2NAve;
   vector<float>   *rjrAX2NGMean;
   vector<float>   *rjrAX2NVSum;
   vector<float>   *rjrAX2VSum;
   vector<float>   *rjrAX2aMass;
   vector<float>   *rjrAX2bMass;
   vector<float>   *rjrDiffPJX2a;
   vector<float>   *rjrDiffPJX2b;
   vector<float>   *rjrDphiMETV;
   vector<float>   *rjrEVa;
   vector<float>   *rjrEVb;
   vector<float>   *rjrMET;
   vector<float>   *rjrMV;
   vector<float>   *rjrMVa;
   vector<float>   *rjrMVb;
   vector<int>     *rjrNJetsJa;
   vector<int>     *rjrNJetsJb;
   vector<int>     *rjrNRjrJets;
   vector<int>     *rjrNRjrPhotons;
   vector<float>   *rjrPJX2a;
   vector<float>   *rjrPJX2b;
   vector<float>   *rjrPTS;
   vector<float>   *rjrPV;
   vector<float>   *rjrPVa;
   vector<float>   *rjrPVb;
   vector<float>   *rjrPVlab;
   vector<float>   *rjrPX1X2a;
   vector<float>   *rjrPX1X2b;
   vector<float>   *rjrPX2Sa;
   vector<float>   *rjrPX2Sb;
   vector<float>   *rjrPZS;
   vector<float>   *rjrSCosA;
   vector<float>   *rjrSMass;
   vector<float>   *rjrSdphiBV;
   vector<float>   *rjrSdphiDA;
   vector<float>   *rjrX1aCosA;
   vector<float>   *rjrX1aMass;
   vector<float>   *rjrX1bCosA;
   vector<float>   *rjrX1bMass;
   vector<float>   *rjrX2Ave;
   vector<float>   *rjrX2Diff;
   vector<float>   *rjrX2GMean;
   vector<float>   *rjrX2NAve;
   vector<float>   *rjrX2NGMean;
   vector<float>   *rjrX2NVSum;
   vector<float>   *rjrX2VSum;
   vector<float>   *rjrX2aCosA;
   vector<float>   *rjrX2aMass;
   vector<float>   *rjrX2aPtS;
   vector<float>   *rjrX2bCosA;
   vector<float>   *rjrX2bMass;
   vector<float>   *rjrX2bPtS;

   // List of branches
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_DataSetKey;   //!
   TBranch        *b_evtGenWgt;   //!
   TBranch        *b_evtXSection;   //!
   TBranch        *b_selCMet;   //!
   TBranch        *b_selCMetPx;   //!
   TBranch        *b_selCMetPy;   //!
   TBranch        *b_selMet;   //!
   TBranch        *b_selMetPx;   //!
   TBranch        *b_selMetPy;   //!
   TBranch        *b_genCharge;   //!
   TBranch        *b_genMass;   //!
   TBranch        *b_genPartEnergy;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPdgId;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartSusId;   //!
   TBranch        *b_genPx;   //!
   TBranch        *b_genPy;   //!
   TBranch        *b_genPz;   //!
   TBranch        *b_genSGMass;   //!
   TBranch        *b_genSQMass;   //!
   TBranch        *b_genSigType;   //!
   TBranch        *b_genStatus;   //!
   TBranch        *b_genVx;   //!
   TBranch        *b_genVy;   //!
   TBranch        *b_genVz;   //!
   TBranch        *b_elePhoMinDr;   //!
   TBranch        *b_epDrMatch;   //!
   TBranch        *b_epSeedIdMatch;   //!
   TBranch        *b_epSeedIdUnique;   //!
   TBranch        *b_nElectrons;   //!
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
   TBranch        *b_selJetLlpDp;   //!
   TBranch        *b_selJetLlpDr;   //!
   TBranch        *b_selJetMass;   //!
   TBranch        *b_selJetMuEF;   //!
   TBranch        *b_selJetNeEmEF;   //!
   TBranch        *b_selJetNeHEF;   //!
   TBranch        *b_selJetNeHM;   //!
   TBranch        *b_selJetPhi;   //!
   TBranch        *b_selJetPt;   //!
   TBranch        *b_selJetQuality;   //!
   TBranch        *b_selJetSusyId;   //!
   TBranch        *b_selJetTime;   //!
   TBranch        *b_selJetchHEF;   //!
   TBranch        *b_rjrASMass;   //!
   TBranch        *b_rjrAX2Ave;   //!
   TBranch        *b_rjrAX2Diff;   //!
   TBranch        *b_rjrAX2GMean;   //!
   TBranch        *b_rjrAX2NAve;   //!
   TBranch        *b_rjrAX2NGMean;   //!
   TBranch        *b_rjrAX2NVSum;   //!
   TBranch        *b_rjrAX2VSum;   //!
   TBranch        *b_rjrAX2aMass;   //!
   TBranch        *b_rjrAX2bMass;   //!
   TBranch        *b_rjrDiffPJX2a;   //!
   TBranch        *b_rjrDiffPJX2b;   //!
   TBranch        *b_rjrDphiMETV;   //!
   TBranch        *b_rjrEVa;   //!
   TBranch        *b_rjrEVb;   //!
   TBranch        *b_rjrMET;   //!
   TBranch        *b_rjrMV;   //!
   TBranch        *b_rjrMVa;   //!
   TBranch        *b_rjrMVb;   //!
   TBranch        *b_rjrNJetsJa;   //!
   TBranch        *b_rjrNJetsJb;   //!
   TBranch        *b_rjrNRjrJets;   //!
   TBranch        *b_rjrNRjrPhotons;   //!
   TBranch        *b_rjrPJX2a;   //!
   TBranch        *b_rjrPJX2b;   //!
   TBranch        *b_rjrPTS;   //!
   TBranch        *b_rjrPV;   //!
   TBranch        *b_rjrPVa;   //!
   TBranch        *b_rjrPVb;   //!
   TBranch        *b_rjrPVlab;   //!
   TBranch        *b_rjrPX1X2a;   //!
   TBranch        *b_rjrPX1X2b;   //!
   TBranch        *b_rjrPX2Sa;   //!
   TBranch        *b_rjrPX2Sb;   //!
   TBranch        *b_rjrPZS;   //!
   TBranch        *b_rjrSCosA;   //!
   TBranch        *b_rjrSMass;   //!
   TBranch        *b_rjrSdphiBV;   //!
   TBranch        *b_rjrSdphiDA;   //!
   TBranch        *b_rjrX1aCosA;   //!
   TBranch        *b_rjrX1aMass;   //!
   TBranch        *b_rjrX1bCosA;   //!
   TBranch        *b_rjrX1bMass;   //!
   TBranch        *b_rjrX2Ave;   //!
   TBranch        *b_rjrX2Diff;   //!
   TBranch        *b_rjrX2GMean;   //!
   TBranch        *b_rjrX2NAve;   //!
   TBranch        *b_rjrX2NGMean;   //!
   TBranch        *b_rjrX2NVSum;   //!
   TBranch        *b_rjrX2VSum;   //!
   TBranch        *b_rjrX2aCosA;   //!
   TBranch        *b_rjrX2aMass;   //!
   TBranch        *b_rjrX2aPtS;   //!
   TBranch        *b_rjrX2bCosA;   //!
   TBranch        *b_rjrX2bMass;   //!
   TBranch        *b_rjrX2bPtS;   //!

   skim_mod2(TTree *tree=0);
   virtual ~skim_mod2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef skim_mod2_cxx
skim_mod2::skim_mod2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rjr_skim_files/kuntuple_GJets_HT-100To200_Met0_R17_v19_Ntuple_v19_LLPgama_Skim_v19_Met0_rjrcombo.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rjr_skim_files/kuntuple_GJets_HT-100To200_Met0_R17_v19_Ntuple_v19_LLPgama_Skim_v19_Met0_rjrcombo.root");
      }
      f->GetObject("kuSkimTree",tree);

   }
   Init(tree);
}

skim_mod2::~skim_mod2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skim_mod2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skim_mod2::LoadTree(Long64_t entry)
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

void skim_mod2::Init(TTree *tree)
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
   genMass = 0;
   genPartEnergy = 0;
   genPartEta = 0;
   genPartPdgId = 0;
   genPartPhi = 0;
   genPartPt = 0;
   genPartSusId = 0;
   genPx = 0;
   genPy = 0;
   genPz = 0;
   genSGMass = 0;
   genSQMass = 0;
   genVx = 0;
   genVy = 0;
   genVz = 0;
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
   selJetLlpDp = 0;
   selJetLlpDr = 0;
   selJetMass = 0;
   selJetMuEF = 0;
   selJetNeEmEF = 0;
   selJetNeHEF = 0;
   selJetNeHM = 0;
   selJetPhi = 0;
   selJetPt = 0;
   selJetQuality = 0;
   selJetSusyId = 0;
   selJetTime = 0;
   selJetchHEF = 0;
   rjrASMass = 0;
   rjrAX2Ave = 0;
   rjrAX2Diff = 0;
   rjrAX2GMean = 0;
   rjrAX2NAve = 0;
   rjrAX2NGMean = 0;
   rjrAX2NVSum = 0;
   rjrAX2VSum = 0;
   rjrAX2aMass = 0;
   rjrAX2bMass = 0;
   rjrDiffPJX2a = 0;
   rjrDiffPJX2b = 0;
   rjrDphiMETV = 0;
   rjrEVa = 0;
   rjrEVb = 0;
   rjrMET = 0;
   rjrMV = 0;
   rjrMVa = 0;
   rjrMVb = 0;
   rjrNJetsJa = 0;
   rjrNJetsJb = 0;
   rjrNRjrJets = 0;
   rjrNRjrPhotons = 0;
   rjrPJX2a = 0;
   rjrPJX2b = 0;
   rjrPTS = 0;
   rjrPV = 0;
   rjrPVa = 0;
   rjrPVb = 0;
   rjrPVlab = 0;
   rjrPX1X2a = 0;
   rjrPX1X2b = 0;
   rjrPX2Sa = 0;
   rjrPX2Sb = 0;
   rjrPZS = 0;
   rjrSCosA = 0;
   rjrSMass = 0;
   rjrSdphiBV = 0;
   rjrSdphiDA = 0;
   rjrX1aCosA = 0;
   rjrX1aMass = 0;
   rjrX1bCosA = 0;
   rjrX1bMass = 0;
   rjrX2Ave = 0;
   rjrX2Diff = 0;
   rjrX2GMean = 0;
   rjrX2NAve = 0;
   rjrX2NGMean = 0;
   rjrX2NVSum = 0;
   rjrX2VSum = 0;
   rjrX2aCosA = 0;
   rjrX2aMass = 0;
   rjrX2aPtS = 0;
   rjrX2bCosA = 0;
   rjrX2bMass = 0;
   rjrX2bPtS = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("DataSetKey", &DataSetKey, &b_DataSetKey);
   fChain->SetBranchAddress("evtGenWgt", &evtGenWgt, &b_evtGenWgt);
   fChain->SetBranchAddress("evtXSection", &evtXSection, &b_evtXSection);
   fChain->SetBranchAddress("selCMet", &selCMet, &b_selCMet);
   fChain->SetBranchAddress("selCMetPx", &selCMetPx, &b_selCMetPx);
   fChain->SetBranchAddress("selCMetPy", &selCMetPy, &b_selCMetPy);
   fChain->SetBranchAddress("selMet", &selMet, &b_selMet);
   fChain->SetBranchAddress("selMetPx", &selMetPx, &b_selMetPx);
   fChain->SetBranchAddress("selMetPy", &selMetPy, &b_selMetPy);
   fChain->SetBranchAddress("genCharge", &genCharge, &b_genCharge);
   fChain->SetBranchAddress("genMass", &genMass, &b_genMass);
   fChain->SetBranchAddress("genPartEnergy", &genPartEnergy, &b_genPartEnergy);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPdgId", &genPartPdgId, &b_genPartPdgId);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartSusId", &genPartSusId, &b_genPartSusId);
   fChain->SetBranchAddress("genPx", &genPx, &b_genPx);
   fChain->SetBranchAddress("genPy", &genPy, &b_genPy);
   fChain->SetBranchAddress("genPz", &genPz, &b_genPz);
   fChain->SetBranchAddress("genSGMass", &genSGMass, &b_genSGMass);
   fChain->SetBranchAddress("genSQMass", &genSQMass, &b_genSQMass);
   fChain->SetBranchAddress("genSigType", &genSigType, &b_genSigType);
   fChain->SetBranchAddress("genStatus", &genStatus, &b_genStatus);
   fChain->SetBranchAddress("genVx", &genVx, &b_genVx);
   fChain->SetBranchAddress("genVy", &genVy, &b_genVy);
   fChain->SetBranchAddress("genVz", &genVz, &b_genVz);
   fChain->SetBranchAddress("elePhoMinDr", &elePhoMinDr, &b_elePhoMinDr);
   fChain->SetBranchAddress("epDrMatch", &epDrMatch, &b_epDrMatch);
   fChain->SetBranchAddress("epSeedIdMatch", &epSeedIdMatch, &b_epSeedIdMatch);
   fChain->SetBranchAddress("epSeedIdUnique", &epSeedIdUnique, &b_epSeedIdUnique);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
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
   fChain->SetBranchAddress("selJetLlpDp", &selJetLlpDp, &b_selJetLlpDp);
   fChain->SetBranchAddress("selJetLlpDr", &selJetLlpDr, &b_selJetLlpDr);
   fChain->SetBranchAddress("selJetMass", &selJetMass, &b_selJetMass);
   fChain->SetBranchAddress("selJetMuEF", &selJetMuEF, &b_selJetMuEF);
   fChain->SetBranchAddress("selJetNeEmEF", &selJetNeEmEF, &b_selJetNeEmEF);
   fChain->SetBranchAddress("selJetNeHEF", &selJetNeHEF, &b_selJetNeHEF);
   fChain->SetBranchAddress("selJetNeHM", &selJetNeHM, &b_selJetNeHM);
   fChain->SetBranchAddress("selJetPhi", &selJetPhi, &b_selJetPhi);
   fChain->SetBranchAddress("selJetPt", &selJetPt, &b_selJetPt);
   fChain->SetBranchAddress("selJetQuality", &selJetQuality, &b_selJetQuality);
   fChain->SetBranchAddress("selJetSusyId", &selJetSusyId, &b_selJetSusyId);
   fChain->SetBranchAddress("selJetTime", &selJetTime, &b_selJetTime);
   fChain->SetBranchAddress("selJetchHEF", &selJetchHEF, &b_selJetchHEF);
   fChain->SetBranchAddress("rjrASMass", &rjrASMass, &b_rjrASMass);
   fChain->SetBranchAddress("rjrAX2Ave", &rjrAX2Ave, &b_rjrAX2Ave);
   fChain->SetBranchAddress("rjrAX2Diff", &rjrAX2Diff, &b_rjrAX2Diff);
   fChain->SetBranchAddress("rjrAX2GMean", &rjrAX2GMean, &b_rjrAX2GMean);
   fChain->SetBranchAddress("rjrAX2NAve", &rjrAX2NAve, &b_rjrAX2NAve);
   fChain->SetBranchAddress("rjrAX2NGMean", &rjrAX2NGMean, &b_rjrAX2NGMean);
   fChain->SetBranchAddress("rjrAX2NVSum", &rjrAX2NVSum, &b_rjrAX2NVSum);
   fChain->SetBranchAddress("rjrAX2VSum", &rjrAX2VSum, &b_rjrAX2VSum);
   fChain->SetBranchAddress("rjrAX2aMass", &rjrAX2aMass, &b_rjrAX2aMass);
   fChain->SetBranchAddress("rjrAX2bMass", &rjrAX2bMass, &b_rjrAX2bMass);
   fChain->SetBranchAddress("rjrDiffPJX2a", &rjrDiffPJX2a, &b_rjrDiffPJX2a);
   fChain->SetBranchAddress("rjrDiffPJX2b", &rjrDiffPJX2b, &b_rjrDiffPJX2b);
   fChain->SetBranchAddress("rjrDphiMETV", &rjrDphiMETV, &b_rjrDphiMETV);
   fChain->SetBranchAddress("rjrEVa", &rjrEVa, &b_rjrEVa);
   fChain->SetBranchAddress("rjrEVb", &rjrEVb, &b_rjrEVb);
   fChain->SetBranchAddress("rjrMET", &rjrMET, &b_rjrMET);
   fChain->SetBranchAddress("rjrMV", &rjrMV, &b_rjrMV);
   fChain->SetBranchAddress("rjrMVa", &rjrMVa, &b_rjrMVa);
   fChain->SetBranchAddress("rjrMVb", &rjrMVb, &b_rjrMVb);
   fChain->SetBranchAddress("rjrNJetsJa", &rjrNJetsJa, &b_rjrNJetsJa);
   fChain->SetBranchAddress("rjrNJetsJb", &rjrNJetsJb, &b_rjrNJetsJb);
   fChain->SetBranchAddress("rjrNRjrJets", &rjrNRjrJets, &b_rjrNRjrJets);
   fChain->SetBranchAddress("rjrNRjrPhotons", &rjrNRjrPhotons, &b_rjrNRjrPhotons);
   fChain->SetBranchAddress("rjrPJX2a", &rjrPJX2a, &b_rjrPJX2a);
   fChain->SetBranchAddress("rjrPJX2b", &rjrPJX2b, &b_rjrPJX2b);
   fChain->SetBranchAddress("rjrPTS", &rjrPTS, &b_rjrPTS);
   fChain->SetBranchAddress("rjrPV", &rjrPV, &b_rjrPV);
   fChain->SetBranchAddress("rjrPVa", &rjrPVa, &b_rjrPVa);
   fChain->SetBranchAddress("rjrPVb", &rjrPVb, &b_rjrPVb);
   fChain->SetBranchAddress("rjrPVlab", &rjrPVlab, &b_rjrPVlab);
   fChain->SetBranchAddress("rjrPX1X2a", &rjrPX1X2a, &b_rjrPX1X2a);
   fChain->SetBranchAddress("rjrPX1X2b", &rjrPX1X2b, &b_rjrPX1X2b);
   fChain->SetBranchAddress("rjrPX2Sa", &rjrPX2Sa, &b_rjrPX2Sa);
   fChain->SetBranchAddress("rjrPX2Sb", &rjrPX2Sb, &b_rjrPX2Sb);
   fChain->SetBranchAddress("rjrPZS", &rjrPZS, &b_rjrPZS);
   fChain->SetBranchAddress("rjrSCosA", &rjrSCosA, &b_rjrSCosA);
   fChain->SetBranchAddress("rjrSMass", &rjrSMass, &b_rjrSMass);
   fChain->SetBranchAddress("rjrSdphiBV", &rjrSdphiBV, &b_rjrSdphiBV);
   fChain->SetBranchAddress("rjrSdphiDA", &rjrSdphiDA, &b_rjrSdphiDA);
   fChain->SetBranchAddress("rjrX1aCosA", &rjrX1aCosA, &b_rjrX1aCosA);
   fChain->SetBranchAddress("rjrX1aMass", &rjrX1aMass, &b_rjrX1aMass);
   fChain->SetBranchAddress("rjrX1bCosA", &rjrX1bCosA, &b_rjrX1bCosA);
   fChain->SetBranchAddress("rjrX1bMass", &rjrX1bMass, &b_rjrX1bMass);
   fChain->SetBranchAddress("rjrX2Ave", &rjrX2Ave, &b_rjrX2Ave);
   fChain->SetBranchAddress("rjrX2Diff", &rjrX2Diff, &b_rjrX2Diff);
   fChain->SetBranchAddress("rjrX2GMean", &rjrX2GMean, &b_rjrX2GMean);
   fChain->SetBranchAddress("rjrX2NAve", &rjrX2NAve, &b_rjrX2NAve);
   fChain->SetBranchAddress("rjrX2NGMean", &rjrX2NGMean, &b_rjrX2NGMean);
   fChain->SetBranchAddress("rjrX2NVSum", &rjrX2NVSum, &b_rjrX2NVSum);
   fChain->SetBranchAddress("rjrX2VSum", &rjrX2VSum, &b_rjrX2VSum);
   fChain->SetBranchAddress("rjrX2aCosA", &rjrX2aCosA, &b_rjrX2aCosA);
   fChain->SetBranchAddress("rjrX2aMass", &rjrX2aMass, &b_rjrX2aMass);
   fChain->SetBranchAddress("rjrX2aPtS", &rjrX2aPtS, &b_rjrX2aPtS);
   fChain->SetBranchAddress("rjrX2bCosA", &rjrX2bCosA, &b_rjrX2bCosA);
   fChain->SetBranchAddress("rjrX2bMass", &rjrX2bMass, &b_rjrX2bMass);
   fChain->SetBranchAddress("rjrX2bPtS", &rjrX2bPtS, &b_rjrX2bPtS);
   Notify();
}

Bool_t skim_mod2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skim_mod2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skim_mod2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skim_mod2_cxx
