//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  8 11:42:53 2023 by ROOT version 6.26/07
// from TTree llpgtree/KUCMSNtuple
// found on file: gmsb_AODSIM_KUCMSNtuplizer_Objectified_v8.root
//////////////////////////////////////////////////////////

#ifndef root_base_h
#define root_base_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

//#define sbDEBUG true
#define sbDEBUG false

class root_base {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<float>   *ECALRecHit_energy;
   std::vector<unsigned int> *ECALRecHit_ID;
   std::vector<float>   *ECALRecHit_swCross;
   std::vector<float>   *ECALRecHit_TOF;
   std::vector<float>   *ECALRecHit_time;
   std::vector<float>   *ECALRecHit_eta;
   std::vector<bool>    *ECALRecHit_isOOT;
   std::vector<float>   *ECALRecHit_phi;
   std::vector<float>   *ECALRecHit_rhx;
   std::vector<float>   *ECALRecHit_rhy;
   std::vector<float>   *ECALRecHit_rhz;

   std::vector<float>   *Electron_energy;
   std::vector<float>   *Electron_eta;
   std::vector<float>   *Electron_genDp;
   std::vector<float>   *Electron_genDr;
   std::vector<int>     *Electron_genIdx;
   std::vector<float>   *Electron_genSDp;
   std::vector<int>     *Electron_genSIdx;
   std::vector<float>   *Electron_phi;
   std::vector<float>   *Electron_pt;
   std::vector<float>   *Electron_px;
   std::vector<float>   *Electron_py;
   std::vector<float>   *Electron_pz;
   std::vector<std::vector<unsigned int> > *Electron_rhIds;
   std::vector<float>   *Electron_seedTOFTime;

   UInt_t          Evt_luminosityBlock;
   UInt_t          Evt_run;
   UInt_t          Evt_event;

   UInt_t          PV_npvs;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;

   std::vector<int>     *Gen_charge;
   std::vector<float>   *Gen_energy;
   std::vector<float>   *Gen_eta;
   std::vector<float>   *Gen_mass;
   std::vector<unsigned int> *Gen_pdgId;
   std::vector<float>   *Gen_phi;
   std::vector<float>   *Gen_pt;
   std::vector<float>   *Gen_px;
   std::vector<float>   *Gen_py;
   std::vector<float>   *Gen_pz;
   std::vector<bool>    *Gen_status;
   std::vector<int>     *Gen_susId;
   std::vector<float>   *Gen_vx;
   std::vector<float>   *Gen_vy;
   std::vector<float>   *Gen_vz;
   Float_t         Evt_genWgt;

   std::vector<float>   *Jet_area;
   std::vector<float>   *Jet_chEmEF;
   std::vector<float>   *Jet_chHEF;
   std::vector<float>   *Jet_chHM;
   std::vector<std::vector<unsigned int> > *Jet_drRhIds;
   std::vector<float>   *Jet_energy;
   std::vector<float>   *Jet_eta;
   std::vector<float>   *Jet_genDptMatch;
   std::vector<float>   *Jet_genDrMatch;
   std::vector<float>   *Jet_genEnergy;
   std::vector<float>   *Jet_genEta;
   std::vector<float>   *Jet_genImpactAngle;
   std::vector<float>   *Jet_genLlpDp;
   std::vector<float>   *Jet_genLlpDr;
   std::vector<float>   *Jet_genLlpId;
   std::vector<float>   *Jet_genPhi;
   std::vector<float>   *Jet_genPt;
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

   std::vector<float>   *Photon_covEtaEta;
   std::vector<float>   *Photon_covEtaPhi;
   std::vector<float>   *Photon_covPhiPhi;
   std::vector<float>   *Photon_ecalRHSumEtConeDR04;
   std::vector<float>   *Photon_energy;
   std::vector<float>   *Photon_energyErr;
   std::vector<float>   *Photon_energyRaw;
   std::vector<float>   *Photon_eta;
   std::vector<bool>    *Photon_excluded;
   std::vector<int>     *Photon_genIdx;
   std::vector<float>   *Photon_genLlpId;
   std::vector<float>   *Photon_hadOverEM;
   std::vector<float>   *Photon_hadTowOverEM;
   std::vector<float>   *Photon_hcalTowerSumEtBcConeDR04;
   std::vector<float>   *Photon_hcalTowerSumEtConeDR04;   //!
   std::vector<float>   *Photon_nTrkHollowConeDR04;   //!
   std::vector<float>   *Photon_nTrkSolidConeDR04;   //!
   std::vector<bool>    *Photon_isOot;
   std::vector<float>   *Photon_phi;
   std::vector<float>   *Photon_pt;
   std::vector<float>   *Photon_px;
   std::vector<float>   *Photon_py;
   std::vector<float>   *Photon_pz;
   std::vector<float>   *Photon_r9;
   std::vector<std::vector<unsigned int> > *Photon_rhIds;
   std::vector<float>   *Photon_s4;
   std::vector<float>   *Photon_salp;
   std::vector<float>   *Photon_smaj;
   std::vector<float>   *Photon_smin;
   std::vector<float>   *Photon_seedTOFTime;
   std::vector<float>   *Photon_SigmaIEtaIEta;
   std::vector<float>   *Photon_trkSumPtHollowConeDR03;
   std::vector<float>   *Photon_trkSumPtHollowConeDR04;
   std::vector<float>   *Photon_trkSumPtSolidConeDR04;
   std::vector<float>   *Photon_ecalPFClusterIso;
   std::vector<bool>    *Photon_electronVeto;
   std::vector<float>   *Photon_etaWidth;
   std::vector<bool>    *hasConversionTracks;
   std::vector<bool>    *Photon_pixelSeed;
   std::vector<float>   *Photon_hcalPFClusterIso;
   std::vector<float>   *Photon_Hoe_PUcorr;
   std::vector<bool>    *Photon_seedIsEB;
   std::vector<bool>    *Photon_isScEtaEB;
   std::vector<bool>    *Photon_isScEtaEE;
   std::vector<float>   *Photon_genSigMomEnergy;
   std::vector<float>   *Photon_genSigMomEta;
   std::vector<float>   *Photon_genSigMomMass;
   std::vector<float>   *Photon_genSigMomPhi;
   std::vector<float>   *Photon_genSigMomPt;
   std::vector<float>   *Photon_genSigMomPx;
   std::vector<float>   *Photon_genSigMomPy;
   std::vector<float>   *Photon_genSigMomPz;
   std::vector<float>   *Photon_genSigMomVx;
   std::vector<float>   *Photon_genSigMomVy;
   std::vector<float>   *Photon_genSigMomVz;
   std::vector<float>   *Photon_pfChargedIso;
   std::vector<float>   *Photon_pfChargedIsoPFPV;
   std::vector<float>   *Photon_pfChargedIsoWorstVtx;
   std::vector<float>   *Photon_pfPhoIso03;
   std::vector<float>   *pfRelIso03_all_quadratic;
   std::vector<float>   *pfRelIso03_chg_quadratic;
   std::vector<float>   *Photon_phiWidth;
   std::vector<int>     *Photon_seediEtaOriX;
   std::vector<int>     *Photon_seediPhiOriY;
   std::vector<float>   *Photon_sieie;
   std::vector<float>   *Photon_sieip;
   std::vector<float>   *Photon_sipip;
   std::vector<float>   *Photon_x_calo;
   std::vector<float>   *Photon_y_calo;
   std::vector<float>   *Photon_z_calo;

   // List of branches
   TBranch        *b_ECALRecHit_energy;   //!
   TBranch        *b_ECALRecHit_ID;   //!
   TBranch        *b_ECALRecHit_swCross;   //!
   TBranch        *b_ECALRecHit_TOF;   //!
   TBranch        *b_ECALRecHit_time;   //!
   TBranch        *b_ECALRecHit_eta;   //!
   TBranch        *b_ECALRecHit_isOOT;   //!
   TBranch        *b_ECALRecHit_phi;   //!
   TBranch        *b_ECALRecHit_rhx;   //!
   TBranch        *b_ECALRecHit_rhy;   //!
   TBranch        *b_ECALRecHit_rhz;   //!

   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_genDp;   //!
   TBranch        *b_Electron_genDr;   //!
   TBranch        *b_Electron_genIdx;   //!
   TBranch        *b_Electron_genSDp;   //!
   TBranch        *b_Electron_genSIdx;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_px;   //!
   TBranch        *b_Electron_py;   //!
   TBranch        *b_Electron_pz;   //!
   TBranch        *b_Electron_rhIds;   //!
   TBranch        *b_Electron_seedTOFTime;   //!

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
   TBranch        *b_Gen_pdgId;   //!
   TBranch        *b_Gen_phi;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_px;   //!
   TBranch        *b_Gen_py;   //!
   TBranch        *b_Gen_pz;   //!
   TBranch        *b_Gen_status;   //!
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

   TBranch        *b_Photon_covEtaEta;   //!
   TBranch        *b_Photon_covEtaPhi;   //!
   TBranch        *b_Photon_covPhiPhi;   //!
   TBranch        *b_Photon_ecalRHSumEtConeDR04;   //!
   TBranch        *b_Photon_energy;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_energyRaw;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_excluded;   //!
   TBranch        *b_Photon_genIdx;   //!
   TBranch        *b_Photon_genLlpId;   //!
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
   TBranch        *b_Photon_rhIds;   //!
   TBranch        *b_Photon_s4;   //!
   TBranch        *b_Photon_salp;   //!
   TBranch        *b_Photon_smaj;   //!
   TBranch        *b_Photon_smin;   //!
   TBranch        *b_Photon_seedTOFTime;   //!
   TBranch        *b_Photon_SigmaIEtaIEta;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_ecalPFClusterIso;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_etaWidth;   //!
   TBranch        *b_hasConversionTracks;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_hcalPFClusterIso;   //!
   TBranch        *b_Photon_Hoe_PUcorr;   //! 
   TBranch        *b_Photon_seedIsEB;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_genSigMomEnergy;   //!
   TBranch        *b_Photon_genSigMomEta;   //!
   TBranch        *b_Photon_genSigMomMass;   //!
   TBranch        *b_Photon_genSigMomPhi;   //!
   TBranch        *b_Photon_genSigMomPt;   //!
   TBranch        *b_Photon_genSigMomPx;   //!
   TBranch        *b_Photon_genSigMomPy;   //!
   TBranch        *b_Photon_genSigMomPz;   //!
   TBranch        *b_Photon_genSigMomVx;   //!
   TBranch        *b_Photon_genSigMomVy;   //!
   TBranch        *b_Photon_genSigMomVz;   //!
   TBranch        *b_Photon_pfChargedIso;   //!
   TBranch        *b_Photon_pfChargedIsoPFPV;   //!
   TBranch        *b_Photon_pfChargedIsoWorstVtx;   //!
   TBranch        *b_Photon_pfPhoIso03;   //!
   TBranch        *b_pfRelIso03_all_quadratic;   //!
   TBranch        *b_pfRelIso03_chg_quadratic;   //!
   TBranch        *b_Photon_phiWidth;   //!
   TBranch        *b_Photon_seediEtaOriX;   //!
   TBranch        *b_Photon_seediPhiOriY;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_sieip;   //!
   TBranch        *b_Photon_sipip;   //!
   TBranch        *b_Photon_x_calo;   //!
   TBranch        *b_Photon_y_calo;   //!
   TBranch        *b_Photon_z_calo;   //!

//   root_base(TTree *tree=0);
//   virtual ~root_base();
//   virtual Int_t    Cut(Long64_t entry);
//   virtual Int_t    GetEntry(Long64_t entry);
//   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init( TTree *tree, bool doGenInfo );
   virtual void     getBranches( Long64_t entry, bool doGenInfo );
//   virtual void     Loop();
//   virtual Bool_t   Notify();
//   virtual void     Show(Long64_t entry = -1);


};

//#endif

//#ifdef root_base_cxx
/*
root_base::root_base(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v8.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v8.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v8.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

root_base::~root_base()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t root_base::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t root_base::LoadTree(Long64_t entry)
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

void root_base::Init(TTree *tree, bool doGenInfo )
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   Evt_genWgt = 1;

   // Set object pointer
   ECALRecHit_energy = 0;
   ECALRecHit_ID = 0;
   ECALRecHit_swCross = 0;
   ECALRecHit_TOF = 0;
   ECALRecHit_time = 0;
   ECALRecHit_eta = 0;
   ECALRecHit_isOOT = 0;
   ECALRecHit_phi = 0;
   ECALRecHit_rhx = 0;
   ECALRecHit_rhy = 0;
   ECALRecHit_rhz = 0;

   Electron_energy = 0;
   Electron_eta = 0;
   Electron_genDp = 0;
   Electron_genDr = 0;
   Electron_genIdx = 0;
   Electron_genSDp = 0;
   Electron_genSIdx = 0;
   Electron_phi = 0;
   Electron_pt = 0;
   Electron_px = 0;
   Electron_py = 0;
   Electron_pz = 0;
   Electron_rhIds = 0;
   Electron_seedTOFTime = 0;

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

   Photon_covEtaEta = 0;
   Photon_covEtaPhi = 0;
   Photon_covPhiPhi = 0;
   Photon_ecalRHSumEtConeDR04 = 0;
   Photon_energy = 0;
   Photon_energyErr = 0;
   Photon_energyRaw = 0;
   Photon_eta = 0;
   Photon_excluded = 0;
   Photon_genIdx = 0;
   Photon_genLlpId = 0;
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
   Photon_rhIds = 0;
   Photon_s4 = 0;
   Photon_salp = 0;
   Photon_smaj = 0;
   Photon_smin = 0;
   Photon_seedTOFTime = 0;
   Photon_SigmaIEtaIEta = 0;
   Photon_trkSumPtHollowConeDR03 = 0;
   Photon_trkSumPtHollowConeDR04 = 0;
   Photon_trkSumPtSolidConeDR04 = 0;
   Photon_ecalPFClusterIso = 0;
   Photon_electronVeto = 0;
   Photon_etaWidth = 0;
   hasConversionTracks = 0;
   Photon_pixelSeed = 0;
   Photon_hcalPFClusterIso = 0;
   Photon_Hoe_PUcorr = 0;
   Photon_seedIsEB = 0;
   Photon_isScEtaEB = 0;
   Photon_isScEtaEE = 0;
   Photon_genSigMomEnergy = 0;
   Photon_genSigMomEta = 0;
   Photon_genSigMomMass = 0;
   Photon_genSigMomPhi = 0;
   Photon_genSigMomPt = 0;
   Photon_genSigMomPx = 0;
   Photon_genSigMomPy = 0;
   Photon_genSigMomPz = 0;
   Photon_genSigMomVx = 0;
   Photon_genSigMomVy = 0;
   Photon_genSigMomVz = 0;
   Photon_pfChargedIso = 0;
   Photon_pfChargedIsoPFPV = 0;
   Photon_pfChargedIsoWorstVtx = 0;
   Photon_pfPhoIso03 = 0;
   pfRelIso03_all_quadratic = 0;
   pfRelIso03_chg_quadratic = 0;
   Photon_phiWidth = 0;
   Photon_seediEtaOriX = 0;
   Photon_seediPhiOriY = 0;
   Photon_sieie = 0;
   Photon_sieip = 0;
   Photon_sipip = 0;
   Photon_x_calo = 0;
   Photon_y_calo = 0;
   Photon_z_calo = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ECALRecHit_energy", &ECALRecHit_energy, &b_ECALRecHit_energy);
   fChain->SetBranchAddress("ECALRecHit_ID", &ECALRecHit_ID, &b_ECALRecHit_ID);
   fChain->SetBranchAddress("ECALRecHit_swCross", &ECALRecHit_swCross, &b_ECALRecHit_swCross);
   fChain->SetBranchAddress("ECALRecHit_TOF", &ECALRecHit_TOF, &b_ECALRecHit_TOF);
   fChain->SetBranchAddress("ECALRecHit_time", &ECALRecHit_time, &b_ECALRecHit_time);
   fChain->SetBranchAddress("ECALRecHit_eta", &ECALRecHit_eta, &b_ECALRecHit_eta);
   fChain->SetBranchAddress("ECALRecHit_isOOT", &ECALRecHit_isOOT, &b_ECALRecHit_isOOT);
   fChain->SetBranchAddress("ECALRecHit_phi", &ECALRecHit_phi, &b_ECALRecHit_phi);
   fChain->SetBranchAddress("ECALRecHit_rhx", &ECALRecHit_rhx, &b_ECALRecHit_rhx);
   fChain->SetBranchAddress("ECALRecHit_rhy", &ECALRecHit_rhy, &b_ECALRecHit_rhy);
   fChain->SetBranchAddress("ECALRecHit_rhz", &ECALRecHit_rhz, &b_ECALRecHit_rhz);

   fChain->SetBranchAddress("Electron_energy", &Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_genDp", &Electron_genDp, &b_Electron_genDp);
   fChain->SetBranchAddress("Electron_genDr", &Electron_genDr, &b_Electron_genDr);
   fChain->SetBranchAddress("Electron_genIdx", &Electron_genIdx, &b_Electron_genIdx);
   fChain->SetBranchAddress("Electron_genSDp", &Electron_genSDp, &b_Electron_genSDp);
   fChain->SetBranchAddress("Electron_genSIdx", &Electron_genSIdx, &b_Electron_genSIdx);
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_px", &Electron_px, &b_Electron_px);
   fChain->SetBranchAddress("Electron_py", &Electron_py, &b_Electron_py);
   fChain->SetBranchAddress("Electron_pz", &Electron_pz, &b_Electron_pz);
   fChain->SetBranchAddress("Electron_rhIds", &Electron_rhIds, &b_Electron_rhIds);
   fChain->SetBranchAddress("Electron_seedTOFTime", &Electron_seedTOFTime, &b_Electron_seedTOFTime);

   fChain->SetBranchAddress("Evt_luminosityBlock", &Evt_luminosityBlock, &b_Evt_luminosityBlock);
   fChain->SetBranchAddress("Evt_run", &Evt_run, &b_Evt_run);
   fChain->SetBranchAddress("Evt_event", &Evt_event, &b_Evt_event);

   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);

   if( doGenInfo ){
   fChain->SetBranchAddress("Gen_charge", &Gen_charge, &b_Gen_charge);
   fChain->SetBranchAddress("Gen_energy", &Gen_energy, &b_Gen_energy);
   fChain->SetBranchAddress("Gen_eta", &Gen_eta, &b_Gen_eta);
   fChain->SetBranchAddress("Gen_mass", &Gen_mass, &b_Gen_mass);
   fChain->SetBranchAddress("Gen_pdgId", &Gen_pdgId, &b_Gen_pdgId);
   fChain->SetBranchAddress("Gen_phi", &Gen_phi, &b_Gen_phi);
   fChain->SetBranchAddress("Gen_pt", &Gen_pt, &b_Gen_pt);
   fChain->SetBranchAddress("Gen_px", &Gen_px, &b_Gen_px);
   fChain->SetBranchAddress("Gen_py", &Gen_py, &b_Gen_py);
   fChain->SetBranchAddress("Gen_pz", &Gen_pz, &b_Gen_pz);
   fChain->SetBranchAddress("Gen_status", &Gen_status, &b_Gen_status);
   fChain->SetBranchAddress("Gen_susId", &Gen_susId, &b_Gen_susId);
   fChain->SetBranchAddress("Gen_vx", &Gen_vx, &b_Gen_vx);
   fChain->SetBranchAddress("Gen_vy", &Gen_vy, &b_Gen_vy);
   fChain->SetBranchAddress("Gen_vz", &Gen_vz, &b_Gen_vz);
   fChain->SetBranchAddress("Evt_genWgt", &Evt_genWgt, &b_Evt_genWgt);
   }//if( doGenInfo )

   fChain->SetBranchAddress("Jet_area", &Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_chEmEF", &Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", &Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_chHM", &Jet_chHM, &b_Jet_chHM);
   fChain->SetBranchAddress("Jet_drRhIds", &Jet_drRhIds, &b_Jet_drRhIds);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);

   if( doGenInfo ){
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
   }//if( doGenInfo )

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

   fChain->SetBranchAddress("Photon_covEtaEta", &Photon_covEtaEta, &b_Photon_covEtaEta);
   fChain->SetBranchAddress("Photon_covEtaPhi", &Photon_covEtaPhi, &b_Photon_covEtaPhi);
   fChain->SetBranchAddress("Photon_covPhiPhi", &Photon_covPhiPhi, &b_Photon_covPhiPhi);
   fChain->SetBranchAddress("Photon_ecalRHSumEtConeDR04", &Photon_ecalRHSumEtConeDR04, &b_Photon_ecalRHSumEtConeDR04);
   fChain->SetBranchAddress("Photon_energy", &Photon_energy, &b_Photon_energy);
   fChain->SetBranchAddress("Photon_energyErr", &Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_energyRaw", &Photon_energyRaw, &b_Photon_energyRaw);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_excluded", &Photon_excluded, &b_Photon_excluded);

   if( doGenInfo ){
   fChain->SetBranchAddress("Photon_genIdx", &Photon_genIdx, &b_Photon_genIdx);
   fChain->SetBranchAddress("Photon_genLlpId", &Photon_genLlpId, &b_Photon_genLlpId);
   fChain->SetBranchAddress("Photon_genSigMomEnergy", &Photon_genSigMomEnergy, &b_Photon_genSigMomEnergy);
   fChain->SetBranchAddress("Photon_genSigMomEta", &Photon_genSigMomEta, &b_Photon_genSigMomEta);
   fChain->SetBranchAddress("Photon_genSigMomMass", &Photon_genSigMomMass, &b_Photon_genSigMomMass);
   fChain->SetBranchAddress("Photon_genSigMomPhi", &Photon_genSigMomPhi, &b_Photon_genSigMomPhi);
   fChain->SetBranchAddress("Photon_genSigMomPt", &Photon_genSigMomPt, &b_Photon_genSigMomPt);
   fChain->SetBranchAddress("Photon_genSigMomPx", &Photon_genSigMomPx, &b_Photon_genSigMomPx);
   fChain->SetBranchAddress("Photon_genSigMomPy", &Photon_genSigMomPy, &b_Photon_genSigMomPy);
   fChain->SetBranchAddress("Photon_genSigMomPz", &Photon_genSigMomPz, &b_Photon_genSigMomPz);
   fChain->SetBranchAddress("Photon_genSigMomVx", &Photon_genSigMomVx, &b_Photon_genSigMomVx);
   fChain->SetBranchAddress("Photon_genSigMomVy", &Photon_genSigMomVy, &b_Photon_genSigMomVy);
   fChain->SetBranchAddress("Photon_genSigMomVz", &Photon_genSigMomVz, &b_Photon_genSigMomVz);  
   }//if( doGenInfo )

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
   fChain->SetBranchAddress("Photon_rhIds", &Photon_rhIds, &b_Photon_rhIds);
   fChain->SetBranchAddress("Photon_s4", &Photon_s4, &b_Photon_s4);
   fChain->SetBranchAddress("Photon_salp", &Photon_salp, &b_Photon_salp);
   fChain->SetBranchAddress("Photon_smaj", &Photon_smaj, &b_Photon_smaj);
   fChain->SetBranchAddress("Photon_smin", &Photon_smin, &b_Photon_smin);
   fChain->SetBranchAddress("Photon_seedTOFTime", &Photon_seedTOFTime, &b_Photon_seedTOFTime);
   fChain->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", &Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", &Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", &Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_ecalPFClusterIso", &Photon_ecalPFClusterIso, &b_Photon_ecalPFClusterIso);
   fChain->SetBranchAddress("Photon_electronVeto", &Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_etaWidth", &Photon_etaWidth, &b_Photon_etaWidth);
   fChain->SetBranchAddress("hasConversionTracks", &hasConversionTracks, &b_hasConversionTracks);
   fChain->SetBranchAddress("Photon_pixelSeed", &Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", &Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_Hoe_PUcorr", &Photon_Hoe_PUcorr, &b_Photon_Hoe_PUcorr);
   fChain->SetBranchAddress("Photon_seedIsEB", &Photon_seedIsEB, &b_Photon_seedIsEB);
   fChain->SetBranchAddress("Photon_isScEtaEB", &Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", &Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_pfChargedIso", &Photon_pfChargedIso, &b_Photon_pfChargedIso);
   fChain->SetBranchAddress("Photon_pfChargedIsoPFPV", &Photon_pfChargedIsoPFPV, &b_Photon_pfChargedIsoPFPV);
   fChain->SetBranchAddress("Photon_pfChargedIsoWorstVtx", &Photon_pfChargedIsoWorstVtx, &b_Photon_pfChargedIsoWorstVtx);
   fChain->SetBranchAddress("Photon_pfPhoIso03", &Photon_pfPhoIso03, &b_Photon_pfPhoIso03);
   fChain->SetBranchAddress("pfRelIso03_all_quadratic", &pfRelIso03_all_quadratic, &b_pfRelIso03_all_quadratic);
   fChain->SetBranchAddress("pfRelIso03_chg_quadratic", &pfRelIso03_chg_quadratic, &b_pfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("Photon_phiWidth", &Photon_phiWidth, &b_Photon_phiWidth);
   fChain->SetBranchAddress("Photon_seediEtaOriX", &Photon_seediEtaOriX, &b_Photon_seediEtaOriX);
   fChain->SetBranchAddress("Photon_seediPhiOriY", &Photon_seediPhiOriY, &b_Photon_seediPhiOriY);
   fChain->SetBranchAddress("Photon_sieie", &Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_sieip", &Photon_sieip, &b_Photon_sieip);
   fChain->SetBranchAddress("Photon_sipip", &Photon_sipip, &b_Photon_sipip);
   fChain->SetBranchAddress("Photon_x_calo", &Photon_x_calo, &b_Photon_x_calo);
   fChain->SetBranchAddress("Photon_y_calo", &Photon_y_calo, &b_Photon_y_calo);
   fChain->SetBranchAddress("Photon_z_calo", &Photon_z_calo, &b_Photon_z_calo);
   //Notify();
}

//   b_selPhoClstrRn->GetEntry(entry);    //!
void root_base::getBranches( Long64_t entry, bool doGenInfo ){

   if( sbDEBUG ) std::cout << "Getting Branches ECALRecHit" << std::endl;	
   b_ECALRecHit_energy->GetEntry(entry);   //!
   b_ECALRecHit_ID->GetEntry(entry);   //!
   b_ECALRecHit_swCross->GetEntry(entry);   //!
   b_ECALRecHit_TOF->GetEntry(entry);   //!
   b_ECALRecHit_time->GetEntry(entry);   //!
   b_ECALRecHit_eta->GetEntry(entry);   //!
   b_ECALRecHit_isOOT->GetEntry(entry);   //!
   b_ECALRecHit_phi->GetEntry(entry);   //!
   b_ECALRecHit_rhx->GetEntry(entry);   //!
   b_ECALRecHit_rhy->GetEntry(entry);   //!
   b_ECALRecHit_rhz->GetEntry(entry);   //!

   if( sbDEBUG ) std::cout << "Getting Branches Electron" << std::endl;
   b_Electron_energy->GetEntry(entry);   //!
   b_Electron_eta->GetEntry(entry);   //!

   if( doGenInfo ){
   b_Electron_genDp->GetEntry(entry);   //!
   b_Electron_genDr->GetEntry(entry);   //!
   b_Electron_genIdx->GetEntry(entry);   //!
   b_Electron_genSDp->GetEntry(entry);   //!
   b_Electron_genSIdx->GetEntry(entry);   //!
   }//if( doGenInfo )

   b_Electron_phi->GetEntry(entry);   //!
   b_Electron_pt->GetEntry(entry);   //!
   b_Electron_px->GetEntry(entry);   //!
   b_Electron_py->GetEntry(entry);   //!
   b_Electron_pz->GetEntry(entry);   //!
   b_Electron_rhIds->GetEntry(entry);   //!
   b_Electron_seedTOFTime->GetEntry(entry);   //!

   if( sbDEBUG ) std::cout << "Getting Branches PV & EVT" << std::endl;
   b_Evt_luminosityBlock->GetEntry(entry);   //!
   b_Evt_run->GetEntry(entry);   //!
   b_Evt_event->GetEntry(entry);   //!

   b_PV_npvs->GetEntry(entry);   //!
   b_PV_x->GetEntry(entry);   //!
   b_PV_y->GetEntry(entry);   //!
   b_PV_z->GetEntry(entry);   //!

   if( sbDEBUG ) std::cout << "Getting Branches Gen" << std::endl;
   if( doGenInfo ){
   b_Gen_energy->GetEntry(entry);   //!
   b_Gen_eta->GetEntry(entry);   //!
   b_Gen_pdgId->GetEntry(entry);   //!
   b_Gen_phi->GetEntry(entry);   //!
   b_Gen_pt->GetEntry(entry);   //!
   b_Gen_px->GetEntry(entry);   //!
   b_Gen_py->GetEntry(entry);   //!
   b_Gen_pz->GetEntry(entry);   //!
   b_Gen_susId->GetEntry(entry);   //!
   b_Evt_genWgt->GetEntry(entry);   //!
   b_Gen_charge->GetEntry(entry);   //!
   b_Gen_mass->GetEntry(entry);   //!   
   b_Gen_status->GetEntry(entry);   //!
   b_Gen_vx->GetEntry(entry);   //!
   b_Gen_vy->GetEntry(entry);   //!
   b_Gen_vz->GetEntry(entry);   //!
   }//if( doGenInfo )

   if( sbDEBUG ) std::cout << "Getting Branches Jet" << std::endl;
   b_Jet_area->GetEntry(entry);   //!
   b_Jet_chEmEF->GetEntry(entry);   //!
   b_Jet_chHEF->GetEntry(entry);   //!
   b_Jet_chHM->GetEntry(entry);   //!
   b_Jet_drRhIds->GetEntry(entry);   //!
   b_Jet_energy->GetEntry(entry);   //!
   b_Jet_eta->GetEntry(entry);   //!

   if( doGenInfo ){
   b_Jet_genDptMatch->GetEntry(entry);   //!
   b_Jet_genDrMatch->GetEntry(entry);   //!
   b_Jet_genEnergy->GetEntry(entry);   //!
   b_Jet_genEta->GetEntry(entry);   //!
   b_Jet_genImpactAngle->GetEntry(entry);   //!
   b_Jet_genLlpDp->GetEntry(entry);   //!
   b_Jet_genLlpDr->GetEntry(entry);   //!
   b_Jet_genLlpId->GetEntry(entry);   //!
   b_Jet_genPhi->GetEntry(entry);   //!
   b_Jet_genPt->GetEntry(entry);   //!
   b_Jet_genTOF->GetEntry(entry);   //!
   b_Jet_genTime->GetEntry(entry);   //!
   b_Jet_genTimeLLP->GetEntry(entry);   //!
   }//if( doGenInfo )

   b_Jet_mass->GetEntry(entry);   //!
   b_Jet_muEF->GetEntry(entry);   //!
   b_Jet_neEmEF->GetEntry(entry);   //!
   b_Jet_neHEF->GetEntry(entry);   //!
   b_Jet_neHM->GetEntry(entry);   //!
   b_Jet_egIndxs->GetEntry(entry);   //!
   b_Jet_phi->GetEntry(entry);   //!
   b_Jet_pt->GetEntry(entry);   //!
   b_Jet_nConstituents->GetEntry(entry);   //!

   if( sbDEBUG ) std::cout << "Getting Branches met" << std::endl;
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

   if( sbDEBUG ) std::cout << "Getting Branches Photon" << std::endl;
   b_Photon_covEtaEta->GetEntry(entry);   //!
   b_Photon_covEtaPhi->GetEntry(entry);   //!
   b_Photon_covPhiPhi->GetEntry(entry);   //!
   b_Photon_ecalRHSumEtConeDR04->GetEntry(entry);   //!
   b_Photon_energy->GetEntry(entry);   //!
   b_Photon_energyErr->GetEntry(entry);   //!
   b_Photon_energyRaw->GetEntry(entry);   //!
   b_Photon_eta->GetEntry(entry);   //!
   b_Photon_excluded->GetEntry(entry);   //!
   b_Photon_hcalTowerSumEtConeDR04->GetEntry(entry);   //!
   b_Photon_nTrkHollowConeDR04->GetEntry(entry);   //!
   b_Photon_nTrkSolidConeDR04->GetEntry(entry);   //!
   b_Photon_SigmaIEtaIEta->GetEntry(entry);   //!
   b_Photon_ecalPFClusterIso->GetEntry(entry);   //!
   b_hasConversionTracks->GetEntry(entry);   //!
   b_Photon_hcalPFClusterIso->GetEntry(entry);   //!
   b_Photon_Hoe_PUcorr->GetEntry(entry);   //!

   if( sbDEBUG ) std::cout << "Getting Branches Photon-GEN" << std::endl;
   if( doGenInfo ){
   b_Photon_genIdx->GetEntry(entry);   //!
   b_Photon_genLlpId->GetEntry(entry);   //!
   b_Photon_genSigMomEnergy->GetEntry(entry);   //!
   b_Photon_genSigMomEta->GetEntry(entry);   //!
   b_Photon_genSigMomMass->GetEntry(entry);   //!
   b_Photon_genSigMomPhi->GetEntry(entry);   //!
   b_Photon_genSigMomPt->GetEntry(entry);   //!
   b_Photon_genSigMomPx->GetEntry(entry);   //!
   b_Photon_genSigMomPy->GetEntry(entry);   //!
   b_Photon_genSigMomPz->GetEntry(entry);   //!
   b_Photon_genSigMomVx->GetEntry(entry);   //!
   b_Photon_genSigMomVy->GetEntry(entry);   //!
   b_Photon_genSigMomVz->GetEntry(entry);   //!

   }//if( doGenInfo )

   if( sbDEBUG ) std::cout << "Getting Branches Photon-shape" << std::endl;
   b_Photon_pfChargedIso->GetEntry(entry);   //!
   b_Photon_pfChargedIsoPFPV->GetEntry(entry);   //!
   b_Photon_pfChargedIsoWorstVtx->GetEntry(entry);   //!
   b_Photon_pfPhoIso03->GetEntry(entry);   //!
   b_pfRelIso03_all_quadratic->GetEntry(entry);   //!
   b_pfRelIso03_chg_quadratic->GetEntry(entry);   //!
   b_Photon_hadOverEM->GetEntry(entry);   //!
   b_Photon_hadTowOverEM->GetEntry(entry);   //!
   b_Photon_hcalTowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_Photon_isOot->GetEntry(entry);   //!
   b_Photon_phi->GetEntry(entry);   //!
   b_Photon_pt->GetEntry(entry);   //!
   b_Photon_px->GetEntry(entry);   //!
   b_Photon_py->GetEntry(entry);   //!
   b_Photon_pz->GetEntry(entry);   //!
   b_Photon_r9->GetEntry(entry);   //!
   b_Photon_rhIds->GetEntry(entry);   //!
   b_Photon_s4->GetEntry(entry);   //!
   b_Photon_salp->GetEntry(entry);   //!
   b_Photon_smaj->GetEntry(entry);   //!
   b_Photon_smin->GetEntry(entry);   //!
   b_Photon_seedTOFTime->GetEntry(entry);   //!
   b_Photon_trkSumPtHollowConeDR03->GetEntry(entry);   //!
   b_Photon_trkSumPtHollowConeDR04->GetEntry(entry);   //!
   b_Photon_trkSumPtSolidConeDR04->GetEntry(entry);   //!
   b_Photon_electronVeto->GetEntry(entry);   //!
   b_Photon_etaWidth->GetEntry(entry);   //!
   b_Photon_pixelSeed->GetEntry(entry);   //!
   b_Photon_seedIsEB->GetEntry(entry);   //!
   b_Photon_isScEtaEB->GetEntry(entry);   //!
   b_Photon_isScEtaEE->GetEntry(entry);   //!
   b_Photon_phiWidth->GetEntry(entry);   //!
   b_Photon_seediEtaOriX->GetEntry(entry);   //!
   b_Photon_seediPhiOriY->GetEntry(entry);   //!
   b_Photon_sieie->GetEntry(entry);   //!
   b_Photon_sieip->GetEntry(entry);   //!
   b_Photon_sipip->GetEntry(entry);   //!
   b_Photon_x_calo->GetEntry(entry);   //!
   b_Photon_y_calo->GetEntry(entry);   //!
   b_Photon_z_calo->GetEntry(entry);   //!

}//<<>>void root_base::getBranches(Long64_t entry)

/*
Bool_t root_base::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void root_base::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t root_base::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/

#endif // #ifdef root_base_cxx
