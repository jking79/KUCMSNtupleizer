//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 13 18:45:20 2024 by ROOT version 6.26/11
// from TTree llpgtree/KUCMSNtuple
// found on file: gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root
//////////////////////////////////////////////////////////

#ifndef llpgtree_h
#define llpgtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class llpgtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *TrackMerger_deltaR;
   vector<int>     *TrackMerger_displacedElectronIndex;
   UInt_t          TrackMerger_nMatches;
   vector<float>   *TrackMerger_ptDiff;
   vector<float>   *TrackMerger_vzDiff;
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
   vector<int>     *Gen_susId;
   vector<float>   *Gen_vx;
   vector<float>   *Gen_vy;
   vector<float>   *Gen_vz;
   Float_t         Evt_genWgt;
   vector<float>   *TimedSV_chi2;
   vector<int>     *TimedSV_displacedElectronIndex;
   vector<float>   *TimedSV_dxy;
   vector<float>   *TimedSV_linearDistance;
   Int_t           TimedSV_nTotal;
   vector<int>     *TimedSV_nTracks;
   vector<float>   *TimedSV_ndof;
   vector<float>   *TimedSV_normalizedChi2;
   vector<float>   *TimedSV_pathLength;
   vector<float>   *TimedSV_trackPt;
   vector<unsigned int> *TimedSV_vertexIndex;
   vector<float>   *TimedSV_x;
   vector<float>   *TimedSV_y;
   vector<float>   *TimedSV_z;

   // List of branches
   TBranch        *b_TrackMerger_deltaR;   //!
   TBranch        *b_TrackMerger_displacedElectronIndex;   //!
   TBranch        *b_TrackMerger_nMatches;   //!
   TBranch        *b_TrackMerger_ptDiff;   //!
   TBranch        *b_TrackMerger_vzDiff;   //!
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
   TBranch        *b_Gen_susId;   //!
   TBranch        *b_Gen_vx;   //!
   TBranch        *b_Gen_vy;   //!
   TBranch        *b_Gen_vz;   //!
   TBranch        *b_Evt_genWgt;   //!
   TBranch        *b_TimedSV_chi2;   //!
   TBranch        *b_TimedSV_displacedElectronIndex;   //!
   TBranch        *b_TimedSV_dxy;   //!
   TBranch        *b_TimedSV_linearDistance;   //!
   TBranch        *b_TimedSV_nTotal;   //!
   TBranch        *b_TimedSV_nTracks;   //!
   TBranch        *b_TimedSV_ndof;   //!
   TBranch        *b_TimedSV_normalizedChi2;   //!
   TBranch        *b_TimedSV_pathLength;   //!
   TBranch        *b_TimedSV_trackPt;   //!
   TBranch        *b_TimedSV_vertexIndex;   //!
   TBranch        *b_TimedSV_x;   //!
   TBranch        *b_TimedSV_y;   //!
   TBranch        *b_TimedSV_z;   //!

   llpgtree(TTree *tree=0);
   virtual ~llpgtree();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
};

llpgtree::llpgtree(TTree *tree) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

llpgtree::~llpgtree() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t llpgtree::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void llpgtree::Init(TTree *tree) {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TrackMerger_deltaR = 0;
   TrackMerger_displacedElectronIndex = 0;
   TrackMerger_ptDiff = 0;
   TrackMerger_vzDiff = 0;
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
   TimedSV_chi2 = 0;
   TimedSV_displacedElectronIndex = 0;
   TimedSV_dxy = 0;
   TimedSV_linearDistance = 0;
   TimedSV_nTracks = 0;
   TimedSV_ndof = 0;
   TimedSV_normalizedChi2 = 0;
   TimedSV_pathLength = 0;
   TimedSV_trackPt = 0;
   TimedSV_vertexIndex = 0;
   TimedSV_x = 0;
   TimedSV_y = 0;
   TimedSV_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("TrackMerger_deltaR", &TrackMerger_deltaR, &b_TrackMerger_deltaR);
   fChain->SetBranchAddress("TrackMerger_displacedElectronIndex", &TrackMerger_displacedElectronIndex, &b_TrackMerger_displacedElectronIndex);
   fChain->SetBranchAddress("TrackMerger_nMatches", &TrackMerger_nMatches, &b_TrackMerger_nMatches);
   fChain->SetBranchAddress("TrackMerger_ptDiff", &TrackMerger_ptDiff, &b_TrackMerger_ptDiff);
   fChain->SetBranchAddress("TrackMerger_vzDiff", &TrackMerger_vzDiff, &b_TrackMerger_vzDiff);
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
   fChain->SetBranchAddress("Gen_susId", &Gen_susId, &b_Gen_susId);
   fChain->SetBranchAddress("Gen_vx", &Gen_vx, &b_Gen_vx);
   fChain->SetBranchAddress("Gen_vy", &Gen_vy, &b_Gen_vy);
   fChain->SetBranchAddress("Gen_vz", &Gen_vz, &b_Gen_vz);
   fChain->SetBranchAddress("Evt_genWgt", &Evt_genWgt, &b_Evt_genWgt);
   fChain->SetBranchAddress("TimedSV_chi2", &TimedSV_chi2, &b_TimedSV_chi2);
   fChain->SetBranchAddress("TimedSV_displacedElectronIndex", &TimedSV_displacedElectronIndex, &b_TimedSV_displacedElectronIndex);
   fChain->SetBranchAddress("TimedSV_dxy", &TimedSV_dxy, &b_TimedSV_dxy);
   fChain->SetBranchAddress("TimedSV_linearDistance", &TimedSV_linearDistance, &b_TimedSV_linearDistance);
   fChain->SetBranchAddress("TimedSV_nTotal", &TimedSV_nTotal, &b_TimedSV_nTotal);
   fChain->SetBranchAddress("TimedSV_nTracks", &TimedSV_nTracks, &b_TimedSV_nTracks);
   fChain->SetBranchAddress("TimedSV_ndof", &TimedSV_ndof, &b_TimedSV_ndof);
   fChain->SetBranchAddress("TimedSV_normalizedChi2", &TimedSV_normalizedChi2, &b_TimedSV_normalizedChi2);
   fChain->SetBranchAddress("TimedSV_pathLength", &TimedSV_pathLength, &b_TimedSV_pathLength);
   fChain->SetBranchAddress("TimedSV_trackPt", &TimedSV_trackPt, &b_TimedSV_trackPt);
   fChain->SetBranchAddress("TimedSV_vertexIndex", &TimedSV_vertexIndex, &b_TimedSV_vertexIndex);
   fChain->SetBranchAddress("TimedSV_x", &TimedSV_x, &b_TimedSV_x);
   fChain->SetBranchAddress("TimedSV_y", &TimedSV_y, &b_TimedSV_y);
   fChain->SetBranchAddress("TimedSV_z", &TimedSV_z, &b_TimedSV_z);
}
#endif // #ifdef llpgtree_cxx
