//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 22 15:06:55 2024 by ROOT version 6.26/11
// from TTree llpgtree/KUCMSNtuple
// found on file: root/GMSB_L-150TeV_Ctau-0_1cm_Fall17.root
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
   Float_t         Beamspot_x;
   Float_t         Beamspot_y;
   Float_t         Beamspot_z;
   vector<int>     *GenParticle_charge;
   vector<float>   *GenParticle_deltaR;
   vector<float>   *GenParticle_dxy;
   vector<float>   *GenParticle_eta;
   vector<bool>    *GenParticle_isSignal;
   vector<int>     *GenParticle_matchedTrackIndex;
   UInt_t          GenParticle_nMatches;
   UInt_t          GenParticle_nSignal;
   vector<float>   *GenParticle_p;
   vector<int>     *GenParticle_pdgId;
   vector<float>   *GenParticle_phi;
   vector<float>   *GenParticle_pt;
   vector<float>   *GenParticle_x;
   vector<float>   *GenParticle_y;
   vector<float>   *GenParticle_z;
   vector<unsigned int> *PVCollection_genIndex;
   vector<unsigned int> *PVCollection_pvIndex;
   vector<float>   *PVCollection_signalWeight;
   vector<unsigned int> *PVCollection_trackIndex;
   vector<float>   *PVCollection_x;
   vector<float>   *PVCollection_y;
   vector<float>   *PVCollection_z;
   Float_t         PV_chi2;
   Float_t         PV_ecalness;
   Bool_t          PV_hasSignal;
   UInt_t          PV_nTracks;
   Float_t         PV_ndof;
   Float_t         PV_normalizedChi2;
   UInt_t          PV_signalCount;
   Float_t         PV_sumPt;
   vector<unsigned int> *PV_trackIndex;
   vector<float>   *PV_trackWeight;
   Float_t         PV_weightedSumPt;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   vector<float>   *SignalSV_chi2;
   vector<float>   *SignalSV_dxy;
   vector<unsigned int> *SignalSV_genIndex;
   Int_t           SignalSV_nTotal;
   vector<int>     *SignalSV_nTracks;
   vector<float>   *SignalSV_ndof;
   vector<float>   *SignalSV_normalizedChi2;
   vector<unsigned int> *SignalSV_trackIndex;
   vector<float>   *SignalSV_trackWeight;
   vector<unsigned int> *SignalSV_vertexIndex;
   vector<float>   *SignalSV_x;
   vector<float>   *SignalSV_y;
   vector<float>   *SignalSV_z;
   vector<int>     *Track_charge;
   vector<float>   *Track_chi2;
   vector<float>   *Track_dxy;
   vector<int>     *Track_ecalTrackIndex;
   vector<float>   *Track_eta;
   vector<float>   *Track_etaError;
   vector<float>   *Track_ip2D;
   vector<float>   *Track_ip3D;
   vector<bool>    *Track_isECAL;
   UInt_t          Track_nECALTracks;
   vector<int>     *Track_nLostHits;
   UInt_t          Track_nTotal;
   vector<int>     *Track_nValidHits;
   vector<float>   *Track_ndof;
   vector<float>   *Track_normalizedChi2;
   vector<float>   *Track_p;
   vector<float>   *Track_phi;
   vector<float>   *Track_phiError;
   vector<float>   *Track_pt;
   vector<float>   *Track_ptError;
   vector<float>   *Track_px;
   vector<float>   *Track_py;
   vector<float>   *Track_pz;
   vector<float>   *Track_qoverp;
   vector<int>     *Track_qualityMask;
   vector<float>   *Track_sip2D;
   vector<float>   *Track_sip3D;
   vector<float>   *Track_vx;
   vector<float>   *Track_vy;
   vector<float>   *Track_vz;
  vector<float>   *Track_SCDR;
  vector<float>   *Vertex_chi2;
   vector<float>   *Vertex_dxy;
   vector<float>   *Vertex_ecalness;
   vector<bool>    *Vertex_hasSignal;
   UInt_t          Vertex_nTotal;
   vector<unsigned int> *Vertex_nTracks;
   vector<float>   *Vertex_ndof;
   vector<float>   *Vertex_normalizedChi2;
   vector<unsigned int> *Vertex_signalCount;
   vector<float>   *Vertex_sumPt;
   vector<unsigned int> *Vertex_trackIndex;
   vector<float>   *Vertex_trackRatioWithPV;
   vector<float>   *Vertex_trackWeight;
   vector<unsigned int> *Vertex_vertexIndex;
   vector<float>   *Vertex_x;
   vector<float>   *Vertex_y;
   vector<float>   *Vertex_z;

   // List of branches
   TBranch        *b_Beamspot_x;   //!
   TBranch        *b_Beamspot_y;   //!
   TBranch        *b_Beamspot_z;   //!
   TBranch        *b_GenParticle_charge;   //!
   TBranch        *b_GenParticle_deltaR;   //!
   TBranch        *b_GenParticle_dxy;   //!
   TBranch        *b_GenParticle_eta;   //!
   TBranch        *b_GenParticle_isSignal;   //!
   TBranch        *b_GenParticle_matchedTrackIndex;   //!
   TBranch        *b_GenParticle_nMatches;   //!
   TBranch        *b_GenParticle_nSignal;   //!
   TBranch        *b_GenParticle_p;   //!
   TBranch        *b_GenParticle_pdgId;   //!
   TBranch        *b_GenParticle_phi;   //!
   TBranch        *b_GenParticle_pt;   //!
   TBranch        *b_GenParticle_x;   //!
   TBranch        *b_GenParticle_y;   //!
   TBranch        *b_GenParticle_z;   //!
   TBranch        *b_PVCollection_genIndex;   //!
   TBranch        *b_PVCollection_pvIndex;   //!
   TBranch        *b_PVCollection_signalWeight;   //!
   TBranch        *b_PVCollection_trackIndex;   //!
   TBranch        *b_PVCollection_x;   //!
   TBranch        *b_PVCollection_y;   //!
   TBranch        *b_PVCollection_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_ecalness;   //!
   TBranch        *b_PV_hasSignal;   //!
   TBranch        *b_PV_nTracks;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_normalizedChi2;   //!
   TBranch        *b_PV_signalCount;   //!
   TBranch        *b_PV_sumPt;   //!
   TBranch        *b_PV_trackIndex;   //!
   TBranch        *b_PV_trackWeight;   //!
   TBranch        *b_PV_weightedSumPt;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_SignalSV_chi2;   //!
   TBranch        *b_SignalSV_dxy;   //!
   TBranch        *b_SignalSV_genIndex;   //!
   TBranch        *b_SignalSV_nTotal;   //!
   TBranch        *b_SignalSV_nTracks;   //!
   TBranch        *b_SignalSV_ndof;   //!
   TBranch        *b_SignalSV_normalizedChi2;   //!
   TBranch        *b_SignalSV_trackIndex;   //!
   TBranch        *b_SignalSV_trackWeight;   //!
   TBranch        *b_SignalSV_vertexIndex;   //!
   TBranch        *b_SignalSV_x;   //!
   TBranch        *b_SignalSV_y;   //!
   TBranch        *b_SignalSV_z;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_Track_dxy;   //!
   TBranch        *b_Track_ecalTrackIndex;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_etaError;   //!
   TBranch        *b_Track_ip2D;   //!
   TBranch        *b_Track_ip3D;   //!
   TBranch        *b_Track_isECAL;   //!
   TBranch        *b_Track_nECALTracks;   //!
   TBranch        *b_Track_nLostHits;   //!
   TBranch        *b_Track_nTotal;   //!
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
   TBranch        *b_Track_qoverp;   //!
   TBranch        *b_Track_qualityMask;   //!
   TBranch        *b_Track_sip2D;   //!
   TBranch        *b_Track_sip3D;   //!
   TBranch        *b_Track_vx;   //!
   TBranch        *b_Track_vy;   //!
   TBranch        *b_Track_vz;   //!
   TBranch        *b_Track_SCDR;   //!  
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_dxy;   //!
   TBranch        *b_Vertex_ecalness;   //!
   TBranch        *b_Vertex_hasSignal;   //!
   TBranch        *b_Vertex_nTotal;   //!
   TBranch        *b_Vertex_nTracks;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_normalizedChi2;   //!
   TBranch        *b_Vertex_signalCount;   //!
   TBranch        *b_Vertex_sumPt;   //!
   TBranch        *b_Vertex_trackIndex;   //!
   TBranch        *b_Vertex_trackRatioWithPV;   //!
   TBranch        *b_Vertex_trackWeight;   //!
   TBranch        *b_Vertex_vertexIndex;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!

   llpgtree(TTree *tree=0);
   virtual ~llpgtree();
   virtual void     Init(TTree *tree);
};

llpgtree::llpgtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root/GMSB_L-150TeV_Ctau-0_1cm_Fall17.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root/GMSB_L-150TeV_Ctau-0_1cm_Fall17.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root/GMSB_L-150TeV_Ctau-0_1cm_Fall17.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

llpgtree::~llpgtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


void llpgtree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   GenParticle_charge = 0;
   GenParticle_deltaR = 0;
   GenParticle_dxy = 0;
   GenParticle_eta = 0;
   GenParticle_isSignal = 0;
   GenParticle_matchedTrackIndex = 0;
   GenParticle_p = 0;
   GenParticle_pdgId = 0;
   GenParticle_phi = 0;
   GenParticle_pt = 0;
   GenParticle_x = 0;
   GenParticle_y = 0;
   GenParticle_z = 0;
   PVCollection_genIndex = 0;
   PVCollection_pvIndex = 0;
   PVCollection_signalWeight = 0;
   PVCollection_trackIndex = 0;
   PVCollection_x = 0;
   PVCollection_y = 0;
   PVCollection_z = 0;
   PV_trackIndex = 0;
   PV_trackWeight = 0;
   SignalSV_chi2 = 0;
   SignalSV_dxy = 0;
   SignalSV_genIndex = 0;
   SignalSV_nTracks = 0;
   SignalSV_ndof = 0;
   SignalSV_normalizedChi2 = 0;
   SignalSV_trackIndex = 0;
   SignalSV_trackWeight = 0;
   SignalSV_vertexIndex = 0;
   SignalSV_x = 0;
   SignalSV_y = 0;
   SignalSV_z = 0;
   Track_charge = 0;
   Track_chi2 = 0;
   Track_dxy = 0;
   Track_ecalTrackIndex = 0;
   Track_eta = 0;
   Track_etaError = 0;
   Track_ip2D = 0;
   Track_ip3D = 0;
   Track_isECAL = 0;
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
   Track_qoverp = 0;
   Track_qualityMask = 0;
   Track_sip2D = 0;
   Track_sip3D = 0;
   Track_vx = 0;
   Track_vy = 0;
   Track_vz = 0;
   Track_SCDR = 0;
   Vertex_chi2 = 0;
   Vertex_dxy = 0;
   Vertex_ecalness = 0;
   Vertex_hasSignal = 0;
   Vertex_nTracks = 0;
   Vertex_ndof = 0;
   Vertex_normalizedChi2 = 0;
   Vertex_signalCount = 0;
   Vertex_sumPt = 0;
   Vertex_trackIndex = 0;
   Vertex_trackRatioWithPV = 0;
   Vertex_trackWeight = 0;
   Vertex_vertexIndex = 0;
   Vertex_x = 0;
   Vertex_y = 0;
   Vertex_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Beamspot_x", &Beamspot_x, &b_Beamspot_x);
   fChain->SetBranchAddress("Beamspot_y", &Beamspot_y, &b_Beamspot_y);
   fChain->SetBranchAddress("Beamspot_z", &Beamspot_z, &b_Beamspot_z);
   fChain->SetBranchAddress("GenParticle_charge", &GenParticle_charge, &b_GenParticle_charge);
   fChain->SetBranchAddress("GenParticle_deltaR", &GenParticle_deltaR, &b_GenParticle_deltaR);
   fChain->SetBranchAddress("GenParticle_dxy", &GenParticle_dxy, &b_GenParticle_dxy);
   fChain->SetBranchAddress("GenParticle_eta", &GenParticle_eta, &b_GenParticle_eta);
   fChain->SetBranchAddress("GenParticle_isSignal", &GenParticle_isSignal, &b_GenParticle_isSignal);
   fChain->SetBranchAddress("GenParticle_matchedTrackIndex", &GenParticle_matchedTrackIndex, &b_GenParticle_matchedTrackIndex);
   fChain->SetBranchAddress("GenParticle_nMatches", &GenParticle_nMatches, &b_GenParticle_nMatches);
   fChain->SetBranchAddress("GenParticle_nSignal", &GenParticle_nSignal, &b_GenParticle_nSignal);
   fChain->SetBranchAddress("GenParticle_p", &GenParticle_p, &b_GenParticle_p);
   fChain->SetBranchAddress("GenParticle_pdgId", &GenParticle_pdgId, &b_GenParticle_pdgId);
   fChain->SetBranchAddress("GenParticle_phi", &GenParticle_phi, &b_GenParticle_phi);
   fChain->SetBranchAddress("GenParticle_pt", &GenParticle_pt, &b_GenParticle_pt);
   fChain->SetBranchAddress("GenParticle_x", &GenParticle_x, &b_GenParticle_x);
   fChain->SetBranchAddress("GenParticle_y", &GenParticle_y, &b_GenParticle_y);
   fChain->SetBranchAddress("GenParticle_z", &GenParticle_z, &b_GenParticle_z);
   fChain->SetBranchAddress("PVCollection_genIndex", &PVCollection_genIndex, &b_PVCollection_genIndex);
   fChain->SetBranchAddress("PVCollection_pvIndex", &PVCollection_pvIndex, &b_PVCollection_pvIndex);
   fChain->SetBranchAddress("PVCollection_signalWeight", &PVCollection_signalWeight, &b_PVCollection_signalWeight);
   fChain->SetBranchAddress("PVCollection_trackIndex", &PVCollection_trackIndex, &b_PVCollection_trackIndex);
   fChain->SetBranchAddress("PVCollection_x", &PVCollection_x, &b_PVCollection_x);
   fChain->SetBranchAddress("PVCollection_y", &PVCollection_y, &b_PVCollection_y);
   fChain->SetBranchAddress("PVCollection_z", &PVCollection_z, &b_PVCollection_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_ecalness", &PV_ecalness, &b_PV_ecalness);
   fChain->SetBranchAddress("PV_hasSignal", &PV_hasSignal, &b_PV_hasSignal);
   fChain->SetBranchAddress("PV_nTracks", &PV_nTracks, &b_PV_nTracks);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_normalizedChi2", &PV_normalizedChi2, &b_PV_normalizedChi2);
   fChain->SetBranchAddress("PV_signalCount", &PV_signalCount, &b_PV_signalCount);
   fChain->SetBranchAddress("PV_sumPt", &PV_sumPt, &b_PV_sumPt);
   fChain->SetBranchAddress("PV_trackIndex", &PV_trackIndex, &b_PV_trackIndex);
   fChain->SetBranchAddress("PV_trackWeight", &PV_trackWeight, &b_PV_trackWeight);
   fChain->SetBranchAddress("PV_weightedSumPt", &PV_weightedSumPt, &b_PV_weightedSumPt);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("SignalSV_chi2", &SignalSV_chi2, &b_SignalSV_chi2);
   fChain->SetBranchAddress("SignalSV_dxy", &SignalSV_dxy, &b_SignalSV_dxy);
   fChain->SetBranchAddress("SignalSV_genIndex", &SignalSV_genIndex, &b_SignalSV_genIndex);
   fChain->SetBranchAddress("SignalSV_nTotal", &SignalSV_nTotal, &b_SignalSV_nTotal);
   fChain->SetBranchAddress("SignalSV_nTracks", &SignalSV_nTracks, &b_SignalSV_nTracks);
   fChain->SetBranchAddress("SignalSV_ndof", &SignalSV_ndof, &b_SignalSV_ndof);
   fChain->SetBranchAddress("SignalSV_normalizedChi2", &SignalSV_normalizedChi2, &b_SignalSV_normalizedChi2);
   fChain->SetBranchAddress("SignalSV_trackIndex", &SignalSV_trackIndex, &b_SignalSV_trackIndex);
   fChain->SetBranchAddress("SignalSV_trackWeight", &SignalSV_trackWeight, &b_SignalSV_trackWeight);
   fChain->SetBranchAddress("SignalSV_vertexIndex", &SignalSV_vertexIndex, &b_SignalSV_vertexIndex);
   fChain->SetBranchAddress("SignalSV_x", &SignalSV_x, &b_SignalSV_x);
   fChain->SetBranchAddress("SignalSV_y", &SignalSV_y, &b_SignalSV_y);
   fChain->SetBranchAddress("SignalSV_z", &SignalSV_z, &b_SignalSV_z);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("Track_dxy", &Track_dxy, &b_Track_dxy);
   fChain->SetBranchAddress("Track_ecalTrackIndex", &Track_ecalTrackIndex, &b_Track_ecalTrackIndex);
   fChain->SetBranchAddress("Track_eta", &Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_etaError", &Track_etaError, &b_Track_etaError);
   fChain->SetBranchAddress("Track_ip2D", &Track_ip2D, &b_Track_ip2D);
   fChain->SetBranchAddress("Track_ip3D", &Track_ip3D, &b_Track_ip3D);
   fChain->SetBranchAddress("Track_isECAL", &Track_isECAL, &b_Track_isECAL);
   fChain->SetBranchAddress("Track_nECALTracks", &Track_nECALTracks, &b_Track_nECALTracks);
   fChain->SetBranchAddress("Track_nLostHits", &Track_nLostHits, &b_Track_nLostHits);
   fChain->SetBranchAddress("Track_nTotal", &Track_nTotal, &b_Track_nTotal);
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
   fChain->SetBranchAddress("Track_qoverp", &Track_qoverp, &b_Track_qoverp);
   fChain->SetBranchAddress("Track_qualityMask", &Track_qualityMask, &b_Track_qualityMask);
   fChain->SetBranchAddress("Track_sip2D", &Track_sip2D, &b_Track_sip2D);
   fChain->SetBranchAddress("Track_sip3D", &Track_sip3D, &b_Track_sip3D);
   fChain->SetBranchAddress("Track_vx", &Track_vx, &b_Track_vx);
   fChain->SetBranchAddress("Track_vy", &Track_vy, &b_Track_vy);
   fChain->SetBranchAddress("Track_vz", &Track_vz, &b_Track_vz);
   fChain->SetBranchAddress("Track_SCDR", &Track_SCDR, &b_Track_SCDR);
   fChain->SetBranchAddress("Vertex_chi2", &Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_dxy", &Vertex_dxy, &b_Vertex_dxy);
   fChain->SetBranchAddress("Vertex_ecalness", &Vertex_ecalness, &b_Vertex_ecalness);
   fChain->SetBranchAddress("Vertex_hasSignal", &Vertex_hasSignal, &b_Vertex_hasSignal);
   fChain->SetBranchAddress("Vertex_nTotal", &Vertex_nTotal, &b_Vertex_nTotal);
   fChain->SetBranchAddress("Vertex_nTracks", &Vertex_nTracks, &b_Vertex_nTracks);
   fChain->SetBranchAddress("Vertex_ndof", &Vertex_ndof, &b_Vertex_ndof);
   fChain->SetBranchAddress("Vertex_normalizedChi2", &Vertex_normalizedChi2, &b_Vertex_normalizedChi2);
   fChain->SetBranchAddress("Vertex_signalCount", &Vertex_signalCount, &b_Vertex_signalCount);
   fChain->SetBranchAddress("Vertex_sumPt", &Vertex_sumPt, &b_Vertex_sumPt);
   fChain->SetBranchAddress("Vertex_trackIndex", &Vertex_trackIndex, &b_Vertex_trackIndex);
   fChain->SetBranchAddress("Vertex_trackRatioWithPV", &Vertex_trackRatioWithPV, &b_Vertex_trackRatioWithPV);
   fChain->SetBranchAddress("Vertex_trackWeight", &Vertex_trackWeight, &b_Vertex_trackWeight);
   fChain->SetBranchAddress("Vertex_vertexIndex", &Vertex_vertexIndex, &b_Vertex_vertexIndex);
   fChain->SetBranchAddress("Vertex_x", &Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", &Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", &Vertex_z, &b_Vertex_z);
}

#endif // #ifdef llpgtree_cxx
