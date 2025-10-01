//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 28 14:33:01 2024 by ROOT version 6.26/11
// Base class for KUCMSNtuple Trees
// found on file: gmsb_AODSIM_KUCMSNtuplizer_justin_test_fullgen_oottrue.root
//////////////////////////////////////////////////////////

#ifndef KUCMSNtuple_Tree_Base_h
#define KUCMSNtuple_Tree_Base_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class KUCMSNtuple_Tree_Base {

public :

   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   llpgtree( TChain *tree=0 );
   virtual void Init( TChain *tree, bool doGenInfo );
   virtual void getBranches( Long64_t entry, bool doGenInfo );
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);

};//<<>>class llpgtree

#endif // #ifdef KUCMSNtuple_Tree_Base_cxx
