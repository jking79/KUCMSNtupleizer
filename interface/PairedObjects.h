#ifndef KUCMSNtupleizer_KUCMSNtupleizer_PairedObjects_h
#define KUCMSNtupleizer_KUCMSNtupleizer_PairedObjects_h

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedPair.h"

//===============================================================================================//
//                                      class: PairedObjects                                     //
//-----------------------------------------------------------------------------------------------//
// Container class for holding any two objects that were matched using either the DeltaRMatch or //
// the DeltaRMatchHungarian classes. This class stores the objects themselves, the indexes w.r.t //
// their original collections and the matched DR calculated from the original algorithm.         //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <class A, class B> class PairedObjects {

 public:

 PairedObjects() : objectA_(A()), objectB_(B()), matchedPair_(MatchedPair()), deltaR_(-999.) {}

 PairedObjects(const A &objectA, const B &objectB, const MatchedPair &matchedPair)
   : objectA_(objectA), objectB_(objectB), matchedPair_(matchedPair) {

    deltaR_ = sqrt(reco::deltaR2(objectA, objectB) );
  }

  A GetObjectA() const { return objectA_;}
  B GetObjectB() const { return objectB_;}
  int GetIndexA() const { return matchedPair_.GetIndexA();}
  int GetIndexB() const { return matchedPair_.GetIndexB();}
  double GetDeltaR() const { return deltaR_;}

  void CompareObjects() const {

    std::cout << "\nComparing objects with deltaR: " << deltaR_ << std::endl;
    std::cout << "Object A Summary: " << std::endl;
    std::cout << "  pT: " << objectA_.pt() << std::endl;
    std::cout << "  eta: " << objectA_.eta() << std::endl;
    std::cout << "  phi: " << objectA_.phi() << std::endl;

    std::cout << "\nObject B Summary: " << std::endl;
    std::cout << "  pT: " << objectB_.pt() << std::endl;
    std::cout << "  eta: " << objectB_.eta() << std::endl;
    std::cout << "  phi: " << objectB_.phi() << std::endl;
  }

 private:

  A objectA_;
  B objectB_;
  MatchedPair matchedPair_;
  double deltaR_;
};

#endif
