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

template <class A, class B> 
class PairedObjectCollection : public std::vector<PairedObjects<A,B> > {

public:

  int FindIndexA(const int indexB) const ;
  int FindIndexB(const A &objectA) const ;
  double FindDeltaR(const A &objectA) const;

  B FindObjectB(const A &objectA) const;
  
  std::vector<A> GetObjectAList() const;
  std::vector<B> GetObjectBList() const;
  
};

template <class A, class B>
int PairedObjectCollection<A,B>::FindIndexA(const int indexB) const {

  int index = -1;
  for(const auto &pair : *this) {
    if(pair.GetIndexB() == indexB)
      index = pair.GetIndexA();
  }
  
  if(index == -1) std::cout << "Warning: Index not found in list of objectB" << std::endl;
  
  return index;
}

template <class A, class B>
int PairedObjectCollection<A,B>::FindIndexB(const A &objectA) const {

  int index = -1;
  for(const auto &pair : *this)	{
    //std::cout << pair.GetIndexA() << std::endl; 
    if(pair.GetObjectA().pt() == objectA.pt() && pair.GetObjectA().eta() == objectA.eta()) {
      index = pair.GetIndexB();
      //std::cout << index << std::endl; 
    }
  }

  if(index == -1) {
    std::cout << "Warning: Index not found in list of objectA" << std::endl;
    
  }
   
  return index;
}

template <class A, class B>
double PairedObjectCollection<A,B>::FindDeltaR(const A &objectA) const {
  double deltaR = -999.;
  for(const auto &pair : *this)
    if(pair.GetObjectA().pt() == objectA.pt() && pair.GetObjectA().eta() == objectA.eta())
      deltaR = pair.GetDeltaR();

  return deltaR;
}

template <class A, class B>
B PairedObjectCollection<A,B>::FindObjectB(const A &objectA) const {

  B objectB;
  for(const auto &pair : *this) 
    if(pair.GetObjectA().pt() == objectA.pt() && pair.GetObjectA().eta() == objectA.eta()) 
      objectB = pair.GetObjectB();
    
  return objectB;
}

template <class A, class B>
std::vector<A> PairedObjectCollection<A,B>::GetObjectAList() const {

  std::vector<A> listA;
  for(const auto &pair : *this) 
    listA.emplace_back(pair.GetObjectA());

  return listA;
}

template <class A, class B>
std::vector<B> PairedObjectCollection<A,B>::GetObjectBList() const {

  std::vector<B> listB;
  for(const auto &pair : *this)	
    listB.emplace_back(pair.GetObjectB());

  return listB;
}

#endif
