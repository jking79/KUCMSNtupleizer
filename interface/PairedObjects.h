#ifndef KUCMSNtupleizer_KUCMSNtupleizer_PairedObjects_h
#define KUCMSNtupleizer_KUCMSNtupleizer_PairedObjects_h

#include <functional>
#include <optional>
#include <type_traits>

#include "DataFormats/Math/interface/deltaR.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchingTools.h"

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
   : objectA_(objectA), objectB_(objectB), matchedPair_(matchedPair), deltaR_(matchedPair.GetMatchCriteria()) {}

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
  B FindObjectB(const A &objectA) const;
  
  std::vector<A> GetObjectAList() const;
  std::vector<B> GetObjectBList() const;

  PairedObjectCollection<reco::Track, reco::GenParticle> ConvertFromTTracks() const {
    return this->Convert<reco::Track, reco::GenParticle>(
        [] (const reco::TransientTrack &ttrack) {
          return ttrack.track();
        },
	[] (const reco::GenParticle &gen) {
	  return gen;
	});
  }
  
private:
  
  template <class C = A, class D = B>
  PairedObjectCollection<C, D> Convert(
    std::optional<std::function<C(const A &)>> convertA = std::nullopt,
    std::optional<std::function<D(const B &)>> convertB = std::nullopt) const;
  
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

template <class A, class B>
template <class C, class D>
PairedObjectCollection<C, D> PairedObjectCollection<A, B>::Convert(
    std::optional<std::function<C(const A &)>> convertA,
    std::optional<std::function<D(const B &)>> convertB) const {

  if constexpr (std::is_same_v<A, C> && std::is_same_v<B, D>) {
    return *this; 
  }
  
  PairedObjectCollection<C, D> newCollection;

  for (const auto &pair : *this) {
    C newA = convertA ? (*convertA)(pair.GetObjectA()) : C();
    D newB = convertB ? (*convertB)(pair.GetObjectB()) : D();

    MatchedPair newMatchedPair(pair.GetIndexA(), pair.GetIndexB(), pair.GetDeltaR()); 

    newCollection.emplace_back(PairedObjects<C, D>(newA, newB, newMatchedPair));
  }

  return newCollection;
}

#endif
