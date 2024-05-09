#ifndef KUCMSNtupleizer_KUCMSNtupleizer_DeltaRMatch_h
#define KUCMSNtupleizer_KUCMSNtupleizer_DeltaRMatch_h

#include <algorithm>

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/Hungarian.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchingTools.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/PairedObjects.h"

template <class A, class B> class DeltaRMatch {

 public:

  DeltaRMatch(const std::vector<A> &objectsA, const std::vector<B> &objectsB, const double minDeltaR, const double minPt = 2.) {
    FindMinDeltaRMatches(objectsA, objectsB, minDeltaR, minPt);
  }

  virtual ~DeltaRMatch() = default;

  int GetNMatches() const {return matchedIndexesA_.size();}
  std::vector<int> GetMatchedIndexesA() const {return matchedIndexesA_;}
  std::vector<int> GetMatchedIndexesB() const {return matchedIndexesB_;}
  std::vector<double> GetMatchedDeltaRs() const {return matchedDeltaRs_;}

 private:
  double deltaR_;
  std::vector<int> matchedIndexesA_;
  std::vector<int> matchedIndexesB_;
  std::vector<double> matchedDeltaRs_;

  bool IsIndexMatched(std::vector<int> matchedIndexes, int index) const {
    return std::find(matchedIndexes.begin(), matchedIndexes.end(), index) != matchedIndexes.end();
  }

  int FindIndex(const std::vector<int>& vec, const int& value) const {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
      return std::distance(vec.begin(), it);
    } else {
      return -1;
    }
  }

  void FindMinDeltaRMatches(const std::vector<A> &objectsA, const std::vector<B> &objectsB, const double minDeltaR, const double minPt) {

    size_t dimensionA = objectsA.size();
    size_t dimensionB = objectsB.size();

    for(size_t ai = 0; ai < dimensionA; ai++) {

      bool isMatched(false);
      int minIndexB(-1);
      float deltaR(-1.), minMatchedDeltaR(minDeltaR);

      const A objectA = objectsA[ai];
      const float etaA = objectA.eta();
      const float phiA = objectA.phi();
      const float ptA = objectA.pt();

      if(ptA < minPt)
	continue;

      for(size_t bi = 0; bi< dimensionB; bi++) {

	const B objectB = objectsB[bi];
	const float etaB = objectB.eta();
	const float phiB = objectB.phi();

	deltaR = sqrt(reco::deltaR2(etaA, phiA, etaB, phiB));

	if(IsIndexMatched(matchedIndexesB_, bi))
	  continue;

	if(deltaR < minMatchedDeltaR) {
	  isMatched = true;

	  minIndexB = bi;
	  minMatchedDeltaR = deltaR;
	}
      }// End of B loop

      if(!isMatched)
	continue;

      const int matchListIndex(FindIndex(matchedIndexesB_, minIndexB));
      if(matchListIndex > 0) {
	if(matchedDeltaRs_[matchListIndex] < minMatchedDeltaR) {
	  matchedIndexesA_[matchListIndex] = ai;
	  continue;
	}
      }

      matchedIndexesA_.push_back(ai);
      matchedIndexesB_.push_back(minIndexB);
      matchedDeltaRs_.push_back(minMatchedDeltaR);

    }// End of A loop
  }

};

template <typename S> using Matrix = std::vector<std::vector<S> >;

//===============================================================================================//
//                               class: DeltaRMatchHungarian                                     //
//-----------------------------------------------------------------------------------------------//
// General delta R matching algorithm between any two lists of objects, as long as they possess  //
// eta() and phi() public methods (this is the only requirement!). Matching is performed using   //
// the Hungarian algorithm. Matches are always guaranteed using this method, which means the net //
// number of matches will always equal the smaller of the two input object vectors given (even   //
// when some of the matches are "bad").                                                          //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <class A, class B> class DeltaRMatchHungarian {
  
 public:

  DeltaRMatchHungarian();

  DeltaRMatchHungarian(Matrix<double> &costMatrix);
    
  DeltaRMatchHungarian(const std::vector<A> &objectsA, const std::vector<B> &objectsB);

  virtual ~DeltaRMatchHungarian() = default;

  std::vector<int> GetIndexListOfMatches() const 
    {return matchedIndexes_;}

  double GetCost() const {return cost_;}
  std::vector<A> GetMatchedObjectsA() const {matchedObjectsA_;}
  std::vector<B> GetMatchedObjectsB() const {matchedObjectsB_;}  
  std::vector<MatchedPair> GetMatchedPairs() const {return matchedPairs_;}
  std::vector<PairedObjects<A,B> > GetPairedObjects() const {return pairedObjects_;}
  std::vector<double> GetDeltaRs() const {
    std::vector<double> deltaRs;
    for(auto const &pair : pairedObjects_)
      deltaRs.push_back(pair.GetDeltaR());
    return deltaRs;
  }  

 protected:
  
  double cost_;
  std::vector<int> matchedIndexes_;
  std::vector<A> matchedObjectsA_;
  std::vector<B> matchedObjectsB_;
  std::vector<MatchedPair> matchedPairs_;
  std::vector<PairedObjects<A,B> > pairedObjects_;

  std::vector<MatchedPair> ConstructMatchedPairs(const Matrix<double> &costMatrix, 
						 const std::vector<A> &objectsA, 
						 const std::vector<B> &objectsB);

  Matrix<double> CalculateCostMatrix(const std::vector<A> &objectsA, const std::vector<B> &objectsB) const;

  void Solve(const std::vector<A> &objectsA, const std::vector<B> &objectsB);

  void PrintCostMatrix(Matrix<double> &costMatrix) const;
};

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                         Default Constructor for DeltaRMatchHungarian                          //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <typename A, typename B>
  DeltaRMatchHungarian<A,B>::DeltaRMatchHungarian() 
  : cost_(-999.),
  matchedIndexes_(std::vector<int>()),
  matchedPairs_(std::vector<MatchedPair>()),
  pairedObjects_(std::vector<PairedObjects<A,B>>()) {}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                            Templated Private Method: Solve<A, B>                              //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================// 
template <typename A, typename B>
  void DeltaRMatchHungarian<A,B>::Solve(const std::vector<A> &objectsA, const std::vector<B> &objectsB) {
  
  Matrix<double> costMatrix = CalculateCostMatrix(objectsA, objectsB);
  HungarianAlgorithm assigner;
  cost_ = assigner.Solve(costMatrix, matchedIndexes_);
  matchedPairs_ = ConstructMatchedPairs(costMatrix, objectsA, objectsB);
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                        Standard Constructor for DeltaRMatchHungarian                          //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <typename A, typename B>
  DeltaRMatchHungarian<A,B>::DeltaRMatchHungarian(const std::vector<A> &objectsA, const std::vector<B> &objectsB) {

  // If either collection is empty, assignment is not possible (obviously)
  // Initialize members to default values and empty vectors
  if(objectsA.size() == 0 || objectsB.size() == 0)// || typeid(B) == typeid(reco::GenParticle)) {
    *this = DeltaRMatchHungarian();
  
  else this->Solve(objectsA, objectsB);

}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//           Alternate constructor for DeltaRMatchHungarian (user has own cost matrix)           //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <typename A, typename B>
  DeltaRMatchHungarian<A,B>::DeltaRMatchHungarian(Matrix<double> &costMatrix) {

  HungarianAlgorithm assigner;
  cost_ = assigner.Solve(costMatrix, matchedIndexes_);
  matchedPairs_ = ConstructMatchedPairs(costMatrix);
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                          Templated Private Method: ConstructMatchedPairs                      //
//-----------------------------------------------------------------------------------------------//
// First step in the book-keeping of paired objects. Saves indexes of each match w.r.t the input //
// collections and the DR of the resulting match.                                                //
//===============================================================================================//
template <typename A, typename B>
  std::vector<MatchedPair> DeltaRMatchHungarian<A,B>::ConstructMatchedPairs(const std::vector<std::vector<double> > &costMatrix,
									    const std::vector<A> &objectsA,
									    const std::vector<B> &objectsB) {

  matchedObjectsA_.clear();
  matchedObjectsB_.clear();

  std::vector<MatchedPair> matchedPairs;
  for(unsigned int mi = 0; mi < costMatrix.size(); mi++) {
    if(matchedIndexes_[mi] == -1)
      continue;

    const double deltaR = costMatrix[mi][matchedIndexes_[mi]];
    const MatchedPair thisMatchedPair(mi, matchedIndexes_[mi], deltaR);

    matchedObjectsA_.emplace_back(objectsA[mi]);
    matchedObjectsB_.emplace_back(objectsB[matchedIndexes_[mi]]);
    matchedPairs.emplace_back(MatchedPair(mi, matchedIndexes_[mi], deltaR));
    pairedObjects_.emplace_back(PairedObjects<A,B>(objectsA[mi], objectsB[matchedIndexes_[mi]], thisMatchedPair) );
  }
  // Sort paired objects from smallest to largest DeltaR
  std::sort(pairedObjects_.begin(), pairedObjects_.end(), [](const PairedObjects<A,B>& a, const PairedObjects<A,B>& b) {
      return a.GetDeltaR() < b.GetDeltaR();
    });
  
  return matchedPairs;
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                         Templated Private Method: CalculateCostMatrix                         //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <typename A, typename B>
  Matrix<double> DeltaRMatchHungarian<A,B>::CalculateCostMatrix(const std::vector<A> &objectsA, const std::vector<B> &objectsB) const {

  const int dimensionA = objectsA.size();
  const int dimensionB = objectsB.size();
  Matrix<double> costMatrix(dimensionA, std::vector<double>(dimensionB) );

  for(int ai = 0; ai < dimensionA; ai++)
    for(int bi = 0; bi < dimensionB; bi++)
      costMatrix[ai][bi] = sqrt(reco::deltaR2(objectsA[ai], objectsB[bi]) );

  return costMatrix;
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                           Templated Private Method: PrintCostMatrix                           //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <typename A, typename B>
  void DeltaRMatchHungarian<A,B>::PrintCostMatrix(Matrix<double> &costMatrix) const {

  std::cout << "Cost Matrix: " << std::endl;
  int indexRow = 0;
  for(const auto &row : costMatrix) {
    std::cout << "track " << indexRow << std::endl;
    for(const auto &element : row) {
      if(element < 0.01)
	std::cout << element << ", ";
    }
    indexRow++;
    std::cout << std::endl;
  }
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                              class: DeltaRGenMatchHungarian                                   //
//-----------------------------------------------------------------------------------------------//
// Specialized version of DeltaRMatchHungarian for gen matching. Makes sure that only status one //
// charged gen particles are used, while preserving the original indexing of the gen collection. //  
//===============================================================================================//
template <class A>
class DeltaRGenMatchHungarian : public DeltaRMatchHungarian<A, reco::GenParticle> {

 public:

 DeltaRGenMatchHungarian() : DeltaRMatchHungarian<A, reco::GenParticle>() {}
   
  DeltaRGenMatchHungarian(const std::vector<A>& objectsA, const std::vector<reco::GenParticle>& objectsB) {
    if(objectsA.size() == 0 || objectsB.size() == 0)
      *this = DeltaRGenMatchHungarian();

    else Solve(objectsA, objectsB);
  }

 private:
  
  void Solve(const std::vector<A> &objectsA, const reco::GenParticleCollection &objectsB) {
    std::map<int, int> indexMap;
    reco::GenParticleCollection statusOne = CleanGenParticles(objectsB, indexMap);

    Matrix<double> costMatrix = this->CalculateCostMatrix(objectsA, statusOne);
    HungarianAlgorithm assigner;
    this->cost_ = assigner.Solve(costMatrix, this->matchedIndexes_);

    for(auto &idx : this->matchedIndexes_) 
      if(idx >= 0) idx = indexMap[idx];
    
    this->matchedPairs_ = this->ConstructMatchedPairs(costMatrix, objectsA, objectsB);
  }
};

#endif
