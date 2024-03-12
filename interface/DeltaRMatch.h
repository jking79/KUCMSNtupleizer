#ifndef KUNTupleFW_KUCMSNtupleizer_DeltaRMatch_h_
#define KUNTupleFW_KUCMSNtupleizer_DeltaRMatch_h_

#include <algorithm>

#include "DataFormats/Math/interface/deltaR.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/Hungarian.h"

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

class MatchedPair {
  
 public: 

 MatchedPair() 
   : indexA_(-1), indexB_(-1), criteria_(999.) {}

 MatchedPair(const int &indexA, const int &indexB, const double &matchCriteria) 
   : indexA_(indexA), indexB_(indexB), criteria_(matchCriteria) {}

 int GetIndexA() const {return indexA_;}
 int GetIndexB() const {return indexB_;}
 double GetMatchCriteria() const {return criteria_;}
  
 private: 
  
  int indexA_;
  int indexB_; 
  double criteria_;

};

template <class A, class B> class PairedObjects {
  
 public:
  
 PairedObjects(const A &objectA, const B &objectB, const MatchedPair &matchedPair) 
   : objectA_(objectA), objectB_(objectB), matchedPair_(matchedPair) {

    deltaR_ = sqrt(reco::deltaR2(objectA, objectB) );
    //if(deltaR_ < 0.5)
    //CompareObjects();
    //std::cout<< "Let's get that deltaR again to be sure: " << deltaR_ << std::endl;
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
  /*
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
  */
};

template <class A, class B> class DeltaRMatchHungarian {
  
 public:

  DeltaRMatchHungarian(std::vector<std::vector<double> > &costMatrix) {
    
    HungarianAlgorithm assigner;
    cost_ = assigner.Solve(costMatrix, matchedIndexes_);

    matchedPairs_ = ConstructMatchedPairs(costMatrix);
  }

  DeltaRMatchHungarian(const std::vector<A> &objectsA, const std::vector<B> &objectsB) {
    
    const int dimensionA = objectsA.size();
    const int dimensionB = objectsB.size();

    // If either collection is empty, assignment is not possible (obviously)
    // Initialize members to default values and empty vectors
    if(dimensionA == 0 || dimensionB == 0) {
      cost_ = 999.;
      matchedIndexes_ = std::vector<int>();
      matchedPairs_ = std::vector<MatchedPair>();
    }

    else {
      //std::cout << dimensionA << " objects in A and " << dimensionB << " objects in B" << std::endl;
      std::vector<std::vector<double> > costMatrix = CalculateCostMatrix(objectsA, objectsB);
      //PrintCostMatrix(costMatrix);
      HungarianAlgorithm assigner;
      cost_ = assigner.Solve(costMatrix, matchedIndexes_);
      //std::cout << "cost: \n" << cost_ << std::endl;
      matchedPairs_ = ConstructMatchedPairs(costMatrix, objectsA, objectsB);
    }
      
  } 

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

 private:

  double cost_;
  std::vector<int> matchedIndexes_;
  std::vector<A> matchedObjectsA_;
  std::vector<B> matchedObjectsB_;
  std::vector<MatchedPair> matchedPairs_;
  std::vector<PairedObjects<A,B> > pairedObjects_;

  std::vector<MatchedPair> ConstructMatchedPairs(const std::vector<std::vector<double> > &costMatrix, 
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
  
  std::vector<std::vector<double>> CalculateCostMatrix(const std::vector<A> &objectsA, const std::vector<B> &objectsB) const {

    const int dimensionA = objectsA.size();
    const int dimensionB = objectsB.size();
    std::vector<std::vector<double>> costMatrix(dimensionA, std::vector<double>(dimensionB) );

    for(int ai = 0; ai < dimensionA; ai++) 
      for(int bi = 0; bi < dimensionB; bi++) 
	costMatrix[ai][bi] = sqrt(reco::deltaR2(objectsA[ai], objectsB[bi]) );

    return costMatrix;
  }

  void PrintCostMatrix( std::vector<std::vector<double>> &costMatrix) {

    std::cout << "Cost Matrix: " << std::endl;
    for(const auto &row : costMatrix) {
      for(const auto &element : row) {

	std::cout << element << ", ";
      }
      std::cout << std::endl;
    }
  }
};

template <typename T>
std::vector<T> EraseIndices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete) {
  if(indicesToDelete.empty())
    return data;

  std::vector<T> ret;
  ret.reserve(data.size() - indicesToDelete.size());

  std::sort(indicesToDelete.begin(), indicesToDelete.end());

  typename std::vector<T>::const_iterator itBlockBegin = data.begin();
  for(std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++ it) {
    typename std::vector<T>::const_iterator itBlockEnd = data.begin() + *it;
    if(itBlockBegin != itBlockEnd) {
      std::copy(itBlockBegin, itBlockEnd, std::back_inserter(ret));
    }
    itBlockBegin = itBlockEnd + 1;
  }

  if(itBlockBegin != data.end()) {
    std::copy(itBlockBegin, data.end(), std::back_inserter(ret));
  }

  return ret;
}

#endif
