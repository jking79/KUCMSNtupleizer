#ifndef KUCMSNtupleizer_KUCMSNtupleizer_MatchedPair_h
#define KUCMSNtupleizer_KUCMSNtupleizer_MatchedPair_h

//===============================================================================================//
//                                  class: MatchedPair                                           //
//-----------------------------------------------------------------------------------------------//
// Class for keeping track of the index and matching criteria (such as deltaR) of a matched pair //
// of any two objects that get matched together. Used in the DeltaRMatch, DeltaRMatchHungarian   //
// and MatchTracksToSC classes.                                                                  //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
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

#endif
