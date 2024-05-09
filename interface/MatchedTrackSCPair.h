#ifndef KUCMSNtupleizer_KUCMSNtupleizer_MatchedTrackSCPair_h
#define KUCMSNtupleizer_KUCMSNtupleizer_MatchedTrackSCPair_h

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenLeptonInfo.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"

//===============================================================================================//                                                                               
//                                  class: MatchedTrackSCPair                                    //                                                                               
//-----------------------------------------------------------------------------------------------//                                                                               
// Matched track-to-supercluster pair class used internally by the MatchTracksToSC class. Aside  //                                                                               
// from the matched track and supercluster it stores the deltaR of the match and the location of //                                                                               
// the extrapolated track at the face of the ECAL. Can take any type of track as input. Gen info //                                                                               
// placeholders have also been added to facilitate gen-matching post track-SC matching.          //                                                                               
//-----------------------------------------------------------------------------------------------//                                                                               
//===============================================================================================//
template <class T> class MatchedTrackSCPair {

 public:
  MatchedTrackSCPair() = default;
  virtual ~MatchedTrackSCPair() = default;

  MatchedTrackSCPair(const MatchedPair &pairInfo,
                     const T &track,
                     const reco::SuperCluster &superCluster,
                     const double deltaR,
                     const GlobalPoint &trackLocationAtECAL)
    : pairInfo_(pairInfo),
    track_(track),
    superCluster_(superCluster),
    deltaR_(deltaR),
    trackLocationAtECAL_(trackLocationAtECAL),
    genInfo_(GenLeptonInfo()),
    isGenMatched_(false),
    isoInfo_(IsolationInfo()){}

  int GetTrackIndex() const { return pairInfo_.GetIndexA(); }
  int GetSCIndex() const { return pairInfo_.GetIndexB(); }
  double GetDeltaR() const { return deltaR_; }
  T GetTrack() const { return track_; }
  reco::SuperCluster GetSuperCluster() const { return superCluster_; }
  GlobalPoint GetTrackLocationAtECAL() const { return trackLocationAtECAL_; }
  GenLeptonInfo GetGenInfo() const { return genInfo_; }
  IsolationInfo GetIsoInfo() const { return isoInfo_; }
  bool isGenMatched() const {return isGenMatched_;}

  void setGenMatchInfo(const GenLeptonInfo &genInfo) {
    if(genInfo.isValid()) {
      isGenMatched_ = true;
      genInfo_ = genInfo;
    }
  }

  //void SetIsoInfo(const IsolationInfo &isoInfo) { isoInfo_ = isoInfo; }
  void AddIsoInfo(const std::string key, const double value) { isoInfo_.insert(std::make_pair(key, value)); }

 private:

  MatchedPair pairInfo_;
  T track_;
  reco::SuperCluster superCluster_;
  double deltaR_;
  GlobalPoint trackLocationAtECAL_;
  GenLeptonInfo genInfo_;
  bool isGenMatched_;
  IsolationInfo isoInfo_;
};

template <class T> class MatchedTrackSCPairs : public std::vector<MatchedTrackSCPair<T>> {
 public:

  MatchedTrackSCPairs()
    : cost_(-1) {}

  MatchedTrackSCPairs(const double &cost)
    : cost_(cost) {}

  virtual ~MatchedTrackSCPairs() = default;

  double GetCost() const {return cost_;}

  reco::SuperClusterCollection GetSCs() const {
    reco::SuperClusterCollection superClusters;
    for(const auto &pair : *this)
      superClusters.emplace_back(pair.GetSuperCluster());
    return superClusters;
  }

 private:

  double cost_;

};

#endif
