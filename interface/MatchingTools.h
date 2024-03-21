#ifndef MatchingTools_h
#define MatchingTools_h

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

class TrackInfo;

typedef std::vector<TrackInfo> TrackInfoCollection;
typedef std::vector<std::vector<GlobalPoint> > TrackRecHitLocations;
template <typename S> using Matrix = std::vector<std::vector<S>>;
enum LepType {kW, kZ, kTau, kConversion, kLight, kHeavy, kSusy, kOther, kUnmatched};

class GenLeptonType {

 public:

  GenLeptonType() 
    : index_(-1), genElectron_(reco::GenParticle()), type_(kUnmatched), deltaR_(-1.) {}

 GenLeptonType(const int index, const reco::GenParticle &genElectron, const LepType &type, const double deltaR = -1)
   : index_(index), genElectron_(genElectron), type_(type), deltaR_(deltaR) {}

  int GetIndex() const { return index_; }
  reco::GenParticle GetGenElectron() const { return genElectron_; }
  LepType GetLepType() const { return type_; }
  double GetDeltaR() const { return deltaR_; }
  bool isValid() const { return (index_ > 0 && deltaR_ > 0); }

 private:
  
  int index_;
  reco::GenParticle genElectron_;
  LepType type_;
  double deltaR_;
};

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
    genInfo_(GenLeptonType()),
    isGenMatched_(false) {}

  int GetTrackIndex() const { return pairInfo_.GetIndexA(); }
  int GetSCIndex() const { return pairInfo_.GetIndexB(); }
  double GetDeltaR() const { return deltaR_; }
  T GetTrack() const { return track_; }
  reco::SuperCluster GetSuperCluster() const { return superCluster_; }
  GlobalPoint GetTrackLocationAtECAL() const { return trackLocationAtECAL_; }
  GenLeptonType GetGenInfo() const { return genInfo_; }
  bool isGenMatched() const {return isGenMatched_;}

  void setGenMatchInfo(const GenLeptonType &genInfo) {
    if(genInfo.isValid()) {
      isGenMatched_ = true;
      genInfo_ = genInfo;
    }
  }
  
 private:

  MatchedPair pairInfo_;
  T track_;
  reco::SuperCluster superCluster_;
  double deltaR_;
  GlobalPoint trackLocationAtECAL_;
  GenLeptonType genInfo_;
  bool isGenMatched_;
};

class TrackInfo {

 public:

  template <typename T> TrackInfo(const T &track, const int index) 
    : pt_(track.pt()), eta_(track.eta()), phi_(track.phi()), index_(index) {
    if(typeid(T) == typeid(reco::Track))
      type_ = kGeneral;
    else if(typeid(T) == typeid(reco::GsfTrack))
      type_ = kGsf;
  }

  TrackInfo(const double &pt, const double &eta, const double &phi)
    : pt_(pt), eta_(eta), phi_(phi), index_(-1), type_(kUnkown) {}

  virtual ~TrackInfo() = default;

  // Getters
  int GetIndex() const {return index_;} 
  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}

  bool isGeneral() const { return type_ == kGeneral; }
  bool isGsf() const { return type_ == kGsf; }
  bool isEmpty() {
    return pt_ < 0 && eta_ < 0 && phi_ < 0;
  }

  bool operator == (const TrackInfo& other) const {
    return pt_ == other.pt_ && eta_ == other.eta_ && phi_ == other.phi_;
  }

 private:
  enum type {kGeneral, kGsf, kUnkown};
  double pt_, eta_, phi_;
  int index_;
  type type_;
};

template <typename T>
inline TrackInfoCollection GetTrackInfo(const std::vector<T> &tracks) {
  
  TrackInfoCollection trackInfo;
  for(size_t t = 0; t < tracks.size(); t++) {
    const T track(tracks[t]);
    trackInfo.emplace_back(TrackInfo(track, int(t)));
  }
  return trackInfo;  
}

template <class T> class MatchedTrackSCPairs : public std::vector<MatchedTrackSCPair<T>> {
 public:

  MatchedTrackSCPairs()
    : cost_(-1) {}

  MatchedTrackSCPairs(const double &cost)
    : cost_(cost) {}

  virtual ~MatchedTrackSCPairs() = default;

  double GetCost() const {return cost_;}

 private:

  double cost_;

};

template <typename T>
std::vector<T> RemoveDataAtIndices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete) {
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

inline reco::GenParticleCollection CleanGenParticles(const reco::GenParticleCollection &genParticles) {
  std::vector<size_t> indecesToRemove;

  for(size_t i = 0; i < genParticles.size(); i++) {
    if(genParticles[i].status() == 1)
      continue;    
    indecesToRemove.push_back(i);
  }
  
  return RemoveDataAtIndices<reco::GenParticle>(genParticles, indecesToRemove);
}

#endif
