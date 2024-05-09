#ifndef KUCMSNtupleizer_KUCMSNtupleizer_TrackInfo_h
#define KUCMSNtupleizer_KUCMSNtupleizer_TrackInfo_h

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"

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

  void SetIsolationInfo(const IsolationInfo &isoInfo) { isoInfo_ = isoInfo; }

 private:
  enum type {kGeneral, kGsf, kUnkown};
  double pt_, eta_, phi_;
  int index_;
  type type_;

  IsolationInfo isoInfo_;

};

typedef std::vector<TrackInfo> TrackInfoCollection;

template <typename T>
inline TrackInfoCollection GetTrackInfo(const std::vector<T> &tracks) {

  TrackInfoCollection trackInfo;
  for(size_t t = 0; t < tracks.size(); t++) {
    const T track(tracks[t]);
    trackInfo.emplace_back(TrackInfo(track, int(t)));
  }
  return trackInfo;
}

#endif
