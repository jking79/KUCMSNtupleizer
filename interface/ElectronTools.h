#ifndef ElectronTools_hh
#define ElectronTools_hh
/*
class TrackInfo {

public:
  
  enum TrackType {kGeneral, kGsf};

  TrackInfo(const size_t index, const reco::Track &track) 
    : index_(index),
      pt_(track.pt()),
      eta_(track.eta()),
      phi_(track.phi()),
      type_(kGeneral),
      track_(track),
      gsfTrack_(reco::GsfTrack()),
      genMatchIndex_(-1),
      genMatchDeltaR_(-1.),
      isGenMatched_(false) {}

  TrackInfo(const size_t index, const reco::GsfTrack &track)
    : index_(index),
      pt_(track.pt()),
      eta_(track.eta()),
      phi_(track.phi()),
      type_(kGeneral),
      track_(track),
      gsfTrack_(reco::GsfTrack()),
      genMatchIndex_(-1),
      genMatchDeltaR_(-1.),
      isGenMatched_(false) {}

  size_t index() const {return index_;}
  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}
  bool isGenMatched() const {return isGenMatched_;}
  TrackType type() const {return type_;}
  reco::Track track() const {return track_;}
  reco::Track gsfTrack() const {return gsfTrack_;}

  void setGenMatchInfo(const int genElectronIndex, const double deltaR) {
    genMatchIndex_ = genElectronIndex;
    genMatchDeltaR_ = deltaR;
    isGenMatched_ = (genMatchIndex_ > 0 && genMatchDeltaR_ > 0);
  }

private:
  
  const size_t index_;
  const double pt_;
  const double eta_;
  const double phi_;
  const TrackType type_;

  reco::Track track_;
  reco::GsfTrack gsfTrack_;

  int genMatchIndex_;
  double genMatchDeltaR_;
  bool isGenMatched_;

};
*/

class TrackInfo {

 public:

  TrackInfo(const double &pt, const double &eta, const double &phi)
    : pt_(pt), eta_(eta), phi_(phi) {}

  // Getters
  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}

  bool isEmpty() {
    return pt_ < 0 && eta_ < 0 && phi_ < 0;
  }

  bool operator == (const TrackInfo& other) const {
    return pt_ == other.pt_ && eta_ == other.eta_ && phi_ == other.phi_;
  }

 private:

  double pt_, eta_, phi_;
};

template <typename T>
double GetDXY(T &object) {
  const double objectDx = object.vx();
  const double objectDy = object.vy();
  const double objectDxy = sqrt(objectDx*objectDx + objectDy*objectDy);
  return objectDxy;
}

#endif
