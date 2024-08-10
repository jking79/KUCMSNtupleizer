#ifndef KUCMSNtupleizer_KUCMSNtupleizer_TrackHelper_h
#define KUCMSNtupleizer_KUCMSNtupleizer_TrackHelper_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class TrackHelper {

 public:

  template <typename T>
    static double GetDXY(const T &object);
  static double GetDxyDiff(const reco::Track &track, const reco::GenParticle &genParticle);
  static double Get3Ddiff(const reco::Track &track, const reco::GenParticle &genParticle);
  static int FindTrackIndex(const reco::Track &track, const reco::TrackCollection &trackCollection);
  template <typename T>
    static int FindIndex(const T &object, const std::vector<T> &objectCollection);
  static bool SameTrack(const reco::Track &track1, const reco::Track &track2);

};

template <typename T> double TrackHelper::GetDXY(const T &object) {
  const double objectDx = object.vx();
  const double objectDy = object.vy();
  const double objectDxy = sqrt(objectDx*objectDx + objectDy*objectDy);
  return objectDxy;
}

inline double TrackHelper::GetDxyDiff(const reco::Track &track, const reco::GenParticle &genParticle) {
  const double diffDx(track.vx()-genParticle.vx());
  const double diffDy(track.vy()-genParticle.vy());
  const double diffDxy = sqrt(diffDx*diffDx + diffDy*diffDy);
  return diffDxy;
}

inline double TrackHelper::Get3Ddiff(const reco::Track &track, const reco::GenParticle &genParticle) {
  const double diffDx(track.vx()-genParticle.vx());
  const double diffDy(track.vy()-genParticle.vy());
  const double diffDz(track.vz()-genParticle.vz());
  const double diff3D = sqrt(diffDx*diffDx + diffDy*diffDy + diffDz*diffDz);
  return diff3D;
}

inline int TrackHelper::FindTrackIndex(const reco::Track &track, const reco::TrackCollection &trackCollection) {

  int index = -1;
  for(size_t i = 0; i < trackCollection.size(); i++) {
    reco::Track trk(trackCollection[i]);
    if(track.pt() == trk.pt() && track.eta() == trk.eta()) {
      index = int(i);
      break;
    }
  }
  return index;
}

template <typename T>
inline int TrackHelper::FindIndex(const T &object, const std::vector<T> &objectCollection) {

  int index = -1;
  for(size_t i = 0; i < objectCollection.size(); i++) {
    T obj(objectCollection[i]);
    if(object.pt() == obj.pt() && object.eta() == obj.eta()) {
      index = int(i);
      break;
    }
  }
  return index;
}

inline bool TrackHelper::SameTrack(const reco::Track &track1, const reco::Track &track2) {
  return track1.pt() == track2.pt() && track1.eta() == track2.eta() && track1.phi() == track2.phi();
}

#endif
