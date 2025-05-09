#ifndef KUCMSNtupleizer_KUCMSNtupleizer_TrackHelper_h
#define KUCMSNtupleizer_KUCMSNtupleizer_TrackHelper_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class TrackHelper {

 public:

  template <typename T>
  static double GetDXY(const T &object);
  static double GetDxyDiff(const reco::Track &track, const reco::GenParticle &genParticle);
  static double Get3Ddiff(const reco::Track &track, const reco::GenParticle &genParticle);
  static double CalculateCosTheta(const reco::Vertex &primaryVertex, const reco::Vertex &trackVertex, const reco::Track &track);
  static double CalculateCosTheta(const reco::Vertex &primaryVertex, const GlobalPoint &position, const reco::GenParticle &track);
  static int FindTrackIndex(const reco::Track &track, const reco::TrackCollection &trackCollection);
  static int FindTrackIndex(const reco::TransientTrack &track, const std::vector<reco::TransientTrack> &trackCollection);
  template <typename T>
  static int FindIndex(const T &object, const std::vector<T> &objectCollection);
  static bool SameTrack(const reco::Track &track1, const reco::Track &track2);
  static bool OverlappingTrack(const reco::Track &track1, const reco::Track &track2, const TransientTrackBuilder* ttBuilder);
  
  static reco::TrackRef GetTrackRef(const reco::Track &track, const edm::Handle<reco::TrackCollection> &trackHandle);
  static bool AreTrackCollectionsIdentical(const reco::TrackCollection &collection1, const reco::TrackCollection &collection2);
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

inline double TrackHelper::CalculateCosTheta(const reco::Vertex &primaryVertex, const reco::Vertex &trackVertex, const reco::Track &track) {

  const double pvX(primaryVertex.x()), pvY(primaryVertex.y()), pvZ(primaryVertex.z());
  const double vX(trackVertex.position().x()), vY(trackVertex.position().y()), vZ(trackVertex.position().z());
  const double xMag( sqrt((pvX-vX)*(pvX-vX) + (pvY-vY)*(pvY-vY) + (pvZ-vZ)*(pvZ-vZ)));
  const double tPx(track.px()), tPy(track.py()), tPz(track.pz());
  const double pMagTrack(sqrt(tPx*tPx + tPy*tPy + tPz*tPz));
  const double cosTheta( ((vX-pvX)*tPx + (vY-pvY)*tPy + (vZ-pvZ)*tPz) / (xMag*pMagTrack));

  return cosTheta;
}

inline double TrackHelper::CalculateCosTheta(const reco::Vertex &primaryVertex, const GlobalPoint &position, const reco::GenParticle &track) {

  const double pvX(primaryVertex.x()), pvY(primaryVertex.y()), pvZ(primaryVertex.z());
  const double vX(position.x()), vY(position.y()), vZ(position.z());
  const double xMag( sqrt((pvX-vX)*(pvX-vX) + (pvY-vY)*(pvY-vY) + (pvZ-vZ)*(pvZ-vZ)));
  const double tPx(track.px()), tPy(track.py()), tPz(track.pz());
  const double pMagTrack(sqrt(tPx*tPx + tPy*tPy + tPz*tPz));
  const double cosTheta( ((vX-pvX)*tPx + (vY-pvY)*tPy + (vZ-pvZ)*tPz) / (xMag*pMagTrack));

  return cosTheta;
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

inline int TrackHelper::FindTrackIndex(const reco::TransientTrack &track, const std::vector<reco::TransientTrack> &trackCollection) {

  int index = -1;
  for(size_t i = 0; i < trackCollection.size(); i++) {
    reco::TransientTrack trk(trackCollection[i]);
    if(track.track().pt() == trk.track().pt() && track.track().eta() == trk.track().eta()) {
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

inline bool TrackHelper::OverlappingTrack(const reco::Track &track1, const reco::Track &track2, const TransientTrackBuilder* ttBuilder) {
  const reco::TransientTrack ttrack1(ttBuilder->build(track1)), ttrack2(ttBuilder->build(track2));
  GlobalTrajectoryParameters par1(ttrack1.initialFreeState().parameters()), par2(ttrack2.initialFreeState().parameters());

  return (par1.position() - par2.position()).mag2() < 1e-7f && (par1.momentum() - par2.momentum()).mag2() < 1e-7f && par1.charge() == par2.charge();
}

inline reco::TrackRef TrackHelper::GetTrackRef(const reco::Track &track, const edm::Handle<reco::TrackCollection> &trackHandle) {

  double pt = track.pt();
  double eta = track.eta();

  reco::TrackRef trackRef;
  for(size_t i(0); i < trackHandle->size(); i++)
    if(pt == trackHandle->at(i).pt() && eta == trackHandle->at(i).eta()) {
      trackRef = reco::TrackRef(trackHandle, i);
      break;
    }

  return trackRef;
}

inline bool TrackHelper::AreTrackCollectionsIdentical(const reco::TrackCollection &collection1, const reco::TrackCollection &collection2) {
  // First, check if the collections have the same size
  if(collection1.size() != collection2.size()) 
    return false;
    
  // Iterate over each track in the collections and compare them
  for (size_t i = 0; i < collection1.size(); ++i) 
    if(FindTrackIndex(collection1[i], collection2) < 0) 
      return false;  
  
  // All tracks are identical
  return true;
}

#endif
