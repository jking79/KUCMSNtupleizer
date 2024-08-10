#ifndef KUCMSNtupleizer_KUCMSNtupleizer_VertexHelper_h
#define KUCMSNtupleizer_KUCMSNtupleizer_VertexHelper_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

class VertexHelper {

 public:
  
  static int CountInstances(const reco::Vertex &vertex, const reco::TrackCollection &tracks);
  static double CalculateTotalPt(const reco::Vertex &vertex);
  static double CalculateTotalPt(const reco::TrackCollection &tracks);
  static double CalculateTotalWeightedPt(const reco::Vertex &vertex);
  static double CalculateTotalTrackWeight(const reco::Vertex &vertex);
  static double CalculateEcalness(const reco::Vertex &vertex, const reco::TrackCollection &ecalTracks);
  static double CalculateTrackOverlap(const reco::Vertex &vertex1, const reco::Vertex &vertex2);
  static double GetDxyDiff(const reco::Vertex &vertex, const reco::GenParticle &genParticle);
  static double Get3Ddiff(const reco::Vertex &vertex, const reco::GenParticle &genParticle);
  static reco::TrackCollection GetTracks(const reco::Vertex &vertex, const bool aboveThreshold = false);
  static reco::Vertex TimeStampVertex(const reco::Vertex &vertex, const double tof); 

  static void SortByPt(reco::VertexCollection &vertices);

  template <typename T>
    static bool isInVertex(const reco::Vertex &vertex, const T &track);

};

template <typename T>
bool VertexHelper::isInVertex(const reco::Vertex &vertex, const T &track) {
  
  bool isInVertex = false;
  for(const auto &vtxTrack : vertex.tracks())
    if(vtxTrack->pt() == track.pt() && vtxTrack->eta() == track.eta())
      isInVertex = true;
  return isInVertex;
}

#endif
