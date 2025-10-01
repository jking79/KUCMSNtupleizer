#ifndef KUCMSNtupleizer_KUCMSNtupleizer_VertexHelper_h
#define KUCMSNtupleizer_KUCMSNtupleizer_VertexHelper_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"
#include <Math/Vector4D.h>
#include <TLorentzVector.h>

typedef ROOT::Math::PtEtaPhiMVector LorentzVec;

class VertexHelper {

 public:

  static int CountInstances(const reco::Vertex &vertex, const reco::TrackCollection &tracks);
  static int CountInstances(const reco::TrackCollection &tracks1, const reco::TrackCollection &tracks2);
  static int CalculateTotalCharge(const reco::Vertex &vertex);
  static double CalculateTotalPt(const reco::Vertex &vertex);
  static double CalculateTotalPt(const reco::TrackCollection &tracks);
  static double CalculateTotalWeightedPt(const reco::Vertex &vertex);
  static double CalculateTotalTrackWeight(const reco::Vertex &vertex);
  static double CalculateEcalness(const reco::Vertex &vertex, const reco::TrackCollection &ecalTracks);
  static double CalculateTrackOverlap(const reco::Vertex &vertex1, const reco::Vertex &vertex2);
  static double CalculateDxy(const reco::Vertex &vertex, const reco::Vertex &primaryVertex);
  static double CalculateDxyError(const reco::Vertex &vertex);
  static double CalculateDxyError(const reco::Vertex &vertex, const reco::Vertex &primaryVertex);
  static double GetDxyDiff(const reco::Vertex &vertex, const reco::GenParticle &genParticle);
  static double Get3Ddiff(const reco::Vertex &vertex, const reco::GenParticle &genParticle);
  static double CalculateCosTheta(const reco::Vertex &primaryVertex, const reco::Vertex &vertex);
  static double CalculateCMCosTheta(const reco::Vertex &vertex, const reco::Track &track);
  static double CalculateDecayAngle(const reco::Vertex &vertex);
  static LorentzVec GetVertex4Vector(const reco::Vertex &vertex);
  static LorentzVec GetVertexWeighted4Vector(const reco::Vertex &vertex);
  static TLorentzVector GetTracks4Vector(const reco::TrackCollection &tracks);
  static reco::TrackCollection GetTracks(const reco::Vertex &vertex, const bool aboveThreshold = false);
  static reco::Vertex TimeStampVertex(const reco::Vertex &vertex, const double tof); 
  static std::vector<reco::TrackCollection> GetAllTrackPairs(const reco::TrackCollection &tracks);  
  static std::pair<reco::TrackCollection, reco::TrackCollection> FindMinMassHemispheres(const reco::Vertex &vertex);
  static double CalculateHadronicDecayAngle(const reco::Vertex &vertex);
  // VertexMerger Helper Functions
  static bool IsSameVertex(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static bool HasDuplicateVertices(const std::vector<CachingVertex<5>> &vertices);
  static bool IsInCollection(const TransientVertex &vertex, const std::vector<CachingVertex<5>> &collection);
  static int CountSharedTracks(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static double GetDxySeparation(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static double GetDxySeparationError(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static double GetDxySeparationSignificance(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static double GetDistance(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static double	GetDistanceError(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static double	GetDistanceSignificance(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static TransientVertex Merge(const TransientVertex &vertex1, const TransientVertex &vertex2);
  static reco::TransientTrack GetCommonTrack(const TransientVertex &vertex1, const TransientVertex &vertex2);
  
  /*
  static double CalculateDistance(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2);
  static double CalculateDistanceUncertainty(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2);
  static double CalculateDistanceSignificance(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2);
  */
  static void SortByPt(reco::VertexCollection &vertices);

  template <typename T>
    static bool isInVertex(const reco::Vertex &vertex, const T &track);

  static bool isInAnyVertex(const reco::VertexCollection &vertices, const reco::Track &track);

  static reco::Vertex ConvertFitVertex(const TransientVertex &vertex);
  
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
