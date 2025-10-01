#ifndef KUCMSNtupleizer_KUCMSNtupleizer_VertexMerger_h
#define KUCMSNtupleizer_KUCMSNtupleizer_VertexMerger_h

#include <boost/math/distributions/chi_squared.hpp>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/ArbitrationInfo.h"

typedef CachingVertex<5>::RefCountedVertexTrack RefCountedVertexTrack;

class VertexMerger {

 public:
  VertexMerger() : iteration_(0), hasGenVertex_(false), debugMode_(false), hasPV_(false), primaryVertex_(reco::Vertex()){}
  
  VertexMerger(const reco::Vertex &primaryVertex) : iteration_(0), hasGenVertex_(false), debugMode_(false), hasPV_(true), primaryVertex_(primaryVertex){}
  
  virtual ~VertexMerger() = default;
  
  std::vector<CachingVertex<5>> mergePairs(const std::vector<CachingVertex<5>>& initialVertices);

  reco::VertexCollection initialSeeds(const std::vector<reco::TransientTrack> &tracks, const edm::Handle<reco::TrackCollection> &trackHandle) const;
  reco::VertexCollection vertices(const std::vector<reco::TransientTrack> &tracks, const edm::Handle<reco::TrackCollection> &trackHandle);

  std::vector<ArbitrationInfo> getArbitrationDecisions() const {return arbitrationDecisions_;}
  std::vector<MergeInfo> getMergeDecisions() const {return mergeDecisions_;}
  
  void setGenVertices(const GenVertices &genVertices, const bool debug = false);

  std::vector<CachingVertex<5>> getRejects() const {return rejects_;}
  
 private:

  // Store working set of vertices
  int iteration_;
  std::vector<CachingVertex<5>> vertices_, rejects_, masterList_;
  GenVertices genVertices_;
  bool hasGenVertex_, debugMode_, hasPV_;
  reco::Vertex primaryVertex_;
  std::vector<ArbitrationInfo> arbitrationDecisions_;
  std::vector<MergeInfo> mergeDecisions_;
  
  double calculateDistance(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const;
  double calculateDistanceUncertainty(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const;
  double calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const;
  double cosThetaWithPV(const TransientVertex &vertex) const;
  double getMass(const TransientVertex &vertex) const;
  
  bool isValidVertex(const CachingVertex<5>& vertex) const;
  bool hasSharedTracks(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const;
  bool hasSharedTracks(const std::vector<CachingVertex<5>> &vertices) const;
  bool tryMergeVertices(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2,
			CachingVertex<5>& merged);
  void arbitrateTracks(CachingVertex<5>& vertex1, CachingVertex<5>& vertex2);

  RefCountedVertexTrack sharedTrack(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const;
  std::vector<RefCountedVertexTrack> getAllUniqueTracks(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const;
  std::vector<RefCountedVertexTrack> allTracksExceptGiven(const CachingVertex<5>& vertex, const RefCountedVertexTrack &toRemove) const;

  CachingVertex<5> refitTracks(const std::vector<RefCountedVertexTrack> &tracks) const;
  std::vector<CachingVertex<5>> getVertexSeeds(const std::vector<reco::TransientTrack> &tracks) const;
  std::vector<CachingVertex<5>> clean(const std::vector<CachingVertex<5>> &vertices) const; 
  
  reco::TrackCollection getOriginalTracks(const CachingVertex<5> &vertex) const;
  
};

#endif
