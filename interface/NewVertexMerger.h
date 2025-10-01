#ifndef KUCMSNtupleizer_KUCMSNtupleizer_NewVertexMerger_h
#define KUCMSNtupleizer_KUCMSNtupleizer_NewVertexMerger_h

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSet.h"

class NewVertexMerger {

 public:
  
  NewVertexMerger() = default;
  NewVertexMerger(const reco::Vertex &primaryVertex, const double seedCosThetaCut = 0.75); 
  
  virtual ~NewVertexMerger() = default;

  std::vector<double> mergeDecisions() const { return mergeDecisions_; }
  std::vector<double> goldMergeDecisions() const { return goldMergeDecisions_; }
  
  reco::VertexCollection seeds(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) const;
  reco::VertexCollection vertices(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder);
  TrackVertexSetCollection trackVertexSets(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder);
  TrackVertexSetCollection trackVertexSetSeeds(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder);
  
  void setGenVertices(const GenVertices &genVertices) { genVertices_ = genVertices; }
    
 private:

  reco::Vertex primaryVertex_;
  double seedCosThetaCut_;
  TrackVertexSetCollection masterList_, vertices_;
  std::vector<double> mergeDecisions_, goldMergeDecisions_;
  GenVertices genVertices_;

  bool debug = false;
  
  TrackVertexSetCollection createSeeds(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) const;
  TrackVertexSetCollection mergePairs(const TrackVertexSetCollection &seeds);
  
  bool isValidVertex(const TrackVertexSet &set) const;

  void filterWithAlignment(TrackVertexSetCollection &collection) const;
};

#endif
