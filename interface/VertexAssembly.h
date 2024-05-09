#ifndef KUCMSNtupleizer_KUCMSNtupleizer_VertexAssembly_h
#define KUCMSNtupleizer_KUCMSNtupleizer_VertexAssembly_h

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

class VertexAssembly {

 public:

  typedef edm::Handle<reco::TrackCollection> TrackHandle;

 VertexAssembly() : ttBuilder_(nullptr), trackHandle_(TrackHandle()), beamSpot_(reco::BeamSpot()) {}

  VertexAssembly(const TransientTrackBuilder* ttBuilder, TrackHandle &trackHandle)
    : ttBuilder_(ttBuilder), trackHandle_(trackHandle), beamSpot_(reco::BeamSpot()) {}

 VertexAssembly(const TransientTrackBuilder* ttBuilder, TrackHandle &trackHandle, reco::BeamSpot &beamSpot) 
   : ttBuilder_(ttBuilder), trackHandle_(trackHandle), beamSpot_(beamSpot) {}

  virtual ~VertexAssembly() = default;

  bool isValid() {return ttBuilder_ == nullptr;}

  reco::VertexCollection CreateVertexCollection(const double ptCut = 0) const;

  template <typename T>
    reco::VertexCollection CreateVertexCollection(const std::vector<T> &tracksSource) const;

 private:
  
  const TransientTrackBuilder* ttBuilder_;
  TrackHandle trackHandle_;
  reco::BeamSpot beamSpot_;

  std::vector<reco::TransientTrack> BuildTransientTracks(const reco::ElectronCollection &electrons) const;
  std::vector<reco::TransientTrack> BuildTransientTracks(const reco::TrackCollection &tracks) const;
  reco::Vertex ConvertFitVertex(const TransientVertex &vertex) const;
  reco::TrackRef GetTrackRef(const reco::Track &track) const;
};

template <typename T>
reco::VertexCollection VertexAssembly::CreateVertexCollection(const std::vector<T> &tracksSource) const {

  reco::VertexCollection vertexCollection;
  if(tracksSource.size() < 2)
    return vertexCollection;

  AdaptiveVertexReconstructor vertexConstructor(2.0, 6.0, 0.5, true);
  std::vector<TransientVertex> transientVertexCollection = vertexConstructor.vertices(BuildTransientTracks(tracksSource), beamSpot_);

  int vtxIndex(0);
  for(const auto &vtx : transientVertexCollection) {
    vertexCollection.emplace_back(ConvertFitVertex(vtx));
    
    for(const auto &track : vtx.originalTracks()) {
      reco::TrackRef trackRef = GetTrackRef(track.track());
      vertexCollection[vtxIndex].add(reco::TrackBaseRef(trackRef), track.track(), vtx.trackWeight(track));
    }
    vtxIndex++;
  }

  //VertexHelper::SortByPt(vertexCollection);

  return vertexCollection;
}

#endif


