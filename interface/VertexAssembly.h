#ifndef KUCMSNtupleizer_KUCMSNtupleizer_VertexAssembly_h
#define KUCMSNtupleizer_KUCMSNtupleizer_VertexAssembly_h

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DisplacedAVR.h"
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
  reco::VertexCollection CreateVertexCollection(const std::vector<T> &tracksSource, const bool useBeamSpot=true) const;

  template <typename T>
  reco::VertexCollection CreateVertexCollection(const std::vector<T> &tracksSource, const edm::ParameterSet &s, const bool useBeamSpot=true) const;
  
  reco::VertexCollection RefineVertices(const reco::VertexCollection &vertices);
  reco::VertexCollection TryRefineVertex(const reco::Vertex &vertex, const bool useBeamSpot = false) const;
 private:
  
  const TransientTrackBuilder* ttBuilder_;
  TrackHandle trackHandle_;
  reco::BeamSpot beamSpot_;

  std::vector<reco::TransientTrack> BuildTransientTracks(const reco::ElectronCollection &electrons) const;
  std::vector<reco::TransientTrack> BuildTransientTracks(const reco::TrackCollection &tracks) const;
  reco::Vertex ConvertFitVertex(const TransientVertex &vertex) const;
  reco::TrackRef GetTrackRef(const reco::Track &track) const;
  //reco::VertexCollection TryRefineVertex(const reco::Vertex &vertex, const bool useBeamSpot = false) const;

  bool IsVertexUnrefined(const reco::Vertex &vertex) const;
  
};

template <typename T>
reco::VertexCollection VertexAssembly::CreateVertexCollection(const std::vector<T> &tracksSource, const bool useBeamSpot) const {

  reco::VertexCollection vertexCollection;
  if(tracksSource.size() < 2)
    return vertexCollection;

  //AdaptiveVertexReconstructor vertexConstructor(2.0, 6.0, 0.5, false);
  DisplacedAVR vertexConstructor(2.0, 6.0, 0.5, false);
  std::vector<TransientVertex> transientVertexCollection;
  if(useBeamSpot)
    transientVertexCollection = vertexConstructor.vertices(BuildTransientTracks(tracksSource), beamSpot_);
  else
    transientVertexCollection = vertexConstructor.vertices(BuildTransientTracks(tracksSource));

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

template <typename T>
reco::VertexCollection VertexAssembly::CreateVertexCollection(const std::vector<T> &tracksSource, const edm::ParameterSet &s, const bool useBeamSpot) const {

  reco::VertexCollection vertexCollection;
  if(tracksSource.size() < 2)
    return vertexCollection;

  //AdaptiveVertexReconstructor vertexConstructor(s);
  DisplacedAVR vertexConstructor(s); 
  std::vector<TransientVertex> transientVertexCollection;
  if(useBeamSpot)
    transientVertexCollection = vertexConstructor.vertices(BuildTransientTracks(tracksSource), beamSpot_);
  else
    transientVertexCollection = vertexConstructor.vertices(BuildTransientTracks(tracksSource));

  int vtxIndex(0);
  for(const auto &vtx : transientVertexCollection) {
    vertexCollection.emplace_back(ConvertFitVertex(vtx));

    for(const auto &track : vtx.originalTracks()) {
      reco::TrackRef trackRef = GetTrackRef(track.track());
      vertexCollection[vtxIndex].add(reco::TrackBaseRef(trackRef), track.track(), vtx.trackWeight(track));
    }
    vtxIndex++;
  }

  return vertexCollection;
}

#endif


