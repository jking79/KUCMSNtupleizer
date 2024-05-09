#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexAssembly.h"

reco::VertexCollection VertexAssembly::CreateVertexCollection(const double ptCut) const {

  if(ptCut < 1e-2)
    return CreateVertexCollection<reco::Track>(*trackHandle_);

  reco::TrackCollection selectedTracks;
  for(const auto &track : *trackHandle_)
    if(track.pt() > ptCut)
      selectedTracks.emplace_back(track);

  return CreateVertexCollection<reco::Track>(selectedTracks);  
}

std::vector<reco::TransientTrack> VertexAssembly::BuildTransientTracks(const reco::ElectronCollection &electrons) const {

  std::vector<reco::TransientTrack> ttracks;
  for(const auto &electron : electrons)
    ttracks.emplace_back(ttBuilder_->build(electron.track()));

  return ttracks;
}

std::vector<reco::TransientTrack> VertexAssembly::BuildTransientTracks(const reco::TrackCollection &tracks) const {

  std::vector<reco::TransientTrack> ttracks;
  for(const auto &track : tracks)
    ttracks.emplace_back(ttBuilder_->build(track));

  return ttracks;
}

reco::Vertex VertexAssembly::ConvertFitVertex(const TransientVertex &vertex) const {

  // Parameters needed for initializing new vertex object
  const math::XYZPoint vertexPosition = math::XYZPoint(vertex.position().x(), vertex.position().y(), vertex.position().z());
  const double vertexChi2 = vertex.totalChiSquared();
  const double vertexNdof = vertex.degreesOfFreedom();
  const size_t nTracks = vertex.originalTracks().size();
  const GlobalError error = vertex.positionError();
  const math::Error<3>::type vertexError(ROOT::Math::SVector<double, 6>(error.cxx(), error.cyx(), error.cyy(),
                                                                        error.czx(), error.czy(), error.czz()));

  // Save the new fitted vertex to a reco::Vertex and reserve space for the tracks
  reco::Vertex outputVertex(vertexPosition, vertexError, vertexChi2, vertexNdof, nTracks);
  outputVertex.reserve(nTracks);

  return outputVertex;
}

reco::TrackRef VertexAssembly::GetTrackRef(const reco::Track &track) const {

  double pt = track.pt();
  double eta = track.eta();

  reco::TrackRef trackRef;
  for(size_t i(0); i < trackHandle_->size(); i++) 
    if(pt == trackHandle_->at(i).pt() && eta == trackHandle_->at(i).eta()) 
      trackRef = reco::TrackRef(trackHandle_, i);
  
  return trackRef;
}
