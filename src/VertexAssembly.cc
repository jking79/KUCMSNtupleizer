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

reco::VertexCollection RefineVertices(const reco::VertexCollection &vertices) {

  return reco::VertexCollection();
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

  //std::cout << "\nStarting ConvertFitVertex" << std::endl;
  // Parameters needed for initializing new vertex object
  const math::XYZPoint vertexPosition = math::XYZPoint(vertex.position().x(), vertex.position().y(), vertex.position().z());
  //std::cout << "vertex position: (" << vertexPosition.x() << ", " << vertexPosition.y() << ", " << vertexPosition.z() << ") " <<  std::endl;
  const double vertexChi2 = vertex.totalChiSquared();
  //std::cout << "Chi2: " << vertexChi2 << std::endl;
  const double vertexNdof = vertex.degreesOfFreedom();
  //std::cout << "ndof: " << vertexNdof << std::endl;
  const size_t nTracks = vertex.originalTracks().size();
  //std::cout << "Total tracks: " << nTracks << std::endl;
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

reco::VertexCollection VertexAssembly::TryRefineVertex(const reco::Vertex &vertex, const bool useBeamSpot) const {

  const int totalTracks(vertex.tracksSize());
  reco::VertexCollection splitVertices;//(CreateVertexCollection(VertexHelper::GetTracks(vertex), useBeamSpot));

  // Check whether there are any tracks below threshold weight, 
  if(IsVertexUnrefined(vertex))
    splitVertices = CreateVertexCollection(VertexHelper::GetTracks(vertex, true), useBeamSpot);
  else {
    splitVertices.emplace_back(vertex);
    return splitVertices;
  }
    
  reco::VertexCollection refinedVertices;
  const double originalWeight(VertexHelper::CalculateTotalTrackWeight(vertex));

  std::cout << "original vertex weight: " << originalWeight << " and weight ratio: " << originalWeight/totalTracks << ", total tracks: " << totalTracks << std::endl; 
  
  for(size_t v = 0; v < splitVertices.size(); v++) {
    const reco::Vertex vtx(splitVertices[v]);
    const double totalWeight(VertexHelper::CalculateTotalTrackWeight(vtx));
    
    if(!IsVertexUnrefined(vtx)) {
      refinedVertices.emplace_back(vtx);
      continue;
    }
    
    std::cout << "  first split vertex weight: " << totalWeight << " and weight ratio: " << totalWeight/vtx.tracksSize() << ", total tracks: " << vtx.tracksSize() << std::endl;

    reco::VertexCollection newSplit;
    if(v == 0 && useBeamSpot)
      newSplit = CreateVertexCollection(VertexHelper::GetTracks(vtx, true), true);
    else
      newSplit = CreateVertexCollection(VertexHelper::GetTracks(vtx, true), false);
    
    for(const auto &newVtx : newSplit)
      //if(IsVertexUnrefined(newVtx))
	std::cout << "    second split vertex weight: " << VertexHelper::CalculateTotalTrackWeight(newVtx)
		  << " and weight ratio: " << VertexHelper::CalculateTotalTrackWeight(newVtx)/newVtx.tracksSize() << ", total tracks: " << newVtx.tracksSize() <<  std::endl;
      
    //}
    
  }

  return refinedVertices;
}

bool VertexAssembly::IsVertexUnrefined(const reco::Vertex &vertex) const {

  bool isUnrefined(false);
  for(const auto &trackRef : vertex.tracks())
    if(vertex.trackWeight(trackRef) < 0.5) {
      isUnrefined = true;
      break;
    }
  return isUnrefined;
}
