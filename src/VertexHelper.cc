#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

int VertexHelper::CountInstances(const reco::Vertex &vertex, const reco::TrackCollection &tracks) {
  
  int instances(0);
  
  for(const auto &vtxTrack : vertex.tracks())
    for(const auto &track: tracks)
      if(vtxTrack->pt() == track.pt() && vtxTrack->eta() == track.eta())
	instances++;
 
  return instances;
}

double VertexHelper::CalculateTotalPt(const reco::TrackCollection &tracks) {

  double pt(0);

  for(const auto &track : tracks)
    pt += track.pt();

  return pt;
}

double VertexHelper::CalculateTotalPt(const reco::Vertex &vertex) {

  double pt(0);

  for(const auto &track : vertex.tracks()) 
    pt += track->pt();
  
  return pt;
}

double VertexHelper::CalculateTotalWeightedPt(const reco::Vertex &vertex) {

  double weightedPt(0);
  for(const auto &trackRef : vertex.tracks())
    weightedPt += trackRef->pt() * vertex.trackWeight(trackRef);

  return weightedPt;
}

void VertexHelper::SortByPt(reco::VertexCollection &vertices) {

  // Lambda function for sorting in descending order based on pt                                                                                                                
  auto compareByPt = [](const reco::Vertex& a, const reco::Vertex& b) {
    return CalculateTotalPt(a) > CalculateTotalPt(b);
  };

  // Sorting the vector using the lambda function                                                                                                                               
  std::sort(vertices.begin(), vertices.end(), compareByPt);
}

double VertexHelper::CalculateTotalTrackWeight(const reco::Vertex &vertex) {

  double totalWeight(0.);
  for(const auto &trackRef : vertex.tracks())
    totalWeight += vertex.trackWeight(trackRef);
  return totalWeight;
}

double VertexHelper::CalculateEcalness(const reco::Vertex &vertex, const reco::TrackCollection &ecalTracks) {
  
  double ecalTrackPt(0);
  for(const auto &track : vertex.tracks())
    for(const auto &ecalTrack: ecalTracks)
      if(track->pt() == ecalTrack.pt() && track->eta() == ecalTrack.eta())
	ecalTrackPt += ecalTrack.pt();

  double ecalness = ecalTrackPt / CalculateTotalPt(vertex);

  return ecalness;
}

double VertexHelper::CalculateTrackOverlap(const reco::Vertex &vertex1, const reco::Vertex &vertex2) {

  reco::TrackCollection tracks(GetTracks(vertex1));

  int sharedTrackCount(0);
  for(const auto &trackRef : vertex2.tracks())
    if(isInVertex(vertex1, *trackRef))
      sharedTrackCount++;
    else
      tracks.emplace_back(*trackRef);

  return double(sharedTrackCount) / double(tracks.size());
}

reco::TrackCollection VertexHelper::GetTracks(const reco::Vertex &vertex, const bool aboveThreshold) {

  reco::TrackCollection tracks;

  for(const auto &trackRef : vertex.tracks())
    if(aboveThreshold && vertex.trackWeight(trackRef) > 0.5)
      tracks.emplace_back(*trackRef);
    else if (!aboveThreshold)
      tracks.emplace_back(*trackRef);

  return tracks;
}

reco::Vertex VertexHelper::TimeStampVertex(const reco::Vertex &vertex, const double tof) {

  const size_t nTracks(vertex.tracksSize());
  const double tof_error(0.1); // error on the time measurement (use 10% for now)

  const GlobalError error = vertex.error();
  const math::Error<4>::type vertexError(ROOT::Math::SVector<double, 10>(error.cxx(), error.cyx(), error.cyy(), error.czx(), error.czy(),
									error.czz(), tof_error, tof_error, tof_error, tof_error));
  
  reco::Vertex outputVertex(vertex.position(), vertexError, tof, vertex.chi2(), vertex.ndof(), nTracks);
  outputVertex.reserve(nTracks);

  return outputVertex;
}

LorentzVec VertexHelper::GetVertex4Vector(const reco::Vertex &vertex) {

  LorentzVec vertex4Vec;

  for(const auto &trackRef : vertex.tracks()) {
    const double pt(trackRef->pt());
    const double eta(trackRef->eta());
    const double phi(trackRef->phi());
    const double mass(0.);

    LorentzVec track4Vec(pt, eta, phi, mass);
    vertex4Vec += track4Vec;
  }

  return vertex4Vec;
}

double VertexHelper::GetDxyDiff(const reco::Vertex &vertex, const reco::GenParticle &genParticle) {
  const double diffDx(vertex.x()-genParticle.vx());
  const double diffDy(vertex.y()-genParticle.vy());
  const double diffDxy = sqrt(diffDx*diffDx + diffDy*diffDy);
  return diffDxy;
}

double VertexHelper::Get3Ddiff(const reco::Vertex &vertex, const reco::GenParticle &genParticle) {
  const double diffDx(vertex.x()-genParticle.vx());
  const double diffDy(vertex.y()-genParticle.vy());
  const double diffDz(vertex.z()-genParticle.vz());
  const double diff3D = sqrt(diffDx*diffDx + diffDy*diffDy + diffDz*diffDz);
  return diff3D;
}
