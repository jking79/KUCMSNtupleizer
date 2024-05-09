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

void VertexHelper::SortByPt(reco::VertexCollection &vertices) {

  // Lambda function for sorting in descending order based on pt                                                                                                                
  auto compareByPt = [](const reco::Vertex& a, const reco::Vertex& b) {
    return CalculateTotalPt(a) > CalculateTotalPt(b);
  };

  // Sorting the vector using the lambda function                                                                                                                               
  std::sort(vertices.begin(), vertices.end(), compareByPt);
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

reco::TrackCollection VertexHelper::GetTracks(const reco::Vertex &vertex) {

  reco::TrackCollection tracks;

  for(const auto &trackRef : vertex.tracks())
    tracks.emplace_back(*trackRef);

  return tracks;
}
