#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

int VertexHelper::CountInstances(const reco::Vertex &vertex, const reco::TrackCollection &tracks) {
  
  int instances(0);
  if(vertex.hasRefittedTracks()) {
    for(const auto &vtxTrack : vertex.refittedTracks())
      for(const auto &track: tracks)
        if(vtxTrack.pt() == track.pt() && vtxTrack.eta() == track.eta())
          instances++;
  }
  else {
    for(const auto &vtxTrack : vertex.tracks())
      for(const auto &track: tracks)
	if(vtxTrack->pt() == track.pt() && vtxTrack->eta() == track.eta())
	  instances++;
  }
  return instances;
}

int VertexHelper::CountInstances(const reco::TrackCollection &tracks1, const reco::TrackCollection &tracks2) {

  int instances(0);

  for(const auto &track1 : tracks1)
    for(const auto &track2: tracks2)
      if(track1.pt() == track2.pt() && track1.eta() == track2.eta())
        instances++;

  return instances;
}

int VertexHelper::CalculateTotalCharge(const reco::Vertex &vertex) {
  int sumCharge(0);
  
  if(vertex.hasRefittedTracks()) {
    for(const auto &track: vertex.refittedTracks())
    sumCharge += track.charge();
  }
  else {
    for(const auto &refTrack: vertex.tracks())
      sumCharge += refTrack->charge();
  }

  return sumCharge;
}

double VertexHelper::CalculateTotalPt(const reco::TrackCollection &tracks) {

  double pt(0);

  for(const auto &track : tracks)
    pt += track.pt();

  return pt;
}

double VertexHelper::CalculateTotalPt(const reco::Vertex &vertex) {

  double pt(0);
  if(vertex.hasRefittedTracks()) {
    for(const auto &track : vertex.refittedTracks())
      pt += track.pt();
  }
  else {
    for(const auto &track : vertex.tracks()) 
      pt += track->pt();
  }
  
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

double VertexHelper::CalculateCosTheta(const reco::Vertex &primaryVertex, const reco::Vertex &vertex) {

  LorentzVec lorentzVec;

  const double pvX(primaryVertex.x()), pvY(primaryVertex.y()), pvZ(primaryVertex.z());
  const double vX(vertex.position().x()), vY(vertex.position().y()), vZ(vertex.position().z());
  const double xMag( sqrt((pvX-vX)*(pvX-vX) + (pvY-vY)*(pvY-vY) + (pvZ-vZ)*(pvZ-vZ)));

  for(const auto &track : vertex.tracks()) {
    const double pt(track->pt()), eta(track->eta()), phi(track->phi()), mass(0.);
    lorentzVec += LorentzVec(pt, eta, phi, mass);
  }

  const double vPx(lorentzVec.px()), vPy(lorentzVec.py()), vPz(lorentzVec.pz());
  const double pMagVertex(sqrt(vPx*vPx + vPy*vPy + vPz*vPz));
  const double cosTheta( ((vX-pvX)*vPx + (vY-pvY)*vPy + (vZ-pvZ)*vPz) / (xMag*pMagVertex));

  return cosTheta;
}

double VertexHelper::CalculateCMCosTheta(const reco::Vertex &vertex, const reco::Track &track) {

  if(!isInVertex(vertex, track))
    return -999;
  
  TLorentzVector lorentzVec;
  for(const auto &trackRef : vertex.tracks()) {
    const double pt(trackRef->pt()), eta(trackRef->eta()), phi(trackRef->phi()), mass(0.);
    TLorentzVector track4Vec;
    track4Vec.SetPtEtaPhiM(pt, eta, phi, mass);
    lorentzVec += track4Vec;
  }

  TVector3 beta = lorentzVec.BoostVector();
  TLorentzVector trackVec;
  trackVec.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), 0.);

  trackVec.Boost(-beta);
  return trackVec.Vect().Unit().Dot(beta.Unit());
}

double VertexHelper::CalculateDecayAngle(const reco::Vertex &vertex) {

  if(vertex.tracksSize() > 2)
    return CalculateHadronicDecayAngle(vertex);
  
  TLorentzVector lorentzVec;
  std::vector<TLorentzVector> track4Vecs;
  
  for(const auto &track : vertex.tracks()) {
    const double pt(track->pt()), eta(track->eta()), phi(track->phi()), mass(0.);
    TLorentzVector track4Vec;
    track4Vec.SetPtEtaPhiM(pt, eta, phi, mass);
    lorentzVec += track4Vec;
    track4Vecs.emplace_back(track4Vec);
  }

  TVector3 beta = lorentzVec.BoostVector();
  TLorentzVector track1 = track4Vecs[0];
  track1.Boost(-beta);
  
  return track1.Vect().Unit().Dot(beta.Unit());
}

double VertexHelper::CalculateDxy(const reco::Vertex &vertex, const reco::Vertex &primaryVertex) {

  const double x(vertex.x()), y(vertex.y());
  const double pvx(primaryVertex.x()), pvy(primaryVertex.y());
  const double dx(x-pvx), dy(y-pvy);
  const double dxy(sqrt(dx*dx + dy*dy));
  
  return dxy;
}

double VertexHelper::CalculateDxyError(const reco::Vertex &vertex) {

  const double C_delta_xx = vertex.error().At(0,0);
  const double C_delta_yy = vertex.error().At(1,1);
  const double C_delta_xy = vertex.error().At(0,1);

  double dxy_error = std::sqrt(C_delta_xx + C_delta_yy + 2 * C_delta_xy);

  return dxy_error;
}

double VertexHelper::CalculateDxyError(const reco::Vertex &vertex, const reco::Vertex &primaryVertex) {
  const double x(vertex.x()), y(vertex.y());
  const double pvx(primaryVertex.x()), pvy(primaryVertex.y());
  const double cxx = vertex.error().At(0,0);
  const double cyy = vertex.error().At(1,1);
  const double cxy = vertex.error().At(0,1);
  const double pvcxx = primaryVertex.error().At(0,0);
  const double pvcyy = primaryVertex.error().At(1,1);
  const double pvcxy = primaryVertex.error().At(0,1);
  const double dx(x-pvx), dy(y-pvy);
  const double dxy(sqrt(dx*dx + dy*dy));
  const double dxyError(sqrt(dx*dx*(cxx+pvcxx) + dy*dy*(cyy+pvcyy) + 2*dx*dy*(cxy+pvcxy))/dxy);

  return dxyError;
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

std::vector<reco::TrackCollection> VertexHelper::GetAllTrackPairs(const reco::TrackCollection &tracks) {

  std::vector<std::vector<reco::Track>> pairs;
  const size_t n = tracks.size();
  
  // Reserve space to avoid reallocations
  pairs.reserve(n * (n - 1) / 2);
  
  for (size_t i = 0; i < n - 1; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      pairs.push_back({tracks[i], tracks[j]});
    }
  }

  return pairs;
}
/*
std::pair<reco::TrackCollection, reco::TrackCollection> VertexHelper::FindMinMassHemispheres(const reco::Vertex &vertex) {

  TLorentzVector lorentzVec;
  std::vector<TLorentzVector> track4Vecs;
  
  for(const auto &track : vertex.tracks()) {
    const double pt(track->pt()), eta(track->eta()), phi(track->phi()), mass(0.);
    TLorentzVector track4Vec;
    track4Vec.SetPtEtaPhiM(pt, eta, phi, mass);
    lorentzVec += track4Vec;
    track4Vecs.emplace_back(track4Vec);
  }

  TVector3 betaVec = lorentzVec.BoostVector();

  for(TLorentzVector &fourVec : track4Vecs)
    fourVec.Boost(-betaVec);

  reco::TrackCollection hemisphere1, hemisphere2;
  
  return std::make_pair(hemisphere1, hemisphere2);
}
*/

std::pair<reco::TrackCollection, reco::TrackCollection> VertexHelper::FindMinMassHemispheres(const reco::Vertex &vertex) {
  
  TLorentzVector lorentzVec;
  std::vector<TLorentzVector> track4Vecs;
  std::vector<const reco::Track*> trackPtrs;

  // Create 4-vectors and compute the total 4-momentum
  for (const auto &track : vertex.tracks()) {
    const double pt = track->pt(), eta = track->eta(), phi = track->phi(), mass = 0.0;
    TLorentzVector track4Vec;
    track4Vec.SetPtEtaPhiM(pt, eta, phi, mass);
    lorentzVec += track4Vec;
    track4Vecs.emplace_back(track4Vec);
    trackPtrs.emplace_back(track.get());
  }

  // Boost all track 4-vectors to the CM frame
  TVector3 betaVec = lorentzVec.BoostVector();
  for (TLorentzVector &vec : track4Vecs) {
    vec.Boost(-betaVec);
  }

  // Initialize best configuration
  reco::TrackCollection bestHem1, bestHem2;
  double bestMetric = -1;

  const int N = track4Vecs.size();

  // Loop over all seed pairs to define candidate hemisphere planes
  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {

      TVector3 refNormal = track4Vecs[i].Vect().Cross(track4Vecs[j].Vect());

      // Skip degenerate pairs (collinear)
      if (refNormal.Mag2() < 1e-6) continue;

      reco::TrackCollection hem1, hem2;
      TLorentzVector hemVec1, hemVec2;

      // Assign tracks to candidate hemispheres
      for (int k = 0; k < N; ++k) {
        if (k == i || k == j) continue;
        const double dot = track4Vecs[k].Vect().Dot(refNormal);
        if (dot > 0) {
          hem1.push_back(*trackPtrs[k]);
          hemVec1 += track4Vecs[k];
        } else {
          hem2.push_back(*trackPtrs[k]);
          hemVec2 += track4Vecs[k];
        }
      }

      // Attempt all track seed assignment configurations
      for(int h1 = 0; h1 < 2; ++h1) {
	for(int h2 = 0; h2 < 2; ++h2) {

	  // Don't assign both tracks to the same hemisphere if the opposite one is empty
	  if(h1 == h2 && (h1 == 0? hem2.empty() : hem1.empty())) continue;

	  reco::TrackCollection tmpHem[2] = {hem1, hem2};
	  TLorentzVector tmpVec[2] = {hemVec1, hemVec2};

	  tmpHem[h1].emplace_back(*trackPtrs[i]);
	  tmpHem[h2].emplace_back(*trackPtrs[j]);

	  tmpVec[h1] += track4Vecs[i];
	  tmpVec[h2] +=	track4Vecs[j];

	  double metric = tmpVec[0].P() + tmpVec[1].P();
	  if (metric > bestMetric) {
	    bestMetric = metric;
	    bestHem1 = tmpHem[0];
	    bestHem2 = tmpHem[1];
	  }
	}
      }
    }
  }

  //if((bestHem1.size() + bestHem2.size()) != vertex.tracksSize())
  //std::cout << "We have a problem" << std::endl;
  
  return std::make_pair(bestHem1, bestHem2);
}

double VertexHelper::CalculateHadronicDecayAngle(const reco::Vertex &vertex) {

  const auto hemispheres = FindMinMassHemispheres(vertex);

  TLorentzVector vertex4Vec(GetTracks4Vector(GetTracks(vertex)));
  TLorentzVector hem4Vec1(GetTracks4Vector(hemispheres.first));

  TVector3 betaVec = vertex4Vec.BoostVector();
  
  hem4Vec1.Boost(-betaVec);

  return hem4Vec1.Vect().Unit().Dot(betaVec.Unit());
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

LorentzVec VertexHelper::GetVertexWeighted4Vector(const reco::Vertex &vertex) {

  LorentzVec vertex4Vec;

  for(const auto &trackRef : vertex.tracks()) {
    const double pt(trackRef->pt());
    const double eta(trackRef->eta());
    const double phi(trackRef->phi());
    const double mass(0.);

    const double trackWeight(vertex.trackWeight(trackRef));
    const double vertexWeight(CalculateTotalTrackWeight(vertex));
    const double ratio(trackWeight / vertexWeight);

    //std::cout << "track weight / vertex weight = " << ratio << std::endl; 
    
    LorentzVec track4Vec(pt, eta, phi, mass);
    if(ratio > 0)
      vertex4Vec += ratio*track4Vec;
  }

  //std::cout << "four vector mass = " << vertex4Vec << std::endl;

  return vertex4Vec;
}

TLorentzVector VertexHelper::GetTracks4Vector(const reco::TrackCollection &tracks) {

  TLorentzVector tracks4Vec;

  for(const auto &track : tracks) {
    const double pt(track.pt());
    const double eta(track.eta());
    const double phi(track.phi());
    const double mass(0.);

    TLorentzVector track4Vec;
    track4Vec.SetPtEtaPhiM(pt, eta, phi, mass);
    tracks4Vec += track4Vec;
  }
  return tracks4Vec;
  
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

bool VertexHelper::isInAnyVertex(const reco::VertexCollection &vertices, const reco::Track &track) {

  for(const auto &vertex : vertices) {
    if(isInVertex<reco::Track>(vertex, track))
      return true;
  }
    
  return false;
}

reco::Vertex VertexHelper::ConvertFitVertex(const TransientVertex &vertex) {

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

// VertexMerger
bool VertexHelper::IsSameVertex(const TransientVertex &vertex1, const TransientVertex &vertex2) {

  bool isSame(false);
  if(vertex1.originalTracks().size() == vertex2.originalTracks().size() && CountSharedTracks(vertex1,vertex2) == int(vertex1.originalTracks().size()))
    isSame = true;

  return isSame;
}

bool VertexHelper::HasDuplicateVertices(const std::vector<CachingVertex<5>> &vertices) {
  
  bool hasDuplicates(false);
  std::vector<CachingVertex<5>>::const_iterator end = vertices.end();
  std::vector<CachingVertex<5>>::const_iterator endm1 = (end - 1);
  for(std::vector<CachingVertex<5>>::const_iterator x = vertices.begin(); x != endm1; ++x) {
    for(std::vector<CachingVertex<5>>::const_iterator y = x + 1; y != end; ++y) {
      if(x->tracks().size() == y->tracks().size() && CountSharedTracks(*x, *y) == int(x->tracks().size())) {
	hasDuplicates = true;
	break;
      }
    }
    if(hasDuplicates) break;
  }

  return hasDuplicates;
}

bool VertexHelper::IsInCollection(const TransientVertex &vertex, const std::vector<CachingVertex<5>> &collection) {

  if(collection.empty())
    return false;
  
  bool isInCollection(false);
  
  std::vector<CachingVertex<5>>::const_iterator end = collection.end();
  std::vector<CachingVertex<5>>::const_iterator endm1 = (end - 1);
  for(std::vector<CachingVertex<5>>::const_iterator x = collection.begin(); x != endm1; ++x) {
    for(std::vector<CachingVertex<5>>::const_iterator y = x + 1; y != end; ++y) {
      if(IsSameVertex(*x,*y)) {
	//std::cout << "Found the duplicate" << std::endl;
	isInCollection = true;
        break;
      }
    }
    if(isInCollection) break;
  }
  //if(isInCollection)
  //std::cout << "Why isn't this working?" << std::endl;
  return isInCollection;
}

int VertexHelper::CountSharedTracks(const TransientVertex &vertex1, const TransientVertex &vertex2) {
  int nCommonTracks(0);
  for(const auto& trackV1 : vertex1.originalTracks()) {
    for(const auto& trackV2 : vertex2.originalTracks()) {
      if(TrackHelper::SameTrack(trackV1.track(), trackV2.track())) {
	nCommonTracks++;
      }
    }
  }
  return nCommonTracks;
}

reco::TransientTrack VertexHelper::GetCommonTrack(const TransientVertex &vertex1, const TransientVertex &vertex2) {

  int nCommonTracks(0);
  reco::TransientTrack commonTrack;
  for(const auto& trackV1 : vertex1.originalTracks()) {
    for(const auto& trackV2 : vertex2.originalTracks()) {
      if(TrackHelper::SameTrack(trackV1.track(), trackV2.track())) {
	commonTrack = trackV1;
	nCommonTracks++;
      }
    }
  }
  
  //if (nCommonTracks > 1)
  //throw std::runtime_error("Found " + std::to_string(nCommonTracks) + " tracks in common (max is 1)!");

  return commonTrack;
}

double VertexHelper::GetDxySeparation(const TransientVertex &vertex1, const TransientVertex &vertex2) {
  const GlobalPoint point1(vertex1.position()), point2(vertex2.position());

  // Calculate the differences in coordinates
  double deltaX = point2.x() - point1.x();
  double deltaY = point2.y() - point1.y();

  // Calculate the distance D between the points
  return std::sqrt(deltaX * deltaX + deltaY * deltaY);
}

double VertexHelper::GetDxySeparationError(const TransientVertex &vertex1, const TransientVertex &vertex2) {

  const GlobalPoint point1(vertex1.position()), point2(vertex2.position());
  const GlobalError totalError(vertex1.positionError() + vertex2.positionError());

  const double distance(GetDxySeparation(vertex1, vertex2));

  // If the distance is zero (the points are the same), return infinite uncertainty
  if(distance == 0) return std::numeric_limits<double>::min();

  // Get difference along each coordinate
  const double deltaX = point2.x() - point1.x();
  const double deltaY = point2.y() - point1.y();

  // Calculate the elements of the combined covariance matrix C_Delta = C1 + C2
  const double C_delta_xx = totalError.cxx();
  const double C_delta_yy = totalError.cyy();
  const double C_delta_xy = totalError.cyx();

  const double sigmaSquared = (deltaX * deltaX * C_delta_xx +
                               deltaY * deltaY * C_delta_yy +
                               2 * deltaX * deltaY * C_delta_xy) / (distance*distance);

  return std::sqrt(sigmaSquared);
}

double VertexHelper::GetDxySeparationSignificance(const TransientVertex &vertex1, const TransientVertex &vertex2) {
  return GetDxySeparation(vertex1, vertex2)/GetDxySeparationError(vertex1, vertex2);
}

double VertexHelper::GetDistance(const TransientVertex &vertex1, const TransientVertex &vertex2) {

  const GlobalPoint point1(vertex1.position()), point2(vertex2.position());
 
  // Calculate the differences in coordinates
  double deltaX = point2.x() - point1.x();
  double deltaY = point2.y() - point1.y();
  double deltaZ = point2.z() - point1.z();

  // Calculate the distance D between the points
  return std::sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);  
}

double VertexHelper::GetDistanceError(const TransientVertex &vertex1, const TransientVertex &vertex2) {

  const GlobalPoint point1(vertex1.position()), point2(vertex2.position());
  const GlobalError totalError(vertex1.positionError() + vertex2.positionError());

  const double distance(GetDistance(vertex1, vertex2));

  // If the distance is zero (the points are the same), return infinite uncertainty
  if(distance == 0) return std::numeric_limits<double>::min();

  // Get difference along each coordinate
  const double deltaX = point2.x() - point1.x();
  const double deltaY = point2.y() - point1.y();
  const double deltaZ = point2.z() - point1.z();

  // Calculate the elements of the combined covariance matrix C_Delta = C1 + C2
  const double C_delta_xx = totalError.cxx();
  const double C_delta_yy = totalError.cyy();
  const double C_delta_zz = totalError.czz();
  const double C_delta_xy = totalError.cyx();
  const double C_delta_xz = totalError.czx();
  const double C_delta_yz = totalError.czy();

  const double sigmaSquared = (deltaX * deltaX * C_delta_xx +
                               deltaY * deltaY * C_delta_yy +
                               deltaZ * deltaZ * C_delta_zz +
                               2 * deltaX * deltaY * C_delta_xy +
                               2 * deltaX * deltaZ * C_delta_xz +
                               2 * deltaY * deltaZ * C_delta_yz) / (distance*distance);

  return std::sqrt(sigmaSquared);
}

double VertexHelper::GetDistanceSignificance(const TransientVertex &vertex1, const TransientVertex &vertex2) {
  return GetDistance(vertex1, vertex2)/GetDistanceError(vertex1, vertex2);
}

TransientVertex VertexHelper::Merge(const TransientVertex &vertex1, const TransientVertex &vertex2) {

  reco::TransientTrack commonTrack(GetCommonTrack(vertex1, vertex2));
  std::vector<reco::TransientTrack> uniqueTracks({commonTrack});

  for(const auto& trackV1 : vertex1.originalTracks()) {
    if(!TrackHelper::SameTrack(trackV1.track(), commonTrack.track()))
      uniqueTracks.emplace_back(trackV1);
  }
  for(const auto& trackV2 : vertex2.originalTracks()) {
    if(!TrackHelper::SameTrack(trackV2.track(), commonTrack.track()))
      uniqueTracks.emplace_back(trackV2);
  }

  const KalmanVertexFitter fitter;
  
  return fitter.vertex(uniqueTracks);
}
