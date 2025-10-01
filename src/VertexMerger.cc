#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexMerger.h"

std::vector<CachingVertex<5>> VertexMerger::mergePairs(const std::vector<CachingVertex<5>>& initialVertices) {

  if(initialVertices.empty())
    return std::vector<CachingVertex<5>>();

  //std::cout << "Starting vertex merger" << std::endl;
  vertices_ = initialVertices;
  //masterList_ = vertices_;
  //std::cout << "Initial vertices: " << initialVertices.size() << std::endl;
  //int counter(0);
  while(hasSharedTracks(vertices_) && iteration_ < 30) {
    std::vector<CachingVertex<5>> newVertices;
    std::vector<bool> used(vertices_.size(), false);
    /*
    std::cout << "\nIteration: " << counter << ", has " << vertices_.size() << " vertices: " << std::endl;
    for(const auto &vertex : vertices_) {
      std::cout << "\tThis vertex has " << vertex.tracks().size() << " tracks: " << std::endl;
      for(const auto &track : TransientVertex(vertex).originalTracks()) {
        std::cout << "\t\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
      }
    }
    */
    //std::cout << "Iteration: " << counter << std::endl;
    // Compare each pair of vertices
    //if(VertexHelper::HasDuplicateVertices(vertices_))
    //std::cout << "Found duplicates in iteration: " << iteration_ << std::endl;
    
    for(size_t i = 0; i < vertices_.size(); i++) {
      if(used[i]) continue;
      
      CachingVertex<5> current_vertex = vertices_[i];
      /*
      if(TransientVertex(current_vertex).normalisedChiSquared() < 1.) {
	used[i] = true;
	continue;
      }	
      */
      // Compare with all other unused vertices
      for(size_t j = i + 1; j < vertices_.size(); j++) {
	if(used[j]) continue;

	if(VertexHelper::IsSameVertex(current_vertex, vertices_[j])) {
	  used[j] = true;
	  continue;
	}
	
	if(!hasSharedTracks(current_vertex, vertices_[j])) continue;
        
	// Try to merge
	//std::cout << "Trying to merge" << std::endl;
	CachingVertex<5> merged_vertex;
	if(tryMergeVertices(current_vertex, vertices_[j], merged_vertex)) {
	  if(hasGenVertex_ && debugMode_) {
	    for(const auto &genVertex : genVertices_) {
	      const int nSignalInVtx = genVertex.CountInstances(getOriginalTracks(merged_vertex));
	      if(nSignalInVtx == 1) 
		std::cout << "Merged vertex is bronze" << std::endl;
	      if(nSignalInVtx == 2 && merged_vertex.tracks().size() > 2)
		std::cout << "Merged vertex is silver" << std::endl;
	      if(nSignalInVtx == 2 && merged_vertex.tracks().size() == 2)
		std::cout << "Merged vertex is gold" << std::endl;
	    }
	  }
	  //std::cout << "Merged successfully" << std::endl;
	  // Successful merge
	  current_vertex = merged_vertex;
	  used[j] = true;
	} else {
	  //std::cout << "Track arbitration" << std::endl;
	  arbitrateTracks(current_vertex, vertices_[j]);

	  if(isValidVertex(vertices_[j])) {
	    used[j] = true;
	    newVertices.push_back(vertices_[j]);
	  }
	}
      }
      
      // Add the final state of current_vertex if it's still valid
      if(isValidVertex(current_vertex)) {
	newVertices.push_back(current_vertex);
      }
      used[i] = true;
    }

    // Replace old vertices with new set
    vertices_ = std::move(newVertices);
    iteration_++;
  }
  //std::cout << "Finished" << std::endl;
  return vertices_;
}

reco::VertexCollection VertexMerger::initialSeeds(const std::vector<reco::TransientTrack> &tracks, const edm::Handle<reco::TrackCollection> &trackHandle) const {
  
  int count(0);
  reco::VertexCollection seeds;

  std::vector<TransientVertex> convertedVertices;
  for(const auto &vertex : getVertexSeeds(tracks)) {
    convertedVertices.emplace_back(vertex);
  }
  
  for(const auto &vertex : convertedVertices) {
    reco::Vertex recoVtx(VertexHelper::ConvertFitVertex(vertex));
    for(const auto &track : vertex.originalTracks()) {
      reco::TrackRef trackRef = TrackHelper::GetTrackRef(track.track(), trackHandle);
      recoVtx.add(reco::TrackBaseRef(trackRef), recoVtx.trackWeight(reco::TrackBaseRef(trackRef)));
    }
    seeds.emplace_back(recoVtx);
    count++;
  }
  return seeds;
}

reco::VertexCollection VertexMerger::vertices(const std::vector<reco::TransientTrack> &tracks, const edm::Handle<reco::TrackCollection> &trackHandle) {

  if(hasGenVertex_ && debugMode_)
    std::cout << "----------------------------------Vertex Merger With " << genVertices_.size() << " Signal Vertices------------------------------------" << std::endl;

  if(tracks.empty())
    return reco::VertexCollection();
  
  const std::vector<CachingVertex<5>> seeds(getVertexSeeds(tracks));
  mergePairs(seeds);

  std::vector<TransientVertex> convertedVertices;
  for(const auto &vertex : vertices_) {
    convertedVertices.emplace_back(vertex);
  }

  reco::VertexCollection recoVertices;
  for(const auto &vertex : convertedVertices) {
    reco::Vertex recoVtx(VertexHelper::ConvertFitVertex(vertex));
    for(const auto &track : vertex.originalTracks()) {
      reco::TrackRef trackRef = TrackHelper::GetTrackRef(track.track(), trackHandle);
      recoVtx.add(reco::TrackBaseRef(trackRef), recoVtx.trackWeight(reco::TrackBaseRef(trackRef)));
    }
    recoVertices.emplace_back(recoVtx);
    if(hasGenVertex_ && debugMode_) {
      for(const auto &genVertex : genVertices_) {
	reco::TrackCollection vertexTracks;
	for(const auto &track : recoVtx.tracks())
	  vertexTracks.emplace_back(*track);
	const int nSignal = genVertex.CountInstances(vertexTracks);
	if(nSignal > 0 && debugMode_) {
	  std::cout << "Found signal tracks in a vertex" << std::endl;
	  if(genVertex.isGold(recoVtx))
	    std::cout << "\tFound a gold vertex" << std::endl;
	  if(genVertex.isSilver(recoVtx))
	    std::cout << "\tFound a silver vertex" << std::endl;
	  if(genVertex.isBronze(recoVtx))
	    std::cout << "\tFound a bronze vertex" << std::endl;
	  if(genVertex.isNoMatch(recoVtx))
	    std::cout << "Didn't find a match" << std::endl;
	}	
      }
    }
  }

  if(hasGenVertex_ && debugMode_) {
    for(const auto &vertex : rejects_) {
      for(const auto &genVertex : genVertices_) {
	const TransientVertex tvertex(vertex);
	const int nSignalInVtx = genVertex.CountInstances(getOriginalTracks(vertex));
	if(nSignalInVtx == 2 && tvertex.originalTracks().size() == 2)
	  std::cout << "Found a gold vertex in the rejects with chi2: " << tvertex.normalisedChiSquared() << std::endl;
	else if (nSignalInVtx == 2)
	  std::cout << "Found a silver vertex in the rejects with chi2: " << tvertex.normalisedChiSquared() << std::endl;
      }
    }
  }
    
  if(hasGenVertex_ && debugMode_)
    std::cout << "-----------------------------------Finished Vertex Merger With Signal-------------------------------------" << std::endl;
  
  return recoVertices;
}

void VertexMerger::setGenVertices(const GenVertices &genVertices, const bool debug) {

  genVertices_ = genVertices;
  if(!genVertices.empty())
    hasGenVertex_ = true;

  debugMode_ = debug;
}

double VertexMerger::calculateDistance(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const {

  const GlobalPoint point1(vertex1.position()), point2(vertex2.position());

  // Calculate the differences in coordinates
  double deltaX = point2.x() - point1.x();
  double deltaY = point2.y() - point1.y();
  double deltaZ = point2.z() - point1.z();

  // Calculate the distance D between the points
  return std::sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

double VertexMerger::calculateDistanceUncertainty(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const {

  const GlobalPoint point1(vertex1.position()), point2(vertex2.position());
  const GlobalError totalError(vertex1.error() + vertex2.error());

  const double distance(calculateDistance(vertex1, vertex2));

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

bool VertexMerger::isValidVertex(const CachingVertex<5>& vertex) const {
  //const bool isInMasterList(VertexHelper::IsInCollection(vertex, masterList_));
  const double normChi2(TransientVertex(vertex).normalisedChiSquared());
  return vertex.isValid() && normChi2 < 9 && vertex.tracks().size() > 1;// && !isInMasterList;
}

bool VertexMerger::hasSharedTracks(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const {

  const auto& tracks1 = vertex1.tracks();
  const auto& tracks2 = vertex2.tracks();

  bool isShared = false;
  for(const auto& track : tracks2) {
    VertexTrackEqual<5> trackComparator(track);

    for (const auto& trackInVertex1 : tracks1) {
      if (trackComparator(trackInVertex1)) {
        isShared = true;
        break;
      }
    }
    if(isShared) break;
  }
  return isShared;
}

bool VertexMerger::hasSharedTracks(const std::vector<CachingVertex<5>> &vertices) const {

  bool hasSharedTrack(false);
  std::vector<CachingVertex<5>>::const_iterator end = vertices.end();
  std::vector<CachingVertex<5>>::const_iterator endm1 = (end - 1);
  for(std::vector<CachingVertex<5>>::const_iterator x = vertices.begin(); x != endm1; ++x) {
    for(std::vector<CachingVertex<5>>::const_iterator y = x + 1; y != end; ++y) {
      if(hasSharedTracks(*x, *y)) {
        hasSharedTrack = true;
        break;
      }
    }
    if(hasSharedTrack) break;
  }
  return hasSharedTrack;
}


bool VertexMerger::tryMergeVertices(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2,
				    CachingVertex<5>& merged) {

  const double distance(calculateDistance(vertex1, vertex2));
  const double distanceUncertainty(calculateDistanceUncertainty(vertex1, vertex2));

  if(hasGenVertex_ && debugMode_) {
    for(const auto &genVertex : genVertices_) {
      const int nSignalInVtx1 = genVertex.CountInstances(getOriginalTracks(vertex1));
      const int nSignalInVtx2 = genVertex.CountInstances(getOriginalTracks(vertex2));
      if(nSignalInVtx1 == 2 || nSignalInVtx2 == 2) {

	std::cout << "\nAttempting to merge gold seed vertex: " << std::endl;
	std::cout << "Found seed vertex with signal. Distance significance = " << distance/distanceUncertainty << std::endl;
	if(distance/distanceUncertainty < 4)
	  std::cout << "  Merging criteria satisfied" << std::endl;
	else
	  std::cout << "  Merging criteria failed" << std::endl;
      }
    }
  }
  /*
  // Save merging info
  if(hasGenVertex_) {
    for(const auto &genVertex : genVertices_) {
      const int nSignalInVtx1 = genVertex.CountInstances(getOriginalTracks(vertex1));
      const int nSignalInVtx2 = genVertex.CountInstances(getOriginalTracks(vertex2));
      
      if(nSignalInVtx1 == 2 && getOriginalTracks(vertex1).size() == 2) {
	MergeInfo info(vertex1, vertex2, true);
	mergeDecisions_.emplace_back(info);
      }
      else if(nSignalInVtx2 == 2 && getOriginalTracks(vertex2).size() == 2) {
	MergeInfo info(vertex2, vertex1, true);
	mergeDecisions_.emplace_back(info);
      }
      else {
	MergeInfo info(vertex1, vertex2);
	mergeDecisions_.emplace_back(info);
      }
    }
  }
  else {
    MergeInfo info(vertex1, vertex2);
    mergeDecisions_.emplace_back(info);
  }
  */
  if(distance < 4*distanceUncertainty) {
    // Save original normChi2 for comparison with merged
    const TransientVertex v1(vertex1), v2(vertex2);
    const double stdRes1(v1.normalisedChiSquared()), stdRes2(v2.normalisedChiSquared());
    const double pValue1(calculateChiSquaredPValue(v1.totalChiSquared(), v1.degreesOfFreedom()));
    const double pValue2(calculateChiSquaredPValue(v2.totalChiSquared(), v2.degreesOfFreedom()));
    
    // Fit merged vertex
    const KalmanVertexFitter fitter;
    const CachingVertex<5> newVertex = fitter.vertex(getAllUniqueTracks(vertex1, vertex2));
    const TransientVertex vNew(newVertex);
    const double newNormChi2 = TransientVertex(newVertex).normalisedChiSquared();
    const double pValueNew(isValidVertex(newVertex)? calculateChiSquaredPValue(vNew.totalChiSquared(), vNew.degreesOfFreedom()) : 999.);

    /*
    if(VertexHelper::CountSharedTracks(vertex1, vertex2) > 1) {
      std::cout << "\nFound two vertices that share " << VertexHelper::CountSharedTracks(vertex1, vertex2) << " tracks pre-merging in iteration " << iteration_ << std::endl;
      std::cout << "vertex 1 has " << vertex1.tracks().size() << " tracks" << std::endl;
        for(const auto &track : TransientVertex(vertex1).originalTracks())
          std::cout << "\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
        std::cout << "vertex 2 has " << vertex2.tracks().size() << " tracks" << std::endl;
        for(const auto &track : TransientVertex(vertex2).originalTracks())
          std::cout << "\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;

	if(vertex1.tracks().size() == 2 && vertex2.tracks().size() == 2) {
	  std::cout << "full list: " << std::endl;
	  for(const auto &vertex : vertices_) {
	    std::cout << "\tThis vertex has " << vertex.tracks().size() << " tracks: " << std::endl;
	    for(const auto &track : TransientVertex(vertex).originalTracks()) {
	      std::cout << "\t\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
	    }
	  }
	  }
    }
    */
    /*
    if(VertexHelper::CountSharedTracks(vertex1, newVertex) > 1 && isValidVertex(newVertex)) {
      std::cout	<< "Found " << VertexHelper::CountSharedTracks(vertex1, newVertex) << " tracks shared between vertex1 and the merged vertex" << std::endl;
    }
    if(VertexHelper::CountSharedTracks(vertex2, newVertex) > 1 && isValidVertex(newVertex)) {
      std::cout << "Found " << VertexHelper::CountSharedTracks(vertex2, newVertex) << " tracks shared between vertex2 and the merged vertex" << std::endl;
    }
    */
    if(isValidVertex(newVertex)) {// && newNormChi2 < stdRes1 && newNormChi2 < stdRes2) {
      /*
      const int trackSize1(vertex1.tracks().size()), trackSize2(vertex2.tracks().size()), newTrackSize(newVertex.tracks().size());
      if(newTrackSize != trackSize1 + trackSize2 - 1) {
	std::cout << "\nvertex 1 has " << trackSize1 << " tracks" << std::endl;
	for(const auto &track : TransientVertex(vertex1).originalTracks())
	  std::cout << "\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl; 
	std::cout << "vertex 2 has " << trackSize2 << " tracks" << std::endl;
	for(const auto &track : TransientVertex(vertex2).originalTracks())
          std::cout << "\ttrack pt = " << track.track().pt() <<	", eta = " << track.track().eta() << std::endl;
	std::cout << "the new vertex has " << newTrackSize << " tracks" << std::endl;
	for(const auto &track : TransientVertex(newVertex).originalTracks())
          std::cout << "\ttrack pt = " << track.track().pt() <<	", eta = " << track.track().eta() << std::endl;
      }
      */
      if(hasGenVertex_ && debugMode_) {
	for(const auto &genVertex : genVertices_) {
	  const int nSignalInVtx1 = genVertex.CountInstances(getOriginalTracks(vertex1));
	  const int nSignalInVtx2 = genVertex.CountInstances(getOriginalTracks(vertex2));
	  const int nSignalInMerged = genVertex.CountInstances(getOriginalTracks(newVertex));

	  if(nSignalInVtx1 > 0 || nSignalInVtx2 > 0)
	    std::cout << "\nVertex with signal merged" << std::endl;
	  if(nSignalInVtx1 > 0)
	    std::cout << "\tvertex 1 has " << nSignalInVtx1 << " signal tracks and a p-value = " << pValue1
		      << " (normChi2 = " << stdRes1 << ")"<< std::endl;
	  if(nSignalInVtx2 > 0)	
            std::cout << "\tvertex 2 has " << nSignalInVtx2 << " signal tracks and a p-value = " << pValue2
		      << " (normChi2 = " << stdRes2 << ")"<< std::endl;
	  if(nSignalInMerged > 0)
            std::cout << "\tmerged vertex has " << nSignalInMerged << " signal tracks and a p-value = "
		      << pValueNew << " (normChi2 = " << newNormChi2 << ")" <<std::endl;
	}
      }
      merged = newVertex;
      //masterList_.emplace_back(newVertex);
      return true;
    }
  }
  return false;
}

void VertexMerger::arbitrateTracks(CachingVertex<5> &vertex1, CachingVertex<5> &vertex2) {

  const CachingVertex<5> vertex1Copy(vertex1), vertex2Copy(vertex2);
  
  const KalmanVertexTrackCompatibilityEstimator<5> estimator;
  const RefCountedVertexTrack shared(sharedTrack(vertex1, vertex2));
  const auto compatibilityPair1(estimator.estimate(vertex1, shared));
  const auto compatibilityPair2(estimator.estimate(vertex2, shared));
  
  if(!compatibilityPair1.first || !compatibilityPair2.first) 
    throw std::runtime_error("Track compatibility estimation failed");

  else if(compatibilityPair1.first && compatibilityPair2.first) {

    const int size1(vertex1.tracks().size()), size2(vertex2.tracks().size());
    const double sigma1(std::sqrt(compatibilityPair1.second));
    const double sigma2(std::sqrt(compatibilityPair2.second));

    if(hasGenVertex_) {

      const TransientVertex v1(vertex1), v2(vertex2);
      const double stdRes1(v1.normalisedChiSquared()), stdRes2(v2.normalisedChiSquared());

      for(const auto &genVertex : genVertices_) {
	const int nSignalInVtx1 = genVertex.CountInstances(getOriginalTracks(vertex1));
	const int nSignalInVtx2 = genVertex.CountInstances(getOriginalTracks(vertex2));

	if(debugMode_ && (nSignalInVtx1 == 2 || nSignalInVtx2 == 2))  
	  std::cout << "\nTrack arbitration with gold seed vertex: " << std::endl;
	
	if(nSignalInVtx1 == 2) {
	  if(debugMode_) {
	    std::cout << "Vertex 1 (normChi2 = " << stdRes1 << ") is signal seed vertex and has size "
		      << size1 << " and arbitrated track has a sigma = " << sigma1 << std::endl;
	    std::cout << "  The compared vertex (normChi2 = " << stdRes2 << ", cosTheta = " << cosThetaWithPV(vertex2)
		      << " and mass = " << getMass(vertex2) << ")  has " <<  nSignalInVtx2
		      << " signal, and " << size2 << " tracks with a sigma = " << sigma2 << std::endl;
	    std::cout << "gold cosTheta = " << cosThetaWithPV(vertex1) << std::endl;
	  }
	  ArbitrationInfo info(v1, v2);
	  if(hasPV_) info.setTrackCosTheta(primaryVertex_);
	  arbitrationDecisions_.emplace_back(info);
	}

	if(nSignalInVtx2 == 2) {
	  if(debugMode_) {
	    std::cout << "Vertex 2 (normChi2 = " << stdRes2 << ") is signal seed vertex and has size "
		      << size2 << " and arbitrated track has a sigma = " << sigma2 << std::endl;
	    std::cout << "  The compared vertex (normChi2 = " << stdRes1 << ", cosTheta = " << cosThetaWithPV(vertex1)
		      << " and mass = " << getMass(vertex1) << ") has " <<  nSignalInVtx1
		      << " signal, and " << size1 << " tracks with a sigma = " << sigma1 << std::endl;
	    std::cout << "gold cosTheta = " << cosThetaWithPV(vertex2) << std::endl;
	  }
	  ArbitrationInfo info(v2, v1);
	  if(hasPV_) info.setTrackCosTheta(primaryVertex_);
          arbitrationDecisions_.emplace_back(info);
	}
      }
    }

    // If the track is within 1.5 standard deviations of both vertices assign to larger vertex
    if(sigma1 < 1.5 && sigma2 < 1.5 && size1 != size2) {
      if(size1 > size2) 
	vertex2 = refitTracks(allTracksExceptGiven(vertex2, shared));
      else
	vertex1 = refitTracks(allTracksExceptGiven(vertex1, shared));
    }
    
    // If the track is more than 5 standard deviations from either vertex remove it
    else if(sigma1 > 5) 
      vertex1 = refitTracks(allTracksExceptGiven(vertex1, shared));
    else if(sigma2 > 5)	
      vertex2 = refitTracks(allTracksExceptGiven(vertex2, shared));

    // otherwise, the shared track is assigned to the vertex with the smallest IP significance
    else {
      if(sigma1 < sigma2)
	vertex2 = refitTracks(allTracksExceptGiven(vertex2, shared));
      else
        vertex1	= refitTracks(allTracksExceptGiven(vertex1, shared));
    }

    // Save vertices rejected by arbitration
    if(!isValidVertex(vertex1))
      rejects_.emplace_back(vertex1Copy);
    if(!isValidVertex(vertex2))
      rejects_.emplace_back(vertex2Copy);
    
  }
  
  if(VertexHelper::CountSharedTracks(vertex1, vertex2) > 1) {
    std::cout << "\n---------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "The two vertices share " << VertexHelper::CountSharedTracks(vertex1, vertex2) << " tracks" << std::endl;
    std::cout << "The original vertices: vertex1 has " << vertex1Copy.tracks().size() << " and vertex 2 has " << vertex2Copy.tracks().size() << " tracks" << std::endl;
    std::cout << "\tvertex1: " << std::endl;
    for(const auto &track : TransientVertex(vertex1Copy).originalTracks())
      std::cout << "\t\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
    std::cout << "\tvertex2: " << std::endl;
    for(const auto &track : TransientVertex(vertex2Copy).originalTracks())
      std::cout << "\t\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
    std::cout << "\nThe arbitrated vertices: vertex1 has " << vertex1.tracks().size() << " and vertex 2 has " << vertex2.tracks().size() << " tracks" << std::endl;
    std::cout << "\tvertex1: " << std::endl;
    for(const auto &track : TransientVertex(vertex1).originalTracks())
      std::cout << "\t\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
    std::cout << "\tvertex2: " << std::endl;
    for(const auto &track : TransientVertex(vertex2).originalTracks())
      std::cout << "\t\ttrack pt = " << track.track().pt() << ", eta = " << track.track().eta() << std::endl;
     std::cout << "---------------------------------------------------------------------------------------------------------" << std::endl;
  }
}

RefCountedVertexTrack VertexMerger::sharedTrack(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const {

  RefCountedVertexTrack sharedTrack;
  const auto& tracks1 = vertex1.tracks();
  const auto& tracks2 = vertex2.tracks();

  int nCommonTracks(0);
  
  bool isShared = false;
  for(const auto& track : tracks2) {
    VertexTrackEqual<5> trackComparator(track);

    for (const auto& trackInVertex1 : tracks1) {
      if(trackComparator(trackInVertex1)) {
        sharedTrack = track;
	nCommonTracks++;
	isShared = true;
	break;
      }
    }
    if(isShared) break;
  }
  
  if(nCommonTracks > 1)
    std::cout << "Found " << nCommonTracks << " tracks in common" << std::endl;

  return sharedTrack;
}

std::vector<RefCountedVertexTrack> VertexMerger::getAllUniqueTracks(const CachingVertex<5>& vertex1, const CachingVertex<5>& vertex2) const {

  std::vector<RefCountedVertexTrack> uniqueTracks;
  bool hasSharedTrack = false;
  
  const auto& tracks1 = vertex1.tracks();
  const auto& tracks2 = vertex2.tracks();
  
  for (const auto& track : tracks1) {
    uniqueTracks.push_back(track);
  }
  
  // Use VertexTrackEqual to check for shared tracks in vertex2
  for (const auto& track : tracks2) {
    VertexTrackEqual<5> trackComparator(track);
    
    bool isShared = false;
    for (const auto& trackInVertex1 : tracks1) {
      if (trackComparator(trackInVertex1)) {
	hasSharedTrack = true;
	isShared = true;
	break;
      }
    }
    
    // If track is unique to vertex2, add it to uniqueTracks
    if (!isShared) {
      uniqueTracks.push_back(track);
    }
  }
  
  // Return all unique tracks if there's a shared track, otherwise return an empty vector
  return hasSharedTrack ? uniqueTracks : std::vector<RefCountedVertexTrack>();
}

std::vector<RefCountedVertexTrack> VertexMerger::allTracksExceptGiven(const CachingVertex<5>& vertex, const RefCountedVertexTrack &toRemove) const {

  const VertexTrackEqual<5> trackComparator(toRemove);
  std::vector<RefCountedVertexTrack> filteredTracks;

  for(const auto& track : vertex.tracks())
    if(!trackComparator(track))
      filteredTracks.emplace_back(track);

  return filteredTracks;
}

CachingVertex<5> VertexMerger::refitTracks(const std::vector<RefCountedVertexTrack> &tracks) const {

  const KalmanVertexFitter fitter;

  if(tracks.size() > 1) 
    return fitter.vertex(tracks);

  return CachingVertex<5>();
}

std::vector<CachingVertex<5>> VertexMerger::getVertexSeeds(const std::vector<reco::TransientTrack> &tracks) const {

  //std::cout << "Total initial tracks: " << tracks.size() << std::endl;
  //std::cout << "Getting vertex seeds" <<std::endl;  
  std::vector<CachingVertex<5>> seeds;
  KalmanVertexFitter fitter;
  std::vector<reco::TransientTrack>::const_iterator end = tracks.end();
  std::vector<reco::TransientTrack>::const_iterator endm1 = (end - 1);
  for (std::vector<reco::TransientTrack>::const_iterator x = tracks.begin(); x != endm1; ++x) {
    for (std::vector<reco::TransientTrack>::const_iterator y = x + 1; y != end; ++y) {
      CachingVertex<5> tempVertex(fitter.vertex({*x, *y}));
      
      const double chi2(tempVertex.totalChiSquared());
      const double ndof(tempVertex.degreesOfFreedom());

      if(hasGenVertex_ && debugMode_) {
	TransientVertex vtx(tempVertex);
	//const double x(vtx.position().x()), y(vtx.position().y());//, z(vtx.position().z());
	for(const auto &genVertex : genVertices_) {
	  const reco::TrackCollection signalTracks(genVertex.tracks());
	  reco::TrackCollection vtxTracks;
	  for(const auto &track : vtx.originalTracks())
	    vtxTracks.emplace_back(track.track());
	  if(TrackHelper::AreTrackCollectionsIdentical(signalTracks, vtxTracks)) {
	    std::cout << "Found the seed with the signal tracks (normChi2 = " << chi2/ndof << ", dxy = " << genVertex.dxy()
		      << ", cosTheta = " << cosThetaWithPV(tempVertex) << " and mass = " << getMass(tempVertex) << ")" << std::endl;
	    for(const auto &pair : genVertex.genMatches())
	      std::cout << "  track gen deltaR: " << pair.GetDeltaR() << std::endl;
	    //std::cout << "Vertex displacement dxy: " << sqrt(x*x + y*y) << std::endl;
	  }
	}
      }
    
      //if(tempVertex.isValid() && chi2/ndof < 9)
      if(isValidVertex(tempVertex) && cosThetaWithPV(tempVertex) > 0.5) {
	if(tempVertex.tracks().size() < 2)
	  std::cout << "A vertex with " << tempVertex.tracks().size() << " shouldn't be in here" << std::endl;
	seeds.emplace_back(tempVertex);
      }
    }
  }
  return seeds;
}

reco::TrackCollection VertexMerger::getOriginalTracks(const CachingVertex<5> &vertex) const {

  const TransientVertex convertedVertex(vertex);
  
  reco::TrackCollection ogTracks;
  for(const auto &track : convertedVertex.originalTracks())
    ogTracks.emplace_back(track.track());

  return ogTracks;
}

std::vector<CachingVertex<5>> VertexMerger::clean(const std::vector<CachingVertex<5>> &vertices) const {

  if(!VertexHelper::HasDuplicateVertices(vertices))
    return vertices;

  std::vector<CachingVertex<5>> cleanVertices;
  for(const auto &vertex : vertices) {
    if(!VertexHelper::IsInCollection(vertex, cleanVertices)) {
      std::cout << "Not in collection" << std::endl;
      cleanVertices.emplace_back(vertex);
    }
    //else std::cout << "In collection" << std::endl;
  }
  if(VertexHelper::HasDuplicateVertices(cleanVertices))
    std::cout << "Vertices not clean" << std::endl;
  return cleanVertices;
}

double VertexMerger::calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const {
    // Define the chi-squared distribution with the specified degrees of freedom
    boost::math::chi_squared chiSqDist(degreesOfFreedom);

    // Compute the p-value as 1 - CDF(chiSquaredValue)
    double pValue = 1.0 - boost::math::cdf(chiSqDist, chiSquaredValue);

    return pValue;
}

double VertexMerger::cosThetaWithPV(const TransientVertex &vertex) const {

  if(!hasPV_) return 1.;
  
  ROOT::Math::PtEtaPhiMVector lorentzVec;
  
  const double pvX(primaryVertex_.x()), pvY(primaryVertex_.y()), pvZ(primaryVertex_.z());
  const	double vX(vertex.position().x()), vY(vertex.position().y()), vZ(vertex.position().z());
  const double xMag( sqrt((pvX-vX)*(pvX-vX) + (pvY-vY)*(pvY-vY) + (pvZ-vZ)*(pvZ-vZ)));
  
  for(const auto &ttrack : vertex.originalTracks()) {
    const reco::Track track(ttrack.track());
    const double pt(track.pt()), eta(track.eta()), phi(track.phi()), mass(0.);
    lorentzVec += ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
  }
  
  const double vPx(lorentzVec.px()), vPy(lorentzVec.py()), vPz(lorentzVec.pz());
  const double pMagVertex(sqrt(vPx*vPx + vPy*vPy + vPz*vPz));
  const double cosTheta( ((vX-pvX)*vPx + (vY-pvY)*vPy + (vZ-pvZ)*vPz) / (xMag*pMagVertex));
  
  return cosTheta;
}

double VertexMerger::getMass(const TransientVertex &vertex) const {

  ROOT::Math::PtEtaPhiMVector lorentzVec;

  for(const auto &ttrack : vertex.originalTracks()) {
    const reco::Track track(ttrack.track());
    const double pt(track.pt()), eta(track.eta()), phi(track.phi()), mass(0.);
    lorentzVec += ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
  }

  return lorentzVec.M();
}
