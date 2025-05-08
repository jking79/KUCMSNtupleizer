#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/NewVertexMerger.h"

NewVertexMerger::NewVertexMerger(const reco::Vertex &primaryVertex, const double seedCosThetaCut) :
  primaryVertex_(primaryVertex), seedCosThetaCut_(seedCosThetaCut) {}

reco::VertexCollection NewVertexMerger::seeds(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) const {
  return createSeeds(tracks, ttBuilder).vertices();
}

reco::VertexCollection NewVertexMerger::vertices(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) {
  return mergePairs(createSeeds(tracks, ttBuilder)).vertices();
}

TrackVertexSetCollection NewVertexMerger::trackVertexSets(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) {
  return mergePairs(createSeeds(tracks, ttBuilder));
}

TrackVertexSetCollection NewVertexMerger::trackVertexSetSeeds(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) {
  return createSeeds(tracks, ttBuilder);
}

TrackVertexSetCollection NewVertexMerger::createSeeds(const std::vector<reco::TrackRef> &tracks, const TransientTrackBuilder* ttBuilder) const {

  TrackVertexSetCollection setCollection;
  
  if(tracks.size() < 2)
    return setCollection;

  std::vector<reco::TrackRef>::const_iterator end = tracks.end();
  std::vector<reco::TrackRef>::const_iterator endm1 = (end - 1);
  for (std::vector<reco::TrackRef>::const_iterator x = tracks.begin(); x != endm1; ++x) {
    for (std::vector<reco::TrackRef>::const_iterator y = x + 1; y != end; ++y) {
      if(TrackHelper::OverlappingTrack(**x, **y, ttBuilder)) continue;
      TrackVertexSet newSet({*x, *y}, ttBuilder);
      if(newSet.isValid() && newSet.normChi2() < 5) {
        setCollection.insert(newSet);
      }
    }
  }

  filterWithAlignment(setCollection);

  return setCollection;
}

/*
TrackVertexSetCollection NewVertexMerger::mergePairs(const TrackVertexSetCollection &seeds) {
  
  if(seeds.empty())
    return seeds;

  if(DEBUG)
    std::cout << "\n----------------------------------Vertex Merger With " << seeds.size() << " initial seeds------------------------------------" << std::endl;
  
  vertices_ = seeds;
  masterList_ =	vertices_;
  TrackVertexSetCollection ignoreList, acceptedList;

  int iteration(0);
  bool hasMerged(true);
  while(hasMerged && iteration < 30) {
    hasMerged = false;
    if(DEBUG)
      std::cout << "\niteration: " << iteration << std::endl;
    for (auto it1 = vertices_.begin(); it1 != vertices_.end(); ++it1) {
      
      if(ignoreList.contains(*it1)) continue;
      TrackVertexSet vertex1(*it1);
      
      for (auto it2 = std::next(it1); it2 != vertices_.end(); ++it2) {
	
	if(ignoreList.contains(*it2) || ignoreList.contains(vertex1)) continue;
	TrackVertexSet vertex2(*it2);
	
	// If the two vertices don't even have one track in common skip to next vertex
	if((vertex1 & vertex2) < 1)
	  continue;
	
	// Attempt to merge the two vertices that share tracks if they are "close" enough
	else if(vertex1.distanceSignificance(vertex2) < 4) {
	  if(DEBUG) {
	    std::cout << "\nMerge criteria satisfied. Attempting to merge the following two vertices: " << std::endl;
	    vertex1.printTrackInfo();
	    vertex2.printTrackInfo();
	  }
	  
	  TrackVertexSet mergedVertex(vertex1 + vertex2);
	  
	  // If the vertex is valid and passes the Chi2 cut replace vertex1 with the merged
	  // vertex, remove vertex2, add merged vertex to the master list and go to next one in line
	  if(isValidVertex(mergedVertex) && masterList_.doesNotContain(mergedVertex)) {
	    if(DEBUG) {
	      std::cout << "Merge successful! New vertex info: " << std::endl; 
	      mergedVertex.printTrackInfo();
	    }
	    
	    ignoreList.add(vertex1);
	    ignoreList.add(vertex2);
	    masterList_.add(mergedVertex);
	    acceptedList.add(mergedVertex);
	    hasMerged = true;
	    continue;
	  }
	  else if (masterList_.contains(mergedVertex)) {
	    if(DEBUG)
	      std::cout << "Merged vertex already exists." << std::endl;
	    ignoreList.add(vertex1);
	    ignoreList.add(vertex2);
	    continue;
	  }
	  else if(DEBUG) {
	    if(mergedVertex.isValid() && masterList_.doesNotContain(mergedVertex))
	      std::cout << "Merge failed due to normChi2 = " << mergedVertex.normChi2() << std::endl;
	    else std::cout << "Merge is invalid" << std::endl;
	  }
	}

	else {
	  std::cout << "Merge criteria failed. Commencing track arbitration: " << std::endl;
	  
	  std::cout << "\nThe vertices have " << (vertex1 & vertex2) << " tracks in common" << std::endl;
	  
	  for(const auto &track : vertex1.commonTracks(vertex2)) {
	    
	    std::cout << "\tCurrent track in track arbitration: " << std::endl;
	    std::cout << "\t\tpt: " << track->pt() << std::endl;
	    std::cout << "\t\teta: " << track->eta() << std::endl;
	    std::cout << "\t\tphi: " << track->phi() << std::endl;
	    
	    const int size1(vertex1.size()), size2(vertex2.size());
	    const double sigma1(vertex1.compatibility(track)), sigma2(vertex2.compatibility(track));
	    
	    std::cout << "\tcompatibility with vertex1 (" << "size = " << size1 << "): " << sigma1 <<std::endl;
	    std::cout << "\tcompatibility with vertex2 (" << "size = " << size2 << "): " << sigma2 <<std::endl;
	    
	    // If the track is compatible with both vertices choose the larger one
	    if(sigma1 < 1.5 && sigma2 < 1.5 && size1 != size2) {
	      
	      std::cout << "\tThe track is compatible with both vertices!" << std::endl;
	      if(vertex1.isLargerThan(vertex2)) {
		std::cout << "\t\tVertex1 is larger and gets to keep the track. Vertex2 is discarded." << std::endl;
		acceptedList.add(vertex1);
		ignoreList.add(vertex2);
		continue;
	      }
	      else {
		std::cout << "\t\tVertex2 is larger and gets to keep the track. Vertex1 is discarded." << std::endl;
		acceptedList.add(vertex2);
		ignoreList.add(vertex1);
		continue;
	      }
	    }
	    
	    // If the track is more than 5 standard deviations from either vertex remove it 
	    else if(sigma1 > 5) {
	      std::cout << "\tTrack is too incompatible (sigma1 = " << sigma1 << ") with vertex1 attempting to remove it and refit:" << std::endl;
	      vertex1.removeTrack(track);
	      if(isValidVertex(vertex1) && masterList_.doesNotContain(vertex1)) {
		std::cout << "\t\tVertex1 is still valid. Keeping it to try again." <<std::endl;
		masterList_.add(vertex1);
		acceptedList.add(vertex1);
	      }
	      else
		ignoreList.add(vertex1);
	      
	    }
	    else if(sigma2 > 5) {
	      std::cout << "\tTrack is too incompatible (sigma2 = " << sigma1 << ") with vertex2 attempting to remove it and refit:" << std::endl;
	      vertex2.removeTrack(track);
	      if(isValidVertex(vertex2) && masterList_.doesNotContain(vertex1)) {
		std::cout << "\t\tVertex1 is still valid. Keeping it to try again." << std::endl;
		masterList_.add(vertex2);
		acceptedList.add(vertex2);
	      }
	      else
		ignoreList.add(vertex2);
	      
	    }
	    else {
	      std::cout << "\t Since they are the same size, choosing the one with better compatibility: " << std::endl;
	      if(sigma1 < sigma2) {
		acceptedList.add(vertex1);
		ignoreList.add(vertex2);
		std::cout << "\t\tvertex1 is better with sigma = " << sigma1 << std::endl;
	      }
	      else {
		acceptedList.add(vertex2);
		ignoreList.add(vertex1);
		std::cout	<< "\t\tvertex2 is better with sigma = " << sigma2 << std::endl;
	      }
	    }
	  }
	}
      }
    }

    vertices_ = acceptedList;
    ignoreList.clear();
    //if(vertices_.hasExclusiveVertices())
    if(hasMerged && DEBUG)
      std::cout << "Found new mergers" << std::endl;
      //std::cout << "All vertices are unique and exclusive." << std::endl;
    else
      if(DEBUG)
	std::cout << "No new mergers found" << std::endl;
    iteration++;
  }
  if(vertices_.hasExclusiveVertices() && DEBUG)
    std::cout << "All vertices are unique and exclusive." << std::endl;
  if(hasMerged && DEBUG)
    std::cout << "Shouldn't have finished" << std::endl;
  if(DEBUG) {
    std::cout << "Finished in " << iteration+1 << " iterations." << std::endl; 
    std::cout << "----------------------------------------Finished Vertex Merger-------------------------------------------" << std::endl;
  }
  return vertices_;
}
*/

TrackVertexSetCollection NewVertexMerger::mergePairs(const TrackVertexSetCollection &seeds) {

  if(seeds.empty())
    return seeds;

  if(debug)
    std::cout << "\n----------------------------------Vertex Merger With " << seeds.size() << " initial seeds------------------------------------" << std::endl;
  
  vertices_ = seeds;
  masterList_ = vertices_;
  TrackVertexSetCollection ignoreList, mergedList;
  
  int iteration(0);
  bool hasMerged(true);

  while(hasMerged && iteration < 30) {
    hasMerged = false;
    if(debug)
      std::cout << "\niteration: " << iteration << std::endl;
    for (auto it1 = vertices_.begin(); it1 != vertices_.end(); ++it1) {

      if(ignoreList.contains(*it1)) continue;
      TrackVertexSet vertex1(*it1);

      for (auto it2 = std::next(it1); it2 != vertices_.end(); ++it2) {

        if(ignoreList.contains(*it2) || ignoreList.contains(vertex1)) continue;
        TrackVertexSet vertex2(*it2);

	if((vertex1 & vertex2) == 0)
          continue;

	mergeDecisions_.push_back(vertex1.distanceSignificance(vertex2));

	for(const auto &genVertex : genVertices_) {

	  if(genVertex.isGold(vertex1) || genVertex.isGold(vertex2))
	    goldMergeDecisions_.push_back(vertex1.distanceSignificance(vertex2));
	}
	
	if(vertex1.distanceSignificance(vertex2) < 4) {
          if(debug) {
            std::cout << "\nMerge criteria satisfied. Attempting to merge the following two vertices: " << std::endl;
            vertex1.printTrackInfo();
            vertex2.printTrackInfo();
          } 
          TrackVertexSet mergedVertex(vertex1 + vertex2);

	  if(isValidVertex(mergedVertex) && masterList_.doesNotContain(mergedVertex)) {

	    if(debug) {
              std::cout << "Merge successful! New vertex info: " << std::endl;
              mergedVertex.printTrackInfo();
            }

            ignoreList.insert(vertex1);
            ignoreList.insert(vertex2);
            masterList_.insert(mergedVertex);
            mergedList.add(mergedVertex);
            hasMerged = true;
	    continue;
          }
	  
	  else if (masterList_.contains(mergedVertex)) {
            if(debug)
              std::cout << "Merged vertex already exists." << std::endl;
            ignoreList.insert(vertex1);
            ignoreList.insert(vertex2);
	    continue;
	  }
	  else
	    if(debug)
	      std::cout << "Merge unsuccessful." << std::endl;
	}
      }
    }
    
    vertices_ -= ignoreList;
    vertices_ += mergedList;
    
    if(hasMerged && debug)
      std::cout << "Found new mergers" << std::endl;
    else if(debug)
        std::cout << "No new mergers found" << std::endl;
    iteration++;
  }

  if(vertices_.hasExclusiveVertices() && debug)
    std::cout << "All vertices are unique and exclusive." << std::endl;
  if(hasMerged && debug)
    std::cout << "Shouldn't have finished" << std::endl;
  if(debug) {
    std::cout << "Finished in " << iteration << " iterations." << std::endl;
    std::cout << "----------------------------------------Finished Vertex Merger-------------------------------------------" << std::endl;
  }

  return vertices_;
}

bool NewVertexMerger::isValidVertex(const TrackVertexSet &set) const {
  return set.isValid() && set.normChi2() < 5;
}

void NewVertexMerger::filterWithAlignment(TrackVertexSetCollection &collection) const {
    if (!primaryVertex_.isValid())
        return;

    for (auto it = collection.begin(); it != collection.end(); ) {
        if (it->cosTheta(primaryVertex_) < seedCosThetaCut_) {
            it = collection.erase(it); // Erase and get the next valid iterator
        } else {
            ++it; // Move to the next element
        }
    }
}
