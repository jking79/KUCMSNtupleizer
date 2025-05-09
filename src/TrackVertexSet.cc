#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSet.h"

TrackVertexSet::TrackVertexSet(const std::vector<reco::TrackRef> &init, const TransientTrackBuilder* ttBuilder) :
  std::set<reco::TrackRef>(init.begin(), init.end()),
  ttBuilder_(ttBuilder),
  fitter_(std::make_unique<KalmanVertexFitter>()) {
  fit();
}

TrackVertexSet::TrackVertexSet(std::initializer_list<reco::TrackRef> init, const TransientTrackBuilder* ttBuilder) :
  std::set<reco::TrackRef>(init),
  ttBuilder_(ttBuilder),
  fitter_(std::make_unique<KalmanVertexFitter>()) {
  fit();
}

TrackVertexSet::TrackVertexSet(const TrackVertexSet& other) :
  std::set<reco::TrackRef>(other),
  ttBuilder_(other.ttBuilder_),
  fitter_(std::make_unique<KalmanVertexFitter>())
{fit();}

// Transverse distance between two vertices
double TrackVertexSet::dxySeparation(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparation(this->vertex_, other.vertex_);
}

// Error in transverse distance between two vertices
double TrackVertexSet::dxySeparationError(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparation(this->vertex_, other.vertex_);
}

// Significance in transverse distance between two vertices
double TrackVertexSet::dxySeparationSignificance(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparationSignificance(this->vertex_, other.vertex_);
}

// Distance between two vertices
double TrackVertexSet::distance(const TrackVertexSet& other) const {
  return VertexHelper::GetDistance(this->vertex_, other.vertex_);
}

// Error in the distance between two vertices
double TrackVertexSet::distanceError(const TrackVertexSet& other) const {
  return VertexHelper::GetDistanceError(this->vertex_, other.vertex_);
}

// Significance of the distance between two vertices: distance/distanceError
double TrackVertexSet::distanceSignificance(const TrackVertexSet& other) const {
  return VertexHelper::GetDistanceSignificance(this->vertex_, other.vertex_);
}

// Uses KalmanVertexTrackCompatibilityEstimator to check compatibility between this vertex and one of its tracks
double TrackVertexSet::compatibility(const reco::TrackRef &track) const {

  if(!contains(track)) {
    std::cout << "Warning in TrackVertexSet::compatibility: Given track is not in vertex!" << std::endl;
    return 999.;
  }
  const KalmanVertexTrackCompatibilityEstimator<5> estimator;
  auto result(estimator.estimate(*this, ttBuilder_->build(*track)));
  
  if(!result.first)
    throw std::runtime_error("Track compatibility estimation failed");
    
  return std::sqrt(result.second);
}

double TrackVertexSet::shiftDzAfterTrackRemoval(const reco::TrackRef &track) const {
  
  TrackVertexSet thisSetMinusTrack = *this;
  thisSetMinusTrack.removeTrack(track);

  return thisSetMinusTrack.isValid()? fabs(this->position().z() - thisSetMinusTrack.position().z()) : -1.;
}

double TrackVertexSet::shift3DAfterTrackRemoval(const reco::TrackRef &track) const {

  TrackVertexSet thisSetMinusTrack = *this;
  thisSetMinusTrack.removeTrack(track);

  return thisSetMinusTrack.isValid()? distance(thisSetMinusTrack) : -1.;
}

std::vector<reco::Track> TrackVertexSet::trackList() const {
  std::vector<reco::TrackRef> returnTracks(this->begin(), this->end());
  std::vector<reco::Track> tracks;
  for(const auto& trackRef : returnTracks)
    tracks.emplace_back(*trackRef);
  return tracks;
}

std::vector<reco::TrackRef> TrackVertexSet::tracks() const {

  std::vector<reco::TrackRef> returnTracks(this->begin(), this->end());
  std::sort(returnTracks.begin(), returnTracks.end(), 
    [](const reco::TrackRef& a, const reco::TrackRef& b) {
      return a->pt() < b->pt();
    });
  return returnTracks;
}

// Function to return common tracks between two sets
std::vector<reco::TrackRef> TrackVertexSet::commonTracks(const TrackVertexSet& other) const {
  std::vector<reco::TrackRef> common;
  for (const auto& track : other) {
    if (this->contains(track)) {
      common.emplace_back(track);
    }
  }
  return common;
}

void TrackVertexSet::printTrackInfo() const {

  int count(0);
  std::cout << "\nThis set (normChi2 = " << normChi2() << ") has " << this->size() << " tracks: " << std::endl;
  for (const auto& track : *this) {
    std::cout << "\ttrack " << count << std::endl;
    std::cout << "\t\tpt: " << track->pt() << std::endl;
    std::cout << "\t\teta: " << track->eta() << std::endl;
    std::cout << "\t\tphi: " << track->phi() << std::endl;
    std::cout << "\t\tcompatibility: " << this->compatibility(track) << std::endl;
    count++;
  }
}

// Merge two TrackVertexSets together
TrackVertexSet TrackVertexSet::merge(const TrackVertexSet &other) const {
  return *this + other;
}

// Add a track to the set
void TrackVertexSet::addTrack(const reco::TrackRef& track) {
  this->insert(track);
  fit();
}

// Check if a track exists in the set
bool TrackVertexSet::contains(const reco::TrackRef& track) const {
  return this->find(track) != this->end();
}

// Remove a track from the set
void TrackVertexSet::removeTrack(const reco::TrackRef& track) {
  this->erase(track);
  if(this->size() > 1)
    fit();
  else {
    this->clear();
    vertex_ = TransientVertex();
  }
}

// Clear all tracks from the set
void TrackVertexSet::clearTracks() {
  this->clear();
}

TrackVertexSet& TrackVertexSet::operator=(const TrackVertexSet& other) {
  if (this != &other) {
    std::set<reco::TrackRef>::operator=(other);
    ttBuilder_ = other.ttBuilder_;
    fitter_ = std::make_unique<KalmanVertexFitter>(*other.fitter_);
    fit();
  }
  return *this;
}

// Overloaded operator to combine two sets into one
TrackVertexSet TrackVertexSet::operator+(const TrackVertexSet& other) const {
  TrackVertexSet result = *this; // Start with the current set
  for(const auto& track : other) {
    result.addTrack(track);
  }

  return result;
}

// Overloaded operator to count overlapping tracks between two sets
size_t TrackVertexSet::operator&(const TrackVertexSet& other) const {
  size_t count = 0;
  for (const auto& track : other) {
    if (this->contains(track)) {
      ++count;
    }
  }
  return count;
}

// Overloaded comparison operator for maintaining uniqueness in a std::set<TrackVertexSet>
bool TrackVertexSet::operator<(const TrackVertexSet& other) const {
    // Compare sizes first
    if (this->size() != other.size()) 
        return this->size() < other.size();

    // Use std::lexicographical_compare on the set's iterators directly
    return std::lexicographical_compare(this->begin(), this->end(),
                                        other.begin(), other.end());
}

// Overloaded equality operator to compare elements
bool TrackVertexSet::operator==(const TrackVertexSet& other) const {
  if(this->size() != other.size()) {
    return false;
  }
  return std::equal(this->begin(), this->end(), other.begin());
}

bool TrackVertexSet::operator!=(const TrackVertexSet& other) const {
  return(!(*this == other));
}

// Overloaded subset operator to check if all tracks in *this are in other
bool TrackVertexSet::operator|=(const TrackVertexSet& other) const {
  if (this->size() > other.size()) {
    return false;
  }
  for (const auto& track : *this) {
    if (!other.contains(track)) {
      return false;
    }
  }
  return true;
}

// private methods
void TrackVertexSet::fit() {
  vertex_ = this->size() < 2? TransientVertex() : fitter_->vertex(convertTracks());
}

std::vector<reco::TransientTrack> TrackVertexSet::convertTracks() const {

  std::vector<reco::TransientTrack> ttracks;
  for(const auto& track : *this) {
    ttracks.emplace_back(ttBuilder_->build(*track));
  }  
  return ttracks;
}

double TrackVertexSet::calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const {

  // Define the chi-squared distribution with the specified degrees of freedom
  boost::math::chi_squared chiSqDist(degreesOfFreedom);

  // Compute the p-value as 1 - CDF(chiSquaredValue)
  double pValue = 1.0 - boost::math::cdf(chiSqDist, chiSquaredValue);

  return pValue;
}

// Direct casting from TrackVertexSet to reco::Vertex
TrackVertexSet::operator reco::Vertex() const {

  if(!vertex_.isValid())
    return reco::Vertex();
  
  reco::Vertex recoVertex(VertexHelper::ConvertFitVertex(vertex_));

  // Check for rare edge case where reco::Vertex can fail matrix inversion
  //leading to an exception when checking track compatibility
  int ifail(0);
  if(this->isValid()) {
    GlobalError err(recoVertex.covariance());
    err.matrix().Inverse(ifail);
  }

  if(!recoVertex.isValid() || ifail != 0)
    return reco::Vertex();

  for(const auto &trackRef : *this) {
    recoVertex.add(reco::TrackBaseRef(trackRef), recoVertex.trackWeight(reco::TrackBaseRef(trackRef)));
  }

  return recoVertex;
}

//---------------------------------------------------------------//
//////////////////// TrackVertexSetCollection /////////////////////
//---------------------------------------------------------------//

bool TrackVertexSetCollection::contains(const TrackVertexSet &set) const {
  return this->find(set) != this->end();
}

bool TrackVertexSetCollection::isVertexUnique(const TrackVertexSet &set) const {

  for (auto it = this->begin(); it != this->end(); ++it) {
    if(*it != set && (*it & set) > 0)
      return false;
  }

  return true;
}

TrackVertexSetCollection TrackVertexSetCollection::testUniqueness(const TrackVertexSet &set) const {

  TrackVertexSetCollection retCollection;
  if(!isVertexUnique(set)) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if(*it != set && (*it & set) > 0)
	retCollection.add(*it);
    }
  }
  
  return retCollection;  
}

// Add a set to the collection, avoiding subsets and replacing subsets with larger sets
void TrackVertexSetCollection::add(const TrackVertexSet &set) {

  if (!set.isValid())
    return;
  
  // Iterate over existing sets to check for subsets
  for (auto it = this->begin(); it != this->end();) {
    if (*it |= set) { // If the new set is a superset of an existing set
      it = this->erase(it); // Remove the smaller subset
    } else if (set |= *it) { // If an existing set is a superset of the new set
      return; // Do not insert the new set
    } else {
      ++it;
    }
  }
  
  // Insert the new set if it is not a subset or superset
  this->insert(set);
}

// Subtraction operator: returns a new collection with unique elements from lhs not in rhs
TrackVertexSetCollection TrackVertexSetCollection::operator-(const TrackVertexSetCollection& rhs) const {
  TrackVertexSetCollection result = *this; // Copy current object
  result -= rhs;                           // Use -= to handle subtraction logic
  return result;                           // Return the result
}

// Subtraction assignment operator as a member function
TrackVertexSetCollection& TrackVertexSetCollection::operator-=(const TrackVertexSetCollection& rhs) {
  TrackVertexSetCollection result;

  for (const auto& set : *this) {
    if (rhs.doesNotContain(set)) { // Keep only elements not in rhs
      result.add(set);           // Use add to ensure subset handling
    }
  }
  
  *this = std::move(result); // Replace the current object with the result
  return *this;
}

// Addition operator: combines both collections, ensuring no subsets
TrackVertexSetCollection TrackVertexSetCollection::operator+(const TrackVertexSetCollection& rhs) const {
  TrackVertexSetCollection result = *this; // Copy current object
  result += rhs;                           // Use += to handle addition logic
  return result;                           // Return the result
}

TrackVertexSetCollection& TrackVertexSetCollection::operator+=(const TrackVertexSetCollection& rhs) {
  for (const auto& set : rhs) {
    this->add(set); // Use add to ensure no subsets are inserted
  }
  return *this;
}

reco::VertexCollection TrackVertexSetCollection::vertices() const {
  
  reco::VertexCollection uniqueVertices;
  for(const auto &set : *this) {
    if(set.isValid())
      uniqueVertices.emplace_back(set);
  }

  return uniqueVertices;
}

reco::TrackCollection TrackVertexSetCollection::tracks() const {

  reco::TrackCollection tracks;
  for(const auto &trackRef : this->completeTrackSet())
    tracks.emplace_back(*trackRef);

  return tracks;
}

std::set<reco::TrackRef> TrackVertexSetCollection::completeTrackSet() const {

  std::set<reco::TrackRef> allTracks;
  for(const auto &vertex : *this) {
    for(const auto &track : vertex)
      allTracks.insert(track);
  }
  
  return allTracks;
}

std::set<reco::TrackRef> TrackVertexSetCollection::overlappingTracks() const {
  // Step 1: Track occurrence count
  std::map<reco::TrackRef, int> trackCount;

  // Step 2: Iterate over all TrackVertexSets to count each track
  for (const auto& vertexSet : *this) {
    for (const auto& track : vertexSet.tracks()) {
      trackCount[track]++;
    }
  }

  // Step 3: Collect tracks that appear more than once
  std::set<reco::TrackRef> overlappingTracks;
  for (const auto& trackEntry : trackCount) {
    if (trackEntry.second > 1) {
      overlappingTracks.insert(trackEntry.first);
    }
  }

  return overlappingTracks;
}

bool TrackVertexSetCollection::hasExclusiveVertices() const {

  std::set<reco::TrackRef> encounteredTracks;

  // Iterate through each TrackVertexSet in the collection
  for (const auto& trackVertexSet : *this) {
    // Iterate through each track in the TrackVertexSet
    for (const auto& track : trackVertexSet) {
      // Check if the track is already in the encountered set
      if (encounteredTracks.find(track) != encounteredTracks.end()) {
	return false; // Duplicate found
      }
      // Add the track to the encountered set
      encounteredTracks.insert(track);
    }
  }
  
  return true;
}

size_t TrackVertexSetCollection::countOverlaps() const {

  size_t overlaps = 0;
  
  // Convert to vector for random access iteration
  std::vector<std::reference_wrapper<const TrackVertexSet>> sets(this->begin(), this->end());
  
  // Use indices for clearer bounds and potentially better optimization
  for (size_t i = 0; i < sets.size() - 1; ++i) {
    for (size_t j = i + 1; j < sets.size(); ++j) {
      if ((sets[i].get() & sets[j].get()) > 0) {
	overlaps++;
      }
    }
  }  
  return overlaps;
}

/**
 * @brief Cleans the vertices in the collection by removing incompatible tracks 
 *        and ensuring only valid vertices remain.
 *
 * This function iterates over the vertices in the current collection and evaluates 
 * each track within the vertices based on its compatibility with the vertex and its 
 * alignment relative to a given primary vertex. Tracks that fail the compatibility 
 * or alignment criteria are removed, provided the vertex still contains more than 
 * two tracks after the removal. Vertices with fewer than three tracks are directly 
 * added to the cleaned collection without modification.
 *
 * @param primaryVertex The primary vertex used to evaluate track alignment (cosTheta).
 * 
 * Function details:
 * - Vertices with fewer than 3 tracks are skipped directly and added to the cleaned 
 *   collection as they are.
 * - For each vertex with 3 or more tracks:
 *   - A working copy of the vertex is created to allow safe modification.
 *   - Tracks are evaluated for removal based on:
 *     1. Compatibility value exceeding 1.5.
 *     2. Alignment (cosTheta) being less than 0.5.
 *   - Tracks marked for removal are collected first, then removed from the vertex 
 *     only if the vertex still contains more than 2 tracks.
 * - After removing incompatible tracks, the modified vertex is added to the cleaned 
 *   collection if it remains valid (e.g., satisfies vertex fitting criteria).
 * 
 * At the end of the process, the current collection (`*this`) is replaced with the 
 * cleaned collection, ensuring only valid and cleaned vertices remain.
 */
void TrackVertexSetCollection::cleanVertices(const reco::Vertex& primaryVertex) {
  
  TrackVertexSetCollection cleanedVertices;
  
  for(const auto &vertex : *this) {
    // Skip small vertices directly
    if (vertex.size() < 3) {
      cleanedVertices.add(vertex);
      continue;
    }
    
    // Create a copy of the current vertex for manipulation
    TrackVertexSet workingVertex(vertex);
    std::vector<reco::TrackRef> tracksToRemove;
    
    // Evaluate track compatibility and collect tracks to remove
    for(const auto &trackRef : workingVertex.tracks()) {
      const double compatibility = vertex.compatibility(trackRef);
      const double cosTheta = vertex.trackCosTheta(primaryVertex, trackRef);
      //const bool condition(vertex.size() > 5? /*compatibility > 1.5*/ false : compatibility > 1.5 || cosTheta < 0.5);
      
      if(compatibility > 1.5 || cosTheta < 0.5) {
	tracksToRemove.push_back(trackRef);
      }
    }
    
    // Remove tracks collected for removal
    for(const auto &trackRef : tracksToRemove) {
      if (workingVertex.size() > 2) {
	workingVertex.removeTrack(trackRef);
      }
    }
    /*
    for(const auto &trackRef : workingVertex.tracks()) {
      double compatibility = vertex.compatibility(trackRef);
      double cosTheta = vertex.trackCosTheta(primaryVertex, trackRef);

      if((compatibility > 1.5 || cosTheta < 0.5) && workingVertex.tracks().size() > 2) {
	std::cout << "Found one" << std::endl;
      }
    }
    */
    // Add copy if it's still valid
    if(workingVertex.isValid())
      cleanedVertices.add(workingVertex);
  }
  
  *this = std::move(cleanedVertices);
}

TrackVertexSetCollection TrackVertexSetCollection::postCleaningCuts(const reco::Vertex &primaryVertex) const {

  TrackVertexSetCollection cleanedVertices;
  for(const auto &vertex : *this) {

    LorentzVec vertex4Vec(VertexHelper::GetVertex4Vector(vertex));
    const double p(vertex4Vec.P());
    const double mass(vertex4Vec.M());
    const double pOverE(p/sqrt(p*p + mass*mass));
    const double cosTheta(vertex.cosTheta(primaryVertex));
    const double normChi2(vertex.vertex().normalizedChi2());

    if(mass < 2 || pOverE < 0.6 || cosTheta < 0.75 || normChi2 > 5)
      continue;

    cleanedVertices.add(vertex);
  }
  
  return cleanedVertices;
}

TrackVertexSetCollection TrackVertexSetCollection::disambiguate() const {

  // Step 1: Map each track to the vertex sets containing it
  std::map<reco::TrackRef, TrackVertexSetCollection> trackToVertexMap;
  for (const auto& vertexSet : *this) {
    for (const auto& track : vertexSet) {
      trackToVertexMap[track].add(vertexSet);
    }
  }

  // Step 2: Calculate the CM decay angles and keep the best vertex for each track
  TrackVertexSetCollection rejectCollection;
  for(const auto& [track, vertexList] : trackToVertexMap) {
    TrackVertexSetCollection tempVertexSet(vertexList-rejectCollection);
    if (tempVertexSet.size() > 1) {
      double minCosTheta = std::numeric_limits<double>::max();
      TrackVertexSet bestVertex;

      for(const auto& vertex : tempVertexSet) {
        double cosTheta = VertexHelper::CalculateCMCosTheta(vertex, *track);
        if (cosTheta < minCosTheta) {
          minCosTheta = cosTheta;
          bestVertex = vertex;
        }
      }

      tempVertexSet.erase(bestVertex);
      rejectCollection += tempVertexSet;
    }
  }

  return *this-rejectCollection;
}

TrackVertexSetCollection TrackVertexSetCollection::candidates(const reco::Vertex &primaryVertex) const {
  TrackVertexSetCollection candidateVertices(TrackVertexSetHelper::applyPostDisambiguationTrackCuts(*this, primaryVertex));
  candidateVertices = TrackVertexSetHelper::applySignificanceCut(candidateVertices, primaryVertex);
  return candidateVertices;
}

TrackVertexSetCollection TrackVertexSetCollection::disambiguatedVerticesNoFinalCuts(const reco::Vertex &primaryVertex) const {

  // Step 1: Clean vertices with 3 or more tracks (converts most silver vertices into gold)
  TrackVertexSetCollection finalCandidates(*this);
  finalCandidates.cleanVertices(primaryVertex);

  // Step 2: Apply vertex-level variable cuts to remove a majority of background vertices.
  // Applying these before disambiguation was found to maximize signal efficiency.
  finalCandidates = finalCandidates.postCleaningCuts(primaryVertex);

  // Step 3: Perform disambiguation for non-unique vertices
  finalCandidates = finalCandidates.disambiguate();

  return finalCandidates;
}

// Full selection for dilepton vertex candidates (only includes two-track vertices)
TrackVertexSetCollection TrackVertexSetCollection::fullSelection(const reco::Vertex &primaryVertex) const {
  return this->disambiguatedVerticesNoFinalCuts(primaryVertex).candidates(primaryVertex);
}

void TrackVertexSetHelper::cleanDuplicates(TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex) {

  TrackVertexSetCollection cleanedCollection;
  std::vector<TrackVertexSet> vertices(collection.begin(), collection.end());
  std::sort(vertices.begin(), vertices.end(), [](const TrackVertexSet& a, const TrackVertexSet& b) {return a.normChi2() < b.normChi2();});

  // Vector to track which vertices are duplicates
  std::vector<bool> isDuplicate(vertices.size(), false);

  // Compare each vertex with all others that come after it
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (isDuplicate[i]) continue; // Skip if already marked as duplicate
    
    for (size_t j = i + 1; j < vertices.size(); ++j) {
      if (isDuplicate[j]) continue; // Skip if already marked as duplicate
      
      const TrackVertexSet& vertexA = vertices[i];
      const TrackVertexSet& vertexB = vertices[j];
      
      // Calculate relative differences
      /*
      double dxyA = vertexA.dxy(primaryVertex);
      double dxyB = vertexB.dxy(primaryVertex);
      //double dxyRelDiff = 2.0 * std::abs(dxyA - dxyB) / (dxyA + dxyB);
      
      */
      // Check if this is a duplicate based on all three criteria
      if(vertices[i].size() == vertices[j].size() && vertexA.dxySeparation(vertexB) < 3*vertexA.distanceError(vertexB)) {
	// Mark the vertex with higher normChi2 as duplicate
	// (Since we sorted by normChi2, j is always higher)
	isDuplicate[j] = true;
      }
    }
  }
  
  // Add non-duplicate vertices to cleanedCollection
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (!isDuplicate[i]) {
      cleanedCollection.insert(vertices[i]);
    }
  }

  collection = cleanedCollection;
}

TrackVertexSetCollection TrackVertexSetHelper::cleanedVertices(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex) {
  TrackVertexSetCollection cleanedSet(collection);
  cleanedSet.cleanVertices(primaryVertex);
  return cleanedSet;
}

TrackVertexSetCollection TrackVertexSetHelper::applyPostDisambiguationTrackCuts(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex) {

  TrackVertexSetCollection filteredVertices;
  for(const auto &vertex : collection) {
    // Skip vertices that don't have exactly 2 tracks
    if(vertex.size() != 2) {
      //filteredVertices.insert(vertex);
      continue;
    }
      
    // Check if the vertex passes the track angular cuts
    bool passCuts = true;
    for(const auto &trackRef : vertex) {
      const double trackCosTheta = vertex.trackCosTheta(primaryVertex, trackRef);
      const double trackCosThetaAtCM = vertex.trackDecayAngleCM(trackRef);

      // Apply track-level cuts
      if(trackCosTheta < 0.5 ||
         fabs(trackCosThetaAtCM) > 1.8 - trackCosTheta ||
         fabs(trackCosThetaAtCM) > 0.95) {
        passCuts = false;
        break;
      }
    }
    if(passCuts)
      filteredVertices.insert(vertex);
  }

  return filteredVertices;
}

TrackVertexSetCollection TrackVertexSetHelper::applySignificanceCut(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex, const double cutVal) {

  TrackVertexSetCollection candidateVertices;

  for(const auto &vertex : collection) {
    // Calculate displacement significance
    const int charge(VertexHelper::CalculateTotalCharge(vertex));
    const double dxy = VertexHelper::CalculateDxy(vertex, primaryVertex);
    const double dxyError = VertexHelper::CalculateDxyError(vertex, primaryVertex);
    const double dxySignificance = dxy / dxyError;

    // Apply vertex-level cuts and add to candidates if passed
    if(dxySignificance > cutVal && (vertex.size() != 2 || charge == 0)) {
      candidateVertices.add(vertex);
    }
  }

  return candidateVertices;
}

TrackVertexSetCollection TrackVertexSetHelper::highMultiplicityFilter(const TrackVertexSetCollection &collection, const size_t threshold) {

  TrackVertexSetCollection outputCollection;

  for(const auto &vertex : collection) {
    if(vertex.size() >= threshold)
      outputCollection.insert(vertex);
  }

  return outputCollection;
}

TrackVertexSetCollection TrackVertexSetHelper::highMultiplicityDisambiguation(const TrackVertexSetCollection &collection) {

  TrackVertexSetCollection outputCollection;

  std::map<reco::TrackRef, TrackVertexSetCollection> trackToVertexMap;
  for(const auto& vertexSet : collection) {
    for(const auto& track : vertexSet) {
      trackToVertexMap[track].add(vertexSet);
    }
  }

  TrackVertexSetCollection rejectCollection;
  for(const auto& [track, vertexList] : trackToVertexMap) {
    TrackVertexSetCollection tempVertexSet(vertexList-rejectCollection);
    if(tempVertexSet.size() > 1) {

      size_t minSize(0);
      TrackVertexSet largestVertex;
      for(const auto& vertex : tempVertexSet) {
	if(vertex.size() > minSize) {
	  minSize = vertex.size();
	  largestVertex = vertex;
	}
      }

      tempVertexSet.erase(largestVertex);
      rejectCollection += tempVertexSet;
      
    }
  }
  return collection-rejectCollection;
}

TrackVertexSetCollection TrackVertexSetHelper::hadronicVertexCuts(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex) {

  TrackVertexSetCollection outputCollection;

  for(const auto &vertex : collection) {

    LorentzVec vertex4Vec(VertexHelper::GetVertex4Vector(vertex));
    const double p(vertex4Vec.P());
    const double mass(vertex4Vec.M());
    const double pOverE(p/sqrt(p*p + mass*mass));
    const double cosTheta(vertex.cosTheta(primaryVertex));
    const double normChi2(vertex.vertex().normalizedChi2());
    const double decayAngle(vertex.decayAngle());

    if(vertex.size() < 5 || mass < 2 || pOverE < 0.6 || cosTheta < 0. || fabs(decayAngle) > 0.9 || normChi2 > 5)
      continue;

    outputCollection.insert(vertex);
  }

  return outputCollection;
}

TrackVertexSetCollection TrackVertexSetHelper::hadronicVertexCandidates(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex) {
  TrackVertexSetCollection finalCandidates = hadronicVertexCuts(collection, primaryVertex);
  finalCandidates = highMultiplicityDisambiguation(finalCandidates);
  finalCandidates = applySignificanceCut(finalCandidates, primaryVertex, 40);
  return finalCandidates;
}
