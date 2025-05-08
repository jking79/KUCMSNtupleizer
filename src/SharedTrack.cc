#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/SharedTrack.h"

std::optional<VertexInfo> SharedTrack::getGoldVertex() const {
  for (const auto &vertex : vertices_) {
    if (vertex.isGold()) {
      return vertex;
    }
  }
  return std::nullopt;
}

std::vector<VertexInfo> SharedTrack::getVerticesByType(VertexType type) const {
  std::vector<VertexInfo> result;
  for (const auto &vertex : vertices_) {
    if (vertex.getType() == type) {
      result.push_back(vertex);
    }
  }
  return result;
}

bool SharedTrack::operator<(const SharedTrack &other) const {
  return this->trackID_ < other.trackID_;
}


void SharedTrack::operator+=(const VertexInfo& vertex) {
  vertices_.insert(vertex);
  flagsCached_ = false;
}

void SharedTrack::merge(const SharedTrack &other) {
  for (const auto &vertex : other.vertices_) {
    vertices_.insert(vertex);
  }
  flagsCached_ = false;
}

void SharedTrack::setVertexType(int vertexID, VertexType type) {
  auto it = std::find_if(vertices_.begin(), vertices_.end(),
			 [vertexID](const VertexInfo& v) { return v.id() == vertexID; });
  
  if (it != vertices_.end()) {
    VertexInfo updatedVertex = *it;
    vertices_.erase(it);
    updatedVertex.setType(type);
    vertices_.insert(updatedVertex);
    flagsCached_ = false;
  }
}

void SharedTrack::removeVertices(const std::vector<int>& verticesToRemove) {
  for (const auto& vertexID : verticesToRemove) {
    auto it = std::find_if(vertices_.begin(), vertices_.end(),
			   [vertexID](const VertexInfo& v) { return v.id() == vertexID; });
    if (it != vertices_.end()) {
      vertices_.erase(it);
    }
  }
  flagsCached_ = false;
}

void SharedTrack::print() const {
  std::cout << "SharedTrack ID: " << trackID_ << " | Vertices:\n";
  for (const auto &vertex : vertices_) {
    std::cout << "  Vertex ID: " << vertex.id()
	      << " -> Type: " << vertexTypeToString(vertex.getType())
	      << " (Association ID: " << vertex.associationID() << ")\n";
  }
}

bool SharedTrack::hasBronze() const {
  updateTypeCache();
  return typeFlags_ & 1;
}

bool SharedTrack::hasSilver() const {
  updateTypeCache();
  return typeFlags_ & 2;
}

bool SharedTrack::hasGold() const {
  updateTypeCache();
  return typeFlags_ & 4;
}

bool SharedTrack::hasBackground() const {
  updateTypeCache();
  return typeFlags_ & 8;
}

size_t SharedTrack::nBronzeWithGold() const {
  if (!hasGold()) return 0;
  return nType(kBronze);
}

size_t SharedTrack::nBronzeWithBronze() const {
  size_t count = nType(kBronze);
  return (count > 1) ? count : 0;
}

size_t SharedTrack::nBronzeWithBackground() const {
  if (!hasBackground()) return 0;
  return nType(kBronze);
}

size_t SharedTrack::nBackgroundWithBackground() const {
  size_t count = nType(kBackground);
  return (count > 1) ? count : 0;
}

// SharedTrack private methods
std::string SharedTrack::vertexTypeToString(VertexType type) {
  switch (type) {
  case kGold: return "gold";
  case kSilver: return "silver";
  case kBronze: return "bronze";
  case kBackground: return "background";
  case kNone: return "none";
  default: return "unknown";
  }
}

void SharedTrack::updateTypeCache() const {
  if (!flagsCached_) {
    typeFlags_ = 0;
    for (const auto& vertex : vertices_) {
      switch (vertex.getType()) {
      case kBronze: typeFlags_ |= 1; break;
      case kSilver: typeFlags_ |= 2; break;
      case kGold: typeFlags_ |= 4; break;
      case kBackground: typeFlags_ |= 8; break;
      default: break;
      }
    }
    flagsCached_ = true;
  }
}

size_t SharedTrack::nType(VertexType type) const {
  size_t count = 0;
  for (const auto& vertex : vertices_) {
    if (vertex.getType() == type) count++;
  }
  return count;
}


// SharedTrackSet
SharedTrackSet::SharedTrackSet(const std::vector<unsigned int> &trackIDs,
			       const std::vector<unsigned int> &vertexIDs) {
  if (trackIDs.size() != vertexIDs.size()) {
    throw std::invalid_argument("trackIDs and vertexIDs must have the same size");
  }
  
  for (size_t i = 0; i < trackIDs.size(); ++i) {
    SharedTrack track(trackIDs[i]);
    track += VertexInfo(vertexIDs[i], i);
    *this += track;
  }
  
  cleanUniqueTracks();
}

void SharedTrackSet::SetVertexType(const std::vector<bool> &isBronze,
				   const std::vector<bool> &isSilver,
				   const std::vector<bool> &isGold) {
  std::set<SharedTrack> updatedTracks;
  
  for (const auto& track : *this) {
    SharedTrack updatedTrack = track;
    
    for (const auto& vertex : track.getVertices()) {
      const size_t vertexID = vertex.id();
      if (vertexID >= isBronze.size() || vertexID >= isSilver.size() || vertexID >= isGold.size()) {
	continue;
      }
      
      if (isGold[vertexID]) {
	updatedTrack.setVertexType(vertexID, SharedTrack::kGold);
      } else if (isSilver[vertexID]) {
	updatedTrack.setVertexType(vertexID, SharedTrack::kSilver);
      } else if (isBronze[vertexID]) {
	updatedTrack.setVertexType(vertexID, SharedTrack::kBronze);
      } else {
	updatedTrack.setVertexType(vertexID, SharedTrack::kBackground);
      }
    }
    
    updatedTracks.insert(updatedTrack);
  }
  
  this->swap(updatedTracks);
}

void SharedTrackSet::operator+=(const SharedTrack &sharedTrack) {
    auto it = this->find(sharedTrack);
    if (it != this->end()) {
      SharedTrack updatedTrack = *it;
      updatedTrack.merge(sharedTrack);
      this->erase(it);
      this->insert(updatedTrack);
    } else {
      this->insert(sharedTrack);
    }
  }

void SharedTrackSet::print() const {
  std::cout << "Shared Tracks:\n";
  for (const auto &track : *this) {
    track.print();
  }
}

std::vector<VertexInfo> SharedTrackSet::disambiguateThrough(const std::vector<float> &vertexVar) const {

  std::vector<VertexInfo> disambiguatedVertices;
  
  for (auto it = this->begin(); it != this->end(); ++it) {

    VertexInfo minVertex;
    double minVal(2.);
    for(const auto& vertex : it->getVertices()) {
      const double val((vertexVar[vertex.associationID()]));
      
      if(minVal > val) {
	minVal = val;
	minVertex = vertex;
      }
    }
    
    if(minVertex.id() >= 0) 
      disambiguatedVertices.push_back(minVertex);

  }

  return disambiguatedVertices;
}

// SharedTrackSet Private methods
void SharedTrackSet::cleanUniqueTracks() {
  for (auto it = this->begin(); it != this->end();) {
    if (it->size() == 1) {
      it = this->erase(it);
    } else {
      ++it;
    }
  }
}
