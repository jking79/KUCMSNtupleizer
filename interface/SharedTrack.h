#ifndef KUCMSNtupleizer_KUCMSNtupleizer_SharedTrack_h
#define KUCMSNtupleizer_KUCMSNtupleizer_SharedTrack_h

#include <set>
#include <map>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <optional>
#include <algorithm>
#include <cmath>

class VertexInfo {

  public:
  enum VertexType { kBronze, kSilver, kGold, kBackground, kNone };

  VertexInfo() : vertexID_(-1), trackAssociationID_(-1), type_(kNone) {}
  
  VertexInfo(const size_t vertexID, const size_t trackAssociationID, const VertexType vertexType=kNone) :
    vertexID_(vertexID),
    trackAssociationID_(trackAssociationID),
    type_(vertexType){}

  bool operator<(const VertexInfo& other) const {
    if (vertexID_ != other.vertexID_) return vertexID_ < other.vertexID_;
    if (trackAssociationID_ != other.trackAssociationID_) return trackAssociationID_ < other.trackAssociationID_;
    return type_ < other.type_; // Ensures unique ordering in std::set
  }

  int id() const { return vertexID_; }
  int associationID() const { return trackAssociationID_; }

  VertexType getType() const {return type_;}

  void setType(const VertexType type) {type_ = type;}

  bool isBronze() const { return type_ == kBronze; }
  bool isSilver() const { return type_ == kSilver; }
  bool isGold() const { return type_ == kGold; }
  bool isBackground() const { return type_ == kBackground; }

private:
  int vertexID_, trackAssociationID_;
  VertexType type_;
};

class SharedTrack {
 public:
  using VertexType = VertexInfo::VertexType;
  static constexpr VertexType kBronze = VertexInfo::kBronze;
  static constexpr VertexType kSilver = VertexInfo::kSilver;
  static constexpr VertexType kGold = VertexInfo::kGold;
  static constexpr VertexType kBackground = VertexInfo::kBackground;
  static constexpr VertexType kNone = VertexInfo::kNone;

  explicit SharedTrack(size_t trackID) : trackID_(trackID) {}

  size_t size() const { return vertices_.size(); }
  size_t id() const { return trackID_; }

  // Getter functions
  const std::set<VertexInfo>& getVertices() const { return vertices_; }
  std::optional<VertexInfo> getGoldVertex() const;
  std::vector<VertexInfo> getVerticesByType(VertexType type) const;

  // operators
  bool operator<(const SharedTrack &other) const;
  void operator+=(const VertexInfo& vertex);
  
  void merge(const SharedTrack &other);
  void setVertexType(int vertexID, VertexType type);
  void removeVertices(const std::vector<int>& verticesToRemove);
  
  bool hasBronze() const;
  bool hasSilver() const;
  bool hasGold() const;
  bool hasBackground() const;

  // Counting methods
  size_t nBronze() const {return nType(kBronze);}
  size_t nSilver() const {return nType(kSilver);}
  size_t nGold() const {return nType(kGold);}
  size_t nBackground() const {return nType(kBackground);}

  size_t nBronzeWithGold() const;
  size_t nBronzeWithBronze() const;
  size_t nBronzeWithBackground() const;
  size_t nBackgroundWithBackground() const;

  void print() const;
  
 private:
  size_t trackID_;
  std::set<VertexInfo> vertices_;
  mutable uint8_t typeFlags_ = 0;
  mutable bool flagsCached_ = false;
  
  static std::string vertexTypeToString(VertexType type);
  void updateTypeCache() const;
  size_t nType(VertexType type) const;
};

class SharedTrackSet : public std::set<SharedTrack> {
 public:
  SharedTrackSet() = default;
  
  SharedTrackSet(const std::vector<unsigned int> &trackIDs,
                 const std::vector<unsigned int> &vertexIDs);

  void SetVertexType(const std::vector<bool> &isBronze,
                     const std::vector<bool> &isSilver,
                     const std::vector<bool> &isGold);

  template <typename T, typename Predicate>
  void cut(const std::vector<T>& vertexValues, Predicate shouldRemove);

  template <typename T>
  void cutBelow(const std::vector<T>& values, T threshold) {
    cut(values, [threshold](T val) { return val < threshold; });
  }
  
  template <typename T>
  void cutAbove(const std::vector<T>& values, T threshold) {
    cut(values, [threshold](T val) { return val > threshold; });
  }
  
  template <typename T>
  void cutAbsoluteValue(const std::vector<T>& values, T threshold) {
    cut(values, [threshold](T val) { return std::abs(val) > threshold; });
  }

  template <typename T>
  void cutIfNotAt(const std::vector<T>& values, T threshold) {
    cut(values, [threshold](T val) { return val != threshold; });
  }
  
  bool contains(const SharedTrack &sharedTrack) const {
    return this->find(sharedTrack) != this->end();
  }
  
  void operator+=(const SharedTrack &sharedTrack);
  
  void print() const;

  std::vector<VertexInfo> disambiguateThrough(const std::vector<float> &vertexVar) const;
  
 private:
  void cleanUniqueTracks();

};

template <typename T, typename Predicate>
void SharedTrackSet::cut(const std::vector<T>& vertexValues, Predicate shouldRemove) {
  std::set<SharedTrack> updatedTracks;
  
  for (const auto& track : *this) {
    std::vector<int> vertexIDsToRemove;
    
    for (const auto& vertex : track.getVertices()) {
      const size_t vertexID = vertex.id();
      if (vertexID < vertexValues.size() && shouldRemove(vertexValues[vertexID])) {
	vertexIDsToRemove.push_back(vertexID);
      }
    }
    
    if (vertexIDsToRemove.empty()) {
      updatedTracks.insert(track);
      continue;
    }
    
    SharedTrack updatedTrack = track;
    updatedTrack.removeVertices(vertexIDsToRemove);
    
    if (!updatedTrack.getVertices().empty()) {
      updatedTracks.insert(updatedTrack);
    }
  }
  
  this->swap(updatedTracks);
  cleanUniqueTracks();
}

#endif
