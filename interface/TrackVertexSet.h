#ifndef KUCMSNtupleizer_KUCMSNtupleizer_TrackVertexSet_h
#define KUCMSNtupleizer_KUCMSNtupleizer_TrackVertexSet_h

#include <unordered_map>
#include <limits>
#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/SharedTrack.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchingTools.h"

class TrackVertexSetHelper;

class TrackVertexSet : public std::set<reco::TrackRef> {
 public:

  // Constructor with initializer list
  TrackVertexSet() = default;
  TrackVertexSet(const std::vector<reco::TrackRef> &init, const TransientTrackBuilder* ttBuilder);
  TrackVertexSet(std::initializer_list<reco::TrackRef> init, const TransientTrackBuilder* ttBuilder);
  
  // Copy constructor - explicitly inherit from base class
  TrackVertexSet(const TrackVertexSet& other);
  
  // Destructor
  ~TrackVertexSet() = default;
  
  TransientVertex tvertex() const { return vertex_; }
  reco::Vertex vertex() const {return reco::Vertex(*this);}
  
  bool isValid() const {return vertex_.isValid();}
  double pt() const {return VertexHelper::GetVertex4Vector(*this).pt();}
  double eta() const {return VertexHelper::GetVertex4Vector(*this).eta();}
  double phi() const {return VertexHelper::GetVertex4Vector(*this).phi();}
  double mass() const {return VertexHelper::GetVertex4Vector(*this).M();}
  double decayAngle() const {return VertexHelper::CalculateDecayAngle(*this);}
  double normChi2() const {return vertex_.normalisedChiSquared();}
  double pValue() const {return calculateChiSquaredPValue(vertex_.totalChiSquared(), vertex_.degreesOfFreedom());}
  double dxySeparation(const TrackVertexSet& other) const;
  double dxySeparationError(const TrackVertexSet& other) const;
  double dxySeparationSignificance(const TrackVertexSet& other) const;
  double distance(const TrackVertexSet& other) const;
  double distanceError(const TrackVertexSet& other) const;
  double distanceSignificance(const TrackVertexSet& other) const;
  double compatibility(const reco::TrackRef &track) const;
  double dxy(const reco::Vertex &primaryVertex)	const {return VertexHelper::CalculateDxy(*this, primaryVertex);}
  double dxyError(const reco::Vertex &primaryVertex) const {return VertexHelper::CalculateDxyError(*this, primaryVertex);}
  double cosTheta(const reco::Vertex &primaryVertex) const { return VertexHelper::CalculateCosTheta(primaryVertex, *this); }
  double trackCosTheta(const reco::Vertex &primaryVertex, const reco::TrackRef &track) const{ return TrackHelper::CalculateCosTheta(primaryVertex, *this, *track); }
  double trackDecayAngleCM(const reco::TrackRef &track) const { return VertexHelper::CalculateCMCosTheta(*this, *track); }
  double shiftDzAfterTrackRemoval(const reco::TrackRef &track) const;
  double shift3DAfterTrackRemoval(const reco::TrackRef &track) const;
  GlobalPoint position() const {return vertex_.position();}
  std::vector<reco::Track> trackList() const;
  std::vector<reco::TrackRef> tracks() const;
  std::vector<reco::TrackRef> commonTracks(const TrackVertexSet& other) const;
  
  void printTrackInfo() const;

  // element manipulation
  void clearTracks();
  void addTrack(const reco::TrackRef& track);
  void removeTrack(const reco::TrackRef& track);
  bool contains(const reco::TrackRef& track) const;
  
  // interaction with other TrackVertexSets
  bool isSmallerThan(const TrackVertexSet &other) const {return this->size() < other.size();}
  bool isLargerThan(const TrackVertexSet &other) const {return this->size() > other.size();}
  TrackVertexSet merge(const TrackVertexSet &other) const;

  TrackVertexSet& operator=(const TrackVertexSet& other); 
  TrackVertexSet operator+(const TrackVertexSet& other) const;
  size_t operator&(const TrackVertexSet& other) const;
  bool operator<(const TrackVertexSet& other) const;
  bool operator==(const TrackVertexSet& other) const;
  bool operator!=(const TrackVertexSet& other) const;
  bool operator|=(const TrackVertexSet& other) const;

  operator reco::Vertex() const;
  
 private:
  const TransientTrackBuilder* ttBuilder_;
  std::unique_ptr<KalmanVertexFitter> fitter_; 
  TransientVertex vertex_;
  
  void fit();
  std::vector<reco::TransientTrack> convertTracks() const;
  double calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const;
};

class TrackVertexSetCollection : public std::set<TrackVertexSet> {

public:
  TrackVertexSetCollection() = default;

  TrackVertexSetCollection(const TrackVertexSetCollection &other) : std::set<TrackVertexSet>(other) {}

  virtual ~TrackVertexSetCollection() = default;

  bool contains(const TrackVertexSet &set) const;
  bool doesNotContain(const TrackVertexSet &set) const {return !contains(set); }
  bool isVertexUnique(const TrackVertexSet &set) const;
  TrackVertexSetCollection testUniqueness(const TrackVertexSet &set) const;
  
  TrackVertexSetCollection operator+(const TrackVertexSetCollection& rhs) const;
  TrackVertexSetCollection operator-(const TrackVertexSetCollection& rhs) const;
  TrackVertexSetCollection& operator+=(const TrackVertexSetCollection& rhs);
  TrackVertexSetCollection& operator-=(const TrackVertexSetCollection& rhs);
  
  void add(const TrackVertexSet &set);
  
  reco::VertexCollection vertices() const;
  reco::TrackCollection tracks() const;
  std::set<reco::TrackRef> completeTrackSet() const;
  std::set<reco::TrackRef> overlappingTracks() const;
  bool hasExclusiveVertices() const;
  size_t countOverlaps() const;

  // Transformations (Cutting and Cleaning)
  void cleanVertices(const reco::Vertex& primaryVertex);
  TrackVertexSetCollection disambiguate() const;
  TrackVertexSetCollection postCleaningCuts(const reco::Vertex &primaryVertex) const;
  TrackVertexSetCollection candidates(const reco::Vertex &primaryVertex) const;
  TrackVertexSetCollection disambiguatedVerticesNoFinalCuts(const reco::Vertex &primaryVertex) const;
  TrackVertexSetCollection fullSelection(const reco::Vertex &primaryVertex) const;

};

class TrackVertexSetHelper {
public:
  static void cleanDuplicates(TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex);
  static TrackVertexSetCollection cleanedVertices(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex);
  static TrackVertexSetCollection applyPostDisambiguationTrackCuts(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex);
  static TrackVertexSetCollection applySignificanceCut(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex, const double cutVal = 25);
  static TrackVertexSetCollection highMultiplicityFilter(const TrackVertexSetCollection &collection, const size_t threshold = 4);
  static TrackVertexSetCollection highMultiplicityDisambiguation(const TrackVertexSetCollection &collection);
  static TrackVertexSetCollection hadronicVertexCuts(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex);
  static TrackVertexSetCollection hadronicVertexCandidates(const TrackVertexSetCollection &collection, const reco::Vertex &primaryVertex);
};

#endif
