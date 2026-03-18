#pragma once

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/HYDDRABase.h"

// Hadronic displaced vertex reconstruction.
// Targets high-multiplicity vertices (size >= minSize).
// Disambiguation keeps the largest vertex for each shared track.
class HadronicHYDDRA : public HYDDRABase<HadronicHYDDRA> {

 public:

  HadronicHYDDRA(const edm::ParameterSet& pset) : HYDDRABase(pset) {
    minSize_       = pset.getParameter<int>("minSize");
    minCosTheta_   = pset.getParameter<double>("minCosTheta");
    maxDecayAngle_ = pset.getParameter<double>("maxDecayAngle");
  }

  // Stage 1: Seeding
  void seedingImpl(const std::vector<reco::TrackRef>& tracks) {
    this->generateSeeds(tracks);
  }

  // Stage 2: Merging
  void mergingImpl() {
    this->mergeVertices();
  }

  // Stage 3: Cleaning
  // Apply vertex-level kinematic cuts. No per-track removal for hadronic vertices.
  void cleaningImpl() {
    if (this->empty() || !primaryVertex_ || !primaryVertex_->isValid()) return;

    HYDDRA_DBG("[Hadronic] Cleaning " << this->size() << " vertices...\n");

    TrackVertexSetCollection validCandidates;
    for (const auto& vertex : *this) {
      if (!failsHadronicCuts(vertex)) {
	validCandidates.add(vertex);
      }
    }

    this->clear();
    for (const auto& v : validCandidates) {
      this->add(v);
    }
  }

  // Stage 4: Disambiguation
  // For each track shared between vertices, keep the largest vertex
  // and reject the rest.
  void disambiguationImpl() {
    if (this->empty()) return;

    // Map each track to the vertices that contain it
    std::map<reco::TrackRef, std::vector<const TrackVertexSet*>> trackToVertexMap;
    for (const auto& vertex : *this) {
      for (const auto& track : vertex) {
	trackToVertexMap[track].push_back(&vertex);
      }
    }

    // Identify vertices to reject
    std::set<const TrackVertexSet*> rejectedPtrs;

    for (const auto& [track, vtxPtrList] : trackToVertexMap) {
      std::vector<const TrackVertexSet*> candidates;
      for (const auto* vtxPtr : vtxPtrList) {
	if (rejectedPtrs.find(vtxPtr) == rejectedPtrs.end()) {
	  candidates.push_back(vtxPtr);
	}
      }

      if (candidates.size() > 1) {
	size_t maxSize = 0;
	const TrackVertexSet* bestVertex = nullptr;

	for (const auto* vtxPtr : candidates) {
	  if (vtxPtr->size() > maxSize) {
	    maxSize = vtxPtr->size();
	    bestVertex = vtxPtr;
	  }
	}

	for (const auto* vtxPtr : candidates) {
	  if (vtxPtr != bestVertex) {
	    rejectedPtrs.insert(vtxPtr);
	  }
	}
      }
    }

    // Rebuild collection without rejected vertices
    if (!rejectedPtrs.empty()) {
      HYDDRA_DBG("[Hadronic] Disambiguation dropping " << rejectedPtrs.size() << " vertices.\n");

      TrackVertexSetCollection survivors;
      for (const auto& vtx : *this) {
	if (rejectedPtrs.find(&vtx) == rejectedPtrs.end()) {
	  survivors.add(vtx);
	}
      }
      this->clear();
      for (const auto& v : survivors) {
	this->add(v);
      }
    }
  }

  // Stage 5: Filtering
  // Apply displacement significance cut.
  void filteringImpl() {
    if (this->empty() || !primaryVertex_ || !primaryVertex_->isValid()) return;

    [[maybe_unused]] const size_t nInput = this->size();
    HYDDRA_DBG("[Hadronic] Final Filtering on " << nInput << " vertices...\n");

    // Cut counters
    size_t nFailDxyError = 0;
    size_t nFailDxySignificance = 0;

    TrackVertexSetCollection finalVertices;

    for (const auto& vertex : *this) {
      const double dxy      = VertexHelper::CalculateDxy(vertex, *primaryVertex_);
      const double dxyError = VertexHelper::CalculateDxyError(vertex, *primaryVertex_);
      if (dxyError <= 0) {
        nFailDxyError++;
        continue;
      }

      if (dxy / dxyError <= minDxySignificance_) {
        nFailDxySignificance++;
        continue;
      }

      finalVertices.add(vertex);
    }

    this->clear();
    for (const auto& v : finalVertices) {
      this->add(v);
    }

    HYDDRA_DBG("[Hadronic] Filtering summary (" << nInput << " -> " << this->size() << "):\n"
               << "  dxyError <= 0:          " << nFailDxyError << " failed\n"
               << "  dxy significance:       " << nFailDxySignificance << " failed (cut: " << minDxySignificance_ << ")\n"
               << "  PASSED:                 " << this->size() << "\n");
  }

 private:

  // Hadronic-specific cut thresholds
  int minSize_;
  double minCosTheta_;
  double maxDecayAngle_;

  // Returns true if the vertex fails hadronic selection.
  // Cuts are ordered cheapest-first.
  bool failsHadronicCuts(const TrackVertexSet& v) const {
    if (!isValidVertex(v)) return true;
    if (static_cast<int>(v.size()) < minSize_) return true;
    if (v.vertex().normalizedChi2() > maxNormChi2_) return true;
    if (v.cosTheta(*primaryVertex_) < minCosTheta_) return true;
    if (std::fabs(v.decayAngle()) > maxDecayAngle_) return true;

    auto p4 = VertexHelper::GetVertex4Vector(v);
    double mass = p4.M();
    if (mass < minMass_) return true;

    double p = p4.P();
    double e = std::sqrt(p * p + mass * mass);
    double pOverE = (e > 1e-6) ? p / e : 0.0;
    if (pOverE < minPOverE_) return true;

    return false;
  }
};
