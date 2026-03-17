#pragma once

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/HYDDRABase.h"

// Leptonic displaced vertex reconstruction.
// Targets 2-track vertices (e.g. displaced dileptons).
// Disambiguation keeps the vertex with the smallest CM cosTheta.
class LeptonicHYDDRA : public HYDDRABase<LeptonicHYDDRA> {

 public:

  LeptonicHYDDRA(const edm::ParameterSet& pset) : HYDDRABase(pset) {
    // Cleaning thresholds
    maxCompatibility_ = pset.getParameter<double>("maxCompatibility");
    minCleanCosTheta_ = pset.getParameter<double>("minCleanCosTheta");
    useDiagonalCut_   = pset.getParameter<bool>("useDiagonalCut");
    cleanCutSlope_    = pset.getParameter<double>("cleanCutSlope");

    // Final filtering cuts (post-disambiguation, 2-track only)
    minTrackCosTheta_             = pset.getParameter<double>("minTrackCosTheta");
    maxTrackCosThetaCM_Limit_     = pset.getParameter<double>("maxTrackCosThetaCM_Limit");
    maxTrackCosThetaCM_Intercept_ = pset.getParameter<double>("maxTrackCosThetaCM_Intercept");
    trackCosThetaCM_Slope_        = pset.getParameter<double>("trackCosThetaCM_Slope");
    requireChargeNeutrality_      = pset.getParameter<bool>("requireChargeNeutrality");
    minVtxCosTheta_               = pset.getParameter<double>("minVtxCosTheta");
    useAbsVtxCosTheta_            = pset.getParameter<bool>("useAbsVtxCosTheta");
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
  // For 2-track vertices, only apply kinematic cuts.
  // For larger vertices, remove tracks with poor compatibility or alignment,
  // then apply kinematic cuts to whatever remains.
  void cleaningImpl() {
    if (this->empty() || !primaryVertex_ || !primaryVertex_->isValid()) return;

    HYDDRA_DBG("[Leptonic] Cleaning " << this->size() << " vertices...\n");

    TrackVertexSetCollection cleanedVertices;

    for (const auto& vertex : *this) {
      // Two-track vertices skip track removal
      if (vertex.size() <= 2) {
	if (!failsKinematics(vertex)) cleanedVertices.add(vertex);
	continue;
      }

      // Scan for bad tracks before modifying
      std::vector<reco::TrackRef> tracksToRemove;
      for (const auto& trackRef : vertex.tracks()) {
	double compatibility = vertex.compatibility(trackRef);
	double cosTheta      = vertex.trackCosTheta(*primaryVertex_, trackRef);

	const bool isBadTrack = useDiagonalCut_
	  ? cosTheta < cleanCutSlope_ * compatibility + minCleanCosTheta_
	  : compatibility > maxCompatibility_ || cosTheta < minCleanCosTheta_;
	if (isBadTrack) {
	  tracksToRemove.push_back(trackRef);
	}
      }

      // No bad tracks — keep original if it passes kinematics
      if (tracksToRemove.empty()) {
	if (!failsKinematics(vertex)) cleanedVertices.add(vertex);
	continue;
      }

      // Remove bad tracks from a working copy
      TrackVertexSet workingVertex(vertex);
      for (const auto& trackRef : tracksToRemove) {
	if (workingVertex.size() > 2) {
	  workingVertex.removeTrack(trackRef);
	}
      }

      if (workingVertex.isValid() && !failsKinematics(workingVertex)) {
	cleanedVertices.add(workingVertex);
      }
    }

    this->clear();
    for (const auto& v : cleanedVertices) {
      this->add(v);
    }
  }

  // Stage 4: Disambiguation
  // For each track shared between vertices, keep the vertex with the
  // smallest CM cosTheta (most isotropic decay) and reject the rest.
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
	double minCosTheta = std::numeric_limits<double>::max();
	const TrackVertexSet* bestVertex = nullptr;

	for (const auto* vtxPtr : candidates) {
	  double cosTheta = VertexHelper::CalculateCMCosTheta(*vtxPtr, *track);
	  if (cosTheta < minCosTheta) {
	    minCosTheta = cosTheta;
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
      HYDDRA_DBG("[Leptonic] Disambiguation dropping " << rejectedPtrs.size() << " vertices.\n");

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
  // Keep only 2-track vertices that pass angular and displacement significance
  // cuts.  Charge neutrality is applied when requireChargeNeutrality_ is true.
  void filteringImpl() {
    if (this->empty() || !primaryVertex_ || !primaryVertex_->isValid()) return;

    [[maybe_unused]] const size_t nInput = this->size();
    HYDDRA_DBG("[Leptonic] Final Filtering on " << nInput << " vertices...\n");

    // Cut counters
    size_t nFailSize = 0;
    size_t nFailTrackCosTheta = 0;
    size_t nFailTrackCosThetaCMSlope = 0;
    size_t nFailTrackCosThetaCMLimit = 0;
    size_t nFailCharge = 0;
    size_t nFailVtxCosTheta = 0;
    size_t nFailDxyError = 0;
    size_t nFailDxySignificance = 0;

    TrackVertexSetCollection finalVertices;

    for (const auto& vertex : *this) {
      if (vertex.size() != 2) {
        nFailSize++;
        continue;
      }

      // Track angular cuts
      bool passAngular = true;
      for (const auto& trackRef : vertex) {
        const double trackCosTheta   = vertex.trackCosTheta(*primaryVertex_, trackRef);
        const double trackCosThetaCM = vertex.trackDecayAngleCM(trackRef);

        if (trackCosTheta < minTrackCosTheta_) {
          nFailTrackCosTheta++;
          passAngular = false;
          break;
        }
        if (std::fabs(trackCosThetaCM) > (trackCosThetaCM_Slope_ * trackCosTheta + maxTrackCosThetaCM_Intercept_)) {
          nFailTrackCosThetaCMSlope++;
          passAngular = false;
          break;
        }
        if (std::fabs(trackCosThetaCM) > maxTrackCosThetaCM_Limit_) {
          nFailTrackCosThetaCMLimit++;
          passAngular = false;
          break;
        }
      }
      if (!passAngular) continue;

      // Charge neutrality
      if (requireChargeNeutrality_ && VertexHelper::CalculateTotalCharge(vertex) != 0) {
        nFailCharge++;
        continue;
      }

      // Vertex cos theta (pointing angle wrt PV)
      {
        const double vtxCosTheta = vertex.cosTheta(*primaryVertex_);
        const bool failsVtxCos = useAbsVtxCosTheta_
          ? std::fabs(vtxCosTheta) < minVtxCosTheta_
          : vtxCosTheta < minVtxCosTheta_;
        if (failsVtxCos) {
          nFailVtxCosTheta++;
          continue;
        }
      }

      // Displacement significance
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

    HYDDRA_DBG("[Leptonic] Filtering summary (" << nInput << " -> " << this->size() << "):\n"
               << "  size != 2:              " << nFailSize << " failed\n"
               << "  trackCosTheta:          " << nFailTrackCosTheta << " failed (cut: " << minTrackCosTheta_ << ")\n"
               << "  trackCosThetaCM slope:  " << nFailTrackCosThetaCMSlope << " failed (cut: slope=" << trackCosThetaCM_Slope_ << " intercept=" << maxTrackCosThetaCM_Intercept_ << ")\n"
               << "  trackCosThetaCM limit:  " << nFailTrackCosThetaCMLimit << " failed (cut: " << maxTrackCosThetaCM_Limit_ << ")\n"
               << "  charge neutrality:      " << nFailCharge << " failed" << (requireChargeNeutrality_ ? "" : " (cut disabled)") << "\n"
               << "  vtxCosTheta:            " << nFailVtxCosTheta << " failed (cut: "
               << (useAbsVtxCosTheta_ ? "|cos#theta|" : "cos#theta") << " > " << minVtxCosTheta_ << ")\n"
               << "  dxyError <= 0:          " << nFailDxyError << " failed\n"
               << "  dxy significance:       " << nFailDxySignificance << " failed (cut: " << minDxySignificance_ << ")\n"
               << "  PASSED:                 " << this->size() << "\n");
  }

 private:

  // Cleaning thresholds
  double maxCompatibility_;
  double minCleanCosTheta_;
  bool   useDiagonalCut_;
  double cleanCutSlope_;

  // Final filtering thresholds
  double minTrackCosTheta_;
  double maxTrackCosThetaCM_Limit_;
  double maxTrackCosThetaCM_Intercept_;
  double trackCosThetaCM_Slope_;
  bool   requireChargeNeutrality_;
  double minVtxCosTheta_;
  bool   useAbsVtxCosTheta_;

  // Returns true if the vertex fails kinematic selection.
  // Checked after cleaning and before disambiguation.
  bool failsKinematics(const TrackVertexSet& v) const {
    if (v.vertex().normalizedChi2() > maxNormChi2_) return true;
    if (v.cosTheta(*primaryVertex_) < seedCosThetaCut_) return true;

    auto p4 = VertexHelper::GetVertex4Vector(v);
    double mass = p4.M();
    if (mass < minMass_) return true;

    double p = p4.P();
    double E = std::sqrt(p * p + mass * mass);
    double pOverE = (E > 1e-6) ? p / E : 0.0;
    if (pOverE < minPOverE_) return true;

    return false;
  }
};
