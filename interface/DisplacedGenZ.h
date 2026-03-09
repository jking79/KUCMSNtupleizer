#pragma once

#include <TVector3.h>
#include <Math/Vector4D.h>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/PairedObjects.h"

#include <vector>
#include <cmath>

// Define the decay modes
enum class ZDecayMode { Hadronic, Electron, Muon, Tau, Invisible, Unknown };

class DisplacedGenZ {
 public:
  // Both reco::GenParticle (pruned, leptonic) and pat::PackedGenParticle (packed, hadronic)
  // inherit from reco::Candidate, so one unified daughter/match type covers both cases.
  using ZMatchPairs = PairedObjectCollection<reco::TransientTrack, const reco::Candidate*>;

 private:
  typedef reco::Candidate::LorentzVector LorentzVector;
  const reco::GenParticle* zBoson_;
  // Unified daughter list: pruned GenParticle* for leptonic Z, PackedGenParticle* for hadronic Z
  std::vector<const reco::Candidate*> daughters_;
  ZDecayMode decayMode_;

  double lxy_;
  double lxyz_;
  LorentzVector p4SumNeutral_, p4SumTrackable_, p4SumAll_;

  ZMatchPairs matchedTracks_;

  void cleanDaughters();
  int nMatchedInVertex(const reco::Vertex& vtx) const;

 public:
  // --- Constructor ---
  DisplacedGenZ(const reco::GenParticle* zBoson,
		const std::vector<const reco::Candidate*>& daughters,
		ZDecayMode mode);

  // --- Public Getters ---
  double x() const { return zBoson_->vx(); }
  double y() const { return zBoson_->vy(); }
  double z() const { return zBoson_->vz(); }
  double p() const { return p4SumTrackable_.P(); }
  double px() const { return p4SumTrackable_.px(); }
  double py() const { return p4SumTrackable_.py(); }
  double pz() const { return p4SumTrackable_.pz(); }
  double pt() const { return p4SumTrackable_.pt(); }
  double eta() const { return p4SumTrackable_.eta(); }
  double phi() const { return p4SumTrackable_.phi(); }
  double mass() const { return p4SumTrackable_.M(); }
  double fullMass() const { return p4SumAll_.M(); }

  double Lxy() const { return lxy_; }
  double Lxyz() const { return lxyz_; }

  ZDecayMode mode() const { return decayMode_; }
  const std::vector<const reco::Candidate*>& getAllDaughters() const { return daughters_; }
  std::vector<const reco::Candidate*> getTrackableDaughters() const;
  std::vector<const reco::Candidate*> getNeutralDaughters() const;

  // --- Match results (populated after per-Z Hungarian matching) ---
  void setMatches(const ZMatchPairs& matches) { matchedTracks_ = matches; }
  bool hasTracks() const { return !matchedTracks_.empty(); }
  const ZMatchPairs& getMatches() const { return matchedTracks_; }
  reco::TrackCollection getTracks() const;

  // --- Vertex quality (mirrors GenVertex::isBronze/isSilver/isGold) ---
  // Gold:   all matched tracks are in the vertex and it has exactly that many tracks
  // Silver: all matched tracks are in the vertex but it has extra tracks
  // Bronze: at least one but not all matched tracks are in the vertex
  bool isBronze(const reco::Vertex& vtx) const;
  bool isSilver(const reco::Vertex& vtx) const;
  bool isGold(const reco::Vertex& vtx) const;

  // Fraction of matched tracks in the vertex passing the deltaR+relPt quality cuts
  double matchRatio(const reco::Vertex& vtx, double threshold = 0.1) const;

  // Raw 3D Euclidean distance between this gen vertex position and a reco vertex
  double distance3D(const reco::Vertex& vtx) const;

  // --- Mode helpers ---
  bool isLeptonic() const { return (decayMode_ == ZDecayMode::Electron || decayMode_ == ZDecayMode::Muon); }
  bool isElectron() const { return decayMode_ == ZDecayMode::Electron; }
  bool isMuon() const { return decayMode_ == ZDecayMode::Muon; }
  bool isHadronic() const { return (decayMode_ == ZDecayMode::Hadronic); }
  bool isTau() const { return decayMode_ == ZDecayMode::Tau; }
  bool isInvisible() const { return decayMode_ == ZDecayMode::Invisible; }

  // Factory: searches both the pruned collection (leptonic daughters) and the packed
  // collection (hadronic daughters) so every signal Z gets the right daughter list.
  static std::vector<DisplacedGenZ> build(
      const edm::Handle<std::vector<reco::GenParticle>>& prunedHandle,
      const edm::Handle<std::vector<pat::PackedGenParticle>>& packedHandle);

};
