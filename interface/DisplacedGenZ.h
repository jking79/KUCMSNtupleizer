#pragma once

#include <TVector3.h>
#include <Math/Vector4D.h>

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <vector>
#include <cmath>

// Define the decay modes
enum class ZDecayMode { Hadronic, Electron, Muon, Tau, Invisible, Unknown };

class DisplacedGenZ {
 private:

  // --- encapsulated Data ---
  typedef reco::Candidate::LorentzVector LorentzVector;
  const reco::GenParticle* zBoson_;
  std::vector<const pat::PackedGenParticle*> daughters_;
  ZDecayMode decayMode_;
  
  // Cached variables to avoid recomputing
  double lxy_;
  double lxyz_;
  LorentzVector p4SumNeutral_, p4SumTrackable_, p4SumAll_;

  void cleanDaughters();
  
 public:
  // --- Constructor ---
  DisplacedGenZ(const reco::GenParticle* zBoson,
		const std::vector<const pat::PackedGenParticle*>& daughters,
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
  const std::vector<const pat::PackedGenParticle*>& getAllDaughters() const { return daughters_; }
  std::vector<const pat::PackedGenParticle*> getTrackableDaughters() const;
  std::vector<const pat::PackedGenParticle*> getNeutralDaughters() const;
  
  // --- Helper Methods for the Analyzer ---
  bool isLeptonic() const { return (decayMode_ == ZDecayMode::Electron || decayMode_ == ZDecayMode::Muon); }
  bool isHadronic() const { return (decayMode_ == ZDecayMode::Hadronic); }
  bool isTau() const { return decayMode_ == ZDecayMode::Tau; }
  bool isInvisible() const { return decayMode_ == ZDecayMode::Invisible; }

  static std::vector<DisplacedGenZ> build(
      const edm::Handle<std::vector<reco::GenParticle>>& prunedHandle,
      const edm::Handle<std::vector<pat::PackedGenParticle>>& packedHandle);
  
};

