#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DisplacedGenZ.h"

void DisplacedGenZ::cleanDaughters() {
    std::vector<const pat::PackedGenParticle*> cleaned;
    cleaned.reserve(daughters_.size());

    for(const auto* d : daughters_) {
	if(d->status() == 1) {
            cleaned.push_back(d);
	}
    }
    daughters_ = cleaned;
}

DisplacedGenZ::DisplacedGenZ(const reco::GenParticle* zBoson,
              const std::vector<const pat::PackedGenParticle*>& daughters,
              ZDecayMode mode)
  : zBoson_(zBoson), daughters_(daughters), decayMode_(mode)
{

  cleanDaughters();
  
  // Cache the displacement at construction
  lxy_ = std::hypot(zBoson_->vx(), zBoson_->vy());
  lxyz_ = std::sqrt(zBoson_->vx()*zBoson_->vx() +
                    zBoson_->vy()*zBoson_->vy() +
                    zBoson_->vz()*zBoson_->vz());

  p4SumNeutral_ = LorentzVector(0,0,0,0);
  p4SumTrackable_ = LorentzVector(0,0,0,0);
  p4SumAll_ = LorentzVector(0,0,0,0);

   for(const auto* d : daughters_) {
    p4SumAll_ += d->p4();
    // Enforce CMS tracking acceptance for the trackable 4-vector
    if(d->charge() != 0 && d->pt() > 0.5 && std::abs(d->eta()) < 2.5)
      p4SumTrackable_ += d->p4();
    else
      p4SumNeutral_ += d->p4();
  }
}

std::vector<const pat::PackedGenParticle*> DisplacedGenZ::getTrackableDaughters() const {
  std::vector<const pat::PackedGenParticle*> trackable;
  for(const auto* d : daughters_) {
    if(d->charge() != 0 && d->pt() > 0.5 && std::abs(d->eta()) < 2.5) {
      trackable.push_back(d);
    }
  }
  return trackable;
}

std::vector<const pat::PackedGenParticle*> DisplacedGenZ::getNeutralDaughters() const {
  std::vector<const pat::PackedGenParticle*> neutrals;
  for(const auto* d : daughters_) {
    if(d->charge() == 0) neutrals.push_back(d);
  }
  return neutrals;
}

// --- Factory Method Implementation ---
std::vector<DisplacedGenZ> DisplacedGenZ::build(
    const edm::Handle<std::vector<reco::GenParticle>>& prunedHandle,
    const edm::Handle<std::vector<pat::PackedGenParticle>>& packedHandle) 
{
  std::vector<DisplacedGenZ> genZs;
  genZs.reserve(2); 
  
  std::vector<const reco::GenParticle*> signalZs;
  
  // 1. Find Signal Zs in Pruned
  for (const auto& gp : *prunedHandle) {
    if (gp.pdgId() == 23 && gp.isLastCopy()) {
      const reco::Candidate* mom = gp.mother(0);
      bool isSignal = false;
      while(mom) {
	if (std::abs(mom->pdgId()) == 1000023) { // Check for N2
	  isSignal = true;
	  break;
	}
	mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
      }
      if (isSignal) signalZs.push_back(&gp);
    }
  }
  
  if (signalZs.empty()) return genZs; 
  
  // 2. Map packed daughters to their parent Z
  std::map<const reco::GenParticle*, std::vector<const pat::PackedGenParticle*>> daughterMap;
  for (const auto* z : signalZs) {
    daughterMap[z] = std::vector<const pat::PackedGenParticle*>();
  }
  
  for (const auto& gp : *packedHandle) {
    const reco::Candidate* mom = gp.mother(0);
    while(mom) {
      if (mom->pdgId() == 23) {
	const reco::GenParticle* zCand = dynamic_cast<const reco::GenParticle*>(mom);
	if (daughterMap.find(zCand) != daughterMap.end()) {
	  daughterMap[zCand].push_back(&gp);
	}
	break;
      }
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
  }
  
  // 3. Classify decay mode and instantiate
  for (const auto* z : signalZs) {
    ZDecayMode mode = ZDecayMode::Unknown;
    if (z->numberOfDaughters() > 0) {
      int firstDId = std::abs(z->daughter(0)->pdgId());
      if (firstDId == 11) mode = ZDecayMode::Electron;
      else if (firstDId == 13) mode = ZDecayMode::Muon;
      else if (firstDId == 15) mode = ZDecayMode::Tau;
      else if (firstDId >= 1 && firstDId <= 5) mode = ZDecayMode::Hadronic;
      else if (firstDId == 12 || firstDId == 14 || firstDId == 16) mode = ZDecayMode::Invisible;
    }
    genZs.emplace_back(z, daughterMap[z], mode);
  }
  
  return genZs;
}
