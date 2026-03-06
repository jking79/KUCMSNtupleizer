#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DisplacedGenZ.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

void DisplacedGenZ::cleanDaughters() {
  std::vector<const pat::PackedGenParticle*> cleaned;
  cleaned.reserve(daughters_.size());
  for(const auto* d : daughters_) {
    if(d->status() == 1)
      cleaned.push_back(d);
  }
  daughters_ = cleaned;
}

DisplacedGenZ::DisplacedGenZ(const reco::GenParticle* zBoson,
			     const std::vector<const pat::PackedGenParticle*>& daughters,
			     ZDecayMode mode)
  : zBoson_(zBoson), daughters_(daughters), decayMode_(mode)
{
  cleanDaughters();

  lxy_  = std::hypot(zBoson_->vx(), zBoson_->vy());
  lxyz_ = std::sqrt(zBoson_->vx()*zBoson_->vx() +
                    zBoson_->vy()*zBoson_->vy() +
                    zBoson_->vz()*zBoson_->vz());

  p4SumNeutral_   = LorentzVector(0,0,0,0);
  p4SumTrackable_ = LorentzVector(0,0,0,0);
  p4SumAll_       = LorentzVector(0,0,0,0);

  for(const auto* d : daughters_) {
    p4SumAll_ += d->p4();
    if(d->charge() != 0)
      p4SumTrackable_ += d->p4();
    else
      p4SumNeutral_ += d->p4();
  }
}

std::vector<const pat::PackedGenParticle*> DisplacedGenZ::getTrackableDaughters() const {
  std::vector<const pat::PackedGenParticle*> trackable;
  for(const auto* d : daughters_) {
    if(d->charge() != 0)
      trackable.push_back(d);
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

reco::TrackCollection DisplacedGenZ::getTracks() const {
  reco::TrackCollection tracks;
  tracks.reserve(matchedTracks_.size());
  for(const auto& pair : matchedTracks_)
    tracks.emplace_back(pair.GetObjectA().track());
  return tracks;
}

int DisplacedGenZ::nMatchedInVertex(const reco::Vertex& vtx) const {
  int count(0);
  for(const auto& pair : matchedTracks_) {
    if(VertexHelper::isInVertex(vtx, pair.GetObjectA().track()))
      count++;
  }
  return count;
}

bool DisplacedGenZ::isBronze(const reco::Vertex& vtx) const {
  if(matchedTracks_.empty()) return false;
  const int nInVtx = nMatchedInVertex(vtx);
  return nInVtx > 0 && nInVtx < (int)matchedTracks_.size();
}

bool DisplacedGenZ::isSilver(const reco::Vertex& vtx) const {
  if(matchedTracks_.empty()) return false;
  const int nInVtx = nMatchedInVertex(vtx);
  return nInVtx == (int)matchedTracks_.size() && (int)vtx.tracksSize() > (int)matchedTracks_.size();
}

bool DisplacedGenZ::isGold(const reco::Vertex& vtx) const {
  if(matchedTracks_.empty()) return false;
  const int nInVtx = nMatchedInVertex(vtx);
  return nInVtx == (int)matchedTracks_.size() && (int)vtx.tracksSize() == (int)matchedTracks_.size();
}

double DisplacedGenZ::matchRatio(const reco::Vertex& vtx, double threshold) const {
  if(matchedTracks_.empty() || vtx.tracksSize() == 0) return -1.;

  int count(0);
  for(const auto& pair : matchedTracks_) {
    const reco::Track track(pair.GetObjectA().track());
    const double genPt(pair.GetObjectB()->pt());
    const double trackPt(track.pt());
    const double relPtDiff((genPt - trackPt) / genPt);
    // Recompute naive deltaR (track direction at beamline PCA vs gen direction).
    // The helical-corrected cost stored in GetDeltaR() can be corrupted by imprecise
    // PackedGenParticle vertex coordinates, producing deltaR values larger than the
    // naive one and causing the quality cut to fail spuriously.
    const pat::PackedGenParticle* gen(pair.GetObjectB());
    const double dEta(track.eta() - gen->eta()), dPhi(track.phi() - gen->phi());
    const double deltaR(std::sqrt(dEta*dEta + dPhi*dPhi));

    if(!VertexHelper::isInVertex(vtx, track)) continue;

    if((deltaR < 0.01 && std::fabs(relPtDiff) < (-10*deltaR + 0.15)) ||
       (deltaR > 0.01 && deltaR < threshold && std::fabs(relPtDiff) < (-deltaR/2 + 0.055))) {
      count++;
    }
  }
  return double(count) / vtx.tracksSize();
}

double DisplacedGenZ::distance3D(const reco::Vertex& vtx) const {
  const double dx(x() - vtx.x()), dy(y() - vtx.y()), dz(z() - vtx.z());
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// --- Factory Method ---
std::vector<DisplacedGenZ> DisplacedGenZ::build(
    const edm::Handle<std::vector<reco::GenParticle>>& prunedHandle,
    const edm::Handle<std::vector<pat::PackedGenParticle>>& packedHandle)
{
  std::vector<DisplacedGenZ> genZs;
  genZs.reserve(2);

  std::vector<const reco::GenParticle*> signalZs;

  // 1. Find signal Zs (last copy, descended from N2) in pruned collection
  for(const auto& gp : *prunedHandle) {
    if(gp.pdgId() == 23 && gp.isLastCopy()) {
      const reco::Candidate* mom = gp.mother(0);
      bool isSignal = false;
      while(mom) {
        if(std::abs(mom->pdgId()) == 1000023) { isSignal = true; break; }
        mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
      }
      if(isSignal) signalZs.push_back(&gp);
    }
  }

  if(signalZs.empty()) return genZs;

  // 2. Map packed daughters to their parent Z
  std::map<const reco::GenParticle*, std::vector<const pat::PackedGenParticle*>> daughterMap;
  for(const auto* z : signalZs)
    daughterMap[z] = {};

  for(const auto& gp : *packedHandle) {
    const reco::Candidate* mom = gp.mother(0);
    while(mom) {
      if(mom->pdgId() == 23) {
        const reco::GenParticle* zCand = dynamic_cast<const reco::GenParticle*>(mom);
        if(daughterMap.find(zCand) != daughterMap.end())
          daughterMap[zCand].push_back(&gp);
        break;
      }
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
  }

  // 3. Classify decay mode and instantiate
  for(const auto* z : signalZs) {
    ZDecayMode mode = ZDecayMode::Unknown;
    if(z->numberOfDaughters() > 0) {
      int firstDId = std::abs(z->daughter(0)->pdgId());
      if(firstDId == 11)                          mode = ZDecayMode::Electron;
      else if(firstDId == 13)                     mode = ZDecayMode::Muon;
      else if(firstDId == 15)                     mode = ZDecayMode::Tau;
      else if(firstDId >= 1 && firstDId <= 5)     mode = ZDecayMode::Hadronic;
      else if(firstDId == 12 || firstDId == 14 || firstDId == 16) mode = ZDecayMode::Invisible;
    }
    genZs.emplace_back(z, daughterMap[z], mode);
  }

  return genZs;
}
