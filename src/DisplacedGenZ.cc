#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DisplacedGenZ.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

void DisplacedGenZ::cleanDaughters() {
  std::vector<const reco::Candidate*> cleaned;
  cleaned.reserve(daughters_.size());
  for(const auto* d : daughters_) {
    if(d->status() == 1)
      cleaned.push_back(d);
  }
  daughters_ = cleaned;
}

DisplacedGenZ::DisplacedGenZ(const reco::GenParticle* zBoson,
			     const std::vector<const reco::Candidate*>& daughters,
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

std::vector<const reco::Candidate*> DisplacedGenZ::getTrackableDaughters() const {
  std::vector<const reco::Candidate*> trackable;
  for(const auto* d : daughters_) {
    if(d->charge() != 0)
      trackable.push_back(d);
  }
  return trackable;
}

std::vector<const reco::Candidate*> DisplacedGenZ::getNeutralDaughters() const {
  std::vector<const reco::Candidate*> neutrals;
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
    const reco::Candidate* gen(pair.GetObjectB());
    const double genPt(gen->pt());
    const double trackPt(track.pt());
    const double relPtDiff((genPt - trackPt) / genPt);
    // Use naive eta/phi deltaR. The helical-corrected cost stored in GetDeltaR()
    // can be corrupted by imprecise gen vertex coordinates (esp. PackedGenParticle),
    // producing deltaR values larger than the naive one and causing spurious failures.
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

  // 2. Determine decay mode for each Z (needed before collecting daughters)
  auto classifyMode = [](const reco::GenParticle* z) -> ZDecayMode {
    if(z->numberOfDaughters() == 0) return ZDecayMode::Unknown;
    int firstDId = std::abs(z->daughter(0)->pdgId());
    if(firstDId == 11)                                    return ZDecayMode::Electron;
    if(firstDId == 13)                                    return ZDecayMode::Muon;
    if(firstDId == 15)                                    return ZDecayMode::Tau;
    if(firstDId >= 1 && firstDId <= 5)                    return ZDecayMode::Hadronic;
    if(firstDId == 12 || firstDId == 14 || firstDId == 16) return ZDecayMode::Invisible;
    return ZDecayMode::Unknown;
  };

  // 3. Collect daughters, routing by collection:
  //    - Leptonic Z (e, mu, tau): search the PRUNED collection for status=1 charged
  //      particles whose mother chain leads back to this signal Z. This correctly finds
  //      the final-state lepton after FSR, which may not be a direct daughter of the Z.
  //    - Hadronic Z: final-state pions/kaons live in the PACKED collection; find them by
  //      tracing each packed particle's mother chain back to the signal Z.
  //    Both cases use the same mother-chain pattern, just in different collections.
  std::map<const reco::GenParticle*, std::vector<const reco::Candidate*>> daughterMap;
  for(const auto* z : signalZs)
    daughterMap[z] = {};

  // Helper: walk a particle's mother chain; if a signal Z is found, return its pointer.
  auto getSignalZAncestor = [&](const reco::Candidate* mom) -> const reco::GenParticle* {
    while(mom) {
      if(mom->pdgId() == 23) {
        const reco::GenParticle* zCand = dynamic_cast<const reco::GenParticle*>(mom);
        return (zCand && daughterMap.count(zCand)) ? zCand : nullptr;
      }
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
    return nullptr;
  };

  // DEBUG: dump every pruned gen particle and whether it passes the leptonic filter
  std::cout << "[DisplacedGenZ::build] DEBUG --- pruned collection scan ---\n";
  for(const auto& gp : *prunedHandle) {
    std::cout << "  pruned pdgId=" << gp.pdgId()
              << " status=" << gp.status()
              << " charge=" << gp.charge()
              << " pt=" << gp.pt()
              << " isLastCopy=" << gp.isLastCopy()
              << "\n";
  }

  // Leptonic daughters: status=1 charged particles in the pruned collection
  for(const auto& gp : *prunedHandle) {
    if(gp.status() != 1 || gp.charge() == 0) continue;
    if(const auto* z = getSignalZAncestor(gp.mother(0)))
      if(classifyMode(z) != ZDecayMode::Hadronic)
        daughterMap[z].push_back(&gp);
  }

  // Hadronic daughters: all particles in the packed collection
  for(const auto& gp : *packedHandle) {
    if(const auto* z = getSignalZAncestor(gp.mother(0)))
      if(classifyMode(z) == ZDecayMode::Hadronic)
        daughterMap[z].push_back(&gp);
  }

  // 4. Instantiate one DisplacedGenZ per signal Z
  for(const auto* z : signalZs) {
    ZDecayMode mode = classifyMode(z);
    const auto& daus = daughterMap[z];
    std::cout << "[DisplacedGenZ::build] DEBUG signal Z pdgId=23"
              << " mode=" << int(mode)
              << " nDaughters=" << daus.size()
              << " Lxy=" << std::hypot(z->vx(), z->vy()) << "\n";
    genZs.emplace_back(z, daus, mode);
  }

  return genZs;
}
