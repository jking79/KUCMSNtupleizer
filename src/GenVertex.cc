#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"

GenVertex::GenVertex(const GenMatches &matchedPairs) :
  isGenElectron_(false),
  isGenMuon_(false),
  isGenHadronic_(false),
  hasMatch_(true),
  size_(matchedPairs.size())
{
 
  genMatches_ = matchedPairs;
  LorentzVec genVertex4Vec;

  vertexCharge_ = 0;
  for(const auto &pair : matchedPairs) {
    reco::GenParticle gen(pair.GetObjectB());

    const double pt(gen.pt());
    const double eta(gen.eta());
    const double phi(gen.phi());
    const double mass(0.);

    LorentzVec gen4Vec(pt, eta, phi, mass);
    genVertex4Vec += gen4Vec;
    vertexCharge_ += gen.charge();
  }

  const reco::GenParticle genParticle0(matchedPairs[0].GetObjectB());
  const reco::GenParticle genParticle1(matchedPairs[1].GetObjectB());

  setMembers(genVertex4Vec, genParticle0);
  genPair_ = std::make_pair(genParticle0, genParticle1);
  
}

GenVertex::GenVertex(const std::pair<reco::GenParticle, reco::GenParticle> &genPair) :
  isGenElectron_(false),
  isGenMuon_(false),
  isGenHadronic_(false),
  hasMatch_(false),
  size_(2){

  reco::GenParticle gen1(genPair.first);
  reco::GenParticle gen2(genPair.second);
  
  if(abs(gen1.pdgId()) != abs(gen2.pdgId()))
    throw std::runtime_error("Attempted to group particles with different PDG IDs: " + 
                           std::to_string(gen1.pdgId()) + " and " + 
                           std::to_string(gen2.pdgId()));
  
  LorentzVec genVertex4Vec = LorentzVec(gen1.pt(), gen1.eta(), gen1.phi(), 0.) + LorentzVec(gen2.pt(), gen2.eta(), gen2.phi(), 0.);
  vertexCharge_ = gen1.charge() + gen2.charge();

  setMembers(genVertex4Vec, gen1);
  genPair_ = genPair;
}

double GenVertex::trackDz() const {
  if(size_ != 2 || !this->hasTracks()) return -1.;

  reco::Track track1(this->tracks()[0]);
  reco::Track track2(this->tracks()[1]);

  return fabs(track1.vz() - track2.vz()); 
}

double GenVertex::decayAngle() const {
  if(size_ != 2 || !this->hasTracks()) return -1.;

  reco::Track track1(this->tracks()[0]);
  reco::Track track2(this->tracks()[1]);
  
  TVector3 direction1(track1.vx(), track1.vy(), track1.vz());
  TVector3 direction2(track2.vx(), track2.vy(), track2.vz());

  return direction1.Angle(direction2);
}

double GenVertex::momentumAngle() const {
  if(size_ != 2 || !this->hasTracks()) return -1.;

  reco::Track track1(this->tracks()[0]);
  reco::Track track2(this->tracks()[1]);

  TVector3 direction1(track1.px(), track1.py(), track1.pz());
  TVector3 direction2(track2.px(), track2.py(), track2.pz());

  return direction1.Angle(direction2);
}

double GenVertex::distance3D(const reco::Vertex &vertex) const {
  
  const double gx(this->x()), gy(this->y()), gz(this->z());
  const double vx(vertex.x()), vy(vertex.y()), vz(vertex.z());
  const double dx(gx-vx), dy(gy-vy), dz(gz-vz);
  return sqrt(dx*dx + dy*dy + dz*dz);
}

int GenVertex::nMatched(const reco::Vertex &vertex) const {
  
  int sharedTracks(0);
  //std::cout << "Tracks in gen vertex: " << std::endl;
  for(const reco::Track &track : this->tracks()) {
    //std::cout << "  pt = " << track.pt() << ", eta = " << track.eta() << std::endl;
    if(VertexHelper::isInVertex(vertex, track)) sharedTracks++;
  }

  return sharedTracks;
}

GenVertex::VertexQuality GenVertex::DetermineQuality(const reco::Vertex &vertex) const {
  if (size_ == 0) return VertexQuality::kEmpty;
  
  int matches = nMatched(vertex);
  int reconstructedSize = static_cast<int>(vertex.tracksSize());
  
  if (matches == 0) return VertexQuality::kNoMatch;
  if (matches == size_ && reconstructedSize == size_) return VertexQuality::kGold;
  if (matches == size_ && reconstructedSize > size_) return VertexQuality::kSilver;
  if (matches > 0 && matches < size_) return VertexQuality::kBronze;
  
  // Fallback case (should not occur if logic is correct)
  std::cerr << "Warning: Unexpected case in determineQuality. "
	    << "Original size: " << size_
	    << ", Reconstructed size: " << reconstructedSize
	    << ", Matches: " << matches << std::endl;
  return VertexQuality::kBronze;
}

int GenVertex::CountInstances(const reco::TrackCollection &tracks) const {

  int instances(0);

  for(const auto &track1 : this->tracks())
    for(const auto &track2: tracks)
      if(track1.pt() == track2.pt() && track1.eta() == track2.eta())
        instances++;

  return instances;
}

TransientVertex GenVertex::fitTracks(const TransientTrackBuilder* ttBuilder) const {

  if(!this->hasTracks())
    return TransientVertex();
  
  std::vector<reco::TransientTrack> ttracks;
  for(const auto &track : this->tracks())
    ttracks.emplace_back(ttBuilder->build(track));

  const KalmanVertexFitter fitter;
  const TransientVertex newVertex(fitter.vertex(ttracks));

  return newVertex;
}

reco::Vertex GenVertex::fitTracks(const TransientTrackBuilder* ttBuilder, const edm::Handle<reco::TrackCollection> &trackHandle) const {

  TransientVertex thisVertex(fitTracks(ttBuilder));
  if(!thisVertex.isValid())
    return reco::Vertex();
  
  reco::Vertex recoVtx(VertexHelper::ConvertFitVertex(thisVertex));

  for(const auto &track : thisVertex.originalTracks()) {
    reco::TrackRef trackRef = TrackHelper::GetTrackRef(track.track(), trackHandle);
    recoVtx.add(reco::TrackBaseRef(trackRef), recoVtx.trackWeight(reco::TrackBaseRef(trackRef)));
  }

  return recoVtx;
}

reco::GenParticleCollection GenVertex::getStableDaughters(const reco::GenParticleCollection &genParticles) const {

  reco::GenParticleCollection stableDaughters;
  
  auto genZ = genPair_.first.mother();
  for(const auto &gen : genParticles) {
    if(gen.status() != 1)
      continue;

    std::vector<int> motherChainIDs(MomIDs(gen));
    if(std::find(motherChainIDs.begin(), motherChainIDs.end(), 23) != motherChainIDs.end()) {
      auto mother(gen.mother());
      while(mother->pdgId() != 23)
	mother = mother->mother();
      
      if(mother == genZ)
	stableDaughters.emplace_back(gen);
    }
  }
  return stableDaughters;
}

reco::GenParticleCollection GenVertex::getStableChargedDaughters(const reco::GenParticleCollection &genParticles) const {

  reco::GenParticleCollection stableChargedDaughters;
  for(const auto &gen : getStableDaughters(genParticles)) {
    if(gen.charge() != 0)
      stableChargedDaughters.emplace_back(gen);
  }
  return stableChargedDaughters;
}

reco::GenParticleCollection GenVertex::getStableNeutralDaughters(const reco::GenParticleCollection &genParticles) const {

  reco::GenParticleCollection stableNeutralDaughters;
  for(const auto &gen : getStableDaughters(genParticles)) {
    if(gen.charge() == 0)
      stableNeutralDaughters.emplace_back(gen);
  }
  return stableNeutralDaughters;
}

void GenVertex::setMembers(const LorentzVec &genVertex4Vec, const reco::GenParticle &gen) {

  vertexP_ = genVertex4Vec.P();
  vertexPx_ = genVertex4Vec.px();
  vertexPy_ = genVertex4Vec.py();
  vertexPz_ = genVertex4Vec.pz();
  vertexPt_ = genVertex4Vec.pt();
  vertexEta_ = genVertex4Vec.eta();
  vertexPhi_ = genVertex4Vec.phi();
  vertexMass_ = genVertex4Vec.M();
  vertexX_ = gen.vx();
  vertexY_ = gen.vy();
  vertexZ_ = gen.vz();
  vertexDxy_ = sqrt(vertexX_*vertexX_ + vertexY_*vertexY_);

  if(abs(gen.pdgId()) == 11) isGenElectron_ = true;
  else if (abs(gen.pdgId()) == 13) isGenMuon_ = true;
  else if (isSignalGenJet(gen)) isGenHadronic_ = true;
}

//GenVertices class
GenVertices::GenVertices(const GenMatches &matchedPairs, const double deltaRCut) {
  std::vector<GenMatches> quasiVertexVector = findGenParticlePairs(matchedPairs, deltaRCut);

  for(const auto &quasiVertex : quasiVertexVector) {

    GenVertex vertex(quasiVertex);
    this->emplace_back(vertex);
  }
}

GenVertices::GenVertices(const std::vector<reco::GenParticle> &genParticles) {

  // Group signal particles by their LLP (SUSY) ancestor pointer rather than
  // by production vertex coordinates.  In MiniAOD (prunedGenParticles) the
  // vx/vy/vz of LLP daughters is unreliable and often collapses to (0,0,0),
  // which caused all four signal leptons in the event to be grouped into one
  // spurious vertex.  The ancestor pointer is stable across pruning and
  // correctly separates decays from different LLPs.
  std::unordered_map<const reco::Candidate*, std::vector<reco::GenParticle>> ancestorMap;

  int nSignalElectrons(0), nSignalMuons(0), nSignalJets(0), nNoAncestor(0);
  for(const auto &gen : genParticles) {
    bool isElectron(isSignalGenElectron(gen));
    bool isMuon(isSignalGenMuon(gen));
    bool isJet(isSignalGenJet(gen));
    bool isSignal(isElectron || isMuon || isJet);
    if(isElectron) nSignalElectrons++;
    if(isMuon) nSignalMuons++;
    if(isJet) nSignalJets++;
    if(!isSignal) continue;

    const reco::Candidate* llp = findSUSYAncestor(gen);
    if(!llp) {
      nNoAncestor++;
      // Print the full mother chain for diagnosis
      std::cerr << "[GenVertices] SKIP: signal particle (pdgId=" << gen.pdgId()
                << ", status=" << gen.status()
                << ", pt=" << gen.pt()
                << ") has no SUSY ancestor. Mother chain: ";
      const reco::Candidate* mom = (gen.numberOfMothers() > 0) ? gen.mother(0) : nullptr;
      while(mom) {
        std::cerr << mom->pdgId() << "(st=" << mom->status() << ") ";
        mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
      }
      std::cerr << std::endl;
      continue;
    }

    // Group by the direct Z parent (if present) rather than the N2 ancestor.
    // This correctly splits cases where multiple Z bosons from the same N2
    // (or from two N2s that share an ancestor in the pruned collection) would
    // otherwise pile 4+ daughters into one map entry.
    const reco::Candidate* groupKey = llp;
    {
      const reco::Candidate* mom = (gen.numberOfMothers() > 0) ? gen.mother(0) : nullptr;
      while(mom) {
        if(abs(mom->pdgId()) == 23) { groupKey = mom; break; }
        if(abs(mom->pdgId()) > 1000000) break; // stop at LLP boundary
        mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
      }
    }
    ancestorMap[groupKey].push_back(gen);
  }

  std::cerr << "[GenVertices] Signal particles found: "
            << nSignalElectrons << " e, "
            << nSignalMuons << " mu, "
            << nSignalJets << " jets, "
            << nNoAncestor << " skipped (no SUSY ancestor)"
            << std::endl;

  for(const auto &entry : ancestorMap) {
    if(entry.second.size() != 2) {
      std::cerr << "[GenVertices] SKIP: LLP (pdgId=" << entry.first->pdgId()
                << ") grouped " << entry.second.size()
                << " particles (expected 2):";
      for(const auto &p : entry.second)
        std::cerr << " pdgId=" << p.pdgId() << " status=" << p.status() << " pt=" << p.pt();
      std::cerr << std::endl;
      continue;
    }

    GenVertex vtx(std::make_pair(entry.second[0], entry.second[1]));

    // In MiniAOD, lepton vx/vy/vz is often (0,0,0) for LLP daughters.
    // Recover the LLP decay vertex by finding any daughter of the ancestor
    // in the collection that carries a valid (non-zero) production vertex.
    if(vtx.x() == 0 && vtx.y() == 0 && vtx.z() == 0) {
      for(const auto &gen : genParticles) {
        bool isDaughter = false;
        for(size_t m = 0; m < gen.numberOfMothers(); ++m) {
          if(gen.mother(m) == entry.first) { isDaughter = true; break; }
        }
        if(isDaughter && (gen.vx() != 0 || gen.vy() != 0 || gen.vz() != 0)) {
          vtx.setVertex(gen.vx(), gen.vy(), gen.vz());
          break;
        }
      }
    }

    this->emplace_back(vtx);
  }

  std::cerr << "[GenVertices] Built " << this->size() << " gen vertices ("
            << ancestorMap.size() << " LLP ancestors found)" << std::endl;
}

bool GenVertices::contains(const GenVertex& genVertex) const {

  for(const auto &vertex : *this) {

    if(genVertex == vertex && vertex.hasTracks())
      return true;
  }

  return false;
}

std::vector<GenMatches> GenVertices::findGenParticlePairs(const GenMatches &matchedPairs, const double deltaRCut) const {

  // Group matched track-gen pairs by LLP ancestor pointer when possible
  // (MiniAOD compatible), falling back to production-vertex string for
  // non-SUSY signals (e.g. Z→ll in AOD where positions are reliable).
  std::unordered_map<const reco::Candidate*, GenMatches> ancestorMap;
  std::unordered_map<std::string, GenMatches>            positionMap;

  for(const auto& pair : FindSignalGenCollection(matchedPairs, deltaRCut)) {
    const reco::GenParticle particle(pair.GetObjectB());
    const reco::Candidate* llp = findSUSYAncestor(particle);
    if(llp) {
      ancestorMap[llp].push_back(pair);
    } else {
      const std::string key = std::to_string(particle.vx()) + "_" +
                              std::to_string(particle.vy()) + "_" +
                              std::to_string(particle.vz());
      positionMap[key].push_back(pair);
    }
  }

  std::vector<GenMatches> pairs;
  for(const auto& entry : ancestorMap)
    if(entry.second.size() == 2) pairs.push_back(entry.second);
  for(const auto& entry : positionMap)
    if(entry.second.size() == 2) pairs.push_back(entry.second);
  return pairs;
}

GenMatches GenVertices::FindSignalGenCollection(const GenMatches &matchedPairs, const double deltaRCut) const {

  GenMatches genCollection;

  for(const auto &pair : matchedPairs) {
    const reco::Track track(pair.GetObjectA());
    const reco::GenParticle gen(pair.GetObjectB());

    const bool isSignalElectron(isSignalGenElectron(gen));
    const bool isSignalMuon(isSignalGenMuon(gen));
    const bool isSignal = isSignalElectron || isSignalMuon ;

    if(isSignal && pair.GetDeltaR() < deltaRCut)
      genCollection.emplace_back(pair);
  }

  return genCollection;
}

std::vector<TransientVertex> GenVertices::getValidVertices(const TransientTrackBuilder* ttBuilder, const double chi2Cut) const {

  std::vector<TransientVertex> fitVertices;
  for(const auto &genVtx : *this) {
    TransientVertex tmpVertex(genVtx.fitTracks(ttBuilder));
    if(tmpVertex.isValid() && tmpVertex.normalisedChiSquared() < chi2Cut)
      fitVertices.emplace_back(tmpVertex);
  }
  return fitVertices;
}

reco::GenParticleCollection GenVertices::getAllGenParticles() const {

  reco::GenParticleCollection genParticles;
  for(const auto &genVertex : *this) {
    if(!genVertex.isGenHadronic()) {
      std::pair<reco::GenParticle, reco::GenParticle> pair(genVertex.genPair());
      genParticles.emplace_back(pair.first);
      genParticles.emplace_back(pair.second);
    }
  }
  return genParticles;
}

// Returns true if two gen pairs represent the same two particles (in either
// ordering).  Exact float comparison is safe here because the particles in
// the track-matched vertices are bitwise copies of those in the gen-only
// vertices, all originating from the same prunedGenParticles collection.
static bool sameGenPair(
    const std::pair<reco::GenParticle,reco::GenParticle>& a,
    const std::pair<reco::GenParticle,reco::GenParticle>& b) {
  auto kinMatch = [](const reco::GenParticle& p, const reco::GenParticle& q) {
    return p.pt() == q.pt() && p.eta() == q.eta() && p.phi() == q.phi();
  };
  return (kinMatch(a.first, b.first)  && kinMatch(a.second, b.second)) ||
         (kinMatch(a.first, b.second) && kinMatch(a.second, b.first));
}

void GenVertices::operator+=(const GenVertices &other) {
  // Merge `other` into this collection.  For each incoming vertex:
  //   • Match by spatial position (operator==, x/y/z) — works in AOD.
  //   • Fall back to genPair kinematics — works in MiniAOD where the
  //     track-matched vertex carries the original (0,0,0) coordinates.
  // When a match is found and the incoming vertex has tracks, it replaces
  // the existing one; the recovered position is copied over if needed.
  for(const auto &incoming : other) {
    bool found = false;
    for(auto &existing : *this) {
      if((existing == incoming) || sameGenPair(existing.genPair(), incoming.genPair())) {
        if(incoming.hasTracks()) {
          GenVertex enriched = incoming;
          // Preserve the gen-only vertex's recovered position when the
          // track-matched vertex still has (0,0,0) coordinates (MiniAOD).
          if(enriched.x() == 0 && enriched.y() == 0 && enriched.z() == 0)
            enriched.setVertex(existing.x(), existing.y(), existing.z());
          existing = enriched;
        }
        found = true;
        break;
      }
    }
    if(!found) this->push_back(incoming);
  }
}

reco::GenParticleCollection getStableChargedDaughtersFromPacked(
    const GenVertex& genVertex,
    const std::vector<pat::PackedGenParticle>& packedGenParticles) {

  // Walk up the first gen particle's mother chain to find the decay parent:
  // either a Z boson (pdgId 23, prompt signal) or a SUSY LLP (|pdgId| > 1e6).
  // Using pointer identity is robust and works for both AOD and MiniAOD.
  const reco::Candidate* parent = nullptr;
  {
    const reco::Candidate* mom =
        (genVertex.genPair().first.numberOfMothers() > 0)
        ? genVertex.genPair().first.mother(0) : nullptr;
    while(mom) {
      const int id = abs(mom->pdgId());
      if(id == 23 || id > 1000000) { parent = mom; break; }
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
  }
  if(!parent) return {};

  reco::GenParticleCollection result;
  for(const auto& packed : packedGenParticles) {
    if(packed.status() != 1 || packed.charge() == 0) continue;
    const reco::Candidate* mom =
        (packed.numberOfMothers() > 0) ? packed.mother(0) : nullptr;
    const reco::Candidate* prev = nullptr;
    while(mom && mom != prev) {
      if(mom == parent) {
        result.emplace_back(reco::GenParticle(
            packed.charge(), packed.p4(), packed.vertex(),
            packed.pdgId(), packed.status(), true));
        break;
      }
      prev = mom;
      mom = (mom->numberOfMothers() > 0) ? mom->mother(0) : nullptr;
    }
  }
  return result;
}

GenVertex GenVertices::getGenVertexFromTrack(const reco::Track &track) const {

  GenVertex genVertex;
  for(const auto &gvertex : *this) {
    if(gvertex.hasTracks()) {
      for(const auto &matchTrack : gvertex.tracks())
	if(TrackHelper::SameTrack(track, matchTrack))
	  genVertex = gvertex;
    }
  }
  return genVertex;
}
