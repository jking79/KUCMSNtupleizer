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

  std::unordered_map<std::string, std::vector<reco::GenParticle>> positionMap;
  
  for(const auto &gen : genParticles) {
    
    //if(gen.status() != 1) continue;

    //isSignalGenJet(gen);
    
    bool isSignal(isSignalGenElectron(gen) || isSignalGenMuon(gen) || isSignalGenJet(gen));
    
    if(!isSignal) continue;
    
    std::string key = std::to_string(gen.vx()) + "_" +
                      std::to_string(gen.vy()) + "_" +
                      std::to_string(gen.vz());

    positionMap[key].push_back(gen);
  }
  
  for (auto it = positionMap.begin(); it != positionMap.end(); ++it) {
    if(it->second.size() != 2)
      throw std::runtime_error("Major Problem: found a gen vertex with " + std::to_string(it->second.size()) + " particles!");

    this->emplace_back(GenVertex(std::make_pair(it->second[0], it->second[1])));
  }
}

bool GenVertices::contains(const GenVertex& genVertex) const {

  for(const auto &vertex : *this) {

    if(genVertex == vertex && vertex.hasTracks())
      return true;
  }

  return false;
}

std::vector<GenMatches> GenVertices::findGenParticlePairs(const GenMatches &matchedPairs, const double deltaRCut) const {
  std::vector<GenMatches> pairs;
  std::unordered_map<std::string, GenMatches> positionMap;

  for (const auto& pair : FindSignalGenCollection(matchedPairs, deltaRCut)) {
    reco::GenParticle particle(pair.GetObjectB());
    
    // Create a unique key for each position
    std::string key = std::to_string(particle.vx()) + "_" +
                      std::to_string(particle.vy()) + "_" +
                      std::to_string(particle.vz());

    positionMap[key].push_back(pair);
  }

  for (const auto& entry : positionMap) {
    if (entry.second.size() == 2 || entry.second.size() == 4) {
      // We found a pair
      pairs.push_back(entry.second);
    }
    if(entry.second.size() == 4)
      std::cout << "Found a case with four signal gen particles at the same vertex!" << std::endl;

    // If size is 1 or 3, we ignore it as it's not a complete pair
  }

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

void GenVertices::operator+=(const GenVertices &other) {
  // Use a set to merge and avoid duplicates
  std::set<GenVertex> mergedVertices(this->begin(), this->end());
  
  // Insert vertices from the other set
  for(const auto &vertex : other) {
    if(vertex.hasTracks()) {
      mergedVertices.erase(vertex);
      mergedVertices.insert(vertex);
    }
    else if(mergedVertices.find(vertex) == mergedVertices.end())
      mergedVertices.insert(vertex);
  }
  
  // Clear the current container and populate it with the merged vertices
  this->clear();
  for(const auto &vertex : mergedVertices) {
    this->emplace_back(vertex);
  }
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
