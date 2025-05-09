#ifndef KUCMSNtupleizer_KUCMSNtupleizer_GenVertex_h
#define KUCMSNtupleizer_KUCMSNtupleizer_GenVertex_h

#include <TVector3.h>
#include <Math/Vector4D.h>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/PairedObjects.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenTools.h"

class VertexHelper;

typedef PairedObjectCollection<reco::Track, reco::GenParticle> GenMatches;
typedef ROOT::Math::PtEtaPhiMVector LorentzVec;

class GenVertex {
 public:
  GenVertex() = default;
  virtual ~GenVertex() = default;

  GenVertex(const GenMatches &matchedPairs);
  GenVertex(const std::pair<reco::GenParticle, reco::GenParticle> &genPair);
  
  bool isGenElectron() const {return isGenElectron_;}
  bool isGenMuon() const {return isGenMuon_;}
  bool isGenHadronic() const {return isGenHadronic_;}
  bool isNoMatch(const reco::Vertex &vertex) const { return DetermineQuality(vertex) == VertexQuality::kNoMatch; }
  bool isEmpty(const reco::Vertex &vertex) const { return DetermineQuality(vertex) == VertexQuality::kEmpty; }
  bool isBronze(const reco::Vertex &vertex) const { return DetermineQuality(vertex) == VertexQuality::kBronze; }
  bool isSilver(const reco::Vertex &vertex) const { return DetermineQuality(vertex) == VertexQuality::kSilver; }
  bool isGold(const reco::Vertex &vertex) const { return DetermineQuality(vertex) == VertexQuality::kGold; }
  bool isInCollection(const reco::TrackCollection &tracks) const { return CountInstances(tracks) == size_; }
  bool isValidPair(const TransientTrackBuilder* ttBuilder) const { return fitTracks(ttBuilder).isValid(); }
  bool hasTracks() const {return hasMatch_;}
  bool isDefault() const {return size_==0;}
  bool isValid() const {return size_==2;}

  bool operator<(const GenVertex &other) const {return this->dxy() < other.dxy(); }
  bool operator==(const GenVertex &other) const { return this->x() == other.x() && this->y() == other.y() && this->z() == other.z(); }
  
  int size() const {return genMatches_.size();}
  int CountInstances(const reco::TrackCollection &tracks) const;
  
  double p() const {return vertexP_;}
  double px() const {return vertexPx_;}
  double py() const {return vertexPy_;}
  double pz() const {return vertexPz_;}
  double pt() const {return vertexPt_;}
  double eta() const {return vertexEta_;}
  double phi() const {return vertexPhi_;}
  double mass() const {return vertexMass_;}
  double x() const {return vertexX_;}
  double y() const {return vertexY_;}
  double z() const {return vertexZ_;}
  double dxy() const {return vertexDxy_;}
  double charge() const {return vertexCharge_;}
  double trackDz() const;
  double decayAngle() const;
  double momentumAngle() const;
  double distance3D(const reco::Vertex &vertex) const;
  
  reco::TrackCollection tracks() const {return genMatches_.GetObjectAList();}
  
  GenMatches genMatches() const {return genMatches_;}
  std::pair<reco::GenParticle, reco::GenParticle> genPair() const {return genPair_;}
  
  TransientVertex fitTracks(const TransientTrackBuilder* ttBuilder) const;
  reco::Vertex fitTracks(const TransientTrackBuilder* ttBuilder, const edm::Handle<reco::TrackCollection> &trackHandle) const;  

  GlobalPoint position() const {return GlobalPoint(this->x(), this->y(), this->z());}
  
  reco::GenParticleCollection getStableDaughters(const reco::GenParticleCollection &genParticles) const;
  reco::GenParticleCollection getStableChargedDaughters(const reco::GenParticleCollection &genParticles) const;
  reco::GenParticleCollection getStableNeutralDaughters(const reco::GenParticleCollection &genParticles) const;
  
private:
  enum class VertexQuality {kNoMatch, kEmpty, kBronze, kSilver, kGold};

  bool isGenElectron_;
  bool isGenMuon_;
  bool isGenHadronic_;
  bool hasMatch_;
  
  int size_;
  
  double vertexP_;
  double vertexPx_;
  double vertexPy_;
  double vertexPz_;
  double vertexPt_;
  double vertexEta_;
  double vertexPhi_;
  double vertexMass_;
  double vertexX_;
  double vertexY_;
  double vertexZ_;
  double vertexDxy_;
  double vertexCharge_;
  
  GenMatches genMatches_;
  std::pair<reco::GenParticle, reco::GenParticle> genPair_;
  
  int nMatched(const reco::Vertex &vertex) const;
  VertexQuality DetermineQuality(const reco::Vertex &vertex) const;
  void setMembers(const LorentzVec &genVertex4Vec, const reco::GenParticle &gen);
};

class GenVertices : public std::vector<GenVertex> {

 public:
  GenVertices() = default;
  virtual ~GenVertices() = default;

  GenVertices(const GenVertices& other) : std::vector<GenVertex>(other) {}
  GenVertices(const GenMatches &matchedPairs, const double deltaRCut = 0.02);
  GenVertices(const std::vector<reco::GenParticle> &genParticles);

  bool contains(const GenVertex& genVertex) const;
  
  std::vector<TransientVertex> getValidVertices(const TransientTrackBuilder* ttBuilder, const double chi2Cut = 5) const;
  reco::GenParticleCollection getAllGenParticles() const;
  
  void operator+=(const GenVertices &other);

  GenVertex getGenVertexFromTrack(const reco::Track &track) const;
  
 private:
  std::vector<GenMatches> findGenParticlePairs(const GenMatches &matchedPairs, const double deltaRCut) const;
  GenMatches FindSignalGenCollection(const GenMatches &matchedPairs, const double deltaRCut) const;

};

#endif
