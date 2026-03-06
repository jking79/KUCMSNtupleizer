#ifndef KUCMSNtupleizer_KUCMSNtupleizer_MatchingTools_h
#define KUCMSNtupleizer_KUCMSNtupleizer_MatchingTools_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"

inline reco::GenParticleCollection CleanGenParticles(const reco::GenParticleCollection &genParticles, std::map<int, int> &indexMap) {

  reco::GenParticleCollection statusOne;
  for(size_t i = 0; i < genParticles.size(); i++) {

    if(genParticles[i].status() == 1 && genParticles[i].charge() != 0) {
      indexMap.emplace(statusOne.size(), i);
      statusOne.emplace_back(genParticles[i]);
    }
  }

  return statusOne;
}

// Overloaded CleanGen dispatchers used by DeltaRGenMatchHungarian
inline reco::GenParticleCollection CleanGen(const reco::GenParticleCollection &gens, std::map<int,int> &indexMap) {
  return CleanGenParticles(gens, indexMap);
}

inline std::vector<const pat::PackedGenParticle*> CleanGen(
    const std::vector<const pat::PackedGenParticle*> &gens, std::map<int,int> &indexMap) {
  std::vector<const pat::PackedGenParticle*> charged;
  for(size_t i = 0; i < gens.size(); i++) {
    if(gens[i]->charge() != 0) {
      indexMap.emplace(charged.size(), i);
      charged.push_back(gens[i]);
    }
  }
  return charged;
}

template <typename A, typename B>
inline double DeltaR2(const A &objectA, const B &objectB) {

  const double dEta(objectA.eta() - objectB.eta()), dPhi(objectA.phi() - objectB.phi());
  //std::cout << "delta R = " << sqrt(dEta*dEta + dPhi*dPhi) << std::endl;
  return dEta*dEta + dPhi*dPhi;
}

// Pointer overload — lets DeltaRMatchHungarian work when B is a pointer type
template <typename A, typename B>
inline double DeltaR2(const A &objectA, const B *objectB) {
  return DeltaR2(objectA, *objectB);
}
/*
template <typename A>
inline double DeltaR2(const A &object1, const A &object2) {
  if(&object1 == &object2) {
    std::cout << "Found the same object twice: pt1 = " << object1.pt() << ", pt2 = " << object2.pt() << std::endl;
    return 100;//std::numeric_limits<double>::max();
    
  }
  return DeltaR2(object1, object2);
}
*/

inline double DeltaR2(const reco::Track &track1, const reco::Track &track2) {
  if(TrackHelper::SameTrack(track1, track2) && &track1 == &track2) {
    //std::cout << "Found the same object twice: pt1 = " << track1.pt() << ", pt2 = " << track2.pt() << std::endl;
    return std::numeric_limits<double>::max();
  }

  const double dEta(track1.eta() - track2.eta()), dPhi(track1.phi() - track2.phi());
  return dEta*dEta + dPhi*dPhi;
}

template <typename B>
inline double DeltaR2(const reco::TransientTrack &ttrack, const B &objectB) {
  return DeltaR2(ttrack.track(), objectB);
}

inline double GetCorrectedPhi(const reco::TransientTrack &ttrack, const reco::Vertex &primaryVertex) {

  const reco::Track track(ttrack.track());
  SteppingHelixPropagator propagator(ttrack.field(), anyDirection);
  const GlobalPoint pvLocation(primaryVertex.x(), primaryVertex.y(), primaryVertex.z());

  const FreeTrajectoryState propLocation = propagator.propagateWithPath(ttrack.initialFreeState(), pvLocation).first;
  const GlobalPoint pca(propLocation.position());
  const double tmin = (pca.z() - track.vz())*track.qoverp()/(87.78*sin(track.lambda()));

  return track.phi()-tmin;
}

inline double DeltaR2(const reco::TransientTrack &ttrack1, const reco::TransientTrack &ttrack2, const reco::Vertex &primaryVertex) {

  const reco::Track track1(ttrack1.track()), track2(ttrack2.track());
  const double pt1(track1.pt()), pt2(track2.pt()), ptRel(2*fabs(pt1-pt2)/(pt1+pt2));
  const double eta1(track1.eta()), eta2(track2.eta()), etaDiff(fabs(eta1-eta2));
  const double phi1(track1.phi()), phi2(track2.phi()), phiDiff(fabs(phi1-phi2));

  if(ptRel > 0.2 || etaDiff > 0.1 || phiDiff > 1.57 || track1.charge() != track2.charge())
    return DeltaR2(track1, track2);

  const double dEta(track1.eta() - track2.eta()), dPhi(GetCorrectedPhi(ttrack1, primaryVertex) - GetCorrectedPhi(ttrack2, primaryVertex));
  
  return dEta*dEta + dPhi*dPhi;
}

// Templated GenDeltaR2: works for any gen type with the reco::Candidate interface
// (reco::GenParticle, pat::PackedGenParticle, etc.)
template <typename GenType>
inline double GenDeltaR2(const reco::TransientTrack &ttrack, const GenType &gen) {

  const reco::Track track(ttrack.track());
  const double pt(track.pt()), genPt(gen.pt()), ptRel(fabs(genPt-pt)/genPt);
  const double eta(track.eta()), genEta(gen.eta()), etaDiff(fabs(genEta-eta));
  const double phi(track.phi()), genPhi(gen.phi()), phiDiff(fabs(genPhi-phi));

  if(ptRel > 0.2 || etaDiff > 0.1 || phiDiff > 1.57 || track.charge() != gen.charge())
    return DeltaR2(track, gen);

  SteppingHelixPropagator propagator(ttrack.field(), anyDirection);

  const GlobalPoint genLocation(gen.vx(), gen.vy(), gen.vz());

  // Propagate from initial track state to PCA of gen vertex location
  //const FreeTrajectoryState propLocation = propagator.propagateWithPath(ttrack.initialFreeState(), genLocation).first;
  const TrajectoryStateClosestToPoint propLocation = ttrack.trajectoryStateClosestToPoint(genLocation);
  //const GlobalPoint pca(propLocation.position());
  const GlobalPoint pca(propLocation.referencePoint());
  // Calculate parameter t at the PCA
  const double tmin = (pca.z() - track.vz())*track.qoverp()/(87.78*sin(track.lambda()));

  // Estimate delta R with corrected phi (phi0 - tmin)
  const double dEta(track.eta() - gen.eta()), dPhi((track.phi() - tmin) - gen.phi());

  //std::cout << "delta R: " << sqrt(dEta*dEta + dPhi*dPhi) << std::endl;

  return dEta*dEta + dPhi*dPhi;
}

// Pointer overload — dereferences and delegates to the above
template <typename GenType>
inline double GenDeltaR2(const reco::TransientTrack &ttrack, const GenType *gen) {
  return GenDeltaR2(ttrack, *gen);
}

// ComputeGenCost: dispatches to GenDeltaR2 (helical correction) for TransientTrack,
// or plain DeltaR2 for all other track types.
template <typename A, typename GenType>
inline double ComputeGenCost(const A &a, const GenType &b) { return sqrt(DeltaR2(a, b)); }

template <typename GenType>
inline double ComputeGenCost(const reco::TransientTrack &a, const GenType &gen) {
  return sqrt(GenDeltaR2(a, gen));
}

template <typename T>
inline std::vector<T> RemoveDataAtIndices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete) {
  if(indicesToDelete.empty())
    return data;
  
  std::vector<T> ret;
  ret.reserve(data.size() - indicesToDelete.size());

  std::sort(indicesToDelete.begin(), indicesToDelete.end());

  typename std::vector<T>::const_iterator itBlockBegin = data.begin();
  for(std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++ it) {
    typename std::vector<T>::const_iterator itBlockEnd = data.begin() + *it;
    if(itBlockBegin != itBlockEnd) {
      std::copy(itBlockBegin, itBlockEnd, std::back_inserter(ret));
    }
    itBlockBegin = itBlockEnd + 1;
  }

  if(itBlockBegin != data.end()) {
    std::copy(itBlockBegin, data.end(), std::back_inserter(ret));
  }

  return ret;
}

#endif
