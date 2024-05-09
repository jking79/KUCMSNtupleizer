#ifndef IsolationTools_h
#define IsolationTools_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "DataFormats/Math/interface/deltaR.h"

class IsolationInfo : public std::map<std::string, double> {
  
 public:
  
  IsolationInfo() = default; 
  virtual ~IsolationInfo() = default;
  
  static std::vector<double> GetTrackIsolation(const reco::TrackCollection &tracks, 
					       const reco::BeamSpot &beamSpot,
					       const double outerRadius = 0.4, 
					       const double innerRadius = 0.015);

  static std::vector<double> GetPFIsolation(const reco::TrackCollection &tracks,
					    const reco::PFCandidateCollection &pfCands,
					    const double outerRadius = 0.4,
					    const double innerRadius = 0.);

  static std::vector<double> GetSCIsolation(const reco::ElectronCollection &electronCandidates,
					    const reco::SuperClusterCollection &superClusters,
					    const double radius = 0.3);

};


#endif
