#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"

std::vector<double> IsolationInfo::GetTrackIsolation(const reco::TrackCollection &tracks, 
						      const reco::BeamSpot &beamSpot,
						      const double outerRadius, 
						      const double innerRadius) {

  std::vector<double> isoValues(tracks.size());
  ElectronTkIsolation myTkIsolation(outerRadius,
                                    innerRadius, //intRadiusBarrel
                                    innerRadius, //intRadiusEndCap
                                    innerRadius, //stripBarrel
                                    innerRadius, //stripEndcap
                                    0.95, //ptMin
                                    0.2, //MaxVertexDist
                                    9999., //MaxVertexDXY
                                    &tracks,
                                    beamSpot.position());

  for(unsigned int t = 0; t < tracks.size(); t++) {
    const double isoValue = myTkIsolation.getIso(&(tracks.at(t))).second;
    isoValues[t] = isoValue;
  }

  return isoValues;
}

std::vector<double> IsolationInfo::GetPFIsolation(const reco::TrackCollection&tracks,
						  const reco::PFCandidateCollection &pfCands,
						  const double outerRadius,
						  const double innerRadius) {


  std::vector<double> isoValues;

  for(const auto &track : tracks) {
    
    double ptSum = 0;
    for(const auto &cand : pfCands) {
      
      if(cand.charge() != 0)
	continue;
      
      const double deltaR = sqrt(reco::deltaR2(track, cand));
      if(deltaR*deltaR < outerRadius*outerRadius && deltaR*deltaR >= innerRadius*innerRadius)
	ptSum += cand.pt();
    }

    isoValues.push_back(ptSum);
  }
  
  return isoValues;
}

std::vector<double> IsolationInfo::GetSCIsolation(const reco::ElectronCollection &electronCandidates,
						  const reco::SuperClusterCollection &superClusters,
						  const double radius) {
  
  std::vector<double> isoValues;
  for(const auto &cand : electronCandidates) {

    reco::SuperCluster electronSC(*(cand.superCluster().get()));
    double ptSum = 0;
    for(const auto &sc : superClusters) {
      if(sc.rawEnergy() == electronSC.rawEnergy() && sc.eta() == electronSC.eta())
	continue;

      // Loop over the constituent calo clusters in the current supercluster 
      // If it's located within the radius add its corrected energy
      for (reco::CaloCluster_iterator it = sc.clustersBegin(); it != sc.clustersEnd(); it++)
	if(reco::deltaR2(electronSC, *(*it)) <= radius*radius)
	  ptSum += (*it)->correctedEnergy();
    }
    isoValues.push_back(ptSum);
  }
  return isoValues;
}
