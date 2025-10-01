#ifndef KUCMSNtupleizer_KUCMSNtupleizer_TrackPropagator_h
#define KUCMSNtupleizer_KUCMSNtupleizer_TrackPropagator_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/PropagatedTrack.h"

//===============================================================================================//
//                                      class: TrackPropagator                                   //
//-----------------------------------------------------------------------------------------------//
// Class for propagating any set of cmssw tracks (tested on those that inherit from reco::Track) //
// to the ECAL and HCAL. Returns the set of propagated tracks that made it to the ECAL. The ones //
// that don't make it are discarded.                                                             //  
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <class T>
class TrackPropagator {
  
 public:
  
  TrackPropagator(const edm::Event &iEvent, 
		  const edm::EventSetup &iSetup,
		  const edm::ESTransientHandle<MagneticField> &magneticField,
		  const TrackAssociatorParameters &parameters,
		  const std::vector<T> &tracks);

  virtual ~TrackPropagator() = default;

  PropagatedTracks<T> GetPropagatedTracks() const { return ecalTracks_; }
  std::vector<T> GetTracks() const; 

 private:

  TrackDetectorAssociator trackAssociator_;
  
  PropagatedTracks<T> ecalTracks_;
  
  TrackDetMatchInfo GetTrackDetMatchInfo(const edm::Event &iEvent,
					 const edm::EventSetup &iSetup,
					 const edm::ESTransientHandle<MagneticField> &magneticField,
					 const TrackAssociatorParameters &parameters,
					 const reco::Track& track);
};

template <class T> 
TrackPropagator<T>::TrackPropagator(const edm::Event &iEvent,
				    const edm::EventSetup &iSetup,
				    const edm::ESTransientHandle<MagneticField> &magneticField,
				    const TrackAssociatorParameters &parameters,
				    const std::vector<T> &tracks) {

  trackAssociator_.useDefaultPropagator();

  ecalTracks_.clear();
  int trackIndex = 0;
  for(const auto &track : tracks) {

    TrackDetMatchInfo detInfo(GetTrackDetMatchInfo(iEvent, iSetup, magneticField, parameters, track));
    if(detInfo.crossedEcalIds.size() != 0)// && track.pt() > 0.95 && track.qualityMask() < 7) 
      ecalTracks_.emplace_back(PropagatedTrack(trackIndex++, track, detInfo));

  }
}

template <class T>
std::vector<T> TrackPropagator<T>::GetTracks() const {
  
  std::vector<T> tracks;
  for(const auto &pair : ecalTracks_) {
    tracks.emplace_back(pair.GetTrack());
  }
  return tracks;
}

template <class T>
TrackDetMatchInfo TrackPropagator<T>::GetTrackDetMatchInfo(const edm::Event &iEvent,
							   const edm::EventSetup &iSetup,
							   const edm::ESTransientHandle<MagneticField> &magneticField,
							   const TrackAssociatorParameters &parameters,
							   const reco::Track& track) {

  FreeTrajectoryState initialState = trajectoryStateTransform::initialFreeState(track, magneticField.product());
  
  const TrackDetMatchInfo detInfo(trackAssociator_.associate(iEvent, iSetup, parameters, &initialState) );

  return detInfo;

}

#endif
