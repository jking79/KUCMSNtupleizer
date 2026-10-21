// -*- C++ -*-
//
// Package:    KUNTupleFW/KUCMSNtupleizer
// Class:      ECALTracksProducer
//
/**\class ECALTracksProducer ECALTracksProducer.cc KUNTupleFW/KUCMSNtupleizer/plugins/ECALTracksProducer.cc

Description: Combine in-time and out-of-time ECAL superclusters after removing
overlaps, and publish the resulting collection for displaced-electron studies.

Implementation: The current producer consumes only the two supercluster
collections.  Earlier track-propagation code is no longer part of this module.
*/
//
// Original Author:  Andres Abreu
//         Created:  Sat, 8 March 2024 10:55:00 GMT
//
//

// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/TrackExtrapolation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/HelixExtrapolatorToLine2Order.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

//Local includes
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"


//
// class declaration
//

class ECALTracksProducer : public edm::stream::EDProducer<> {
  public:
  explicit ECALTracksProducer(const edm::ParameterSet&);
  ~ECALTracksProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
  typedef edm::Handle<reco::TrackCollection> TracksHandle;
  typedef edm::Handle<reco::GsfTrackCollection> GsfTracksHandle;
  typedef edm::Handle<reco::SuperClusterCollection> SCsHandle;

  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input supercluster collections used to build the displaced-electron view.
  edm::EDGetTokenT<reco::SuperClusterCollection> superClusterToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> ootSuperClusterToken_;
};

ECALTracksProducer::ECALTracksProducer(const edm::ParameterSet& iConfig) :
superClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters")) ),
ootSuperClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters")) )
{
  produces<reco::SuperClusterCollection>("displacedElectronSCs").setBranchAlias("displacedElectronSCs");
}

ECALTracksProducer::~ECALTracksProducer() {}

//
// member functions
//
// ------------ method called to produce the data  ------------
void ECALTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  edm::Handle<reco::SuperClusterCollection> superClusterHandle;
  iEvent.getByToken(superClusterToken_, superClusterHandle);

  edm::Handle<reco::SuperClusterCollection> ootSuperClusterHandle;
  iEvent.getByToken(ootSuperClusterToken_, ootSuperClusterHandle);
  std::unique_ptr<reco::SuperClusterCollection> displacedElectronSCs = std::make_unique<reco::SuperClusterCollection>();
  // Remove overlap between superClusters and out-of-time superClusters
  reco::SuperClusterCollection superClusters(*superClusterHandle);
  reco::SuperClusterCollection ootSuperClusters(*ootSuperClusterHandle);

  superClusters = RemoveOverlap<reco::SuperCluster, reco::SuperCluster>(superClusters, ootSuperClusters, 0.1);
  displacedElectronSCs->insert(displacedElectronSCs->end(), superClusters.begin(), superClusters.end());
  displacedElectronSCs->insert(displacedElectronSCs->end(), ootSuperClusters.begin(), ootSuperClusters.end());

  iEvent.put(std::move(displacedElectronSCs), "displacedElectronSCs");
}//<<>>void ECALTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ECALTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}//<<>>void ECALTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)

//define this as a plug-in
DEFINE_FWK_MODULE(ECALTracksProducer);
