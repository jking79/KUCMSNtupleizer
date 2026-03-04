// -*- C++ -*-
//
// Package:    KUCMSNtupleizer
// Class:      HyddraSVsProducer
//
// Description: Produces leptonic and hadronic displaced secondary vertices
//              using the HYDDRA reconstruction framework.
//
// Original Author:  Andres Abreu
//

#include <memory>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// HYDDRA
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/LeptonicHYDDRA.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/HadronicHYDDRA.h"

class HyddraSVsProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraSVsProducer(const edm::ParameterSet&);
  ~HyddraSVsProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // HYDDRA reconstruction objects
  LeptonicHYDDRA leptonic_;
  HadronicHYDDRA hadronic_;
};

HyddraSVsProducer::HyddraSVsProducer(const edm::ParameterSet& iConfig) :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  leptonic_(iConfig.getParameter<edm::ParameterSet>("leptonic")),
  hadronic_(iConfig.getParameter<edm::ParameterSet>("hadronic"))
{
  produces<reco::VertexCollection>("leptonicVertices");
  produces<reco::VertexCollection>("hadronicVertices");
}

void HyddraSVsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get inputs
  edm::Handle<reco::TrackCollection> tracksHandle;
  iEvent.getByToken(tracksToken_, tracksHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);

  // Build TrackRef vector from the input collection
  std::vector<reco::TrackRef> trackRefs;
  trackRefs.reserve(tracksHandle->size());
  for (size_t i = 0; i < tracksHandle->size(); ++i) {
    trackRefs.emplace_back(tracksHandle, i);
  }

  if(trackRefs.size() > 500) {
    edm::LogWarning("HyddraSVsProducer") << "Large input track collection (" << trackRefs.size()
      << " tracks). Vertex reconstruction may be very slow. "
      << "Consider enabling track cuts (applyCuts=True).";
  }

  // Use leading PV
  const reco::Vertex& pv = pvHandle->at(0);

  // Run both reconstruction paths
  leptonic_.run_reconstruction(trackRefs, ttBuilder, pv);
  hadronic_.run_reconstruction(trackRefs, ttBuilder, pv);

  // Convert to reco::VertexCollection and put into event
  auto leptonicVertices = std::make_unique<reco::VertexCollection>(leptonic_.vertices());
  auto hadronicVertices = std::make_unique<reco::VertexCollection>(hadronic_.vertices());

  iEvent.put(std::move(leptonicVertices), "leptonicVertices");
  iEvent.put(std::move(hadronicVertices), "hadronicVertices");
}

void HyddraSVsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collections
  desc.add<edm::InputTag>("tracks", edm::InputTag("miniAODMuonEnhancedTracks", "sip2DMuonEnhancedTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlineSlimmedPrimaryVertices"));

  // Leptonic PSet
  edm::ParameterSetDescription leptonicDesc;
  leptonicDesc.add<double>("seedCosThetaCut",              0.75);
  leptonicDesc.add<double>("minMass",                      2.0);
  leptonicDesc.add<double>("minPOverE",                    0.6);
  leptonicDesc.add<double>("maxNormChi2",                  5.0);
  leptonicDesc.add<double>("minDxySignificance",           25.0);
  leptonicDesc.add<double>("maxCompatibility",             1.5);
  leptonicDesc.add<double>("minCleanCosTheta",             0.5);
  leptonicDesc.add<bool>  ("useDiagonalCut",               false);
  leptonicDesc.add<double>("cleanCutSlope",                0.0);
  leptonicDesc.add<double>("minTrackCosTheta",             0.5);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Limit",     0.95);
  leptonicDesc.add<double>("maxTrackCosThetaCM_Intercept", 1.8);
  leptonicDesc.add<double>("trackCosThetaCM_Slope",        -1.0);
  leptonicDesc.add<bool>  ("requireChargeNeutrality",      true);
  leptonicDesc.add<double>("minVtxCosTheta",               -1.0);
  leptonicDesc.add<bool>  ("useAbsVtxCosTheta",            false);
  leptonicDesc.add<bool>  ("doMerging",                    true);
  leptonicDesc.add<bool>  ("doCleaning",                   true);
  leptonicDesc.add<bool>  ("doDisambiguation",             true);
  leptonicDesc.add<bool>  ("doFiltering",                  true);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  // Hadronic PSet
  edm::ParameterSetDescription hadronicDesc;
  hadronicDesc.add<double>("seedCosThetaCut", 0.0);
  hadronicDesc.add<double>("minMass", 2.0);
  hadronicDesc.add<double>("minPOverE", 0.6);
  hadronicDesc.add<double>("maxNormChi2", 5.0);
  hadronicDesc.add<double>("minDxySignificance", 40.0);
  hadronicDesc.add<int>("minSize", 5);
  hadronicDesc.add<double>("minCosTheta", 0.0);
  hadronicDesc.add<double>("maxDecayAngle", 0.9);
  hadronicDesc.add<bool>("doMerging", true);
  hadronicDesc.add<bool>("doCleaning", true);
  hadronicDesc.add<bool>("doDisambiguation", true);
  hadronicDesc.add<bool>("doFiltering", true);
  desc.add<edm::ParameterSetDescription>("hadronic", hadronicDesc);

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HyddraSVsProducer);
