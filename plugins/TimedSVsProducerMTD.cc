// -*- C++ -*-
//
// Package:    TimingWithSVs/SVsTimingProducer
// Class:      SVsTimingProducer
//
/**\class SVsTimingProducer SVsTimingProducer.cc TimingWithSVs/SVsTimingProducer/plugins/SVsTimingProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Sat, 21 Oct 2023 20:06:20 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/HelixExtrapolatorToLine2Order.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/PropagatedTrack.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TimingHelper.h"

//
// class declaration
//

class TimedSVsProducerMTD : public edm::stream::EDProducer<> {
public:
  explicit TimedSVsProducerMTD(const edm::ParameterSet&);
  ~TimedSVsProducerMTD() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  std::string vertexLabel_;

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> svToken_;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;
};

TimedSVsProducerMTD::TimedSVsProducerMTD(const edm::ParameterSet& iConfig) 
  : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
    svToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("svSrc"))),
    transientTrackBuilder_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()) {
  
  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters_.loadParameters(parameters, iC);
  trackAssociator_.useDefaultPropagator();
  
  vertexLabel_ = iConfig.getParameter<std::string>("vertexLabel");
  produces<reco::VertexCollection>(vertexLabel_).setBranchAlias(vertexLabel_);

}

TimedSVsProducerMTD::~TimedSVsProducerMTD() {}

//
// member functions
//
// ------------ method called to produce the data  ------------
void TimedSVsProducerMTD::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::unique_ptr<reco::VertexCollection> timedVertices = std::make_unique<reco::VertexCollection>();

  const edm::Handle<reco::TrackCollection> tracksHandle = iEvent.getHandle(tracksToken_);
  const edm::Handle<reco::VertexCollection> svHandle = iEvent.getHandle(svToken_);
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);

  auto magfield = iSetup.getTransientHandle(magneticFieldToken_);
  
  for(const reco::Vertex &sv : *svHandle) {

    const GlobalPoint svPosition(sv.x(), sv.y(), sv.z());
    TrackPropagator svTracksPropagator(iEvent, iSetup, magfield, trackAssocParameters_, VertexHelper::GetTracks(sv));
    PropagatedTracks<reco::Track> svPropTracks(svTracksPropagator.GetPropagatedTracks());

    double tof(0.);
    std::vector<std::pair<reco::Track, reco::Track> > timedTrackToTrackMap;
    for(const auto &propTrack : svPropTracks) {

      const reco::Track track(propTrack.GetTrack());
      const reco::TransientTrack ttrack(ttBuilder->build(track));
      const TrackDetMatchInfo detInfo(propTrack.GetDetInfo());
      const GlobalPoint trackGlobalPosAtEcal(detInfo.trkGlobPosAtEcal.x(), detInfo.trkGlobPosAtEcal.y(), detInfo.trkGlobPosAtEcal.z());
      const double pathLength(TimingHelper::PathLength(ttrack, svPosition, trackGlobalPosAtEcal));

      timedTrackToTrackMap.emplace_back(std::make_pair(track, TimingHelper::TimeStampTrack(track, pathLength)));
      tof += TimingHelper::Time(track, 0.13957018, pathLength);
    } 

    const double avgTof(tof/svPropTracks.size());
    reco::Vertex timedSV = VertexHelper::TimeStampVertex(sv, avgTof);

    for(const auto &trackPairs : timedTrackToTrackMap) {
      const reco::TransientTrack ttrackOg = ttBuilder->build(trackPairs.first);
      const reco::TransientTrack ttrack = ttBuilder->build(trackPairs.second);
      timedSV.add(ttrack.trackBaseRef(), ttrack.track(), sv.trackWeight(ttrackOg.trackBaseRef()));
    }
    
    timedVertices->emplace_back(VertexHelper::TimeStampVertex(sv, avgTof));
  }

  iEvent.put(std::move(timedVertices), "timedVertices");
}// Producer end

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TimedSVsProducerMTD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TimedSVsProducerMTD);
