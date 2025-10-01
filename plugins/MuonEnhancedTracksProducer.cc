#ifndef MuonEnhancedTracksProducer_cc
#define MuonEnhancedTracksProducer_cc
// -*- C++ -*-
//
// Package:    KUCMSNtupleizer/KUCMSNtupleizer
// Class:      MuonEnhancedTracksProducer
//
/**\class MuonEnhancedTracksProducer MuonEnhancedTracksProducer.cc KUCMSNtupleizer/KUCMSNtupleizer/plugins/MuonEnhancedTracksProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Sat, 1 March 2025 12:11:43 GMT
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//Local includes
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/NewVertexMerger.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSet.h"

//
// class declaration
//

class MuonEnhancedTracksProducer : public edm::stream::EDProducer<> {

public:
  explicit MuonEnhancedTracksProducer(const edm::ParameterSet&);
  virtual ~MuonEnhancedTracksProducer() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------  
  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> dsaMuonTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> displacedGlobalMuonTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> displacedTracksToken_;
  edm::EDGetTokenT<reco::MuonCollection> displacedMuonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;

  edm::Handle<reco::TrackCollection> generalTracksHandle_;
  edm::Handle<reco::TrackCollection> dsaMuonTracksHandle_;
  edm::Handle<reco::TrackCollection> displacedGlobalMuonTracksHandle_;
  edm::Handle<reco::TrackCollection> displacedTracksHandle_;
  edm::Handle<reco::MuonCollection> displacedMuonsHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

  std::set<reco::TrackRef> dsaTrackSet_;
  
  std::vector<reco::Track> toTrackVec(const std::set<reco::TrackRef> &trackSet) const;
  std::set<reco::TrackRef> toTrackSet(const edm::Handle<reco::TrackCollection> &trackHandle) const;
  reco::TrackRef findRef(const reco::Track &track, const std::set<reco::TrackRef> &trackSet) const;
  void subMuonTracks(std::set<reco::TrackRef> &tracks, const std::set<reco::TrackRef> &muonTracks) const;
};

MuonEnhancedTracksProducer::MuonEnhancedTracksProducer(const edm::ParameterSet& iConfig) :
  generalTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTracks")) ),
  dsaMuonTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("dsaMuonTracks")) ),
  displacedGlobalMuonTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("dgMuonTracks")) ),
  displacedTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("displacedTracks")) ),
  displacedMuonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("displacedMuons"))),
  pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection")) ),
  transientTrackBuilder_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder")))
{
  produces<reco::TrackCollection>("combinedMuonTracks").setBranchAlias("combinedMuonTracks");
  produces<reco::TrackCollection>("muonEnhancedTracks").setBranchAlias("muonEnhancedTracks");
  produces<reco::TrackCollection>("sip2DTracks").setBranchAlias("sip2DTracks");
  produces<reco::TrackCollection>("sip2DMuonEnhancedTracks").setBranchAlias("sip2DMuonEnhancedTracks");
  produces<reco::TrackCollection>("selectedTracks").setBranchAlias("selectedTracks");
}

//
// member functions
//
// ------------ method called to produce the data  ------------
void MuonEnhancedTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;

  iEvent.getByToken( generalTracksToken_, generalTracksHandle_);
  iEvent.getByToken( dsaMuonTracksToken_, dsaMuonTracksHandle_);
  iEvent.getByToken( displacedGlobalMuonTracksToken_, displacedGlobalMuonTracksHandle_);
  iEvent.getByToken( displacedTracksToken_, displacedTracksHandle_);
  iEvent.getByToken( displacedMuonsToken_, displacedMuonsHandle_);
  iEvent.getByToken( pvToken_, pvHandle_);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  dsaTrackSet_ = toTrackSet(dsaMuonTracksHandle_);
  
  // Smart pointers to containers of output collections
  unique_ptr<reco::TrackCollection> combinedMuonTracks = make_unique<reco::TrackCollection>();
  unique_ptr<reco::TrackCollection> muonEnhancedTracks = make_unique<reco::TrackCollection>();
  unique_ptr<reco::TrackCollection> sip2DTracks = make_unique<reco::TrackCollection>();
  unique_ptr<reco::TrackCollection> sip2DMuonEnhancedTracks = make_unique<reco::TrackCollection>();
  unique_ptr<reco::TrackCollection> selectedTracks = make_unique<reco::TrackCollection>();

  std::vector<std::set<reco::TrackRef>> allMuonTracks({toTrackSet(displacedGlobalMuonTracksHandle_),
						       toTrackSet(displacedTracksHandle_),
						       toTrackSet(dsaMuonTracksHandle_)});

  std::vector<reco::TrackRef> allMuonTrackVector;
  for(const auto &collection : allMuonTracks) {
    allMuonTrackVector.insert(allMuonTrackVector.end(), collection.begin(), collection.end());
  }

  NewVertexMerger merger(pvHandle_->at(0));
  TrackVertexSetCollection vertexSeeds(merger.trackVertexSetSeeds(allMuonTrackVector, ttBuilder));
  vertexSeeds = TrackVertexSetHelper::applySignificanceCut(vertexSeeds.postCleaningCuts(pvHandle_->at(0)), pvHandle_->at(0));
  TrackVertexSetHelper::cleanDuplicates(vertexSeeds, pvHandle_->at(0));
  set<reco::TrackRef> muonTracks(vertexSeeds.completeTrackSet());
  
  for(const auto &trackRef : muonTracks)
    combinedMuonTracks->emplace_back(*trackRef);

  //cout << "Combining muon tracks with general tracks" << endl;
  set<reco::TrackRef> generalTrackSet(toTrackSet(generalTracksHandle_));  
  // Apply preliminary track cuts to all collections before attempting to merge
  // Use erase-remove idiom with a predicate
  auto it = generalTrackSet.begin();
  while (it != generalTrackSet.end()) {
    
    const auto& track = *it;
    const double pt(track->pt());
    const double ptError(track->ptError());
    const double ptRes(ptError/pt);
    const double normChi2(track->normalizedChi2());

    const bool foundInPV = std::any_of(pvHandle_->begin(), pvHandle_->end(), [&](const auto& pv) {
        return std::any_of(pv.tracks().begin(), pv.tracks().end(), [&](const auto& pvTrackRef) {
            return TrackHelper::SameTrack(*track, *pvTrackRef);
        });
    });
    
    if(pt < 1 || ptRes > 0.5 || normChi2 > 5 || foundInPV) {
      // Erase the current element and get the next valid iterator
      it = generalTrackSet.erase(it);
    } else {
      // Only increment if we didn't erase
      selectedTracks->emplace_back(*track);
      ++it;
    }
  }

  // Add all muon tracks to the selected general tracks. If matching criteria is satisfied erase general track from collection.
  subMuonTracks(generalTrackSet, muonTracks);

  // Fill output tracks with merged collection
  for(const auto &trackRef : generalTrackSet)
    muonEnhancedTracks->emplace_back(*trackRef);

  // Fill Sip2D tracks
  set<reco::TrackRef> sip2DTrackSet;
  for(const auto &track : *generalTracksHandle_) {

    reco::TransientTrack ttrack(ttBuilder->build(track));
    GlobalVector direction(track.px(), track.py(), track.pz());
    std::pair<bool, Measurement1D> ip2D = IPTools::signedTransverseImpactParameter(ttrack, direction, pvHandle_->at(0));

    const double pt(track.pt());
    const double ptError(track.ptError());
    const double ptRes(ptError/pt);
    const double normChi2(track.normalizedChi2());
    const double sip2D(ip2D.second.value() / ip2D.second.error());

    if(pt > 1 && ptRes < 0.5 && normChi2 < 5 && fabs(sip2D) > 4)
      sip2DTrackSet.insert(TrackHelper::GetTrackRef(track, generalTracksHandle_));
  }

  set<reco::TrackRef> sip2DMuonEnhancedTrackSet(sip2DTrackSet);
  subMuonTracks(sip2DMuonEnhancedTrackSet, muonTracks);
    
  for(const auto &trackRef : sip2DTrackSet)
    sip2DTracks->emplace_back(*trackRef);

  for(const auto &trackRef : sip2DMuonEnhancedTrackSet)
    sip2DMuonEnhancedTracks->emplace_back(*trackRef);
  
  //cout << "Finished" << endl;
  
  // Put in EDM
  iEvent.put(std::move(combinedMuonTracks), "combinedMuonTracks");
  iEvent.put(std::move(muonEnhancedTracks), "muonEnhancedTracks");
  iEvent.put(std::move(sip2DTracks), "sip2DTracks");
  iEvent.put(std::move(sip2DMuonEnhancedTracks), "sip2DMuonEnhancedTracks");
  iEvent.put(std::move(selectedTracks), "selectedTracks"); 
}// Producer end

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonEnhancedTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("generalTracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("dsaMuonTracks", edm::InputTag("displacedStandAloneMuons"));
  desc.add<edm::InputTag>("dgMuonTracks", edm::InputTag("displacedGlobalMuons"));
  desc.add<edm::InputTag>("displacedTracks", edm::InputTag("displacedTracks"));
  desc.add<edm::InputTag>("displacedMuons", edm::InputTag("displacedMuons"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"));
  
  descriptions.addDefault(desc);
}

std::vector<reco::Track> MuonEnhancedTracksProducer::toTrackVec(const std::set<reco::TrackRef> &trackSet) const {
  std::vector<reco::Track> tracks;
  for(const auto &track : trackSet)
    tracks.emplace_back(*track);
  return tracks;
}

std::set<reco::TrackRef> MuonEnhancedTracksProducer::toTrackSet(const edm::Handle<reco::TrackCollection> &trackHandle) const {
  std::set<reco::TrackRef> trackSet;
  for(size_t i = 0; i < trackHandle->size(); i++)
    trackSet.insert(reco::TrackRef(trackHandle, i));
  return trackSet;
}

reco::TrackRef MuonEnhancedTracksProducer::findRef(const reco::Track &track, const std::set<reco::TrackRef> &trackSet) const {
  reco::TrackRef trackRef;
  for(const auto &ref : trackSet)
    if(TrackHelper::SameTrack(track, *ref))
      trackRef = ref;

  return trackRef;
}

void MuonEnhancedTracksProducer::subMuonTracks(std::set<reco::TrackRef> &tracks, const std::set<reco::TrackRef> &muonTracks) const {
  DeltaRMatchHungarian<reco::Track, reco::Track> trackAssigner(toTrackVec(tracks), toTrackVec(muonTracks));
  tracks.insert(muonTracks.begin(), muonTracks.end());
  for(const auto &pair : trackAssigner.GetPairedObjects()) {
    const reco::Track track(pair.GetObjectA()), muonTrack(pair.GetObjectB());
    const reco::TrackRef muonTrackRef(findRef(muonTrack, muonTracks));
     
    // Matching criteria
    const double deltaR(pair.GetDeltaR());
    const double relPtDiff = 2*(track.pt() - muonTrack.pt()) / (track.pt() + muonTrack.pt());
    
    if(deltaR < 0.04 && fabs(relPtDiff) < 0.5 && track.charge() == muonTrack.charge())
      tracks.erase(findRef(track, tracks));
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonEnhancedTracksProducer);

#endif
