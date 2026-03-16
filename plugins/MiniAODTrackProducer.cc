// -*- C++ -*-
//
// Package:    HyddraSVProducer
// Class:      MiniAODTrackProducer
//
// Description: Extracts reco::Track collections from MINIAOD packed candidates.
//              Produces separate and merged collections from packedPFCandidates
//              and lostTracks for use with downstream vertex reconstruction.
//
// Original Author:  Andres Abreu
//

#include <memory>
#include <vector>
#include <cmath>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class MiniAODTrackProducer : public edm::stream::EDProducer<> {

public:
  explicit MiniAODTrackProducer(const edm::ParameterSet&);
  ~MiniAODTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Helper to extract tracks from a PackedCandidateCollection (with optional cuts)
  void extractTracks(const pat::PackedCandidateCollection& candidates,
                     reco::TrackCollection& outputTracks,
                     const TransientTrackBuilder* ttBuilder,
                     const reco::Vertex* pv) const;

  // Helper to extract global tracks from pat::Muon collection
  void extractMuonTracks(const std::vector<pat::Muon>& muons,
                         reco::TrackCollection& outputTracks) const;

  // Helper to extract muon tracks using the LLPNanoAOD PatMuonVertex track selection:
  //   isGlobalMuon()     -> combinedMuon()
  //   isStandAloneMuon() -> standAloneMuon()
  //   else               -> tunePMuonBestTrack()
  void extractMuonTracksLLPNano(const std::vector<pat::Muon>& muons,
                                reco::TrackCollection& outputTracks) const;

  // Helper to add tracks to a merged collection, dropping any that fall within
  // deltaR < 0.01 of a PF electron track in the reference collection.
  void addTracksDeduped(const reco::TrackCollection& inputTracks,
                        const reco::TrackCollection& electronTracks,
                        reco::TrackCollection& outputTracks) const;

  // Build the muon-merged sip2D collection:
  //   PF (|pdgId|!=13) + Lost (deduped vs electrons, |pdgId|!=13),
  //   with any track matching a slimmedMuon (deltaR<0.001, |relPt|<0.1) replaced
  //   by the muon bestTrack, then sip2D cuts applied uniformly to all tracks.
  void buildMuonMergedSip2D(const pat::PackedCandidateCollection& pfCands,
                             const pat::PackedCandidateCollection& lostCands,
                             const std::vector<pat::Muon>& muons,
                             const TransientTrackBuilder& ttBuilder,
                             const reco::Vertex& pv,
                             reco::TrackCollection& outputTracks) const;

  // Input tokens
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> eleLostTracksToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> displacedMuonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // Cut parameters
  bool applyCuts_;
  double minPt_;
  double minAbsSip2D_;
  double maxNormalizedChi2_;
};

MiniAODTrackProducer::MiniAODTrackProducer(const edm::ParameterSet& iConfig)
    : pfCandidatesToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      lostTracksToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("lostTracks"))),
      eleLostTracksToken_(consumes<pat::PackedCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("eleLostTracks"))),
      muonsToken_(consumes<std::vector<pat::Muon>>(
          iConfig.getParameter<edm::InputTag>("muons"))),
      displacedMuonsToken_(consumes<std::vector<pat::Muon>>(
          iConfig.getParameter<edm::InputTag>("displacedMuons"))),
      pvToken_(consumes<reco::VertexCollection>(
          iConfig.getParameter<edm::InputTag>("pvCollection"))),
      ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      applyCuts_(iConfig.getParameter<bool>("applyCuts")),
      minPt_(iConfig.getParameter<double>("minPt")),
      minAbsSip2D_(iConfig.getParameter<double>("minAbsSip2D")),
      maxNormalizedChi2_(iConfig.getParameter<double>("maxNormalizedChi2")) {

  // Individual collections from packed candidates
  produces<reco::TrackCollection>("pfCandidateTracks");
  produces<reco::TrackCollection>("lostTracks");
  produces<reco::TrackCollection>("eleLostTracks");

  // Merged collections from packed candidates
  produces<reco::TrackCollection>("merged");        // pfCandidates + lostTracks
  produces<reco::TrackCollection>("mergedWithEle"); // pfCandidates + eleLostTracks
  produces<reco::TrackCollection>("mergedAll");     // pfCandidates + lostTracks + eleLostTracks

  // Muon global tracks
  produces<reco::TrackCollection>("muonGlobalTracks");
  produces<reco::TrackCollection>("displacedMuonGlobalTracks");

  // Muon tracks using LLPNanoAOD PatMuonVertex selection
  produces<reco::TrackCollection>("muonLLPNanoTracks");
  produces<reco::TrackCollection>("displacedMuonLLPNanoTracks");

  // PF+Lost (no muon-ID tracks) deduped vs slimmedMuons, union with all muon bestTracks,
  // with sip2D cuts applied uniformly
  produces<reco::TrackCollection>("mergedMuonSip2D");
}

void MiniAODTrackProducer::extractTracks(
    const pat::PackedCandidateCollection& candidates,
    reco::TrackCollection& outputTracks,
    const TransientTrackBuilder* ttBuilder,
    const reco::Vertex* pv) const {

  for (const auto& cand : candidates) {
    // Must have a track
    if (!cand.hasTrackDetails()) {
      continue;
    }

    // Get the pseudo-track
    const reco::Track track = cand.pseudoTrack();

    // Apply cuts if enabled
    if (applyCuts_) {
      // pT cut
      if (track.pt() <= minPt_) {
        continue;
      }

      // Normalized chi2 cut
      if (track.normalizedChi2() >= maxNormalizedChi2_) {
        continue;
      }

      // sip2D cut (requires transient track and PV)
      if (ttBuilder != nullptr && pv != nullptr) {
        reco::TransientTrack ttrack = ttBuilder->build(track);
        auto ip2dResult = IPTools::signedTransverseImpactParameter(
            ttrack, GlobalVector(track.px(), track.py(), track.pz()), *pv);

        if (ip2dResult.first) {
          double sip2D = ip2dResult.second.significance();
          // Keep tracks with |sip2D| >= minAbsSip2D (displaced tracks)
          if (std::fabs(sip2D) < minAbsSip2D_) {
            continue;
          }
        }
      }
    }

    // Track passed all cuts (or cuts disabled)
    outputTracks.push_back(track);
  }
}

void MiniAODTrackProducer::extractMuonTracks(
    const std::vector<pat::Muon>& muons,
    reco::TrackCollection& outputTracks) const {

  for (const auto& muon : muons) {
    if (muon.muonBestTrack().isNonnull()) {
      outputTracks.push_back(*muon.muonBestTrack());
    }
  }
}

void MiniAODTrackProducer::extractMuonTracksLLPNano(
    const std::vector<pat::Muon>& muons,
    reco::TrackCollection& outputTracks) const {

  // Replicates the LLPNanoAOD PatMuonVertex track selection:
  //   isGlobalMuon()     -> combinedMuon()
  //   isStandAloneMuon() -> standAloneMuon()
  //   else               -> tunePMuonBestTrack()
  for (const auto& muon : muons) {
    reco::TrackRef trackRef;
    if (muon.isGlobalMuon())
      trackRef = muon.combinedMuon();
    else if (muon.isStandAloneMuon())
      trackRef = muon.standAloneMuon();
    else
      trackRef = muon.tunePMuonBestTrack();

    if (trackRef.isNonnull()) {
      outputTracks.push_back(*trackRef);
    }
  }
}

void MiniAODTrackProducer::addTracksDeduped(
    const reco::TrackCollection& inputTracks,
    const reco::TrackCollection& electronTracks,
    reco::TrackCollection& outputTracks) const {

  for (const auto& track : inputTracks) {
    bool overlaps = false;
    for (const auto& eTrack : electronTracks) {
      double dEta = track.eta() - eTrack.eta();
      double dPhi = std::fabs(track.phi() - eTrack.phi());
      if (dPhi > M_PI) dPhi = 2.0 * M_PI - dPhi;
      if (std::sqrt(dEta * dEta + dPhi * dPhi) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (!overlaps) outputTracks.push_back(track);
  }
}

void MiniAODTrackProducer::buildMuonMergedSip2D(
    const pat::PackedCandidateCollection& pfCands,
    const pat::PackedCandidateCollection& lostCands,
    const std::vector<pat::Muon>& muons,
    const TransientTrackBuilder& ttBuilder,
    const reco::Vertex& pv,
    reco::TrackCollection& outputTracks) const {

  // Cache muon (eta, phi, pt) for fast duplicate checking
  struct MuonKin { double eta, phi, pt; };
  std::vector<MuonKin> muonCache;
  muonCache.reserve(muons.size());
  for (const auto& mu : muons)
    muonCache.push_back({mu.eta(), mu.phi(), mu.pt()});

  // Returns true if track is a duplicate of any slimmedMuon
  auto matchesMuon = [&](const reco::Track& track) -> bool {
    for (const auto& mu : muonCache) {
      double dEta = track.eta() - mu.eta;
      double dPhi = std::fabs(track.phi() - mu.phi);
      if (dPhi > M_PI) dPhi = 2.0 * M_PI - dPhi;
      if (std::sqrt(dEta * dEta + dPhi * dPhi) >= 0.001) continue;
      if (mu.pt <= 0) continue;
      if (std::fabs(track.pt() - mu.pt) / mu.pt < 0.1) return true;
    }
    return false;
  };

  // Cache PF electron (eta, phi) for lost-track dedup (no quality cuts here)
  std::vector<std::pair<double, double>> eleCache;
  for (const auto& cand : pfCands) {
    if (!cand.hasTrackDetails() || std::abs(cand.pdgId()) != 11) continue;
    const reco::Track tr = cand.pseudoTrack();
    eleCache.emplace_back(tr.eta(), tr.phi());
  }

  auto overlapsElectron = [&](const reco::Track& track) -> bool {
    for (const auto& [eEta, ePhi] : eleCache) {
      double dEta = track.eta() - eEta;
      double dPhi = std::fabs(track.phi() - ePhi);
      if (dPhi > M_PI) dPhi = 2.0 * M_PI - dPhi;
      if (std::sqrt(dEta * dEta + dPhi * dPhi) < 0.01) return true;
    }
    return false;
  };

  // Build raw merged collection (no sip2D cuts yet)
  reco::TrackCollection raw;

  // PF candidates: skip muon-ID'd tracks and any duplicate of a slimmedMuon
  for (const auto& cand : pfCands) {
    if (!cand.hasTrackDetails()) continue;
    if (std::abs(cand.pdgId()) == 13) continue;
    const reco::Track track = cand.pseudoTrack();
    if (matchesMuon(track)) continue;
    raw.push_back(track);
  }

  // Lost tracks: skip muon-ID'd, electron dedup, and muon duplicates
  for (const auto& cand : lostCands) {
    if (!cand.hasTrackDetails()) continue;
    if (std::abs(cand.pdgId()) == 13) continue;
    const reco::Track track = cand.pseudoTrack();
    if (overlapsElectron(track)) continue;
    if (matchesMuon(track)) continue;
    raw.push_back(track);
  }

  // Add all slimmedMuon bestTracks unconditionally
  for (const auto& mu : muons) {
    if (mu.muonBestTrack().isNonnull())
      raw.push_back(*mu.muonBestTrack());
  }

  // Apply sip2D cuts uniformly to every track in the merged collection
  for (const auto& track : raw) {
    if (track.pt() <= minPt_) continue;
    if (track.normalizedChi2() >= maxNormalizedChi2_) continue;
    reco::TransientTrack ttrack = ttBuilder.build(track);
    auto ip2dResult = IPTools::signedTransverseImpactParameter(
        ttrack, GlobalVector(track.px(), track.py(), track.pz()), pv);
    if (!ip2dResult.first) continue;
    if (std::fabs(ip2dResult.second.significance()) < minAbsSip2D_) continue;
    outputTracks.push_back(track);
  }
}

void MiniAODTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get input collections
  edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken_, pfCandidatesHandle);

  edm::Handle<pat::PackedCandidateCollection> lostTracksHandle;
  iEvent.getByToken(lostTracksToken_, lostTracksHandle);

  edm::Handle<pat::PackedCandidateCollection> eleLostTracksHandle;
  iEvent.getByToken(eleLostTracksToken_, eleLostTracksHandle);

  edm::Handle<std::vector<pat::Muon>> muonsHandle;
  iEvent.getByToken(muonsToken_, muonsHandle);

  edm::Handle<std::vector<pat::Muon>> displacedMuonsHandle;
  iEvent.getByToken(displacedMuonsToken_, displacedMuonsHandle);

  // Always fetch TransientTrackBuilder and PV (needed for mergedMuonSip2D regardless of applyCuts_)
  const TransientTrackBuilder& ttBuilderRef = iSetup.getData(ttBuilderToken_);
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  const reco::Vertex* pv = (pvHandle.isValid() && !pvHandle->empty()) ? &pvHandle->at(0) : nullptr;

  // For the applyCuts_ path: pass pointer (null when disabled preserves existing behavior)
  const TransientTrackBuilder* ttBuilder = applyCuts_ ? &ttBuilderRef : nullptr;

  // Create output collections
  auto pfCandidateTracks = std::make_unique<reco::TrackCollection>();
  auto lostTracks = std::make_unique<reco::TrackCollection>();
  auto eleLostTracks = std::make_unique<reco::TrackCollection>();
  auto mergedTracks = std::make_unique<reco::TrackCollection>();
  auto mergedTracksWithEle = std::make_unique<reco::TrackCollection>();
  auto mergedTracksAll = std::make_unique<reco::TrackCollection>();
  auto muonGlobalTracks = std::make_unique<reco::TrackCollection>();
  auto displacedMuonGlobalTracks = std::make_unique<reco::TrackCollection>();
  auto muonLLPNanoTracks = std::make_unique<reco::TrackCollection>();
  auto displacedMuonLLPNanoTracks = std::make_unique<reco::TrackCollection>();
  auto mergedMuonSip2D = std::make_unique<reco::TrackCollection>();

  // Extract tracks from packed candidates
  if (pfCandidatesHandle.isValid()) {
    extractTracks(*pfCandidatesHandle, *pfCandidateTracks, ttBuilder, pv);
  }

  if (lostTracksHandle.isValid()) {
    extractTracks(*lostTracksHandle, *lostTracks, ttBuilder, pv);
  }

  if (eleLostTracksHandle.isValid()) {
    extractTracks(*eleLostTracksHandle, *eleLostTracks, ttBuilder, pv);
  }

  // Extract global tracks from muons
  if (muonsHandle.isValid()) {
    extractMuonTracks(*muonsHandle, *muonGlobalTracks);
    extractMuonTracksLLPNano(*muonsHandle, *muonLLPNanoTracks);
  }

  if (displacedMuonsHandle.isValid()) {
    extractMuonTracks(*displacedMuonsHandle, *displacedMuonGlobalTracks);
    extractMuonTracksLLPNano(*displacedMuonsHandle, *displacedMuonLLPNanoTracks);
  }

  // Build PF electron reference tracks for deduplication (same quality cuts as extractTracks)
  reco::TrackCollection pfElectronTracks;
  if (pfCandidatesHandle.isValid()) {
    for (const auto& cand : *pfCandidatesHandle) {
      if (!cand.hasTrackDetails() || std::abs(cand.pdgId()) != 11) continue;
      const reco::Track track = cand.pseudoTrack();
      if (applyCuts_) {
        if (track.pt() <= minPt_) continue;
        if (track.normalizedChi2() >= maxNormalizedChi2_) continue;
        if (ttBuilder != nullptr && pv != nullptr) {
          reco::TransientTrack ttrack = ttBuilder->build(track);
          auto ip2dResult = IPTools::signedTransverseImpactParameter(
              ttrack, GlobalVector(track.px(), track.py(), track.pz()), *pv);
          if (ip2dResult.first && std::fabs(ip2dResult.second.significance()) < minAbsSip2D_) continue;
        }
      }
      pfElectronTracks.push_back(track);
    }
  }

  // Build merged collections: keep all PF candidates, deduplicate lost tracks vs PF electrons.
  // Lost tracks within deltaR < 0.01 of any PF electron (abs(pdgId)==11) are dropped.

  // merged = pfCandidateTracks + lostTracks (deduped vs PF electrons)
  mergedTracks->insert(mergedTracks->end(), pfCandidateTracks->begin(), pfCandidateTracks->end());
  addTracksDeduped(*lostTracks, pfElectronTracks, *mergedTracks);

  // mergedWithEle = pfCandidateTracks + eleLostTracks (deduped vs PF electrons)
  mergedTracksWithEle->insert(mergedTracksWithEle->end(), pfCandidateTracks->begin(), pfCandidateTracks->end());
  addTracksDeduped(*eleLostTracks, pfElectronTracks, *mergedTracksWithEle);

  // mergedAll = pfCandidateTracks + lostTracks (deduped) + eleLostTracks (deduped)
  mergedTracksAll->insert(mergedTracksAll->end(), pfCandidateTracks->begin(), pfCandidateTracks->end());
  addTracksDeduped(*lostTracks, pfElectronTracks, *mergedTracksAll);
  addTracksDeduped(*eleLostTracks, pfElectronTracks, *mergedTracksAll);

  // mergedMuonSip2D: PF+Lost (no muon-ID'd tracks) deduped vs slimmedMuons,
  // union with all slimmedMuon bestTracks, sip2D cuts applied to everything.
  if (pfCandidatesHandle.isValid() && lostTracksHandle.isValid() &&
      muonsHandle.isValid() && pv != nullptr) {
    buildMuonMergedSip2D(*pfCandidatesHandle, *lostTracksHandle,
                         *muonsHandle, ttBuilderRef, *pv, *mergedMuonSip2D);
  }

  // Log some stats
  edm::LogInfo("MiniAODTrackProducer")
      << "Extracted tracks - PF: " << pfCandidateTracks->size()
      << ", Lost: " << lostTracks->size()
      << ", EleLost: " << eleLostTracks->size()
      << ", Merged: " << mergedTracks->size()
      << ", MergedWithEle: " << mergedTracksWithEle->size()
      << ", MergedAll: " << mergedTracksAll->size()
      << ", MuonGlobal: " << muonGlobalTracks->size()
      << ", DisplacedMuonGlobal: " << displacedMuonGlobalTracks->size()
      << ", MuonLLPNano: " << muonLLPNanoTracks->size()
      << ", DisplacedMuonLLPNano: " << displacedMuonLLPNanoTracks->size()
      << ", MergedMuonSip2D: " << mergedMuonSip2D->size();

  // Put collections into the event
  iEvent.put(std::move(pfCandidateTracks), "pfCandidateTracks");
  iEvent.put(std::move(lostTracks), "lostTracks");
  iEvent.put(std::move(eleLostTracks), "eleLostTracks");
  iEvent.put(std::move(mergedTracks), "merged");
  iEvent.put(std::move(mergedTracksWithEle), "mergedWithEle");
  iEvent.put(std::move(mergedTracksAll), "mergedAll");
  iEvent.put(std::move(muonGlobalTracks), "muonGlobalTracks");
  iEvent.put(std::move(displacedMuonGlobalTracks), "displacedMuonGlobalTracks");
  iEvent.put(std::move(muonLLPNanoTracks), "muonLLPNanoTracks");
  iEvent.put(std::move(displacedMuonLLPNanoTracks), "displacedMuonLLPNanoTracks");
  iEvent.put(std::move(mergedMuonSip2D), "mergedMuonSip2D");
}

void MiniAODTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collections (MINIAOD defaults)
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("eleLostTracks", edm::InputTag("lostTracks", "eleTracks"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("displacedMuons", edm::InputTag("slimmedDisplacedMuons"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlineSlimmedPrimaryVertices"));

  // Cut configuration
  desc.add<bool>("applyCuts", false);  // Disabled by default for backwards compatibility
  desc.add<double>("minPt", 1.0);
  desc.add<double>("minAbsSip2D", 4.0);  // Keep |sip2D| >= this value (displaced tracks)
  desc.add<double>("maxNormalizedChi2", 5.0);

  descriptions.add("miniAODTrackProducer", desc);
}

DEFINE_FWK_MODULE(MiniAODTrackProducer);
