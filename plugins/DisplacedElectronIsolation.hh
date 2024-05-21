#ifndef DisplacedElectronIsolation_hh
#define DisplacedElectronIsolation_hh
// -*- C++ -*-
//
//
// Original Author:  Andres Abreu
// 
// DisplacedElectronIsolation class
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

// Includes for all object types

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

// Add includes for interface collections
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedTrackSCPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexAssembly.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

// KUCMS Object includes
#include "KUCMSObjectBase.hh"
#include "KUCMSBranchManager.hh"

#define DEBUG false

class DisplacedElectronIsolation : public KUCMSObjectBase {

public:

  // use base class constructor
  DisplacedElectronIsolation( const edm::ParameterSet& iConfig );
  virtual ~DisplacedElectronIsolation() = default; 

  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadGeneralTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ generalTracksToken_ = token; } 
  void LoadGsfTracksToken( edm::EDGetTokenT<reco::GsfTrackCollection> token ){ gsfTracksToken_ = token; }
  void LoadECALTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ ecalTracksToken_ = token; }
  void LoadDisplacedElectrons( edm::EDGetTokenT<reco::ElectronCollection> token ) {electronToken_ = token;}
  void LoadPrimaryVertex( edm::EDGetTokenT<reco::VertexCollection> token ) {pvToken_ = token;}
  void LoadPFCandidatesToken( edm::EDGetTokenT<reco::PFCandidateCollection> token ){ pfCandidatesToken_ = token; }
  void LoadBeamSpot(edm::EDGetTokenT<reco::BeamSpot> token) { beamspotToken_ = token; }
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }  
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> token ){ genToken_ = token; }

  //void LoadDisplacedElectronObj(const KUCMSDisplacedElectron* electrons) { displacedElectronObj_ = electrons; }
  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree );

  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar ) {}
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs() {} // do any jobs that need to be done after main event loop
  
  void CalculateIsolation(MatchedTrackSCPairs<reco::TransientTrack> &candidates);

private:

  reco::TrackCollection mergedTracks_;

  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
  edm::Handle<reco::TrackCollection> generalTracksHandle_;

  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  edm::Handle<reco::GsfTrackCollection> gsfTracksHandle_;

  edm::EDGetTokenT<reco::TrackCollection> ecalTracksToken_;
  edm::Handle<reco::TrackCollection> ecalTracksHandle_;

  edm::EDGetTokenT<reco::ElectronCollection> electronToken_;
  edm::Handle<reco::ElectronCollection> electronHandle_;

  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
  edm::Handle<reco::PFCandidateCollection> pfCandidatesHandle_;

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::Handle<reco::GenParticleCollection> genHandle_;

  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  reco::BeamSpot beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  reco::Vertex primaryVertex_;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  reco::TrackCollection MergeTracks(const reco::TrackCollection &generalTracks, 
				    const reco::GsfTrackCollection &gsfTracks, 
				    double deltaR = 0.1) const;

};//<<>>class DisplacedElectronIsolation : public KUCMSObjectBase

DisplacedElectronIsolation::DisplacedElectronIsolation( const edm::ParameterSet& iConfig ) 
  : caloGeometryToken_(esConsumes()) {   
// ---- end constructor initilizations  --------------------------
  trackAssociator_.useDefaultPropagator();

  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
  cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );

}//<<>>DisplacedElectronIsolation::DisplacedElectronIsolation( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void DisplacedElectronIsolation::InitObject( TTree* fOutTree ) {
  
  Branches.makeBranch("TrackMerger_nMatches","TrackMerger_nMatches", UINT);
  Branches.makeBranch("TrackMerger_deltaR", "TrackMerger_deltaR", VFLOAT);
  Branches.makeBranch("TrackMergerGeneral_pt", "TrackMergerGeneral_pt", VFLOAT);
  Branches.makeBranch("TrackMergerGsf_pt", "TrackMergerGsf_pt", VFLOAT);
  Branches.makeBranch("TrackMerger_ptDiff", "TrackMerger_ptDiff", VFLOAT);
  Branches.makeBranch("TrackMerger_vzDiff", "TrackMerger_vzDiff", VFLOAT);
  Branches.makeBranch("TrackMerger_displacedElectronIndex","TrackMerger_displacedElectronIndex", VINT);
  //Branches.makeBranch("TrackMerger_","TrackMerger_", INT);

  Branches.attachBranches(fOutTree);
}

void DisplacedElectronIsolation::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  // Get handles on the input edm collections
  iEvent.getByToken( generalTracksToken_, generalTracksHandle_);
  iEvent.getByToken( gsfTracksToken_, gsfTracksHandle_);
  iEvent.getByToken( ecalTracksToken_, ecalTracksHandle_);
  iEvent.getByToken( electronToken_, electronHandle_);
  iEvent.getByToken( pfCandidatesToken_, pfCandidatesHandle_);
  iEvent.getByToken(genToken_, genHandle_);
  beamSpot_ = iEvent.get(beamspotToken_);
  primaryVertex_ = iEvent.get(pvToken_).at(0);
  /*
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  VertexAssembly vertexBuilder(ttBuilder, generalTracksHandle_);
  reco::VertexCollection generalVertices = vertexBuilder.CreateVertexCollection(0.95);

  reco::TrackCollection tracks;
  for(const auto &track : *generalTracksHandle_)
    if(track.pt() > 0.95)
      tracks.emplace_back(track);

  TrackPropagator<reco::Track> generalTrackPropagator(iEvent, iSetup, magfield, trackAssocParameters_, tracks);
  reco::TrackCollection ecalTracks = generalTrackPropagator.GetTracks();  

  int sigCount(0);
  for(const auto &gen : *genHandle_) {
    if(abs(gen.pdgId()) == 11 && gen.status() == 1)
      if(isSignalGenElectron(gen))
	sigCount++;
  }
  
  //DeltaRGenMatchHungarian<reco::Track> assigner(ecalTracks, *genHandle_);
  //const auto genMatchedPairs = assigner.GetPairedObjects();

  //reco::TrackCollection signalMatches = GetSignalMatches<reco::Track>(ecalTracks, *genHandle_);
  
  for(const auto &track : signalMatches) {
    GlobalVector direction(track.px(), track.py(), track.pz());
    std::pair<bool, Measurement1D> ip2D = IPTools::signedTransverseImpactParameter(ttBuilder->build(track), direction, primaryVertex_);
    double sip2D = ip2D.second.value() / ip2D.second.error();

    cout << "The signal candidate track has a dxy = " << ip2D.second.value() << " and sip2D = " << sip2D << endl;
  }
  */
  /*
  if(sigCount>0) 
    cout << "This event has " << sigCount << " signal electrons!" << endl;
  if(signalMatches.size() > 0)
    cout << "Found " << signalMatches.size() << " signal tracks in the ecal tracks!" << endl;
  */
  /*
  if(signalMatches.size() > 0) {
    cout << "This event has " << generalVertices.size() << " vertices" << std::endl;
    cout << "The pv has " << primaryVertex_.tracksSize() << " tracks and a Pt: " << VertexHelper::CalculateTotalPt(primaryVertex_) << endl;
    cout << "This event has " << sigCount << " signal electrons!" << endl;
    for(const auto &vertex : generalVertices) {
      double ecalness = VertexHelper::CalculateEcalness(vertex, ecalTracks);
      int signalInstances = VertexHelper::CountInstances(vertex, signalMatches);   
      if(ecalness > 0 && signalInstances > 0) {
	cout << "  This vertex has " << vertex.tracksSize() << " tracks" << std::endl;
	cout << "  Total pt =  " << VertexHelper::CalculateTotalPt(vertex) << std::endl;
	cout << "  Ecalness = " << ecalness << endl;
	cout << "  Found " << signalInstances << " signal tracks in this event!" << endl; 

	reco::TrackCollection vtxTracks;
	for(const auto &track : vertex.tracks()) 
	  vtxTracks.emplace_back(*track);
	  
	signalMatches = GetSignalMatches<reco::Track>(vtxTracks, *genHandle_);

	for(const auto &track : signalMatches) {
	  GlobalVector direction(track.px(), track.py(), track.pz());
	  std::pair<bool, Measurement1D> ip2D = IPTools::signedTransverseImpactParameter(ttBuilder->build(track), direction, primaryVertex_);
	  double sip2D = ip2D.second.value() / ip2D.second.error();
	  cout << "The signal candidate track has a dxy = " << ip2D.second.value() << " and sip2D = " << sip2D << endl;
	}
      }
    }
  }
  */
  if( DEBUG ) std::cout << "Collecting Examples" << std::endl;

  mergedTracks_ = MergeTracks(*generalTracksHandle_, *gsfTracksHandle_);

}//<<>>void DisplacedElectronIsolation::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void DisplacedElectronIsolation::PostProcessEvent( ItemManager<float>& geVar ) {

  Branches.clearBranches();

  const auto matches = GetOverlap<reco::Track, reco::GsfTrack>(*generalTracksHandle_, *gsfTracksHandle_, 0.1);

  Branches.fillBranch("TrackMerger_nMatches", unsigned(matches.size()) );
  for(const auto &match : matches) {
    
    const float deltaR(match.GetDeltaR());
    const reco::Track track(match.GetObjectA());
    const reco::GsfTrack gsfTrack(match.GetObjectB());

    int index = -1;
    for(size_t e(0); e < electronHandle_->size(); e++) {
      reco::Electron candidate(electronHandle_->at(e));
      if(candidate.pt() == gsfTrack.pt() && candidate.eta() == gsfTrack.eta())
	index = e;
    }

    Branches.fillBranch("TrackMerger_deltaR", deltaR);
    Branches.fillBranch("TrackMergerGeneral_pt", float(track.pt()) );
    Branches.fillBranch("TrackMergerGsf_pt", float(gsfTrack.pt()) );
    Branches.fillBranch("TrackMerger_ptDiff", float(gsfTrack.pt() - track.pt()) );
    Branches.fillBranch("TrackMerger_vzDiff", float(gsfTrack.vz() - track.vz()) );
    Branches.fillBranch("TrackMerger_displacedElectronIndex", int(index));
  }
}

void DisplacedElectronIsolation::CalculateIsolation(MatchedTrackSCPairs<reco::TransientTrack> &candidates) {

  std::vector<double> trackIso0p4 = IsolationInfo::GetTrackIsolation(*generalTracksHandle_, beamSpot_);
  std::vector<double> pfIso0p4 = IsolationInfo::GetPFIsolation(*generalTracksHandle_, *pfCandidatesHandle_);
  std::vector<double> scIso0p4 = IsolationInfo::GetSCIsolation(*electronHandle_, candidates.GetSCs());

  for(unsigned int t(0); t < generalTracksHandle_->size(); t++) {
    const reco::Track track(mergedTracks_[t]);
    for(auto &pair : candidates) {
      const auto candidateTrack(pair.GetTrack().track());
      if(track.pt() == candidateTrack.pt() && track.eta() == candidateTrack.eta()) {
        pair.AddIsoInfo("trackIso0p4", trackIso0p4[t]);
	pair.AddIsoInfo("pfIso0p4", pfIso0p4[t]);
      }
    }
  }
  
  int count(0);
  for(auto &pair : candidates) {
    reco::Track track(pair.GetTrack().track());
    GlobalVector direction(track.px(), track.py(), track.pz());
    reco::Vertex beamspotvertex(beamSpot_.position(), beamSpot_.covariance3D());
    std::pair<bool, Measurement1D> ip3D = IPTools::signedImpactParameter3D(pair.GetTrack(), direction, primaryVertex_);
    std::pair<bool, Measurement1D> ip2D = IPTools::signedTransverseImpactParameter(pair.GetTrack(), direction, primaryVertex_);
    std::pair<bool, Measurement1D> ip3Dbs = IPTools::signedImpactParameter3D(pair.GetTrack(), direction, beamspotvertex);
    std::pair<bool, Measurement1D> ip2Dbs = IPTools::signedTransverseImpactParameter(pair.GetTrack(), direction, beamspotvertex);
    pair.AddIsoInfo("ip3D", ip3D.second.value());
    pair.AddIsoInfo("ip2D", ip2D.second.value());
    pair.AddIsoInfo("sip3D", ip3D.second.value() / ip3D.second.error());
    pair.AddIsoInfo("sip2D", ip2D.second.value() / ip2D.second.error());
    pair.AddIsoInfo("ip3Dbs", ip3Dbs.second.value());
    pair.AddIsoInfo("ip2Dbs", ip2Dbs.second.value());
    pair.AddIsoInfo("sip3Dbs", ip3Dbs.second.value() / ip3Dbs.second.error());
    pair.AddIsoInfo("sip2Dbs", ip2Dbs.second.value() / ip2Dbs.second.error());
    pair.AddIsoInfo("scIso0p4", scIso0p4[count]);
    count++;
  }
}

reco::TrackCollection DisplacedElectronIsolation::MergeTracks(const reco::TrackCollection &generalTracks, 
							      const reco::GsfTrackCollection &gsfTracks, 
							      double deltaR) const {

  reco::TrackCollection outputTracks = RemoveOverlap<reco::Track, reco::GsfTrack>(generalTracks, gsfTracks, deltaR);

  for (const auto &track : gsfTracks) {
    outputTracks.emplace_back(track);
  }

  return outputTracks;
}


#endif
