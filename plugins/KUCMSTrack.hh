// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Prompt Track Object
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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

// Add includes for specfic collections
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
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

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSGenObjects.hh"
#include "KUCMSDisplacedVertex.hh"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#ifndef KUCMSTrackObjectHeader
#define KUCMSTrackObjectHeader

//#define TrackDEBUG true
#define TrackDEBUG false

using namespace edm; 

class KUCMSTrackObject : public KUCMSObjectBase {

    public:

	typedef pat::IsolatedTrack::PolarLorentzVector PolarLorentzVector;

    // use base class constructor
    KUCMSTrackObject( const edm::ParameterSet& iConfig );
    //KUCMSTrackObject( const edm::ParameterSet& iConfig );
    ~KUCMSTrackObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    //void LoadConversionTokens( edm::EDGetTokenT<reco::ConversionCollection> conversionsToken ){ conversionsToken_ = conversionsToken; };
    void LoadBeamSpotTokens( edm::EDGetTokenT<reco::BeamSpot> beamLineToken ){ beamLineToken_ = beamLineToken; };
    //void LoadRhoTokens( edm::EDGetTokenT<double> RhoToken ){ RhoToken_ = RhoToken; };
    void LoadGeneralTrackTokens( edm::EDGetTokenT<edm::View<reco::Track>> token ){ generalTracksToken_ = token; }
    void LoadGsfTrackTokens( edm::EDGetTokenT<edm::View<reco::GsfTrack>> token ){ gsfTracksToken_ = token; }
    void LoadVertexTokens( edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_ ){ verticesToken = verticesToken_; };
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
	//void LoadDisplacedVertexObject( KUCMSDisplacedVertex* svObj_ ){ svObj = svObj_; };
    void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }
    void LoadPfcandTokens( edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken ){ pfcandToken_ = pfcandToken; };
  	void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters = parameters;}
  	void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }
    void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
    // do cross talk jobs with other objects, do event processing, and load branches
    void ProcessEvent( ItemManager<float>& geVar );
    void PostProcessEvent( ItemManager<float>& geVar );

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...
	void getTrkIso(const PolarLorentzVector& p4, int pc_idx, pat::PFIsolation& iso, pat::PFIsolation& miniiso) const;

    private:

  	// General Tracks
  	edm::EDGetTokenT<edm::View<reco::Track>> generalTracksToken_;
  	edm::Handle<edm::View<reco::Track> > generalTracksHandle_;

    // GsfTracks
    edm::EDGetTokenT<edm::View<reco::GsfTrack>> gsfTracksToken_;
    edm::Handle<edm::View<reco::GsfTrack> > gsfTracksHandle_;

    // Tracks
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
    TTBuilderWrapper ttBuilder;

	// Track Associator
  	TrackDetectorAssociator trackAssociator;
  	TrackAssociatorParameters trackAssocParameters;

	// Magnetic field
  	edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
	edm::ESTransientHandle<MagneticField> magfield_;

    // BeamSpot ----------------------------------------------------------
    //const edm::InputTag beamSpotTag;
    edm::EDGetTokenT<reco::BeamSpot> beamLineToken_;
    edm::Handle<reco::BeamSpot> beamSpot_;

	// PV
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken;
    edm::Handle<std::vector<reco::Vertex>> vertices_;

    // PF Candidates View 
    //const edm::InputTag pfcandTag;
    edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken_;
    edm::Handle<edm::View<reco::Candidate>> pfcands_;

	//Track Collection
	std::vector<reco::Track> prmtTracks;
    std::vector<bool> isGenTrack;// true general, false gsf 
    std::vector<TrackDetMatchInfo> detIdInfo;

	// other object collections 
    KUCMSEcalRecHitObject* rhObj;

};//<<>>class KUCMSTrack : public KUCMSObjectBase

KUCMSTrackObject::KUCMSTrackObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "minEleE", iConfig.existsAs<double>("minEleE") ? iConfig.getParameter<double>("minEleE") : 2.0 );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );
    cfPrm.set( "minEleE", iConfig.existsAs<double>("minEleE") ? iConfig.getParameter<double>("minEleE") : 2.0 );
    cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );
    prmStr.set( "RunEra", iConfig.existsAs<std::string>("RunEra") ? iConfig.getParameter<std::string>("RunEra") : "Run2" );

}//<<>>KUCMSTrack::KUCMSTrack( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSTrackObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("mpdg","Track_pdgId",VINT);
    Branches.makeBranch("pt","Track_pt",VFLOAT);
    Branches.makeBranch("eta","Track_eta",VFLOAT);
    Branches.makeBranch("phi","Track_phi",VFLOAT);
    Branches.makeBranch("vx","Track_vx",VFLOAT);
    Branches.makeBranch("vy","Track_vy",VFLOAT);
    Branches.makeBranch("vz","Track_vz",VFLOAT);
    //Branches.makeBranch("beta","Track_beta",VFLOAT);
    Branches.makeBranch("p","Track_p",VFLOAT);
    Branches.makeBranch("charge","Track_charge",VINT);

    Branches.makeBranch("px","Track_px",VFLOAT);//float(track.px()) );
    Branches.makeBranch("py","Track_py",VFLOAT);//float(track.py()) );
    Branches.makeBranch("pz","Track_pz",VFLOAT);//float(track.pz()) );
    Branches.makeBranch("theta","Track_theta",VFLOAT);//float(track.theta()) );
    Branches.makeBranch("lambda","Track_lambda",VFLOAT);//float(track.lambda()) );
    Branches.makeBranch("dxy","Track_dxy",VFLOAT);//float(track.dxy()) );
    Branches.makeBranch("ndof","Track_ndof",VFLOAT);//float(track.ndof()) );
    Branches.makeBranch("chi2","Track_chi2",VFLOAT);//float(track.chi2()) );
    Branches.makeBranch("normalizedChi2","Track_normalizedChi2",VFLOAT);//float(track.normalizedChi2()) );
    Branches.makeBranch("qoverp","Track_qoverp",VFLOAT);//float(track.qoverp()) );
    Branches.makeBranch("ptError","Track_ptError",VFLOAT);//float(track.ptError()) );
    Branches.makeBranch("etaError","Track_etaError",VFLOAT);//float(track.etaError()) );
    Branches.makeBranch("phiError","Track_phiError",VFLOAT);//float(track.phiError()) );
    Branches.makeBranch("qualityMask","Track_qualityMask",VINT);//int(track.qualityMask()) );
    Branches.makeBranch("nValidHits","Track_nValidHits",VINT);//int(track.numberOfValidHits()) );
    Branches.makeBranch("nLostHits","Track_nLostHits",VINT);//int(track.numberOfLostHits()) );

    Branches.makeBranch("absiso","Track_absIso",VFLOAT);
    Branches.makeBranch("reliso","Track_relIso",VFLOAT);
    Branches.makeBranch("minireliso","Track_miniRelIso",VFLOAT);

    Branches.makeBranch("relisolep","Track_relisolep",VFLOAT);
    Branches.makeBranch("mrelisolep","Track_mrelisolep",VFLOAT);

    Branches.makeBranch("ip2D","Track_ip2D",VFLOAT);
    Branches.makeBranch("sip2D","Track_sip2D",VFLOAT);
    Branches.makeBranch("ip3D","Track_ip3D",VFLOAT);
    Branches.makeBranch("sip3D","Track_sip3D",VFLOAT);

    Branches.makeBranch("isGen","Track_isGenTrack",VBOOL);
    Branches.makeBranch("nIsoTracks","Track_nIsoTracks",INT);

    Branches.makeBranch("ecalDetIds","Track_ecalDetIds",VVUINT);
    Branches.makeBranch("scIndex","Track_scIndexs",VVINT);
    Branches.makeBranch("hcalDetIds","Track_hcalDetIds",VVUINT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSTrack::InitObject( TTree* fOutTree )

void KUCMSTrackObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){


    if( TrackDEBUG ) std::cout << "Getting Tokens in Tracks" << std::endl;

    iEvent.getByToken( beamLineToken_, beamSpot_ );
    iEvent.getByToken( pfcandToken_, pfcands_ );
  	iEvent.getByToken( generalTracksToken_, generalTracksHandle_ );
    iEvent.getByToken( gsfTracksToken_, gsfTracksHandle_ );
    iEvent.getByToken( verticesToken, vertices_ );

    ttBuilder = TTBuilderWrapper( &iSetup.getData(transientTrackBuilder_) );
  	magfield_ = iSetup.getTransientHandle(magneticFieldToken_);
	trackAssociator.useDefaultPropagator();
	//trackAssocParameters

    if( TrackDEBUG ) std::cout << "Collecting Tracks" << std::endl;

	float ptmin = 0.0;
	prmtTracks.clear();
    isGenTrack.clear();
	detIdInfo.clear();
	for( const auto &track : *generalTracksHandle_ ){ 

        if( track.pt() < ptmin ) continue;
		prmtTracks.push_back( track );
		isGenTrack.push_back(true);
		FreeTrajectoryState initialState = trajectoryStateTransform::initialFreeState( track, magfield_.product() );
		TrackDetMatchInfo detInfo = trackAssociator.associate( iEvent, iSetup, trackAssocParameters, &initialState );
		detIdInfo.push_back( detInfo );

	}//<<>>for(const auto &track : *generalTracksHandle_ )

    for( const auto &track : *gsfTracksHandle_ ){

        if( track.pt() < ptmin ) continue;
		bool dup = false;
		for( const auto &gtrack : prmtTracks ){ 
			if( TrackHelper::SameTrack( track, gtrack ) ){ dup = true; break; }
		}//<<>>for( auto & gtrack : prmtTracks )
		if( dup ) continue;
        prmtTracks.push_back( track );
		isGenTrack.push_back(false);
        FreeTrajectoryState initialState = trajectoryStateTransform::initialFreeState( track, magfield_.product() );
        TrackDetMatchInfo detInfo = trackAssociator.associate( iEvent, iSetup, trackAssocParameters, &initialState );
        detIdInfo.push_back( detInfo );

    }//<<>>for(const auto &track : *generalTracksHandle_ )


}//<<>>void KUCMSTrack::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSTrackObject::PostProcessEvent( ItemManager<float>& geVar ){

}//<<>>void KUCMSTrackObject::PostProcessEvent( ItemManager<float>& geVar )

void KUCMSTrackObject::ProcessEvent( ItemManager<float>& geVar ){

    if( TrackDEBUG ) std::cout << "Processing Tracks" << std::endl;

    Branches.clearBranches();

    int nIsoTracks = 0;
	int trackIndx = -1;

    PolarLorentzVector polarP4;
	double zero(0);
	//float rimax = 0.2;//0.2
    //float minmax = 0.2;//0.2
    //float absval = 5;//5
	float rilepmax = 0.05;
    //float rilepmax = 0.1;
    //float rilepmax = 0.2;
    float trkptmin = 5.0;
	float isoptmin = 20.0;
	for( const auto &track : prmtTracks ){

		trackIndx++;

		int idx = -1;	
        int matched = -1;
		int mpdg = 0;
        for(const auto &pfcan : *pfcands_ ){
			idx++;
            auto pftrackPtr = pfcan.bestTrack();
            if( pftrackPtr == nullptr ) continue;
			int id = std::abs( pfcan.pdgId() );
            if( TrackHelper::SameTrack( track, *pftrackPtr ) ){ matched = idx; mpdg = id; break; }
        }//<<>>for(const auto &pfcan : *pfcands_ )
        //if( mpdg < 0 ) continue;

        float m = 0;  //assume ~ massless
        float pt = track.pt();
		float eta = track.eta();
		float phi = track.phi();
		float vx = track.vx();
        float vy = track.vy();
        float vz = track.vz();
		//float beta = track.beta();
		float p = track.p();
		int charge = track.charge();

        std::vector<DetId> ecalDetIDs = detIdInfo[trackIndx].crossedEcalIds;
        std::vector<uInt> ecalRawIds;
        std::vector<int> scIndex{-9};
        if( ecalDetIDs.size() > 0 ){
            for( auto & id : ecalDetIDs ){ ecalRawIds.push_back( id.rawId() ); }
            scIndex = rhObj->getSuperClusterIndex( ecalRawIds, trackIndx );
        }//<<>>if( ecalDetIDs.size() > 0 )

		if( ( scIndex[0] < 0 ) && ( pt < trkptmin ) ) continue;
        Branches.fillBranch("ecalDetIds", ecalRawIds );
        Branches.fillBranch("scIndex", scIndex );

        std::vector<DetId> hcalDetIDs = detIdInfo[trackIndx].crossedHcalIds;
        std::vector<uInt> hcalRawIds;
        if( hcalDetIDs.size() > 0 ){
            for( auto & id : hcalDetIDs ){ hcalRawIds.push_back( id.rawId() ); }
        }//<<>>if( ecalDetIDs.size() > 0 )
        Branches.fillBranch("hcalDetIds", hcalRawIds );


        polarP4.SetCoordinates( pt, eta, phi, m);
        pat::PFIsolation pfIso;
        pat::PFIsolation miniIso;
        getTrkIso( polarP4, matched, pfIso, miniIso );
        float absiso = pfIso.chargedHadronIso();
        float reliso = pfIso.chargedHadronIso() / polarP4.pt();
        float minireliso = miniIso.chargedHadronIso() / polarP4.pt();
        float rilpu = std::max( zero, pfIso.neutralHadronIso() + pfIso.photonIso() - 0.5*pfIso.puChargedHadronIso() );
        float relisolep =  ( pfIso.chargedHadronIso() + rilpu ) / polarP4.pt();
        float mrilpu = std::max( zero, miniIso.neutralHadronIso() + miniIso.photonIso() - 0.5*miniIso.puChargedHadronIso() );
        float mrelisolep = ( miniIso.chargedHadronIso() + mrilpu ) / polarP4.pt();

        bool ptcut = polarP4.pt() > isoptmin;
        //bool isocut = relisolep < rilepmax;
        bool isocut = mrelisolep < rilepmax;
        bool trkiso = ptcut && isocut;
        if( trkiso ) nIsoTracks++;
		//if( trkiso ) std::cout << " trk iso: : " << absiso << " " << reliso << " " << mrelisolep << " " << relisolep << std::endl;

    	const auto & primevtx = vertices_->front();
        auto ttrack = ttBuilder.build( track );
        GlobalVector direction(track.px(), track.py(), track.pz());

        std::pair<bool, Measurement1D> ip2Dtool = IPTools::signedTransverseImpactParameter(ttrack, direction, primevtx );
		float ip2d = ip2Dtool.second.value();
        float sip2D = std::abs( ip2d / ip2Dtool.second.error() );
        Branches.fillBranch("ip2D",ip2d);
		Branches.fillBranch("sip2D",sip2D);

        std::pair<bool, Measurement1D> ip3Dtool = IPTools::signedImpactParameter3D(ttrack, direction, primevtx );
        float ip3d = std::abs( ip3Dtool.second.value() );
        float sip3D = std::abs( ip3d / ip3Dtool.second.error() );
        Branches.fillBranch("ip3D",ip3d);
        Branches.fillBranch("sip3D",sip3D);

        Branches.fillBranch("mpdg",mpdg);
        Branches.fillBranch("pt",pt);
        Branches.fillBranch("eta",eta);
        Branches.fillBranch("phi",phi);
        Branches.fillBranch("vx",vx);
        Branches.fillBranch("vy",vy);
        Branches.fillBranch("vz",vz);
        //Branches.fillBranch("beta",beta);
        Branches.fillBranch("p",p);
        Branches.fillBranch("charge",charge);

		Branches.fillBranch("isGen", isGenTrack[trackIndx] );
        Branches.fillBranch("absiso",absiso);
        Branches.fillBranch("reliso",reliso);
        Branches.fillBranch("minireliso",minireliso);

		Branches.fillBranch("relisolep",relisolep);
        Branches.fillBranch("mrelisolep",mrelisolep);

        Branches.fillBranch("px", float(track.px()) );
        Branches.fillBranch("py", float(track.py()) );
        Branches.fillBranch("pz", float(track.pz()) );
        Branches.fillBranch("theta", float(track.theta()) );
        Branches.fillBranch("lambda", float(track.lambda()) );
        Branches.fillBranch("dxy", float(track.dxy()) );
        Branches.fillBranch("ndof", float(track.ndof()) );
        Branches.fillBranch("chi2", float(track.chi2()) );
        Branches.fillBranch("normalizedChi2", float(track.normalizedChi2()) );
        Branches.fillBranch("qoverp", float(track.qoverp()) );
        Branches.fillBranch("ptError", float(track.ptError()) );
        Branches.fillBranch("etaError", float(track.etaError()) );
        Branches.fillBranch("phiError", float(track.phiError()) );
        Branches.fillBranch("qualityMask", int(track.qualityMask()) );
        Branches.fillBranch("nValidHits", int(track.numberOfValidHits()) );
        Branches.fillBranch("nLostHits", int(track.numberOfLostHits()) );

    }//<<>>for(const auto &track : *generalTracksHandle_ )

  	//std::cout << std::endl << " ----------------------------------"  << std::endl;
   	Branches.fillBranch("nIsoTracks",nIsoTracks);

}//<<>>void KUCMSTrack::ProcessEvent()

void KUCMSTrackObject::EndJobs(){}

// paramters form PhysicsTools/PatAlgos/pythin/slimming/isolatedTracks_cfi.py
void KUCMSTrackObject::getTrkIso(const PolarLorentzVector& p4, int pc_idx, pat::PFIsolation& iso, pat::PFIsolation& miniiso) const {

  	const auto & primevtx = vertices_->front();
  	float chiso = 0, nhiso = 0, phiso = 0, puiso = 0;      // standard isolation
  	float chmiso = 0, nhmiso = 0, phmiso = 0, pumiso = 0;  // mini isolation
  	// miniIsoParams = cms.vdouble(0.05, 0.2, 10.0),
  	//float miniDR = std::max(miniIsoParams_[0], std::min(miniIsoParams_[1], miniIsoParams_[2] / p4.pt()));
  	float miniDR = std::max( 0.05, std::min( 0.2, 10.0/p4.pt() ) );
  	//for (reco::PFCandidateCollection::const_iterator pf_it = pc->begin(); pf_it != pc->end(); pf_it++) {
  	int nPFCand = -1;
  	for(const auto &pfcan : *pfcands_ ){

		nPFCand++;
    	//if (int(pf_it - pc->begin()) == pc_idx) continue;  //don't count itself
		if( nPFCand == pc_idx ) continue;  //don't count itself
    	int id = std::abs(pfcan.pdgId());
		// pfIsolation_DZ = cms.double(0.1),
		// bool fromPV = (pf_it->fromPV() > 1 || fabs(pf_it->dz()) < pfIsolation_DZ_ );
    	// bool fromPV = (pfcan.fromPV() > 1 || fabs(pfcan.dz()) < 0.1 );
    	// fromPV only checked in cases where pfcand is a charged pion   
		// if pfcan is an electron or a muon then fromPV() will return 2 or 3 in all cases
		// if not -> if pfcan from PV then always true, else if pfcan passes close to PV ( dz < 0.1 ) then returns true 
        auto track = pfcan.bestTrack();
		//float sip2d(0), sip3d(0);
        float ip2d(0), ip3d(0);
        if( track != nullptr ){
			auto ttrack = ttBuilder.build( *track );
			GlobalVector direction( track->px(), track->py(), track->pz() );
			std::pair<bool, Measurement1D> ip2Dtool = IPTools::signedTransverseImpactParameter(ttrack, direction, primevtx );
			ip2d = std::abs( ip2Dtool.second.value() );
			//sip2d = std::abs( ip2Dtool.second.value() / ip2Dtool.second.error() );
            std::pair<bool, Measurement1D> ip3Dtool = IPTools::signedImpactParameter3D(ttrack, direction, primevtx );
			ip3d = std::abs( ip3Dtool.second.value() );
            //sip3d = std::abs( ip3Dtool.second.value() / ip3Dtool.second.error() );
		}//<<>>if( pftrackPtr != nullptr )
		//bool fromPV = sip2d < 3; 
        bool fromPV = ip2d < 0.05 && ip3d < 0.1;
    	float pt = pfcan.p4().pt();
    	float dr = deltaR(p4, pfcan);
    	//if (dr < pfIsolation_DR_) {
    	if (dr < 0.3 ) {
      		if (id == 211 && fromPV) chiso += pt; // charged cands from PV get added to trackIso
      		else if (id == 211) puiso += pt; // charged cands not from PV get added to pileup iso
      		if (id == 130) nhiso += pt; // neutral hadron iso
      		if (id == 22) phiso += pt; // photon iso
    	}//<<>>if (dr < pfIsolation_DR_) 
    	// same for mini isolation
    	if (dr < miniDR) {
      		if (id == 211 && fromPV) chmiso += pt;
      		else if (id == 211) pumiso += pt;
      		if (id == 130) nhmiso += pt;
      		if (id == 22) phmiso += pt;
    	}//<<>>if (dr < miniDR) 

  	}//<<>>for (reco::PFCandidateCollection::const_iterator pf_it = pc->begin(); pf_it != pc->end(); pf_it++)
	
  	iso = pat::PFIsolation(chiso, nhiso, phiso, puiso);
  	miniiso = pat::PFIsolation(chmiso, nhmiso, phmiso, pumiso);

}//<<>>void KUCMSTrackObject::getTrkIso(const PolarLorentzVector& .....

#endif
