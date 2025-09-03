// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Gen Object
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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchingTools.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenLeptonInfo.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSDisplacedVertex.hh"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSGenObjectHeader
#define KUCMSGenObjectHeader

//#define GenDEBUG true
#define GenDEBUG false

using namespace edm; 

class KUCMSGenObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSGenObject( const edm::ParameterSet& iConfig );
    ~KUCMSGenObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadGenEvtInfoTokens( edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken ){ genEvtInfoToken_ = genEvtInfoToken;};
    void LoadGenT0Tokens( edm::EDGetTokenT<float> gent0Token ){ gent0Token_ = gent0Token; };
    void LoadGenXYZ0Tokens( edm::EDGetTokenT<Point3D> genxyz0Token ){ genxyz0Token_ = genxyz0Token; };
    void LoadGenPUInfosTokens( edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfosToken ){ pileupInfosToken_ = pileupInfosToken; };
    void LoadgTokens( edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken ){ genParticlesToken_ = genPartToken; };
    void LoadGenJetsTokens( edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken ){ genJetsToken_ = genJetsToken; }; 
 
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 

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

    void keepMothers( const reco::Candidate* kid );
    void keepMothers( const reco::GenParticle kid );

    string bigKidChase( std::vector<reco::CandidatePtr> kids, float vx );
    //void motherChase( const reco::Candidate* kid, string & depth );

    void motherChase( const reco::Candidate* kid, string depth );
    void motherChase( const reco::GenParticle kid, string  depth );
    int typeChase( const reco::Candidate* kid );
    int typeChase( const reco::GenParticle kid );
	int getEvtType( int genPdgId );

    int llpGenChase( const reco::Candidate* kid );
    int llpGenChase( const reco::GenParticle kid );
    int llpGenChaseP( const reco::Candidate* kid, int depth );
    //int llpGenChaseJ( const reco::Candidate* kid, int depth );
    int llpGenChaseP( const reco::GenParticle & kid, int depth );
    int llpGenChaseP( const reco::GenJet & kid, int depth );
    const reco::Candidate* llpGenSearch( const reco::Candidate* kid, int target );
    void kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz );
    std::vector<float> kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  );
    // old single reco part gen matching
    std::vector<float> getGenPartMatch( const reco::SuperCluster* scptr, float pt );
  	//std::map<std::string, float> getGenPartMatch( const reco::SuperCluster &scptr, float pt ) const; 
  	//std::map<std::string, bool> MotherID(const int genIndex) const;
    std::vector<float> getGenJetInfo( float jetEta, float jetPhi, float jetPt );
    GlobalPoint GenVertex() {return GlobalPoint(genxyz0_->x(), genxyz0_->y(), genxyz0_->z()); };

    // new exclusive gen matching  v3fPoint
    std::vector<int> getGenMatch( const std::vector<v3fPoint> sc, std::vector<float> energies );
    std::vector<int> getGenPhoMatch( const std::vector<v3fPoint> sc, std::vector<float> energies );
    std::vector<int> getGenEleMatch( const std::vector<v3fPoint> sc, std::vector<float> energies, std::vector<reco::TransientTrack> trakcs );
    std::vector<int> getGenMatch( const std::vector<v3fPoint> sc, std::vector<float> energies, int select, std::vector<reco::TransientTrack> trakcs );
	std::vector<int> getGenMuonMatch( const std::vector<v3fPoint> sc, std::vector<float> reco_e, std::vector<reco::TransientTrack> trakcs );
    int getGenSigPhoXMother( uInt genIndex );
    int getGenSigEleXMother( uInt genIndex );
	int getGenSigPhoXMother( uInt genIndex, int loopcnt );

	// basic gen matching
	std::pair<int,std::pair<float,float>> getGenLeptonMatch( const reco::RecoCandidate & parton, int type  );

    // Gen electrons
  /*
  void GenElectronContent() const;
  std::vector<reco::GenParticle> GetSignalGenElectrons() const {return genSignalElectrons_;}
  std::vector<reco::GenParticle> GetGenElectrons() const {return genElectrons_;}
  std::vector<reco::GenParticle> GetGenParticles() const {return fgenparts;}
  LepMomType ClassifyGenElectron(const reco::GenParticle &genElectron) const;
  std::vector<int> MomIDs(const reco::GenParticle &genElectron) const; 
  //template <typename T>
  //GenClassifiedElectrons<T> GetGenClassifiedElectrons(const std::vector<T> &candidateTracks) const;
  //std::map<std::string, float> GenMatchElectrons(const Lorentz4Vec &electronTrack, const int charge) const;
  //std::map<std::string, float> getGenPartMatch( const reco::SuperCluster &scptr, float pt ) const;
  //std::map<std::string, bool> MotherID(const int genIndex) const;
  bool IsMotherZ(const reco::GenParticle &genElectron) const;
  */

    private:

    std::vector<reco::GenParticle> fgenparts;
    std::vector<int> fgpLlp;
    std::vector<int> fgpMomIdx;
    std::vector<int> fgpGMomIdx;

    std::vector<reco::GenParticle> tgenparts;
    std::vector<int> tgpLlp;
    std::vector<bool> tgpGetMom;
    std::vector<bool> tgpKeep;

    std::vector<reco::GenJet> fgenjets;
    std::vector<int> fgenjetllp;

    std::vector<int> nueEvntId;
    std::vector<int> chrEvntId;

    std::vector<float> matdr;
    std::vector<float> matde;

    //std::vector<reco::GenParticle> genElectrons_;
    //std::vector<reco::GenParticle> genSignalElectrons_;

    // genEvtInfo
    //const edm::InputTag genEvtInfoTag;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;
    edm::Handle<GenEventInfoProduct> genEvtInfo_;

    // gen time
    //const edm::InputTag gent0Tag;
    edm::EDGetTokenT<float> gent0Token_;
    edm::Handle<float> gent0_;

    // gen vertex
    //const edm::InputTag genxyz0Tag;
    edm::EDGetTokenT<Point3D> genxyz0Token_;
    edm::Handle<Point3D> genxyz0_;

    // pileups
    //const edm::InputTag pileupInfosTag;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfosToken_;
    edm::Handle<std::vector<PileupSummaryInfo>> pileupInfos_;

    // g ----------------------------------------------------------

    const edm::InputTag genParticlesTag;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;
    edm::Handle<std::vector<reco::GenParticle>> genParticles_;

    // GenJets ----------------------------------------------------------

    const edm::InputTag genJetsTag;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::Handle<std::vector<reco::GenJet>> genJets_;

  //std::vector<int> MomIDs(const reco::GenParticle &genElectron) const;
  //LepMomType AssignLeptonMomType(const int motherID) const;
  //LepType ClassifyGenElectron(const std::vector<int> &motherIDs) const;
  //LepType ClassifyGenElectron(const reco::GenParticle &genElectron) const;
  //	bool isSignalGenElectron(const reco::GenParticle &genElectron) const;
  //	void PrintMother(const LepMomType &momType) const;

	public:

	// gen particle info accessors
    int getGenLlpId( uInt genIndex ) const { return ( genIndex < fgpLlp.size() ) ? fgpLlp[genIndex] : -10; };
    int getGenMomIdx( uInt genIndex ) const { return ( genIndex < fgpMomIdx.size() ) ? fgpMomIdx[genIndex] : -10; };
    int getGenPdgId( uInt genIndex ) const { return ( genIndex < fgpLlp.size() ) ? fgenparts[genIndex].pdgId() : -10; };


    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

};//<<>>class KUCMSg : public KUCMSObjectBase

KUCMSGenObject::KUCMSGenObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : false );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfFlag.set( "motherChase", iConfig.existsAs<bool>("doGenMotherChase") ? iConfig.getParameter<bool>("doGenMotherChase") : false );
    cfPrm.set( "minGenE", iConfig.existsAs<double>("minGenE") ? iConfig.getParameter<double>("minGenE") : 0.0 );

}//<<>>KUCMSg::KUCMSGen( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSGenObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("genNtotal","Gen_nTotal",UINT);

    Branches.makeBranch("genPt","Gen_pt",VFLOAT);
    Branches.makeBranch("genEnergy","Gen_energy",VFLOAT);
    Branches.makeBranch("genPhi","Gen_phi",VFLOAT);
    Branches.makeBranch("genEta","Gen_eta",VFLOAT);
    Branches.makeBranch("genPx","Gen_px",VFLOAT);
    Branches.makeBranch("genPy","Gen_py",VFLOAT);
    Branches.makeBranch("genPz","Gen_pz",VFLOAT);
    Branches.makeBranch("genVx","Gen_vx",VFLOAT);
    Branches.makeBranch("genVy","Gen_vy",VFLOAT);
    Branches.makeBranch("genVz","Gen_vz",VFLOAT);
    Branches.makeBranch("genPdgId","Gen_pdgId",VUINT);
    Branches.makeBranch("genSusId","Gen_susId",VINT,"22=SignalPhoton");
    Branches.makeBranch("genStatus","Gen_status",VINT);
    Branches.makeBranch("genCharge","Gen_charge",VINT);
    Branches.makeBranch("genMass","Gen_mass",VFLOAT);
    Branches.makeBranch("genMomIdx","Gen_motherIdx",VINT);
    Branches.makeBranch("genMomDisplacment","Gen_momDisplacment",VFLOAT);
    Branches.makeBranch("genMomPdgId","Gen_MomPdgId",VINT);
    Branches.makeBranch("genGMomPdgId","Gen_GMomPdgId",VINT);
    Branches.makeBranch("genLWZX","Gen_isLWZX",VBOOL);
    Branches.makeBranch("genLWXQ","Gen_isLWZQ",VBOOL);
	Branches.makeBranch("genMatchDr","Gen_matchDr",VFLOAT);
    Branches.makeBranch("genMatchDe","Gen_matchDe",VFLOAT);

    Branches.makeBranch("genWgt","Evt_genWgt",FLOAT);
    Branches.makeBranch("genSusEvtType","Gen_susEvtType",VINT,"1=squark, 2=gluino, 3=xino, 4=slepton" );
    Branches.makeBranch("hasLWZX","Evt_hasLWZX",BOOL);
    Branches.makeBranch("nLWZX","Evt_nLWZX",INT);
    Branches.makeBranch("hasLWZQ","Evt_hasLWZQ",BOOL);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSGen::InitObject( TTree* fOutTree )

void KUCMSGenObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    if( GenDEBUG ) std::cout << "Getting gen tokens" << std::endl;
	//std::cout << "Getting gen tokens" << std::endl;

    iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);
    iEvent.getByToken(gent0Token_, gent0_);
    iEvent.getByToken(genxyz0Token_, genxyz0_);
    iEvent.getByToken(pileupInfosToken_, pileupInfos_);
    iEvent.getByToken(genParticlesToken_, genParticles_);
    iEvent.getByToken(genJetsToken_, genJets_);

	fgenparts.clear();
    fgpLlp.clear();
    fgpMomIdx.clear();
    fgpGMomIdx.clear();

	matdr.clear();
    matde.clear();

    tgenparts.clear();
    tgpLlp.clear();
    tgpGetMom.clear();
    tgpKeep.clear();

    fgenjets.clear();
    fgenjetllp.clear();
    //genElectrons_.clear();
    //genSignalElectrons_.clear();

    if( GenDEBUG ) std::cout << "Collecting Gen Particles" << std::endl;
	if( cfFlag("motherChase") ){
        std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
		std::cout << "Collecting Gen Particles for New Event --------------------------------------------------------------" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
	}//<<>>if( cfFlag("motherChase") ){
	nueEvntId.clear();
 	chrEvntId.clear();
	int oIndx = 0;
    for( const auto &genPart : *genParticles_ ){

        int genPdgId = std::abs( genPart.pdgId() );
        auto genPartSUSID = llpGenChaseP(genPart,0);

		tgenparts.push_back(genPart);
		tgpLlp.push_back(genPartSUSID);
        tgpGetMom.push_back(false);
        tgpKeep.push_back(false);

        if( genPart.pt() < cfPrm("minGenE") ){ oIndx++; continue; }
        //if( cfFlag("motherChase") && genPart.status() < 40 ){
        bool motherMeson( false );
        if( genPart.numberOfMothers() > 0 ) motherMeson = ( genPart.mother(0)->pdgId() > 100 ) ? true : false;
        if( cfFlag("motherChase") && ( genPart.pdgId() == 22 || genPart.pdgId() == 11 )  && not motherMeson ){
        //if( cfFlag("motherChase") ){
            std::cout << "- New Gen Particles ------------------------------------------------------------------" << std::endl;
			motherChase( genPart, "" );
		}//<<>>if( cfFlag("motherChase") )
        //std::cout << "Select Gen Particles to keep ----------------------------------------------" << std::endl;
		bool goodPdgId( abs(genPart.pdgId()) < 100 || abs(genPart.pdgId()) == 1000022 || abs(genPart.pdgId()) == 1000023 );
        if( genPart.status() == 1 && goodPdgId ){
			tgpKeep[oIndx] = true;
            tgpGetMom[oIndx] = true; // <<<<<<<<  turn on
			//if( genPart.pdgId() == 22 && genPartSUSID == 22 )
			//	{ std::cout << "Stable Gen Photon : " << std::endl; motherChase( genPart, "" ); }
			//if( genPdgId < 100 && genPart.pt() > 10.0 && genPartSUSID < 40 )
			//	{ std::cout << "Stable PO: " << genPartSUSID << " : " << std::endl; motherChase( genPart, "" );}
			//motherChase( genPart, "" );
	    	// Gen electron collections
	    	//if(abs(genPart.pdgId()) == 11) {
	      	//	genElectrons_.push_back(genPart);
	      	//	//std::cout << "\nmatched as electron!" << std::endl;
	      	//	//for(const auto &id : MomIDs(genPart))
	      	//	//std::cout << "  id = " << id << std::endl;
	      	//	//if( isSignalGenElectron(genPart) ) genSignalElectrons_.push_back(genPart);
	    	//}//if(abs(genPart.pdgId()) == 11)
            if( GenDEBUG && genPartSUSID == 22 ) std::cout << "Found genParton signal photon" << std::endl;
        } else {
			bool isTree( genPart.status() < 40 );
            bool isQuark( genPdgId < 9 );
			//bool isGluon( genPdgId == 21 );
			bool isWZ( genPdgId == 23 || genPdgId == 24 );
			bool isSusy( genPdgId > 1000000 && genPdgId < 3000000 );
            bool isNtlino1( genPdgId == 1000022 );
            bool isNtlino2( genPdgId == 1000023 );
            bool isChrino1( genPdgId == 1000024 );
            //bool isGrvtino( genPdgId == 1000039 );
			if( isQuark && isTree ){
                if( genPartSUSID == 41 || genPartSUSID == 42 || genPartSUSID == 43 || genPartSUSID == 33 || genPartSUSID == 34 ){
            		tgpKeep[oIndx] = true;
            		//tgpGetMom[oIndx] = true; //<<<<<<<<<<<<< turn on
				}//<<>>if( genPartSUSID < 50 )
			}//<<>>if( ( isSQuark || isGlino ) )
            if( isWZ && isTree ){
                //std::cout << "Gen Virt Q/WZ Particle Status : " << genPart.status() << std::endl;
                if( genPartSUSID == 23 || genPartSUSID == 24 ){ 
                    tgpKeep[oIndx] = true;
                    tgpGetMom[oIndx] = true; //<<<<<<<<<<<<< turn on
                }//<<>>if( genPartSUSID < 50 )
            }//<<>>if( ( isSQuark || isGlino ) )
			if( isSusy && isTree ){
				tgpKeep[oIndx] = true;
				tgpGetMom[oIndx] = true;
			}//<<>>if( isSusy && isTree ( 
            if( isNtlino1 || isNtlino2 || isChrino1 ){
            	int kPdgId = genPart.pdgId();
            	int mPdgId = (genPart.mother(0))->pdgId();
				bool lastInChian = ( kPdgId == mPdgId ) ? 0 : 1;
				if( lastInChian ){ 
					( isNtlino1 || isNtlino2 ) ? nueEvntId.push_back( typeChase(genPart) ) : chrEvntId.push_back( typeChase(genPart) ); 
				}//<<>>if( lastInChian )  
                //std::cout << "Gen Susy Particle : " << std::endl; motherChase( genPart.mother(0), "" );
			}//<<>>if( ( isSQuark || isGlino ) ) 
		}//<<>> else if( genPart.status() == 1 )

		oIndx++;
    }//<<>>for(const auto& genPart : *genParticles_ )

    //int gpidx = 0;
	//std::cout << "Flagging Moms -------------------------------------------------------------" << std::endl;
    //for (const auto & gp : tgenparts ){ 
    int nTgenparts = tgenparts.size();
    for( int gpit = 0; gpit < nTgenparts; gpit++ ){

		//auto gp = fgenparts[gpit];
		if( tgpGetMom[gpit] ){
			//std::cout << "Saving Moms for llpid: " << tgpLlp[gpidx] << " pdgid : " << gp.pdgId() << std::endl; 
			//motherChase( gp, "" );
			//std::cout << " - Getting Moms : "<< std::endl; 
			keepMothers( tgenparts[gpit] );
		}//<<>>if( tgpGetMom[gpidx] ) 
		//gpidx++; 

	}//<<>>for (const auto & gp : tgenparts )

    //std::cout << "Saving Genparitcle to keep -------------------------------------------------------------" << std::endl;
	//gpidx = 0;
    //for (const auto & gp : tgenparts ){ 
    for( int gpit = 0; gpit < nTgenparts; gpit++ ){

		if( tgpKeep[gpit] ){
            //if( abs(tgenparts[gpit].pdgId()) == 1000022 ) std::cout << "Keeping 1000022" << std::endl;
			//if( abs(tgenparts[gpit].pdgId()) == 1000023 ) std::cout << "Keeping 1000023" << std::endl;
			//if( abs(tgenparts[gpit].pdgId()) > 100 && abs(tgenparts[gpit].pdgId()) < 1000000 ) continue;
            //if( abs(tgenparts[gpit].pdgId()) > 3000000 ) continue;
		    fgenparts.push_back( tgenparts[gpit] );
    		fgpLlp.push_back( tgpLlp[gpit] );
			fgpMomIdx.push_back( -1 );
            fgpGMomIdx.push_back( -1 );
			matdr.push_back( -10.f );
            matde.push_back( -10.f );
		}//<<>>if( tgpKeep[gpidx] )
		//gpidx++;

	}//<<>>for (const auto & gp : tgenparts )

    //std::cout << "Getting GenPart mom indexes -------------------------------------------------------------" << std::endl;
    //gpidx = 0;
	int nFgenparts = fgenparts.size();
    //for( const auto & gp : fgenparts ){
	for( int gpit = 0; gpit < nFgenparts; gpit++ ){

		auto gp = fgenparts[gpit];
		//int gppdg = abs(gp.pdgId());
		if( gp.numberOfMothers() < 1 ) continue;;

		auto gpm = gp.mother(0);
		bool done = false;
		while( not done ){
    		bool notLastInChain( ( gpm->numberOfMothers() == 1 ) && ( gp.pdgId() == gpm->pdgId() ) );
    		bool validStatus( gpm->status() < 40 );
    		if( notLastInChain || not validStatus ){ gpm = gpm->mother(0); }
			else done = true;
		}//<<>>while( not done )

        float mr = hypo(gpm->eta(),gpm->phi());
        int mpdg = gpm->pdgId();
		//int status = gp.status();
		//int mstatus = gpm->status();
		//int canidx = 0;
        //for( const auto &  cangp : fgenparts ){
        //std::cout << "Finding mom for : " << gppdg << " <- " << mpdg << " status : " << status << " <- " << mstatus << std::endl;
		//if( gppdg == 1000022 ) std::cout << "Finding mom for : " << gppdg << " type : " << mpdg << std::endl;
        //if( gppdg == 22 && mpdg == 1000023 ) std::cout << "Finding mom for : " << gppdg << " type : " << mpdg << std::endl;
		bool matched( false );
		for( int cangpit = 0; cangpit < nFgenparts; cangpit++ ){
			if( gpit == cangpit ) continue;
			auto cangp = fgenparts[cangpit];
            float canr = hypo(cangp.eta(),cangp.phi());
            int canpdg = cangp.pdgId();
			//if( canpdg == mpdg ) std::cout << " Matching: " << mpdg << " dr: " << std::abs( canr - mr ) << std::endl;
			//if( canpdg == 1000023 && mpdg == 1000023 ) std::cout << " Matching: 1000023 dr: " << std::abs( canr - mr ) << std::endl;
            if( canpdg == mpdg && std::abs( canr - mr ) == 0 ){
				//std::cout << " ------- Matched" << std::endl;
                fgpMomIdx[gpit] = cangpit;
				matched = true;
                break;
            }//<<>>if( canpt == mpt && canpdg == mpdg )
			//canidx++;
		}//<<>>for( auto cangp : fgenparts )
		if( fgpMomIdx[gpit] == -1 && matched ) std::cout << " gen part mom index bad !!!!" << std::endl;
		if( gp.pdgId() != fgenparts[gpit].pdgId() ) std::cout << " gen part index bad !!!!" << std::endl;
		//if( fgpMomIdx[gpidx] == -1 && gp.pdgId() == 22 ){ // && fgpLlp[gpidx] == 22 ){ 
		//	std::cout << " Pho Mom Not Found : " << gp.numberOfMothers() << std::endl;
		//}//<<>> no match search

		// last step !!!!!!!!!!!!
        //gpidx++;
	}//<<>>for (const auto & gp : tgenparts )

    for( int gpit = 0; gpit < nFgenparts; gpit++ ){

		auto mgpit = fgpMomIdx[gpit];
		if( mgpit < 0 ) continue;
        auto gp = fgenparts[mgpit];
        //int gppdg = abs(gp.pdgId());
        if( gp.numberOfMothers() < 1 ) continue;;

        auto gpm = gp.mother(0);
        bool done = false;
        while( not done ){
            bool notLastInChain( ( gpm->numberOfMothers() == 1 ) && ( gp.pdgId() == gpm->pdgId() ) );
            bool validStatus( gpm->status() < 40 );
            if( notLastInChain || not validStatus ){ gpm = gpm->mother(0); }
            else done = true;
        }//<<>>while( not done )

        float mr = hypo(gpm->eta(),gpm->phi());
        int mpdg = gpm->pdgId();
        bool matched( false );
        for( int cangpit = 0; cangpit < nFgenparts; cangpit++ ){
            if( mgpit == cangpit ) continue;
            auto cangp = fgenparts[cangpit];
            float canr = hypo(cangp.eta(),cangp.phi());
            int canpdg = cangp.pdgId();
            if( canpdg == mpdg && std::abs( canr - mr ) == 0 ){
                fgpGMomIdx[gpit] = cangpit;
                matched = true;
                break;
            }//<<>>if( canpt == mpt && canpdg == mpdg )
        }//<<>>for( auto cangp : fgenparts )
        if( fgpGMomIdx[gpit] == -1 && matched ) std::cout << " gen part gmom index bad !!!!" << std::endl;
        if( gp.pdgId() != fgenparts[mgpit].pdgId() ) std::cout << " gen part index bad !!!!" << std::endl;

    }//<<>>for (const auto & gp : tgenparts )

/*
    gpidx = 0;
    for (const auto & gp : fgenparts ){

		if( fgpMomIdx[gpidx] == -1 && ( ( gp.pdgId() < 100 && gp.pdgId() > 10 ) || ( gp.pdgId() > 1000000 && gp.pdgId() < 3000000  ) ) ){
		std::cout << " Bad Match Investigation: pdgid : " << gp.pdgId() << std::endl;
        float mr = hypo((gp.mother(0))->eta(),(gp.mother(0))->phi());
        int mpdg = (gp.mother(0))->pdgId();
        int canidx = 0;
        for( auto cangp : fgenparts ){
            float canr = hypo(cangp.eta(),cangp.phi());
            int canpdg = cangp.pdgId();
			//std::cout << " No Mom Match -pdg: " << mpdg << " = " << canpdg << " && " << mpt << " = "  << canpt << std::endl;
            if( canpdg == mpdg ){
				std::cout << " No Mom Match -pdg: " << mpdg << " = " << canpdg << " && " << mr << " = "  << canr << std::endl;
                //fgpMomIdx[gpidx] = canidx;
                //break;
            }//<<>>if( canpt == mpt && canpdg == mpdg )
            canidx++;
        }//<<>>for( auto cangp : fgenparts )
        //if( gp.pdgId() != fgenparts[gpidx].pdgId() ) std::cout << " gen part index bad !!!!1 : " << std::endl;
        //if( fgpMomIdx[gpidx] == -1 && gp.pdgId() == 22 ){ // && fgpLlp[gpidx] == 22 ){ 
        //    std::cout << " Pho Mom Not Found : " << gp.numberOfMothers() << std::endl;
        //}//<<>> no match search
		}//<<>> only look at unmatched moms
	
        // last step !!!!!!!!!!!!
        gpidx++;
    }//<<>>for (const auto & gp : tgenparts )
*/

	//std::cout << " -- # Nuetral - # Charge : " << nNuet << " - " << nChar << std::endl;

    if( GenDEBUG ) std::cout << "Collecting GenJets" << std::endl;
    for(const auto &genjet : *genJets_ ){ 
		fgenjets.push_back(genjet);
		auto genJetSUSID = llpGenChaseP(genjet,0); 
		fgenjetllp.push_back(genJetSUSID); 
	}//<<>>for(const auto &genjet : *genJets_ )
	
}//<<>>void KUCMSGen::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSGenObject::PostProcessEvent( ItemManager<float>& geVar ){

	int ngen = matdr.size();
	if( ngen != int(fgenparts.size()) ) std::cout << " !!!! -- MDR v fGen mismatch !!!! " << std::endl;
	for( int idx = 0; idx < ngen; idx++ ){
		if( matde[idx] < 1000 ){
    		Branches.fillBranch("genMatchDr",matdr[idx]);
    		Branches.fillBranch("genMatchDe",matde[idx]);
		} else {
            Branches.fillBranch("genMatchDr",-1.f);
            Branches.fillBranch("genMatchDe",-10.f);
		}//<<>>if( matde[idx] < 1000 )
	}//<<>>for( int idx = 0; idx < ngen; idx++ )

}//<<>>void KUCMSGenObject::PostProcessEvent( ItemManager<float>& geVar )

void KUCMSGenObject::ProcessEvent( ItemManager<float>& geVar ){

    //if( GenDEBUG ) std::cout << "Processing Gens" << std::endl;

    Branches.clearBranches();

    //if( GenDEBUG ) std::cout << " - enetering Gen loop" << std::endl;
	//std::cout << " - enetering Gen loop" << std::endl;	

	bool hasLWZX = false;
    bool hasLWZQ = false;
	int nLWZX = 0;
	int nGenParts = 0;
    for (const auto & genpart : fgenparts ){

        const float genPt = genpart.pt();
        const float genEnergy = genpart.energy();
        const float genPhi = genpart.phi();
        const float genEta = genpart.eta();
        const uInt genPdgId = std::abs(genpart.pdgId());
        const float genPx = genpart.px();
        const float genPy = genpart.py();
        const float genPz = genpart.pz();
		//if( GenDEBUG ) std::cout << "GenPart : 1/2 WAY " << std::endl;
        const int genStatus = genpart.status();
		//const int genSusId = llpGenChaseP( genpart, 0 );
        const int genSusId = fgpLlp[nGenParts];
		const int genMomIdx = fgpMomIdx[nGenParts];
        const int genGMomIdx = fgpGMomIdx[nGenParts];
        const int genCharge = genpart.charge();
        const float genVx = genpart.vx();
        const float genVy = genpart.vy();
        const float genVz = genpart.vz();
        const float genMass = genpart.mass();		

		const float momVx = ( genMomIdx > -1 ) ? fgenparts[genMomIdx].vx() : -999;
        const float momVy = ( genMomIdx > -1 ) ? fgenparts[genMomIdx].vy() : -999;
        const float momVz = ( genMomIdx > -1 ) ? fgenparts[genMomIdx].vz() : -999;
		const float displacment = ( genMomIdx > -1 ) ? hypo( genVx-momVx, genVy-momVy, genVz-momVz ) : -10;		
		
		const int genMomPdgId = ( genMomIdx > -1 ) ? std::abs(fgenparts[genMomIdx].pdgId()) : -9;
        const int genGMomPdgId = ( genGMomIdx > -1 ) ? std::abs(fgenparts[genGMomIdx].pdgId()) : -9;

		bool isLWZX = ( genPdgId == 11 || genPdgId == 13 ) && ( genMomPdgId == 23 || genMomPdgId == 24 ) && genGMomPdgId == 1000023;
		bool isLWZQ = ( genPdgId == 11 || genPdgId == 13 ) && ( genMomPdgId == 23 || genMomPdgId == 24 ) && genGMomPdgId < 10;
		if( isLWZX ){ hasLWZX = true; nLWZX++; }
        if( isLWZQ ) hasLWZQ = true;

		//if( GenDEBUG ) std::cout << "GenPart : genSusId = " << genSusId << std::endl;
		Branches.fillBranch("genNtotal", unsigned(fgenparts.size()) );
        Branches.fillBranch("genPt",genPt);
        Branches.fillBranch("genEnergy",genEnergy);
        Branches.fillBranch("genPhi",genPhi);
        Branches.fillBranch("genEta",genEta);
        Branches.fillBranch("genPdgId",genPdgId);
        Branches.fillBranch("genPx",genPx);
        Branches.fillBranch("genPy",genPy);
        Branches.fillBranch("genPz",genPz);
        Branches.fillBranch("genVx",genVx);
        Branches.fillBranch("genVy",genVy);
        Branches.fillBranch("genVz",genVz);
        Branches.fillBranch("genMomIdx",genMomIdx);
        Branches.fillBranch("genSusId",genSusId);
        Branches.fillBranch("genStatus",genStatus);
        Branches.fillBranch("genCharge",genCharge);
        Branches.fillBranch("genMass",genMass);
        Branches.fillBranch("genMomDisplacment",displacment);
        Branches.fillBranch("genMomPdgId",genMomPdgId);
        Branches.fillBranch("genGMomPdgId",genGMomPdgId);
        Branches.fillBranch("genLWZX",isLWZX);
        Branches.fillBranch("genLWXQ",isLWZQ);

		nGenParts++;
    }//<<>> for (const auto genpart : fgenparts )

	//"1=squark, 2=gluino, 3=Xino" 
	int nNue = nueEvntId.size();
	int nChr = chrEvntId.size();
	int first = 8;
	int second = 8;
	if( nNue == 2 ){
		first  = getEvtType( nueEvntId[0] );
		second = getEvtType( nueEvntId[1] );
	} else if( nNue == 0 ){
		if( nChr > 0 ) first  = getEvtType( chrEvntId[0] );
		else { first = 5; second = 5; }
		if( nChr > 1 ) second = getEvtType( chrEvntId[1] );
		else { second = 5; }
	} else if( nNue == 1 ){
		first  = getEvtType( nueEvntId[0] );
		if( nChr > 0 ){
			for( auto cType : chrEvntId ){ if( cType != nueEvntId[0] ){ second = getEvtType( cType ); break; } }
			if( second == 8 ) { second = 7; }
		} else { second = 6; } 
	} else if( nNue > 2 ){ first = 6; second = 6; } 
	
	int genEvtType = ( first < second ) ? 10*first+second : 10*second+first;
    Branches.fillBranch("genSusEvtType",genEvtType);
    //if( GenDEBUG ) std::cout << "GenPart : LOADING GEN WT " << std::endl;
    float wgt = genEvtInfo_->weight();
    Branches.fillBranch("genWgt",wgt);
    Branches.fillBranch("hasLWZX",hasLWZX);
    Branches.fillBranch("nLWZX",nLWZX);
    Branches.fillBranch("hasLWZQ",hasLWZQ);
	geVar.fill("genWgt",wgt);
    //if( GenDEBUG ) std::cout << "GenPart : Done " << std::endl;

}//<<>>void KUCMSGen::ProcessEvent()

void KUCMSGenObject::EndJobs(){}

std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<v3fPoint> sc, std::vector<float> reco_e ){

	// gen matching to supcer clusters froom reeco phtons to gen photons
	std::vector<reco::TransientTrack> null;
	return getGenMatch( sc, reco_e, 22, null );

}//<<>>std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenEleMatch( const std::vector<v3fPoint> sc, std::vector<float> reco_e, std::vector<reco::TransientTrack> trakcs ){

    // gen matching to supcer clusters froom reco electrons to gen electrons
    return getGenMatch( sc, reco_e, 11, trakcs );

}//<<>>std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenMuonMatch( const std::vector<v3fPoint> sc, std::vector<float> reco_e, std::vector<reco::TransientTrack> trakcs ){

    // gen matching to supcer clusters froom reco electrons to gen electrons
    return getGenMatch( sc, reco_e, 13, trakcs );

}//<<>>std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenMatch( const std::vector<v3fPoint> sc, std::vector<float> reco_e ){

	// matchoing any particle :: matching valid for nutral particles only
	std::vector<reco::TransientTrack> null;
	return getGenMatch( sc, reco_e, 0, null );

}//<<>>std::vector<int> KUCMSGenObject::getGenMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenMatch( const std::vector<v3fPoint> sc, std::vector<float> en, int select, std::vector<reco::TransientTrack> tracks ){

    // GenParticle Info   -------------------------------------------------------------------
    if( GenDEBUG ) std::cout << "Getting phoGenParton Match Information" << std::endl;
	// code caclutes eta & phi of the momentum vector requird for a nuetral particle 
	// produced at a gen particles production vertix to cause the given supercluster and compares
	// this to the actual gen particles mometum vectors eta & phi to find best, exclusive, match
	// between a reco particle with that supercluster and a gen particle

    std::vector<int> genPartIndx;
    std::vector<float> genMDR;
    std::vector<float> genMDE;
	std::vector<int> genRecoIndx;
	float drthres = 0.1;
	float demax = 1.0;
	int nGenPart = fgenparts.size();
    for( int idx = 0; idx < nGenPart; idx++ ){
		if( fgenparts[idx].status() != 1 ) continue;
		if( select != 0 && std::abs(fgenparts[idx].pdgId()) != select ) continue;
        genPartIndx.push_back(idx);
        genMDR.push_back(drthres);
		genMDE.push_back(1000);
        genRecoIndx.push_back(-9);
    }//<<>>for( int idx = 0; idx < nGenPart; idx++ )

    int nSCPtrs = sc.size();
	int nGPIs = genPartIndx.size();
    std::vector<int> recoPartIndx;
    for( int it = 0; it < nSCPtrs; it++ ){ recoPartIndx.push_back(-1); }
    bool done = false;
    int advance = 0;
    while( not done ){

		if( GenDEBUG ) std::cout << " - Finding  Matches " << std::endl;
        for( int it = advance; it < nSCPtrs; it++ ){
            for( int pgidx = 0; pgidx < nGPIs; pgidx++ ){
				//if( GenDEBUG ) std::cout << " -- Getting sc info " << std::endl;
				auto partidx = genPartIndx[pgidx];
        		const auto rhX = sc[it].x();
        		const auto rhY = sc[it].y();
        		const auto rhZ = sc[it].z();
                //const auto rhe = reco_e[it];
                //if( GenDEBUG ) std::cout << " -- Getting genpart info " << std::endl;
                const auto gnEta = fgenparts[partidx].eta();
                const auto gnPhi = fgenparts[partidx].phi();
        		const auto gnX = fgenparts[partidx].vx();
        		const auto gnY = fgenparts[partidx].vy();
        		const auto gnZ = fgenparts[partidx].vz();
        		//const auto ge = fgenparts[partidx].energy();
                //if( GenDEBUG ) std::cout << " -- calc dr match " << std::endl;
                //float cmEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
				//float cmPhi = std::atan2(rhY-gnY,rhX-gnX);
				double cmEta = rhX;
				double cmPhi = rhY;
				if( select == 22 ){
        			cmEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
        			cmPhi = std::atan2(rhY-gnY,rhX-gnX);
				}//<<>>if( select != 13 )
				if( select == 11 || select == 13 ){
        			SteppingHelixPropagator propagator( tracks[it].field(), anyDirection );
        			const GlobalPoint genLocation( gnX, gnY, gnZ );
        			//const TrajectoryStateClosestToPoint propLocation = tracks[it].trajectoryStateClosestToPoint( genLocation );
                    const FreeTrajectoryState propLocation = propagator.propagateWithPath(tracks[it].initialFreeState(), genLocation).first;
        			//const GlobalPoint pca( propLocation.referencePoint() );// with TrajectoryStateClosestToPoint
					const GlobalPoint pca( propLocation.position() );// with FreeTrajectoryState
					const reco::Track track(tracks[it].track());
        			const double tmin = (pca.z() - track.vz())*track.qoverp()/(87.78*sin(track.lambda()));
					//std::cout << " -- GenMatch: cmPhi: " <<  cmPhi << " - " << tmin;
					cmPhi -= tmin;
					//std::cout << " = " << cmPhi << std::endl;
				}//<<>>else if( select == 11 || select == 13 )
        		float dr = std::sqrt(reco::deltaR2( gnEta, gnPhi, cmEta, cmPhi ));
				float de = std::abs( en[it] - fgenparts[partidx].energy() )/fgenparts[partidx].energy();
				bool passDr = dr < genMDR[pgidx];
				bool passDe = de < demax;
                if( passDr && passDe ){ genMDR[pgidx] = dr; genMDE[pgidx] = de; genRecoIndx[pgidx] = it; }
            }//<<>>for( int pgidx = 0; pgidx < genSigPhoIdx.size(); pgidx++ )
        }//<<>>for( int it = 0; it < nSCPtrs; it++ )
        bool unassigned = false;
		if( GenDEBUG ) std::cout << " - Checking Dup Matches " << std::endl;
        for( long unsigned int pgidx = 0; pgidx < genMDR.size(); pgidx++ ){
            for( long unsigned int pgidx2 = pgidx+1; pgidx2 < genMDR.size(); pgidx2++ ){
                if( ( genRecoIndx[pgidx] != -9 ) && ( genRecoIndx[pgidx] == genRecoIndx[pgidx2] ) ){
                    unassigned = true;
                    if( genMDR[pgidx2] >= genMDR[pgidx] ){ genRecoIndx[pgidx2] = -9; genMDR[pgidx2] = drthres; genMDE[pgidx2] = 1000; }
                    else { genRecoIndx[pgidx] = -9; genMDR[pgidx] = drthres; genMDE[pgidx] = 1000; }
                }//<<<>>if( genPhoIndx[pgidx] == genPhoIndx[pgidx2] )
            }//<<>>for( int pgidx2 = pgidx; pgidx2 < genGenPhoMDR.size(); pgidx2++ )
        }//<<>>for( int pgidx = 0; pgidx < genGenPhoMDR.size(); pgidx++ )
		if( GenDEBUG ) std::cout << " - Checking for advancment and done " << std::endl;
        if( unassigned ) advance++; else done = true;
		if( advance == nSCPtrs ) done = true;

    }//<<>>while()
    if( GenDEBUG ) std::cout << " - Creating Reco matched gen part index list " << std::endl;
    for( long unsigned int iter = 0; iter < genRecoIndx.size(); iter++ ){
        if( genRecoIndx[iter] != -9 ) recoPartIndx[genRecoIndx[iter]] = genPartIndx[iter];
		matdr[genPartIndx[iter]] = genMDR[iter]; matde[genPartIndx[iter]] = genMDE[iter];
    }//<<>>for( int iter = 0; iter < genRecoIndx.size(); iter++ )
    if( GenDEBUG ){ std::cout << " -- dRLIST: #" << nGPIs << " "; for( auto it : genMDR ) std::cout << it << " "; std::cout << std::endl; }
	if( GenDEBUG ){ std::cout << " -- IdLIST: #" << nSCPtrs << " "; for( auto it : recoPartIndx ) std::cout << it << " "; std::cout << std::endl; }
	if( GenDEBUG ) std::cout << " - Done with Gen Matching" << std::endl;
    return recoPartIndx;

}//<<>>getGenPhoMatch( reco::SuperClusterCollection *scptr, std::vector<reco::GenParticle> fgenparts )

std::vector<float> KUCMSGenObject::getGenPartMatch( const reco::SuperCluster* scptr, float energy ){

    // GenParticle Info   -------------------------------------------------------------------
    // best match by dr for given SC to a gen particle ( exclusive not guaranteed )
    if( GenDEBUG ) std::cout << "Getting phoGenParton Information Single" << std::endl;
    // set defaults for no match here
    vector<float> results;
    float minDr(0.3);
    float minSDr(0.4);
    float minRe(-1.0);
    float minSRe(-1.0);
    float llpId(-1.0);
    float llpSId(-1.0);
    int matchedIdx(-10);// 1
    int matchedSIdx(-10);// 1
    int index(0);
    for(const auto & genPart : fgenparts ){

		if( genPart.status() != 1 ) continue;

        const auto rhX = scptr->x();
        const auto rhY = scptr->y();
        const auto rhZ = scptr->z();
        const auto gnX = genPart.vx();
        const auto gnY = genPart.vy();
        const auto gnZ = genPart.vz();
        const auto ge = genPart.energy();
		//if( GenDEBUG ) std::cout << " GenPart: " << gnX << " " << gnY << " " << gnZ << std::endl;
		//if( ge < 1.0 ) continue;
        auto cphoEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
        auto cphoPhi = std::atan2(rhY-gnY,rhX-gnX);
        auto dr = std::sqrt(reco::deltaR2(genPart.eta(), genPart.phi(), cphoEta, cphoPhi ));
        auto re = energy/ge;
        if( dr < minDr ){
            minDr = dr;
            minRe = re;
            matchedIdx = index;
            llpId = fgpLlp[index];
        } else if( dr < minSDr ){
            minSDr = dr;
            minSRe = re;
            matchedSIdx = index;
            llpSId = fgpLlp[index];
        }//<<>>if( dr < minDr && dp < minDp )
        index++;

    }//<<>>for(const auto& genPart : fgenparts  )

    float matDr(-1.0), matRe(-1.0), matSDr(-1.0), matSRe(-1.0);
    if( matchedIdx >=0 ){ matDr = minDr; matRe = minRe; matSDr = minSDr; matSRe = minSRe; }

    results.push_back( matchedIdx );
    results.push_back( matDr );
    results.push_back( matRe );
    results.push_back( matchedSIdx );
    results.push_back( matSDr );
    results.push_back( matSRe );
    results.push_back( llpId );
    results.push_back( llpSId );
    return results;

}//<<>>getGenPartMatch( reco::SuperClusterCollection *scptr, std::vector<reco::GenParticle> fgenparts )

std::pair<int,std::pair<float,float>> KUCMSGenObject::getGenLeptonMatch( const reco::RecoCandidate & parton, int type ){

	// basic dr match for best match to a single gen particle   - thrshold 0.8

	float eta = parton.eta();
    float phi = parton.phi();
	float en = parton.energy();
	int parPdgId = parton.pdgId();
	//std::cout << " -- type/pdgid: " << type << " ?= " << parPdgId << std::endl;
    float minRe(100.0);
    float minDr(0.8);
	int matchedIdx(-9);
	//std::cout << " -- Getting fgenparts size " << std::endl;

    int nGenPart = fgenparts.size();
    //std::cout << " - Starting fgenparts loop for : " << nGenPart << std::endl;
    for( int idx = 0; idx < nGenPart; idx++ ){

        //std::cout << " -- in loop with :: " << std::endl;
        if( fgenparts[idx].status() != 1 ) continue;

		int pdgId = std::abs( fgenparts[idx].pdgId() );
		if( pdgId != parPdgId ) continue;
		//std::cout << " -- pdgID: " << pdgId << std::endl;		

		int genMomIdx = fgpMomIdx[idx];
		int genMomPdgId = ( genMomIdx > -1 ) ? std::abs(fgenparts[genMomIdx].pdgId()) : -9;
		if( genMomPdgId != 23 && genMomPdgId != 24 && genMomPdgId != 1000023  ) continue;
		//std::cout << " -- momPdgID: " << genMomPdgId << std::endl;

		float geta = fgenparts[idx].eta();
		float gphi = fgenparts[idx].phi();
		float gen = fgenparts[idx].energy();
        float dr = std::sqrt(reco::deltaR2(geta,gphi,eta,phi));
        float re = std::abs( (en/gen) - 1 );
		bool dRmatch = ( dr < minDr );
		//bool rEmatch = ( re <= ( minRe + 0.25 ) );	
        //if( dRmatch && rEmatch ){
        if( dRmatch ){
            minDr = dr;
            minRe = re;
            matchedIdx = idx;
        }//<<>>if( dr < minDr && dp < minDp )
	
	}//<<>>for(const auto & genPart : fgenparts ){

	std::pair< float, float > matchStats = std::make_pair( minDr, minRe );
    std::pair<int,std::pair<float,float>> result = std::make_pair( matchedIdx, matchStats );
	return result;

}//<<>> int KUCMSGenObject::getGenPartMatch( const T &lepton )

int KUCMSGenObject::getGenSigPhoXMother( uInt genIndex ){

    if( GenDEBUG ) std::cout << "Getting phoSigGenPhoXMother for " << genIndex << std::endl;
    int result( -1 );
    if( genIndex >= fgpLlp.size() ){ std::cout << " --- INDEX Blown !!!!  getGenSigPhoInfo " << std::endl; return result; }
	else { 

		if( genIndex >= fgpMomIdx.size() ) return -10;
		int momidx = fgpMomIdx[genIndex];
		if( momidx < 0 ) return -11;
        if( momidx >= int(fgenparts.size()) ) return -12;
		if( fgenparts[momidx].pdgId() > 1000000 ) result = momidx;
		else result = getGenSigPhoXMother( momidx, 0 );

	}//<<>>else
	//if( result < -10 ) std::cout << " - Found " << result << " of " << fgenparts[genIndex].pdgId() << std::endl;
	//if( result >= 0 ) std::cout << " - Found phoSigGenPhoXMother " << result << " of " << fgenparts[result].pdgId() << std::endl;
	//else std::cout << " - Found phoSigGenPhoXMother result of " << result << " on Index Bounds Error !" << std::endl;
	return result;

}//<<>>int KUCMSGenObject::getGenSigPhoInfo( uInt genIndex )

int KUCMSGenObject::getGenSigPhoXMother( uInt genIndex, int loopcnt ){

    if( GenDEBUG ) std::cout << "Looping phoSigGenPhoXMother with " << genIndex << " on loop " << loopcnt << std::endl;
    int result( -1 );
    if( genIndex >= fgpLlp.size() ){ std::cout << " --- INDEX Blown !!!! phoSigGenPhoXMother " << std::endl; return result; }
    else {

        if( genIndex >= fgpMomIdx.size() ) return -10;
        int momidx = fgpMomIdx[genIndex];
        if( momidx < 0 ) return -11;
        if( momidx >= int(fgenparts.size()) ) return -12;
        if( fgenparts[momidx].pdgId() > 1000000 ) result = momidx;
        else {  
			if( momidx == int(genIndex) ){
				std::cout << " -- CJ!!! On pgdid: " << fgenparts[momidx].pdgId() << " for " << momidx << std::endl;
				return -99;
			}//<<>>if( momidx == int(genIndex) )
			loopcnt++; 
			result = getGenSigPhoXMother( momidx, loopcnt );
		}//<<>> else 

    }//<<>>else
    //if( result < -10 ) std::cout << " - Found " << result << " of " << fgenparts[genIndex].pdgId() << std::endl;
    //if( result >= 0 ) std::cout << " - Found phoSigGenPhoXMother of " << result << " for pdgId " << fgenparts[result].pdgId() << std::endl;
    //else std::cout << " - Found phoSigGenPhoXMother result of " << result << " on Index Bounds Error " << std::endl;
    return result;

}//<<>>int KUCMSGenObject::getGenSigPhoInfo( uInt genIndex )

int KUCMSGenObject::getGenSigEleXMother( uInt genIndex ){

	return getGenSigPhoXMother( genIndex );

}//<<>>std::vector<float> KUCMSGenObject::getGenSigEleInfo( uInt genIndex )


std::vector<float> KUCMSGenObject::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  ){
// redo this function to give tof and impact angle for input gen particle

    vector<float> result;
    vector<float> kidtime;
    vector<float> kide;
    vector<float> impact;
    bool llpcut(false);
    float nextBX(0.0);
    float llpkid(0.0);
    for( auto kid : kids ){

        bool done(false);
        bool first(true);
        bool llp(false);
        bool bquark(false);
        auto mom = kid->mother();
        vector<double> stepp{kid->p()};
        vector<double> stepe{kid->energy()};
        vector<double> orignvx;
        vector<double> orignvy;
        vector<double> orignvz;
        bool top(true);
        int steps(0);
        while( not done ){
            auto pdgid = llpGenChaseP( kid.get(), 0 );
            if( pdgid < 25 && pdgid > 21 ) llp = true;
            if( pdgid < 22 ) bquark = true;
            if( mom->pt() == 0 ) done = true;
            else {
                auto gmom = mom->mother();
                if( top ){
                    if( first ) first = false;
                    else {
                        stepp.push_back(mom->p());
                        stepe.push_back(mom->energy());
                    }//<<>>if( first )
                    steps++;
                    //if( llp || bquark ) llpstep.push_back(1.0); else llpstep.push_back(0.0);
                    orignvx.push_back(mom->vx());
                    orignvy.push_back(mom->vy());
                    orignvz.push_back(mom->vz());
                }//<<>>if( top )
                if( gmom->pt() == 0 ) done = true;
                else {
                    if( mom->vx() == gmom->vx() ) top = false;
                    else top = true;
                    mom = gmom;
                }//<<>>if( gmom->pt() == 0 )
            }//<<>>if( mom->pt() == 0 )                 
        }//<<>>while( not done )
        first = true;

        if( GenDEBUG && false ) std::cout << " ---- jetGenTime Calc Steps : " << steps << std::endl;
        bool stepcut( steps == 0 );
        if( GenDEBUG && false ) std::cout << " ---- jetGenTime Calc llpcut : " << llpcut << std::endl;
        if( not stepcut ){
            float maxe(0.0);
            float totaltime(0.0);
            auto destvx(cx);
            auto destvy(cy);
            auto destvz(cz);
            for( int it(0); it < steps; it++ ){
                if( stepe[it] > maxe ) maxe = stepe[it];
                auto beta = stepp[it]/stepe[it];
                auto legtime = hypo( (destvx-orignvx[it]), (destvy-orignvy[it]), (destvz-orignvz[it]) )/(SOL*beta);
                totaltime += legtime;
                //std::cout << "Sum path legs beta: " << beta << " x: " << (destvx-orignvx[it]) << " y: " << (destvy-orignvy[it]);
                //std::cout << " z: " << (destvz-orignvz[it]) << " legtime: " << legtime << " totaltime: " << totaltime << std::endl;
                //std::cout << " -- dest : " << destvx << ", " << destvy << ", " << destvz; 
                //std::cout << " orig : " << orignvx[it] << ", " << orignvy[it] << ", " << orignvz[it] << std::endl;
                destvx = orignvx[it];
                destvy = orignvy[it];
                destvz = orignvz[it];
            }//<<>>for( int it(0); it < steps; it++ )
            auto r = hypo((cx-orignvx[0]), (cy-orignvy[0]));
            auto iangle = getAngle((cz-orignvz[0]),r);
            //std::cout << "getAngle: " << iangle << " r: " << r << " z: " << cz-orignvz[0] << " x: " << cx-orignvx[0];
            //std::cout << " y: " << cy-orignvy[0] << std::endl; 
            if( totaltime > 25.0 ) {
                nextBX = 1.0;
            } else {
                impact.push_back(iangle);
                kidtime.push_back(totaltime);
                kide.push_back(maxe);
                llpcut = llpcut || llp || bquark;
                if( llp ) llpkid++;
            }//<<>>if( totaltime > 25.0 )
        }//<<>>if( steps == 0 )

    }//<<>>for( auto kid : kids )

    auto mtime(-30.0);
    auto mangle(7.0);
    auto tvar(-2.0);
    auto ktsize = kidtime.size();
    //if( llpcut && ( ktsize > 0 ) ){
    if( ktsize > 0 ){
        if( ktsize > 1 ){
            mtime = mean(kidtime,kide);
            mangle = mean(impact,kide);
            tvar = var(kidtime,mtime,kide);
        } else {
            mtime = kidtime[0];
            mangle = impact[0];
            tvar = -1.0;
        }//<<>>if( kidtime.size() > 1 )
    }//<<>>if( llpcut )
    result.push_back(mtime);//0
    result.push_back(mangle);//1
    result.push_back(tvar);//2
    result.push_back(nextBX);//3
    if( llpcut ) result.push_back(1.0); else result.push_back(0.0);//4
    result.push_back(llpkid/ktsize);//5
    result.push_back(ktsize);//6
    return result;

}//>>>>vector<float> KUCMSGenObject::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  )

void KUCMSGenObject::kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz ){

    for( auto kid : kids ){

        bool done(false);
        auto mom = kid->mother();
        std::cout << " -- Kid pdgID : " << kid->pdgId() << " p : " << kid->p() << " e : " << kid->energy() << std::endl;
        bool top(true);
        while( not done ){
            if( mom->pt() == 0 ) done = true;
            else {
                auto gmom = mom->mother();
                if( top ){
                    std::cout << " ---- KidStep pdgID : " << mom->pdgId() << " p : " << mom->p() << " e : " << mom->energy();
                    std::cout << " vert : (" << mom->vx() << "," << mom->vy() << "," << mom->vz() << ")" << std::endl;
                }//<<>>if( top )
                if( gmom->pt() == 0 ) done = true;
                else {
                    if( mom->vx() == gmom->vx() ) top = false;
                    else top = true;
                    mom = gmom;
                }//<<>>if( gmom->pt() == 0 )
            }//<<>>if( mom->pt() == 0 )                 
        }//<<>>while( not done )
    }//<<>>for( auto kid : kids )

}//<<>>string kidChase( std::vector<reco::CandidatePtr> kids, float vx )

int KUCMSGenObject::getEvtType( int genPdgId ){

	int type = 9;
	genPdgId = std::abs( genPdgId );
	//"1=squark, 2=gluino, 3=Xino, 4=slepton"
	if( genPdgId == 1000021 ) type = 2;
	else if( genPdgId > 1000000 && genPdgId < 1000007 ) type = 1;
	else if( genPdgId > 2000000 && genPdgId < 2000007 ) type = 1;
    else if( genPdgId > 1000021 && genPdgId < 1000038 ) type = 3;
    else if( genPdgId > 1000010 && genPdgId < 1000017 ) type = 4;
    else if( genPdgId > 2000010 && genPdgId < 2000017 ) type = 4;
	//if( type == 9 ) std::cout << " -- found a : " << genPdgId << std::endl;
	return type;

}//<<>>int KUCMSGenObject::getEvtType( int genPdgId )

const reco::Candidate* KUCMSGenObject::llpGenSearch( const reco::Candidate* kid, int target ){

	int depth = 0;
	const reco::Candidate* notfound = nullptr;
    const reco::Candidate* result = nullptr;
    auto kidPdgID = std::abs(kid->pdgId());
    if( GenDEBUG ) std::cout << " --- llpgen chase genpart: " << kidPdgID << " depth: " << depth << std::endl;
    bool wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? true : false;
    bool phFlag = ( kidPdgID == 22 ) ? true : false;
    //bool eleFlag = ( kidPdgID == 11 ) ? true : false;
    bool qFlag = ( kidPdgID < 7 ) ? true : false;
    bool pFlag = ( kidPdgID == 2212 ) ? true : false;
    //bool n1Flag = ( kidPdgID > 1000021 && kidPdgID < 1000026 ) ? true : false;
    //bool n2Flag = ( kidPdgID == 1000035 || kidPdgID == 1000037 ) ? true : false;

    if( pFlag ) return notfound;
    int nMoms = kid->numberOfMothers();
    if( nMoms == 0 ) return notfound;
    for( int gmit(0); gmit < nMoms; gmit++ ){

        auto genmom = kid->mother(gmit);
        auto momPdgID = std::abs(genmom->pdgId());
        bool ni1 = ( momPdgID == 1000022 ) ? true : false;
        bool ni2 = ( momPdgID == 1000023 ) ? true : false;
        bool ni3 = ( momPdgID == 1000025 ) ? true : false;
        bool ni4 = ( momPdgID == 1000035 ) ? true : false;
        bool ci1 = ( momPdgID == 1000024 ) ? true : false;
        bool ci2 = ( momPdgID == 1000037 ) ? true : false;
        bool sq1 = ( momPdgID > 1000000 && momPdgID < 1000007 ) ? true : false;
        bool sq2 = ( momPdgID > 2000000 && momPdgID < 2000007 ) ? true : false;
        //bool qrk = ( momPdgID < 7 || momPdgID == 21 ) ? true : false;

        bool pro = ( momPdgID == 2212 ) ? true : false;
        bool gli = ( momPdgID == 1000021 ) ? true : false;
        bool nti = ni1 || ni2 || ni3 || ni4;
        bool chi = ci1 || ci2;
        bool sqk = sq1 || sq2;

		int type = 100;
        if( phFlag && ni1 ){ if( depth == 0 ) type = 22; else type = 32; }
        else if( phFlag && nti ){ if( depth == 0 ) type = 25; else type = 35; }
        else if( wzFlag && nti ){ if( depth == 0 ) type = 23; else type = 33; }
        else if( wzFlag && chi ){ if( depth == 0 ) type = 24; else type = 34; }
        else if( qFlag && gli ){ if( depth == 0 ) type = 21; else type = 31; }
        else if( qFlag && sqk ){ if( depth == 0 ) type = 20; else type = 30; }
        else if( chi || nti ) type = 36;
        //else if( ( n1Flag || n2Flag ) && ( sqk || gli || pro || qrk ) ) type = 36;
        else if( gli || sqk ) type = 37;
        //else if( nti || chi ) type = 38;
        else if( pro ) type = 97;

		if( type == target ){ result = kid; return result; }
		else if( type == 97 || type == 37 ){ return result; }
        result = llpGenSearch(kid->mother(gmit),target);
		if( result ) return result;
        //if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    //if( GenDEBUG ) std::cout << " --- Found llpgen chase genpart: " << kidPdgID << " susid: " << genPartID << std::endl;
    return result;

}//<<>> void llpChase( Candidate* kid )

int KUCMSGenObject::llpGenChaseP( const reco::GenJet & kid, int depth ){


    int genPartID = 100;
    int nMoms = kid.numberOfMothers();
    if( (nMoms == 1) && (kid.pdgId() == (kid.mother(0))->pdgId()) ) return 101;
    int genPartId = std::abs(kid.pdgId());
    if( GenDEBUG ) std::cout << "Genjet motherCase for : " << genPartId << std::endl;
    for( int gmit(0); gmit < nMoms; gmit++ ){
        //motherChase(kid.mother(gmit),0);
        auto mGenPartID = llpGenChaseP(kid.mother(gmit),genPartId);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;
    }//<<>>for( int gmit(0); gmit < nMoms; gmit++ )
    return genPartID;

}//<<>>int KUCMSGenObject::llpGenChaseP( const reco::GenJet & kid, int depth )

int KUCMSGenObject::llpGenChaseP( const reco::GenParticle & kid, int depth ){

    int genPartSusID = 100;
    int nMoms = kid.numberOfMothers();
	if( nMoms == 1 && kid.pdgId() == kid.mother(0)->pdgId() ) return 101;
    int genPartId = std::abs(kid.pdgId());
    //if( GenDEBUG ) std::cout << "Genpart motherCase for : " << genPartId << std::endl;
    for( int gmit(0); gmit < nMoms; gmit++ ){
        auto mGenPartSusID = llpGenChaseP(kid.mother(gmit),genPartId);
        if( mGenPartSusID < genPartSusID ) genPartSusID = mGenPartSusID;
    }//<<>>for( int gmit(0); gmit < nMoms; gmit++ )
    return genPartSusID;

}//<<>>int KUCMSGenObject::llpGenChaseP( const reco::GenJet & kid, int depth )

int KUCMSGenObject::llpGenChaseP( const reco::Candidate* kid, int childPdgId ){

    int nMoms = kid->numberOfMothers();
    if( nMoms == 0 ) return 99;
	if( nMoms == 1 && kid->pdgId() == kid->mother(0)->pdgId() ) return llpGenChaseP( kid->mother(0), childPdgId );

    auto kidPdgID = std::abs(kid->pdgId());
	//if( GenDEBUG ) std::cout << " --- llpgen chase canadate: " << kidPdgID << " depth: " << depth << std::endl;
    bool wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? true : false;
    bool phFlag = ( kidPdgID == 22 ) ? true : false;
    //bool eleFlag = ( kidPdgID == 11 ) ? true : false;
    bool qFlag = ( kidPdgID < 7 ) ? true : false;
    bool gFlag = ( kidPdgID == 21 ) ? true : false;
    bool pFlag = ( kidPdgID == 2212 ) ? true : false;

    bool gliFlag = ( kidPdgID == 1000021 ) ? true : false;
    bool sq1Flag = ( kidPdgID > 1000000 && kidPdgID < 1000007 ) ? true : false;
    bool sq2Flag = ( kidPdgID > 2000000 && kidPdgID < 2000007 ) ? true : false;
    bool sqkFlag = sq1Flag || sq2Flag;

    //bool n0Flag = ( kidPdgID  == 1000022 ) ? true : false;
    bool n0Flag = ( kidPdgID  == 1000022 || kidPdgID  == 1000023 ) ? true : false;
	bool ntiFlag = ( kidPdgID == 1000022 || kidPdgID == 1000023 || kidPdgID == 1000025 || kidPdgID == 1000035 ) ? true : false;
    //bool ntiFlag = ( kidPdgID == 1000025 || kidPdgID == 1000035 ) ? true : false;
    bool chiFlag = ( kidPdgID == 1000024 || kidPdgID == 1000037 ) ? true : false;

    bool slpFlag = ( kidPdgID > 1000010 && kidPdgID < 1000017  ) ? true : false;

    bool cphFlag = ( childPdgId == 22 ) ? true : false;
    bool cwzFlag = ( childPdgId == 23 || childPdgId == 24 ) ? true : false;
    bool cqFlag = ( childPdgId > 0 && childPdgId < 7 ) ? true : false;
    bool cgFlag = ( childPdgId == 21 ) ? true : false;

    if( pFlag ) return 98;

    if( cphFlag && n0Flag ) return 22;
    else if( cphFlag && ( ntiFlag || chiFlag ) ) return 25;
    else if( cwzFlag && ntiFlag ) return 23;
    else if( cwzFlag && chiFlag ) return 24;

    else if( cqFlag && gliFlag ) return 41;
    else if( cqFlag && sqkFlag ) return 42;
    else if( cqFlag && slpFlag ) return 43;
    else if( cgFlag && gliFlag ) return 44;
    else if( cgFlag && sqkFlag ) return 45;
    else if( cgFlag && slpFlag ) return 46;

    int genPartID = 100;
    for( int gmit(0); gmit < nMoms; gmit++ ){

        auto genmom = kid->mother(gmit);
        auto momPdgID = std::abs(genmom->pdgId());
        bool ni1 = ( momPdgID == 1000022 ) ? true : false;
        bool ni2 = ( momPdgID == 1000023 ) ? true : false;
        bool ni3 = ( momPdgID == 1000025 ) ? true : false;
        bool ni4 = ( momPdgID == 1000035 ) ? true : false;
        bool ci1 = ( momPdgID == 1000024 ) ? true : false;
        bool ci2 = ( momPdgID == 1000037 ) ? true : false;

        bool sq1 = ( momPdgID > 1000000 && momPdgID < 1000007 ) ? true : false;
        bool sq2 = ( momPdgID > 2000000 && momPdgID < 2000007 ) ? true : false;
        bool qrk = ( momPdgID < 7 || momPdgID == 21 ) ? true : false;
        bool slp = ( momPdgID > 1000010 && momPdgID < 1000017 ) ? true : false;

        bool pro = ( momPdgID == 2212 ) ? true : false;
        bool gli = ( momPdgID == 1000021 ) ? true : false;
        bool nti = ni1 || ni2 || ni3 || ni4;
        bool chi = ci1 || ci2;
        bool sqk = sq1 || sq2;

        if( phFlag && ni1 ) return 32;
        else if( phFlag && ( nti || chi ) ) return 35;
        else if( wzFlag && nti ) return 33;
        else if( wzFlag && chi ) return 34;

        else if( qFlag && gli ) return 51;
        else if( qFlag && sqk ) return 52;
        else if( qFlag && slp ) return 53;
        else if( gFlag && gli ) return 54;
        else if( gFlag && sqk ) return 55;
        else if( gFlag && slp ) return 56;

        else if( ( ntiFlag || chiFlag ) && ( slp || sqk || gli || pro || qrk ) ) return 67;
        else if( ( gliFlag || sqkFlag ) && pro ) return 68;
        else if( slpFlag && pro ) return 69;
		//else if( gli || sqk || slp ) return 96;
        else if( pro ) return 97;

        auto mGenPartID = llpGenChaseP(kid->mother(gmit),0);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    //if( GenDEBUG ) std::cout << " --- Found llpgen chase canadate: " << kidPdgID << " susid: " << genPartID << std::endl;
    return genPartID;

}//<<>> void llpChase( const reco::Candidate* kid, int wzFlag ){)

void KUCMSGenObject::keepMothers( const reco::GenParticle kid ){

    if( kid.numberOfMothers() < 1 ) return;
	bool notLastInChain( ( kid.numberOfMothers() == 1 ) && ( kid.pdgId() == (kid.mother(0))->pdgId() ) );
	bool validStatus( kid.status() < 40 );
    if( notLastInChain || not validStatus ){ keepMothers( kid.mother(0) ); }
    else{
		// std::cout << " ------ Getting Moms for : " << kid.pdgId() << std::endl;
        for( long unsigned int gmit(0); gmit < kid.numberOfMothers(); gmit++ ){
            auto mother = kid.mother(gmit);
			if( std::abs(mother->pdgId()) == 2212 ) continue;
            if( std::abs(mother->pdgId()) < 10 || std::abs(mother->pdgId()) == 21 ) continue;
            keepMothers( mother );
        }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    }//else

}//<<>>void KUCMSGenObject::keepMothers( const reco::GenParticle kid )

void KUCMSGenObject::keepMothers( const reco::Candidate* kid ){

	if( kid->numberOfMothers() < 1 ) return;
	bool notLastInChain( ( kid->numberOfMothers() == 1 ) && ( kid->pdgId() == (kid->mother(0))->pdgId() ) );
    bool validStatus( kid->status() < 40 );
    if( notLastInChain || not validStatus ){ keepMothers( kid->mother(0) ); }
    else{
		//  -------   match mother and mark keep  --------------
		float mpt = kid->pt();
		int mpdg = kid->pdgId();
		int tgidx = 0;
		bool match = false;
		//if( mpdg == 1000023 ) std::cout << " -- Gen Part 1000023 found !!!!!! " << std::endl;
		for( auto gp : tgenparts ){
			float canpt = gp.pt();
        	int canpdg = gp.pdgId();
			if( canpt == mpt && canpdg == mpdg ){
				//if( canpdg == 1000023 ) std::cout << " --- Match: " << canpt << " : " << mpt << std::endl; 
				//std::cout << " -- Gen Mom match found !!!!!! " << std::endl;
        		//std::cout << " --- kidID : " << kid->pdgId() << " pt : " << kid->pt();
        		//std::cout << " Vertix (" << kid->vx() << "," << kid->vy() << "," << kid->vz() << ")";
        		//std::cout << " nMothers " << kid->numberOfMothers() << std::endl;
				tgpKeep[tgidx] = true; 
				match = true;
				break;
			}//<<>>if( canpt == mpt && canpdg == mpdg )
			tgidx++;
		}//<<>>for( auto gp : tgenparts ) 
		if( not match ) std::cout << " Failed Mother Match !!!!! " << std::endl;		
		//  ----------------------------------------------------
		// std::cout << " ------ Getting Moms for : " << kid->pdgId() << std::endl;
        for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
            auto mother = kid->mother(gmit);
			if( std::abs(mother->pdgId()) == 2212 ) continue;
            if( std::abs(mother->pdgId()) < 10 || std::abs(mother->pdgId()) == 21 ) continue;
            keepMothers( mother );
        }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    }//else

}//<<>> void keepMothers( Candidate* kid )

void KUCMSGenObject::motherChase( const reco::Candidate* kid, string  depth ){

    //if( ( kid->numberOfMothers() == 1 ) && ( kid->pdgId() == (kid->mother(0))->pdgId() ) ){ motherChase( kid->mother(0), depth ); }
    if( false );
    else{
	    if( true ){ //status()
		//if( kid->status() < 40 ){
    		std::cout <<  depth  << " kidID : " << kid->pdgId() << " pt : " << kid->pt();
    		std::cout << " Vertix (" << kid->vx() << "," << kid->vy() << "," << kid->vz() << ")";
			if( kid->numberOfMothers() > 0 ) std::cout << " mother0 " << kid->mother(0)->pdgId();
    		std::cout << " nMothers " << kid->numberOfMothers() << " status " << kid->status() << std::endl;
    		depth += "-";
		}//<<>>if( kid.status() < 30 )
        if( std::abs(kid->pdgId()) > 100 && std::abs(kid->pdgId()) < 1000000 ) return;
    	for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
        	auto mother = kid->mother(gmit);
            //if( std::abs(mother->pdgId()) < 10 || std::abs(mother->pdgId()) == 21 ) continue;
			//if( std::abs(kid->pdgId()) > 1000000 && ( std::abs(mother->pdgId()) < 10 || std::abs(mother->pdgId()) == 21 ) ) continue;
        //;//std::cout <<  depth  << " gMomID : " << kid->mother(gmit)->pdgId() << " pt : " << kid->mother(gmit)->pt();
        //;//std::cout << " Vertix (" << kid->mother(gmit)->vx() << "," << kid->mother(gmit)->vy() << "," << kid->mother(gmit)->vz() << ")";
        //;//std::cout << " nGMothers " << kid->mother(gmit)->numberOfMothers() << std::endl;
        	motherChase( mother, depth );
    	}//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
	}//<<>>else
    //;//if( not depth.empty() ) depth = depth.substr (0,depth.length()-1);

}//<<>> void MotherChase( Candidate* kid, string depth  )

void KUCMSGenObject::motherChase( const reco::GenParticle kid, string  depth ){
 
	//if( ( kid.numberOfMothers() == 1 ) && ( kid.pdgId() == (kid.mother(0))->pdgId() ) ){ motherChase( kid.mother(0), depth ); }
	if( false );
	else{
	    if( true ){ //status()
	    //if( kid.status() < 40 ){
    		std::cout <<  depth  << " kidID : " << kid.pdgId() << " pt : " << kid.pt();
    		std::cout << " Vertix (" << kid.vx() << "," << kid.vy() << "," << kid.vz() << ")";
            if( kid.numberOfMothers() > 0 ) std::cout << " mother0 " << kid.mother(0)->pdgId();
    		std::cout << " nMothers " << kid.numberOfMothers() << " status " << kid.status() << std::endl;
    		depth += "-";
		}//<<>>if( kid.status() < 30 ){
		if( std::abs(kid.pdgId()) > 100 && std::abs(kid.pdgId()) < 1000000 ) return;
    	for( long unsigned int gmit(0); gmit < kid.numberOfMothers(); gmit++ ){
        	auto mother = kid.mother(gmit);
            //if( std::abs(mother->pdgId()) < 10 || std::abs(mother->pdgId()) == 21 ) continue;
			//if( std::abs(kid.pdgId()) > 1000000 && ( std::abs(mother->pdgId()) < 10 || std::abs(mother->pdgId()) == 21 ) ) continue;
        //std::cout <<  depth  << " gMomID : " << kid->mother(gmit)->pdgId() << " pt : " << kid->mother(gmit)->pt();
        //std::cout << " Vertix (" << kid->mother(gmit)->vx() << "," << kid->mother(gmit)->vy() << "," << kid->mother(gmit)->vz() << ")";
        //std::cout << " nGMothers " << kid->mother(gmit)->numberOfMothers() << std::endl;
        	motherChase( mother, depth );
   		}//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
	}//else
    //;//if( not depth.empty() ) depth = depth.substr (0,depth.length()-1);

}//<<>> void MotherChase( Candidate* kid, string depth  )

int KUCMSGenObject::typeChase( const reco::Candidate* kid ){

    int type = -9999;
    if( kid->numberOfMothers() == 2 ){ type = kid->pdgId(); } //std::cout << " kidID : " << kid->pdgId() << " pt : " << kid->pt() << std::endl; }
    else if( kid->numberOfMothers() == 1 ){ type = typeChase( kid->mother(0) ); }
    else { type = 0; } //std::cout << " nMothers " << kid->numberOfMothers() << std::endl; }
    return std::abs(type);

}//<<>> void MotherChase( Candidate* kid, string depth  )

int KUCMSGenObject::typeChase( const reco::GenParticle kid ){

	int type = -9999;
    if( kid.numberOfMothers() == 2 ){ type = kid.pdgId(); }//std::cout << " kidID : " << kid.pdgId() << " pt : " << kid.pt() << std::endl; }
    else if( kid.numberOfMothers() == 1 ){ type = typeChase( kid.mother(0) ); }
    else { type = 0; }//std::cout << " nMothers " << kid.numberOfMothers() << std::endl; }
    return std::abs(type);

}//<<>> void MotherChase( Candidate* kid, string depth  )

int KUCMSGenObject::llpGenChase( const reco::Candidate* kid ){

    auto match = ( std::abs(kid->pdgId()) > 1000020 ) && ( std::abs(kid->pdgId()) < 1000040 );
    if( match ) return kid->pdgId();
    for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
        if( kid->mother(gmit)->numberOfMothers() > 0 ){ auto r = llpGenChase(kid->mother(gmit)); if(r != 0) return r; }
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return 0;

}//<<>> void llpChase( Candidate* kid )

int KUCMSGenObject::llpGenChase( const reco::GenParticle kid ){

    auto match = ( std::abs(kid.pdgId()) > 1000020 ) && ( std::abs(kid.pdgId()) < 1000040 );
    auto hard = kid.isLastCopy();
    if( match && hard ) return kid.pdgId();
    for( long unsigned int gmit(0); gmit < kid.numberOfMothers(); gmit++ ){
        //if( kid->pdgId() == 6000113 ) return true;
        if( kid.mother(gmit)->numberOfMothers() > 0 ){ auto r = llpGenChase(kid.mother(gmit)); if(r != 0) return r; }
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return 0;

}//<<>> void llpChase( Candidate* kid )

string KUCMSGenObject::bigKidChase( std::vector<reco::CandidatePtr> kids, float vx ){

    std::string result("");
    float maxpt(0.0);
    const reco::Candidate* bigKid(NULL);
    for( auto &kid : kids ){
        if( kid->pt() > maxpt ){
            maxpt = kid->pt();
            bigKid = kid.get();
        }//<<>>if( kid->pt() > maxpt )
    }//<<>>for( auto kid : kids )
    result += std::string(" -- BigKid : pdgID : ") + bigKid->pdgId() + std::string(" pt : ") + bigKid->pt();
    result += std::string(" vtx : ") + bigKid->mother(0)->vx() + std::string(" / ") + vx ;
    return result;

}//<<>>string bigKidChase( std::vector<reco::CandidatePtr> kids, float vx )

std::vector<float> KUCMSGenObject::getGenJetInfo( float jetEta, float jetPhi, float jetEnergy ){

    float genImpactAngle(-10.0);
    float genTime(-99.9);
    float genPt(-10.0);
    float genEta(-10.0);
    float genPhi(-10.0);
    float genEnergy(-10.0);
    float genEMFrac(-1.0);
    float genDrMatch(-1.0);
    float genReMatch(-1.0);
    float genTimeLLP(-99.9);
    float genTOF(-99.9);
    float genPartLlpId(101);
    float genJetLlpId(101);

	//std::cout << " --- Jet-GenJet matching " << std::endl;

    float goodDr(0.3);
    float goodRe(-1.0);
    int matchedIdx(-1);
    int index(0);
    for(const auto &genJet : fgenjets ){

        auto gjeta = genJet.eta();
        auto gjphi = genJet.phi();
        auto gje = genJet.energy();
		if( gje < 1.0 ){ index++; continue; }
        auto jtgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jetEta, jetPhi ));
        auto jtgjre = jetEnergy/gje;
        if( jtgjdr < goodDr ){
            goodDr = jtgjdr;
            goodRe = jtgjre;
            matchedIdx = index;
        }//<<>>if( jtgjdr <= goodDr )
        index++;

    }//<<>>for(const auto& genJet : fgenjets ) 
	

	//std::cout << " --- Jet-Genpart quark matching" << std::endl;
    float genPartDr(0.3);
    float genPartRe(-1.0);
    int mGenPartIdx(-1);
    int idx(0);
    for(const auto &genPart : fgenparts ){

		if( std::abs(genPart.pdgId()) > 9 ){ idx++; continue; }
        auto gjeta = genPart.eta();
        auto gjphi = genPart.phi();
        auto gje = genPart.energy();
        if( gje < 1.0 ){ idx++; continue; }
        auto jgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jetEta, jetPhi ));
        auto jtgjre = jetEnergy/gje;
        if( jgjdr < genPartDr ){
            genPartDr = jgjdr;
            genPartRe = jtgjre;
            mGenPartIdx = idx;
        }//<<>>if( jtgjdr <= goodDr )
        idx++;

    }//<<>>for(const auto& genJet : fgenjets ) 

    if( mGenPartIdx > -1 ){ 
		genPartLlpId = float(fgpLlp[mGenPartIdx]);
    	//if( genPartLlpId == 51 ) motherChase( fgjquarks[mGenPartIdx], "" );
	} else { genPartDr = -1.0; }

	//std::cout << " --- Jet-GenJet genjet info " << std::endl;
    if( matchedIdx > -1 ){

        auto genJet = fgenjets[matchedIdx];
		genJetLlpId = fgenjetllp[matchedIdx];
        //if( GenDEBUG ) std::cout << " --- Jet-GenJet dR match : " << goodDr << std::endl;
        genEta = genJet.eta();
        genPhi = genJet.phi();
        //if( GenDEBUG ) std::cout << " - genJet GenTime : " << genTime << " Angle: " << genImpactAngle << std::endl;
        genPt = genJet.pt();
        genEnergy = genJet.energy();
        genEMFrac = (genJet.chargedEmEnergy() + genJet.neutralEmEnergy())/genEnergy;
        genDrMatch = goodDr; //std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), genJet.eta(), genJet.phi()));
        genReMatch = goodRe;
        //if( GenDEBUG ) std::cout << " -- Energy : " << genEnergy << " Pt : " << genPt << " EMfrac : " << genEMFrac << std::endl;

    }//<<>>if( matchedIdx >= 0 )
	else { goodDr = -1.0; }

	std::vector<float> results;

    results.push_back(genImpactAngle);
    results.push_back(genTime);
    results.push_back(genPt);//2
    results.push_back(genEta);
    results.push_back(genPhi);
    results.push_back(genEnergy);//5
    results.push_back(genEMFrac);
    results.push_back(genDrMatch);
    results.push_back(genReMatch);//8
    results.push_back(genTimeLLP);
    results.push_back(genTOF);

    results.push_back(goodDr);//11
    results.push_back(goodRe);	
    results.push_back(genJetLlpId);

    results.push_back(genPartDr);//14
    results.push_back(genPartRe);
    results.push_back(genPartLlpId);

    return results;

}//<<>>std::vector<float> KUCMSGenObject::getGenJetInfo( float jetEta, float jetPhi, float jetPt )

// Gen Electron methods


/*
std::vector<int>KUCMSGenObject:: MomIDs(const reco::GenParticle &genElectron) const {

  auto mother = genElectron.mother();


  std::vector<int> motherIDs;

  //if( genElectron.numberOfMothers() < 1 ) return motherIDs; 
/*
  auto mother = genElectron.mother(0);
  while(mother->pt() > 0) {
    const int motherID = mother->pdgId();

    mother = mother->mother();

    if(motherID == mother->pdgId())
      continue;

    motherIDs.push_back(motherID);
  }
*/

  return motherIDs;

}

bool KUCMSGenObject::IsMotherZ(const reco::GenParticle &genElectron) const {
  LepMomType momType = ClassifyGenElectron(genElectron);
  return (momType == kZ);
}

bool KUCMSGenObject::isSignalGenElectron(const reco::GenParticle &genElectron) const {

  LepMomType momType = ClassifyGenElectron(genElectron);
  return (momType == kZ || momType == kSusy);

}

LepMomType KUCMSGenObject::ClassifyGenElectron(const reco::GenParticle &genElectron) const {

  std::vector<int> motherIDs(MomIDs(genElectron));

  LepMomType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignLeptonMomType(id);

    if(momType != kUnmatched)
      break;
  }
  return momType;
}

void KUCMSGenObject::PrintMother(const LepMomType &momType) const {

  std::cout << ": ";
  if (momType == kW) std::cout << "W boson" << std::endl;
  else if (momType == kZ) std::cout << "Z boson" << std::endl;
  else if (momType == kTau)  std::cout << "tau lepton" << std::endl;
  else if (momType == kConversion) std::cout << "photon (conversion)" << std::endl;
  else if (momType == kLight) std::cout << "light quark" << std::endl;
  else if (momType == kHeavy) std::cout << "heavy quark" << std::endl;
  else if (momType == kSusy) std::cout << "prompt from Susy particle" << std::endl;
  else if (momType == kUnmatched) std::cout << "unmatched (this shouldn't happen!)" << std::endl;

}

/*
void KUCMSGenObject::GenElectronContent() const {

  int index(0);
  for(const auto & genElectron : genElectrons_ ) {

    const int genCharge(genElectron.charge());
    const float gnX(genElectron.vx() );
    const float gnY(genElectron.vy() );
    const float gnZ(genElectron.vz() );
    const float distXY(sqrt(gnX*gnX + gnY*gnY));
    const float genPt(genElectron.pt());
    const float genEta(genElectron.eta());
    const float genPhi(genElectron.phi());
    const float genEnergy(genElectron.energy());

    std::vector<int> motherIDs = MomIDs(genElectron);

    LepMomType momType = ClassifyGenElectron(genElectron);

    if (genSignalElectrons_.size() > 2 && index == 0)
      std::cout << "There are " << genElectrons_.size() << " gen electrons in this event." << std::endl;

    if ( (momType == kZ || momType == kSusy ) && genSignalElectrons_.size() > 2) {
      std::cout << "gen electron " << index << ": " << std::endl;
      std::cout << "  vertex position: (" << gnX << ", " << gnY <<  ", " << gnZ << ")" << std::endl;
      std::cout << "  Total transverse displacement: " << distXY << std::endl;
      std::cout << "  pT: " << genPt << ", eta: " << genEta << ", phi: " << genPhi << ", energy: " << genEnergy << std::endl;
      std::cout << "  gen charge: " << genCharge <<  std::endl;
      for(auto const& id : motherIDs)
	std::cout << "  motherID: " << id << std::endl;
      PrintMother(momType);
      std::cout << std::endl;

    }

    index++;
  }
}
*/

LepMomType KUCMSGenObject::AssignLeptonMomType(const int motherID) const {
  
  LepMomType type = kUnmatched;

  if(abs(motherID) == 24)
    type = kW;
  else if(motherID == 23)
    type = kZ;
  else if(abs(motherID) == 15)
    type = kTau;
  else if((abs(motherID%1000) > 100 && abs(motherID%1000) < 400)
          || (abs(motherID%1000) > 1000 && abs(motherID%1000) < 4000)
          || (abs(motherID) > 0 && abs(motherID) < 4)
          || motherID == 21)
    type = kLight;
  else if((abs(motherID%1000) > 400 && abs(motherID%1000) < 600)
          || (abs(motherID%1000) > 4000 && abs(motherID%1000) < 6000)
          || (abs(motherID) > 3 && abs(motherID) < 7))
    type = kHeavy;
  else if(motherID == 22)
    type = kConversion;
  else if(motherID == 1000022 || abs(motherID) == 9000006)
    type = kSusy;
  else if(abs(motherID) > 1000000)
    type = kOther;

  return type;
  
}
 
template <typename T>
GenClassifiedElectrons<T> KUCMSGenObject::GetGenClassifiedElectrons(const std::vector<T> &candidateTracks) const {
  
  // Create output collection
  GenClassifiedElectrons<T> classifiedElectrons;

  // Match input tracks to Gen electrons using hungarian algorithm 
  DeltaRMatchHungarian<T, reco::GenParticle> assigner(candidateTracks, genElectrons_);
  
  std::vector<int> matchedIndexes;
  // Loop over matches
  for(const auto &pair : assigner.GetPairedObjects()) {
    const T track(pair.GetObjectA());
    const reco::GenParticle genElectron(pair.GetObjectB());
    const double deltaR(pair.GetDeltaR());
    const LepType type(ClassifyGenElectron(genElectron));
    
    ElectronGenType genType(track, type, deltaR);

    matchedIndexes.push_back(pair.GetIndexA());
    classifiedElectrons.emplace_back(genType);
  }

  // Classify any track that did not get matched in the event as unmatched
  for(int i = 0; i < candidateTracks.size(); i++) {
   
    if (std::find(matchedIndexes.begin(), matchedIndexes.end(), i) != matchedIndexes.end())
      continue;
    
    ElectronGenType genType(candidateTracks[i], kUnmatched);
    classifiedElectrons.emplace_back(genType);
  }

  return classifiedElectrons;
}
*/
#endif
