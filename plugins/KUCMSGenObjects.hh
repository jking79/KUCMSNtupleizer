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

#include "KUNTupleFW/KUCMSNtupleizer/interface/MatchingTools.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/DeltaRMatch.h"

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

    string bigKidChase( std::vector<reco::CandidatePtr> kids, float vx );
    void motherChase( const reco::Candidate* kid, string & depth );
    int llpGenChase( const reco::Candidate* kid );
    int llpGenChase( const reco::GenParticle kid );
    int llpGenChaseP( const reco::Candidate* kid, int depth );
    int llpGenChaseP( const reco::GenParticle & kid, int depth );
    const reco::Candidate* llpGenSearch( const reco::Candidate* kid, int target );
    void kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz );
    std::vector<float> kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  );
    // old single reco part gen matching
    std::vector<float> getGenPartMatch( const reco::SuperCluster* scptr, float pt );
    std::map<std::string, float> getGenPartMatch( const reco::SuperCluster &scptr, float pt ) const; 
    std::map<std::string, bool> MotherID(const int genIndex) const;
    std::vector<float> getGenJetInfo( float jetEta, float jetPhi, float jetPt );
    GlobalPoint GenVertex() {return GlobalPoint(genxyz0_->x(), genxyz0_->y(), genxyz0_->z()); };

    // new exclusive gen matching
    std::vector<int> getGenMatch( const std::vector<reco::SuperCluster> scptrs, std::vector<float> energies );
    std::vector<int> getGenPhoMatch( const std::vector<reco::SuperCluster> scptrs, std::vector<float> energies );
    std::vector<int> getGenEleMatch( const std::vector<reco::SuperCluster> scptrs, std::vector<float> energies );
    std::vector<int> getGenMatch( const std::vector<reco::SuperCluster> scptrs, std::vector<float> energies, int select );
    int getGenSigPhoInfo( uInt genIndex );
    int getGenSigEleInfo( uInt genIndex );

    // Gen electrons
    void GenElectronContent() const;
    std::vector<reco::GenParticle> GetSignalGenElectrons() const {return genSignalElectrons_;}
    std::vector<reco::GenParticle> GetGenElectrons() const {return genElectrons_;}
    std::vector<reco::GenParticle> GetGenParticles() const {return fgenparts;}
    LepType ClassifyGenElectron(const reco::GenParticle &genElectron) const;
    std::vector<int> MomIDs(const reco::GenParticle &genElectron) const; 
 //template <typename T>
    //GenClassifiedElectrons<T> GetGenClassifiedElectrons(const std::vector<T> &candidateTracks) const;
  //std::map<std::string, float> GenMatchElectrons(const Lorentz4Vec &electronTrack, const int charge) const;
  //std::map<std::string, float> getGenPartMatch( const reco::SuperCluster &scptr, float pt ) const;
  //std::map<std::string, bool> MotherID(const int genIndex) const;
    bool IsMotherZ(const reco::GenParticle &genElectron) const;

    private:

    std::vector<reco::GenParticle> fgenparts;
    std::vector<int> fgenpartllp;
    std::vector<reco::GenJet> fgenjets;
    //std::vector<reco::GenParticle> fgjquarks;
    //std::vector<int> fgjquarksllp;
    //std::vector<reco::GenParticle> fgenallpho;
    //std::vector<int> fgenallphollp;
    std::vector<reco::GenParticle> genElectrons_;
    std::vector<reco::GenParticle> genSignalElectrons_;

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
  LepType AssignLeptonMomType(const int motherID) const;
  //LepType ClassifyGenElectron(const std::vector<int> &motherIDs) const;
  //LepType ClassifyGenElectron(const reco::GenParticle &genElectron) const;
  bool isSignalGenElectron(const reco::GenParticle &genElectron) const;
  void PrintMother(const LepType &momType) const;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

};//<<>>class KUCMSg : public KUCMSObjectBase

KUCMSGenObject::KUCMSGenObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : false );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "minGenE", iConfig.existsAs<double>("minGenE") ? iConfig.getParameter<double>("minGenE") : 0.0 );

}//<<>>KUCMSg::KUCMSGen( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSGenObject::InitObject( TTree* fOutTree ){

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
    Branches.makeBranch("genWgt","Evt_genWgt",FLOAT);
    Branches.makeBranch("genSusId","Gen_susId",VINT);
    Branches.makeBranch("genStatus","Gen_status",VBOOL);
    Branches.makeBranch("genCharge","Gen_charge",VINT);
    Branches.makeBranch("genMass","Gen_mass",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSGen::InitObject( TTree* fOutTree )

void KUCMSGenObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    if( GenDEBUG ) std::cout << "Getting gen tokens" << std::endl;

    iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);
    iEvent.getByToken(gent0Token_, gent0_);
    iEvent.getByToken(genxyz0Token_, genxyz0_);
    iEvent.getByToken(pileupInfosToken_, pileupInfos_);
    iEvent.getByToken(genParticlesToken_, genParticles_);
    iEvent.getByToken(genJetsToken_, genJets_);

    fgenparts.clear();
    fgenpartllp.clear();
    fgenjets.clear();
    //fgjquarks.clear();
    //fgjquarksllp.clear();
    //fgenallpho.clear();
    //fgenallphollp.clear();
    genElectrons_.clear();
    genSignalElectrons_.clear();
    int nGenPart(0);

    if( GenDEBUG ) std::cout << "Collecting Gen Particles" << std::endl;
    for( const auto &genPart : *genParticles_ ){

        if( genPart.energy() < cfPrm("minGenE") ) continue;
        auto genPartSUSID = llpGenChaseP(genPart,0);
        nGenPart++;
        if( genPart.status() == 1 ){
            fgenparts.push_back(genPart);
            fgenpartllp.push_back(genPartSUSID);
	    // Gen electron collections
	    if(abs(genPart.pdgId()) == 11) {
	      genElectrons_.push_back(genPart);
	      //std::cout << "\nmatched as electron!" << std::endl;
	      //for(const auto &id : MomIDs(genPart))
	      //std::cout << "  id = " << id << std::endl;
	      if(isSignalGenElectron(genPart))
		genSignalElectrons_.push_back(genPart);
	    }
            if( GenDEBUG && genPartSUSID == 22 ) std::cout << "Found genParton signal photon" << std::endl;
        } else {
			//bool isPhoton( genPart.pdgId() == 22 );
			bool isZW( genPart.pdgId() == 23 || genPart.pdgId() == 24 );
            bool isQuark( genPart.pdgId() < 9 );
            bool isSQuark( genPart.pdgId() > 1000000 && genPart.pdgId() < 1000007 );
            bool isGlino( genPart.pdgId() == 1000021 );
            bool isNtlino( genPart.pdgId() == 1000022 || genPart.pdgId() == 1000023 || genPart.pdgId() == 1000025 || genPart.pdgId() == 1000035 );
            bool isChrino( genPart.pdgId() == 1000024 || genPart.pdgId() == 1000037 );
            bool isGrvtino( genPart.pdgId() == 1000039 );
			//if( GenDEBUG && isNtlino ) std::cout << "Found genParton isNtlino" << std::endl;
			//bool keep( isPhoton || isQuark );
        	if( isQuark || isZW || isNtlino || isChrino || isGrvtino || isSQuark || isGlino ){
            	fgenparts.push_back(genPart);
            	fgenpartllp.push_back(genPartSUSID);
        	}//<<>>if( genPart.pdgId() == 22 )
		}//<<>> else if( genPart.status() == 1 )

    }//<<>>for(const auto& genPart : *genParticles_ )

    if( GenDEBUG ) std::cout << "Collecting GenJets" << std::endl;
    for(const auto &genjet : *genJets_ ){ fgenjets.push_back(genjet); }

}//<<>>void KUCMSGen::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSGenObject::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSGenObject::ProcessEvent( ItemManager<float>& geVar ){

    //if( GenDEBUG ) std::cout << "Processing Gens" << std::endl;

    Branches.clearBranches();

    //if( GenDEBUG ) std::cout << " - enetering Gen loop" << std::endl;

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
        const bool genStatus = genpart.status();
		//const int genSusId = llpGenChaseP( genpart, 0 );
        const int genSusId = fgenpartllp[nGenParts];
        const int genCharge = genpart.charge();
        const float genVx = genpart.vx();
        const float genVy = genpart.vy();
        const float genVz = genpart.vz();
        const float genMass = genpart.mass();		

		//if( GenDEBUG ) std::cout << "GenPart : genSusId = " << genSusId << std::endl;

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
        Branches.fillBranch("genSusId",genSusId);
        Branches.fillBranch("genStatus",genStatus);
        Branches.fillBranch("genCharge",genCharge);
        Branches.fillBranch("genMass",genMass);

		nGenParts++;
    }//<<>> for (const auto genpart : fgenparts )

    //if( GenDEBUG ) std::cout << "GenPart : LOADING GEN WT " << std::endl;
    float wgt = genEvtInfo_->weight();
    Branches.fillBranch("genWgt",wgt);
    //if( GenDEBUG ) std::cout << "GenPart : Done " << std::endl;

}//<<>>void KUCMSGen::ProcessEvent()

void KUCMSGenObject::EndJobs(){}

std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<reco::SuperCluster> sc, std::vector<float> reco_e ){

	// gen matching to supcer clusters froom reeco phtons to gen photons
	return getGenMatch( sc, reco_e, 22 );

}//<<>>std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenEleMatch( const std::vector<reco::SuperCluster> sc, std::vector<float> reco_e ){

    // gen matching to supcer clusters froom reco electrons to gen electrons
    return getGenMatch( sc, reco_e, 11 );

}//<<>>std::vector<int> KUCMSGenObject::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenMatch( const std::vector<reco::SuperCluster> sc, std::vector<float> reco_e ){

	// matchoing any particle :: matching valid for nutral particles only
	return getGenMatch( sc, reco_e, 0 );

}//<<>>std::vector<int> KUCMSGenObject::getGenMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObject::getGenMatch( const std::vector<reco::SuperCluster> sc, std::vector<float> reco_e, int select ){

    // GenParticle Info   -------------------------------------------------------------------
    if( GenDEBUG ) std::cout << "Getting phoGenParton Match Information" << std::endl;
	// code caclutes eta & phi of the momentum vector requird for a nuetral particle 
	// produced at a gen particles production vertix to cause the given supercluster and compares
	// this to the actual gen particles mometum vectors eta & phi to find best, exclusive, match
	// between a reco particle with that supercluster and a gen particle

    std::vector<int> genPartIndx;
    std::vector<float> genMDR;
	std::vector<int> genRecoIndx;
	float drthres = 0.4;
	int nGenPart = fgenparts.size();
    for( int idx = 0; idx < nGenPart; idx++ ){
		if( select != 0 && std::abs(fgenparts[idx].pdgId()) != select ) continue;
        genPartIndx.push_back(idx);
        genMDR.push_back(drthres);
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
        		float cmEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
        		float cmPhi = std::atan2(rhY-gnY,rhX-gnX);
        		float dr = std::sqrt(reco::deltaR2( gnEta, gnPhi, cmEta, cmPhi ));
                if( dr < genMDR[pgidx] ){ genMDR[pgidx] = dr; genRecoIndx[pgidx] = it; }
            }//<<>>for( int pgidx = 0; pgidx < genSigPhoIdx.size(); pgidx++ )
        }//<<>>for( int it = 0; it < nSCPtrs; it++ )
        bool unassigned = false;
		if( GenDEBUG ) std::cout << " - Checking Dup Matches " << std::endl;
        for( long unsigned int pgidx = 0; pgidx < genMDR.size(); pgidx++ ){
            for( long unsigned int pgidx2 = pgidx+1; pgidx2 < genMDR.size(); pgidx2++ ){
                if( ( genRecoIndx[pgidx] != -9 ) && ( genRecoIndx[pgidx] == genRecoIndx[pgidx2] ) ){
                    unassigned = true;
                    if( genMDR[pgidx2] >= genMDR[pgidx] ){ genRecoIndx[pgidx2] = -9; genMDR[pgidx2] = drthres; }
                    else { genRecoIndx[pgidx] = -9; genMDR[pgidx] = drthres; }
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
            llpId = fgenpartllp[index];
        } else if( dr < minSDr ){
            minSDr = dr;
            minSRe = re;
            matchedSIdx = index;
            llpSId = fgenpartllp[index];
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

/*
<<<<<<< HEAD
std::map<std::string, float> KUCMSGenObject::getGenPartMatch( const reco::SuperCluster &scptr, float pt ) const {

  std::map<std::string, float> genInfo;
  float minDr(1.0);
  float minSDr(5.0);
  float minDp(1.0);
  float minSDp(5.0);
  int matchedIdx(-999);// 1
  int matchedSIdx(-999);// 1
  int index(0);

  for(const auto & genPart : fgenparts ) {

    const auto rhX = scptr.x();
    const auto rhY = scptr.y();
    const auto rhZ = scptr.z();
    const auto gpt = genPart.pt();
    const float gnX(genPart.vx() );
    const float gnY(genPart.vy() );
    const float gnZ(genPart.vz() );
    auto cphoEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
    auto cphoPhi = std::atan2(rhY-gnY,rhX-gnX);
    auto dr = std::sqrt(reco::deltaR2(genPart.eta(), genPart.phi(), cphoEta, cphoPhi ));
    auto dp = std::abs(gpt-pt)/gpt;

    //std::cout << "gen coordinates: (" << genPart.vertex().x() << ", " << genPart.vertex().y() << ", " << genPart.vertex().z() << ")" << std::endl;

    if( dr < minDr && dp < minDp ) {
      minDr = dr;
      minDp = dp;
      matchedIdx = index;
    } else if( dr < minSDr && dp < minSDp ) {
      minSDr = dr;
      minSDp = dp;
      matchedSIdx = index;
    }//<<>>if( dr < minDr && dp < minDp )
    index++;

  }//<<>>for(const auto& genPart : fgenparts  )

  // PDG ID defaults to 0 for unmatched
  int pdgId(0), pdgSId(0);//, momID(-1);
  float matDr(-999.0), matDp(-999.0), matSDr(-999.0), matSDp(-999.0);
  float genVtxX(-999.), genVtxY(-999.), genVtxZ(-999.);
  if( matchedIdx >=0 ) { 
    matDr = minDr; 
    matDp = minDp; 
    matSDr = minSDr; 
    matSDp = minDp; 
    pdgId = fgenparts[matchedIdx].pdgId();
    pdgSId = fgenparts[matchedSIdx].pdgId();
    genVtxX = fgenparts[matchedIdx].vx();
    genVtxY = fgenparts[matchedIdx].vy();
    genVtxZ = fgenparts[matchedIdx].vz();
  }
  
  genInfo["matchedIdx"] = matchedIdx;
  genInfo["matchedSIdx"] = matchedSIdx;
  genInfo["matchedId"] = pdgId;
  genInfo["matchedSId"] = pdgSId;
  genInfo["GenDR"] = matDr;
  genInfo["GenDpt"] = matDp;
  genInfo["GenSDR"] = matSDr;
  genInfo["GenSDpt"] = matSDp;
  genInfo["GenVertex_x"] = genVtxX;
  genInfo["GenVertex_y"] = genVtxY;
  genInfo["GenVertex_z"] = genVtxZ;

  return genInfo;
}

std::map<std::string, bool> KUCMSGenObject::MotherID(const int genIndex) const {
  
  std::map<std::string, bool> IDmap = {{"isW", false},
				       {"isZ", false},
				       {"isTau", false},
				       {"isLight", false},
				       {"isHeavy", false},
				       {"isConversion", false},
				       {"isUnmatched", false}};

  if(genIndex < 0) {
    IDmap["isUnmatched"] = true;
    return IDmap;
  }

  reco::GenParticle genCandidate = fgenparts[genIndex];

  LepType momType = kUnmatched;
  auto mother = genCandidate.mother();

  while(mother->pt() > 0) {
    momType = AssignLeptonMomType(mother->pdgId());

    if(momType != kUnmatched)
      break;

    mother = mother->mother();
  }

  if(momType == kUnmatched)
    momType = AssignLeptonMomType(mother->pdgId());

  switch (momType) {
  case kW:
    IDmap["isW"] = true;
    break;
  case kZ:
    IDmap["isZ"] = true;
    break;
  case kTau:
    IDmap["isTau"] = true;
    break;
  case kConversion:
    IDmap["isConversion"] = true;
    break;
  case kLight:
    IDmap["isLight"] = true;
    break;
  case kHeavy:
    IDmap["isHeavy"] = true;
    break;
  default:
    IDmap["isUnmatched"] = true;
    break;
  }

  if(IDmap["isUnmatched"])
    std::cout << "index: " << genIndex << ", momID: " << mother->pdgId() << std::endl;

  return IDmap;
}

KUCMSGenObject::LepType KUCMSGenObject::AssignLeptonMomType(const int motherID) const {

  KUCMSGenObject::LepType type = kUnmatched;

  if(abs(motherID) == 24) 
    type = kW;
  else if(motherID == 23 || motherID == 33)
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
  
return type;
}
=======*/
int KUCMSGenObject::getGenSigPhoInfo( uInt genIndex ){

	int result( -1 );
    if( genIndex > fgenpartllp.size() ){ std::cout << " --- INDEX Blown !!!!  getGenSigPhoInfo " << std::endl; return result; }
    auto gensusid = fgenpartllp[genIndex];
    //if( GenDEBUG ) std::cout << " -- getGenSigPhoInfo: " << genIndex << " gensusid: " << gensusid << std::endl;

	if( gensusid == 22 || gensusid == 25 ){
		
		//if( fgenparts[genIndex].numberOfMothers() == 0 ) return -1;
		//auto nino = fgenparts[genIndex].mother();

        const reco::Candidate* foundpart = nullptr;
        int searchtype = 1000100;
        int nMothers = fgenparts[genIndex].numberOfMothers();
        for( int it = 0; it < nMothers; it++ ){
            auto parent = fgenparts[genIndex].mother(it);
            //foundpart = llpGenSearch( parent, searchtype );
            auto parentPdgId = std::abs(parent->pdgId());
            if( ( parentPdgId > 1000021 && parentPdgId < 1000038 ) && ( parentPdgId < searchtype ) ){ searchtype = parentPdgId; foundpart = parent; }
            //if( foundpart ) break;
        }//<<>>for( it = 0; it < nMothers; it++ )

		if( foundpart ){
		auto ninoeta = foundpart->eta();
        auto ninophi = foundpart->phi();
		int nFgen = fgenparts.size();
		for( int it = 0; it < nFgen; it++  ){
			
			auto fgeneta = fgenparts[it].eta();
            auto fgenphi = fgenparts[it].phi();
			//if( fgenp4 == ninop4 ){ result = it; break; }
			if( fgeneta == ninoeta && fgenphi == ninophi ){ 
	//std::cout << " getGnSigPhoInf2225: " << genIndex << " to " << it << " llpid: " << fgenpartllp[it] << " pdgid: "  << fgenparts[it].pdgId() << std::endl;
				result = it;
				//break;
			}//>><<if( fgeneta == ninoeta && fgenphi == ninophi )
		}//<<>>for( auto fgen : fgenparts )
        }//<<>>if( foundpart )

		if( result == -1 ) result = -10;

	}//<<>>if( gensusid == 22 )
    else if( gensusid == 32 || gensusid == 35 ){

        const reco::Candidate* foundpart = nullptr;
        int searchtype = 36;
        int nMothers = fgenparts[genIndex].numberOfMothers();
        for( int it = 0; it < nMothers; it++ ){
            auto parent = fgenparts[genIndex].mother(it);
			//std::cout << " searching 32/35 with " << parent << " target " << searchtype << " for " << it << " of " << nMothers << std::endl;
            foundpart = llpGenSearch( parent, searchtype );
			//std::cout << " searching 32/35 Done : " << foundpart << std::endl;
            if( foundpart ) break;
        }//<<>>for( it = 0; it < nMothers; it++ )

        //auto nino = parent;
        if( foundpart ){
        auto ninoeta = foundpart->eta();
        auto ninophi = foundpart->phi();
        int nFgen = fgenparts.size();
		//std::cout << " Found Part: eta " << ninoeta << " phi " << ninophi << " pdgID " << foundpart->pdgId() << std::endl;
        for( int it = 0; it < nFgen; it++  ){

            auto fgeneta = fgenparts[it].eta();
            auto fgenphi = fgenparts[it].phi();
            //if( fgenp4 == ninop4 ){ result = it; break; }
            if( fgeneta == ninoeta && fgenphi == ninophi ){
   // std::cout << " getGnSigPhoIfo3235: " << genIndex << " to " << it << " llpid: " << fgenpartllp[it] << " pdgid: "  << fgenparts[it].pdgId() << std::endl;
                result = it;
                //break;
            }//>><<if( fgeneta == ninoeta && fgenphi == ninophi )

        }//<<>>for( auto fgen : fgenparts )
		}//<<>>if( foundpart )
        if( result == -1 ) result = -10;

	}//<<>>else if( gensusid == 32 )
	
	//std::cout << " -- getGenSigPhoInfo Done !! " << std::endl;
	return result;

}//<<>>std::vector<float> KUCMSGenObject::getGenSigPhoInfo( uInt genIndex )

int KUCMSGenObject::getGenSigEleInfo( uInt genIndex ){

    int result( -1 );
	if( genIndex > fgenpartllp.size() ){ std::cout << " --- INDEX Blown !!!!  getGenSigEleInfo " << std::endl; return result; }
    auto gensusid = fgenpartllp[genIndex];
    if( GenDEBUG ) std::cout << " -- getGenSigEleInfo: " << genIndex << " gensusid: " << gensusid << std::endl;
    if( gensusid == 33 || gensusid == 34 ){

		const reco::Candidate* foundpart = nullptr;
		int searchtype = ( gensusid == 33 ) ? 23 : 24;
		int nMothers = fgenparts[genIndex].numberOfMothers();
		for( int it = 0; it < nMothers; it++ ){
			auto parent = fgenparts[genIndex].mother(it);
			foundpart = llpGenSearch( parent, searchtype );
			if( foundpart ) break;
		}//<<>>for( it = 0; it < nMothers; it++ )

        if( foundpart ){
        auto ninoeta = foundpart->eta();
        auto ninophi = foundpart->phi();
        //auto ninoeta = parent->eta();
        //auto ninophi = parent->phi();
        int nFgen = fgenparts.size();
        for( int it = 0; it < nFgen; it++  ){

            auto fgeneta = fgenparts[it].eta();
            auto fgenphi = fgenparts[it].phi();
            //if( fgenp4 == ninop4 ){ result = it; break; }
            if( fgeneta == ninoeta && fgenphi == ninophi ){
    //std::cout << " getGnSigEleInf2324: " << genIndex << " to " << it << " llpid: " << fgenpartllp[it] << " pdgid: "  << fgenparts[it].pdgId() << std::endl;
                result = it;
                //break;
            }//>><<if( fgeneta == ninoeta && fgenphi == ninophi )

        }//<<>>for( auto fgen : fgenparts )
        }//<<>>if( foundpart )
        if( result == -1 ) result = -10;

    }//<<>>if( gensusid == 33 || gensusid == 34 )
	else if( gensusid == 23 || gensusid == 24 ){

/*
        const reco::Candidate* foundpart = nullptr;
        int searchtype = 36;
        int nMothers = fgenparts[genIndex].numberOfMothers();
        for( int it = 0; it < nMothers; it++ ){
            auto parent = fgenparts[genIndex].mother(it);
            foundpart = llpGenSearch( parent, searchtype );
            if( foundpart ) break;
        }//<<>>for( it = 0; it < nMothers; it++ )
*/
        const reco::Candidate* foundpart = nullptr;
        int searchtype = 1000100;
        int nMothers = fgenparts[genIndex].numberOfMothers();
        for( int it = 0; it < nMothers; it++ ){
            auto parent = fgenparts[genIndex].mother(it);
            //foundpart = llpGenSearch( parent, searchtype );
            auto parentPdgId = std::abs(parent->pdgId());
            if( ( parentPdgId > 1000021 && parentPdgId < 1000038 ) && ( parentPdgId < searchtype ) ){ searchtype = parentPdgId; foundpart = parent; }
            //if( foundpart ) break;
        }//<<>>for( it = 0; it < nMothers; it++ )

        if( foundpart ){
        auto ninoeta = foundpart->eta();
        auto ninophi = foundpart->phi();
        int nFgen = fgenparts.size();
        for( int it = 0; it < nFgen; it++  ){

            auto fgeneta = fgenparts[it].eta();
            auto fgenphi = fgenparts[it].phi();
            //if( fgenp4 == ninop4 ){ result = it; break; }
            if( fgeneta == ninoeta && fgenphi == ninophi ){
    //std::cout << " getGnSigEleInf3334: " << genIndex << " to " << it << " llpid: " << fgenpartllp[it] << " pdgid: "  << fgenparts[it].pdgId() << std::endl;
                result = it;
                //break;
            }//>><<if( fgeneta == ninoeta && fgenphi == ninophi )

        }//<<>>for( auto fgen : fgenparts )
        }//<<>>if( foundpart )
        if( result == -1 ) result = -10;

	}//<<>>else if( gensusid == 33 || gensusid == 34 )

    return result;

}//<<>>std::vector<float> KUCMSGenObject::getGenSigEleInfo( uInt genIndex )
>>>>>>> 2f0e277fb0524b1d969b7ed35559585a2c17f0ba

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

int KUCMSGenObject::llpGenChaseP( const reco::GenParticle & kid, int depth ){

    auto kidPdgID = std::abs(kid.pdgId());
    //if( GenDEBUG ) std::cout << " --- llpgen chase genpart: " << kidPdgID << " depth: " << depth << std::endl;
    bool wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? true : false;
    bool phFlag = ( kidPdgID == 22 ) ? true : false;
    //bool eleFlag = ( kidPdgID == 11 ) ? true : false;
    bool qFlag = ( kidPdgID < 7 ) ? true : false;
    bool pFlag = ( kidPdgID == 2212 ) ? true : false;
	bool n1Flag = ( kidPdgID > 1000021 && kidPdgID < 1000026 ) ? true : false;
    bool n2Flag = ( kidPdgID == 1000035 || kidPdgID == 1000037 ) ? true : false;

    if( pFlag ) return 98;
    int nMoms = kid.numberOfMothers();
    if( nMoms == 0 ) return 99;
    int genPartID = 100;
    for( int gmit(0); gmit < nMoms; gmit++ ){

        auto genmom = kid.mother(gmit);
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

        bool pro = ( momPdgID == 2212 ) ? true : false;
        bool gli = ( momPdgID == 1000021 ) ? true : false;
		bool nti = ni1 || ni2 || ni3 || ni4;
        bool chi = ci1 || ci2;
		bool sqk = sq1 || sq2;

        if( phFlag && ni1 ){ if( depth == 0 ) return 22; else return 32; }
		else if( phFlag && nti ){ if( depth == 0 ) return 25; else return 35; }
		else if( wzFlag && nti ){ if( depth == 0 ) return 23; else return 33; }
        else if( wzFlag && chi ){ if( depth == 0 ) return 24; else return 34; }
        else if( qFlag && gli ){ if( depth == 0 ) return 21; else return 31; }
		else if( qFlag && sqk ){ if( depth == 0 ) return 20; else return 30; }
        //else if( chi || nti ) return 36;
        else if( ( n1Flag || n2Flag ) && ( sqk || gli || pro || qrk ) ) return 36;
        else if( gli || sqk ) return 37;
        //else if( nti || chi ) return 38;
        else if( pro ) return 97;

        auto mGenPartID = llpGenChaseP(kid.mother(gmit),depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    //if( GenDEBUG ) std::cout << " --- Found llpgen chase genpart: " << kidPdgID << " susid: " << genPartID << std::endl;
    return genPartID;

}//<<>> void llpChase( Candidate* kid )

const reco::Candidate* KUCMSGenObject::llpGenSearch( const reco::Candidate* kid, int target ){

	int depth = 0;
	const reco::Candidate* notfound = nullptr;
    const reco::Candidate* result = nullptr;
    auto kidPdgID = std::abs(kid->pdgId());
    //if( GenDEBUG ) std::cout << " --- llpgen chase genpart: " << kidPdgID << " depth: " << depth << std::endl;
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

/*
        if( phFlag && ni1 ){ if( depth == 0 ) return 22; else return 32; }
        else if( phFlag && nti ){ if( depth == 0 ) return 25; else return 35; }
        else if( wzFlag && ( ni1 || nti ) ) return 23;
        else if( wzFlag && chi ) return 24;
        else if( qFlag && gli ){ if( depth == 0 ) return 21; else return 31; }
        else if( qFlag && sqk ){ if( depth == 0 ) return 20; else return 30; }
        else if( nti || chi ) return 26;
        else if( gli || sqk ) return 27;
*/

int KUCMSGenObject::llpGenChaseP( const reco::Candidate* kid, int depth ){

    auto kidPdgID = std::abs(kid->pdgId());
	//if( GenDEBUG ) std::cout << " --- llpgen chase canadate: " << kidPdgID << " depth: " << depth << std::endl;
    bool wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? true : false;
    bool phFlag = ( kidPdgID == 22 ) ? true : false;
    //bool eleFlag = ( kidPdgID == 11 ) ? true : false;
    bool qFlag = ( kidPdgID < 7 ) ? true : false;
    bool pFlag = ( kidPdgID == 2212 ) ? true : false;
    bool n1Flag = ( kidPdgID > 1000021 && kidPdgID < 1000026 ) ? true : false;
    bool n2Flag = ( kidPdgID == 1000035 || kidPdgID == 1000037 ) ? true : false;

    if( pFlag ) return 98;
    int nMoms = kid->numberOfMothers();
    if( nMoms == 0 ) return 99;
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

        bool pro = ( momPdgID == 2212 ) ? true : false;
        bool gli = ( momPdgID == 1000021 ) ? true : false;
        bool nti = ni1 || ni2 || ni3 || ni4;
        bool chi = ci1 || ci2;
        bool sqk = sq1 || sq2;

        if( phFlag && ni1 ){ if( depth == 0 ) return 22; else return 32; }
        else if( phFlag && nti ){ if( depth == 0 ) return 25; else return 35; }
        else if( wzFlag && nti ){ if( depth == 0 ) return 23; else return 33; }
        else if( wzFlag && chi ){ if( depth == 0 ) return 24; else return 34; }
        else if( qFlag && gli ){ if( depth == 0 ) return 21; else return 31; }
        else if( qFlag && sqk ){ if( depth == 0 ) return 20; else return 30; }
        //else if( chi || nti ) return 36;
        else if( ( n1Flag || n2Flag ) && ( sqk || gli || pro || qrk ) ) return 36;
        else if( gli || sqk ) return 37;
        //else if( nti || chi ) return 38;
        else if( pro ) return 97;

        auto mGenPartID = llpGenChaseP(kid->mother(gmit),depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    //if( GenDEBUG ) std::cout << " --- Found llpgen chase canadate: " << kidPdgID << " susid: " << genPartID << std::endl;
    return genPartID;

}//<<>> void llpChase( const reco::Candidate* kid, int wzFlag ){)

void KUCMSGenObject::motherChase( const reco::Candidate* kid, string & depth ){

    depth += "-";
    for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
        std::cout <<  depth  << " gMomID : " << kid->mother(gmit)->pdgId() << " pt : " << kid->mother(gmit)->pt();
        std::cout << " Vertix (" << kid->mother(gmit)->vx() << "," << kid->mother(gmit)->vy() << "," << kid->mother(gmit)->vz() << ")";
        std::cout << " nGMothers " << kid->mother(gmit)->numberOfMothers() << std::endl;
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    if( not depth.empty() ) depth = depth.substr (0,depth.length()-1);

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
    float genJetLLP(-1.0);

	//std::cout << " --- Jet-GenJet matching " << std::endl;

    float goodDr(0.3);
    float goodRe(-1.0);
    int matchedIdx(-1);
    int index(0);
    for(const auto &genJet : fgenjets ){

        auto gjeta = genJet.eta();
        auto gjphi = genJet.phi();
        auto gje = genJet.energy();
		if( gje < 5.0 ) continue;
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
    float gDr(0.3);
    float gRe(-1.0);
    int mIdx(-1);
    int idx(0);
    for(const auto &genJet : fgenparts ){

		if( genJet.pdgId() > 9 ){ idx++; continue; }
        auto gjeta = genJet.eta();
        auto gjphi = genJet.phi();
        auto gje = genJet.energy();
        if( gje < 5.0 ) continue;
        auto jgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jetEta, jetPhi ));
        auto jtgjre = jetEnergy/gje;
        if( jgjdr < gDr ){
            gDr = jgjdr;
            gRe = jtgjre;
            mIdx = idx;
        }//<<>>if( jtgjdr <= goodDr )
        idx++;

    }//<<>>for(const auto& genJet : fgenjets ) 

    if( mIdx >= 0 ){ genJetLLP = fgenpartllp[mIdx]; }
    else{ genJetLLP = 35; }

	//std::cout << " --- Jet-GenJet genjet info " << std::endl;
    if( matchedIdx >= 0 ){

        auto genJet = fgenjets[matchedIdx];
        //if( GenDEBUG ) std::cout << " --- Jet-GenJet dR match : " << goodDr << std::endl;

        //auto nSources = genJet.numberOfSourceCandidatePtrs();
//        auto kids = genJet.daughterPtrVector();
//        auto theta = 2*std::atan(std::exp(-1*genJet.eta()));
//        auto cx = 120*std::sin(genJet.phi());
//        auto cy = 120*std::cos(genJet.phi());
//        auto cz = 120/std::tan(theta);
//        auto tofcor = hypo( cx, cy, cz )/SOL;
        //if( GenDEBUG ) kidChase( kids, vtxX, vtxY, vtxZ );
//        auto genKidInfo = kidTOFChain( kids, cx, cy, cz );
        //if( GenDEBUG ) std::cout << " - genJet GenTime noTOF : " << genKidInfo[0] << " rhPos: " << cx;
        //if( GenDEBUG ) std::cout << "," << cy << "," << cz << std::endl;
        genEta = genJet.eta();
        genPhi = genJet.phi();
//        if( genKidInfo[0] > 25.0 ) genTime = -28.0;
//        else if( genKidInfo[0] > -25.0 ) genTime = genKidInfo[0]-tofcor;
//        else genTime = -27.0;
//        genImpactAngle = genKidInfo[1];
        //if( GenDEBUG ) std::cout << " - genJet GenTime : " << genTime << " Angle: " << genImpactAngle << std::endl;
        genPt = genJet.pt();
        genEnergy = genJet.energy();
        genEMFrac = (genJet.chargedEmEnergy() + genJet.neutralEmEnergy())/genEnergy;
        genDrMatch = goodDr; //std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), genJet.eta(), genJet.phi()));
        genReMatch = goodRe;
        //genTimeVar = genKidInfo[2];
        //genNextBX = genKidInfo[3];
//        genTimeLLP = genKidInfo[4];
        //genLLPPurity = genKidInfo[5];
        //genNKids = genKidInfo[6];
        //genTOF = tofcor;
        //if( GenDEBUG ) std::cout << " -- Energy : " << genEnergy << " Pt : " << genPt << " EMfrac : " << genEMFrac << std::endl;

    }//<<>>if( matchedIdx >= 0 )i

	std::vector<float> results;

    results.push_back(genImpactAngle);
    results.push_back(genTime);
    results.push_back(genPt);
    results.push_back(genEta);
    results.push_back(genPhi);
    results.push_back(genEnergy);
    results.push_back(genEMFrac);
    results.push_back(genDrMatch);
    results.push_back(genReMatch);
    results.push_back(genTimeLLP);
    results.push_back(genTOF);

    results.push_back(goodDr);
    results.push_back(goodRe);	

    results.push_back(genJetLLP);
    results.push_back(gDr);
    results.push_back(gRe);

    return results;

}//<<>>std::vector<float> KUCMSGenObject::getGenJetInfo( float jetEta, float jetPhi, float jetPt )

// Gen Electron methods

std::vector<int>KUCMSGenObject:: MomIDs(const reco::GenParticle &genElectron) const {

  auto mother = genElectron.mother();

  std::vector<int> motherIDs;
  while(mother->pt() > 0) {
    const int motherID = mother->pdgId();

    mother = mother->mother();

    if(motherID == mother->pdgId())
      continue;

    motherIDs.push_back(motherID);
  }

  return motherIDs;

}

bool KUCMSGenObject::IsMotherZ(const reco::GenParticle &genElectron) const {
  LepType momType = ClassifyGenElectron(genElectron);
  return (momType == kZ);
}

bool KUCMSGenObject::isSignalGenElectron(const reco::GenParticle &genElectron) const {

  LepType momType = ClassifyGenElectron(genElectron);
  return (momType == kZ || momType == kSusy);

}
/*
LepType KUCMSGenObject::ClassifyGenElectron(const std::vector<int> &motherIDs) const {

  LepType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignLeptonMomType(id);

    if(momType != kUnmatched)
      break;
  }
  return momType;
}
*/
LepType KUCMSGenObject::ClassifyGenElectron(const reco::GenParticle &genElectron) const {

  std::vector<int> motherIDs(MomIDs(genElectron));

  LepType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignLeptonMomType(id);

    if(momType != kUnmatched)
      break;
  }
  return momType;
}

void KUCMSGenObject::PrintMother(const LepType &momType) const {

  std::cout << "mother: ";
  if (momType == kW) std::cout << "W boson" << std::endl;
  else if (momType == kZ) std::cout << "Z boson" << std::endl;
  else if (momType == kTau)  std::cout << "tau lepton" << std::endl;
  else if (momType == kConversion) std::cout << "photon (conversion)" << std::endl;
  else if (momType == kLight) std::cout << "light quark" << std::endl;
  else if (momType == kHeavy) std::cout << "heavy quark" << std::endl;
  else if (momType == kSusy) std::cout << "prompt from Susy particle" << std::endl;
  else if (momType == kUnmatched) std::cout << "unmatched (this shouldn't happen!)" << std::endl;

}

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

    LepType momType = ClassifyGenElectron(genElectron);

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

LepType KUCMSGenObject::AssignLeptonMomType(const int motherID) const {
  /*
 LepType type = kUnmatched;
  int absMotherId = abs(motherID);
  int motherIdMod1000 = absMotherId % 1000;
 
  switch (true) {
  case absMotherId == 24:
    type = kW;
    break;
  case motherID == 23:
    type = kZ;
    break;
  case absMotherId == 15:
    type = kTau;
    break;
  case (motherIdMod1000 > 100 && motherIdMod1000 < 400) ||
    (motherIdMod1000 > 1000 && motherIdMod1000 < 4000) ||
    (absMotherId > 0 && absMotherId < 4) ||
    motherID == 21:
    type = kLight;
    break;
  case (motherIdMod1000 > 400 && motherIdMod1000 < 600) ||
    (motherIdMod1000 > 4000 && motherIdMod1000 < 6000) ||
    (absMotherId > 3 && absMotherId < 7):
    type = kHeavy;
    break;
  case motherID == 22 || absMotherId == 11:
    type = kConversion;
    break;
  case motherID == 1000022:
    type = kSusy;
    break;
  case absMotherId > 1000000:
    type = kOther;
    break;
  }

  return type;
  }*/
  
  LepType type = kUnmatched;

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
  else if(motherID == 22)// || abs(motherID) == 11)
    type = kConversion;
  else if(motherID == 1000022)
    type = kSusy;
  else if(abs(motherID) > 1000000)
    type = kOther;

  return type;
  
}
 /*
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
