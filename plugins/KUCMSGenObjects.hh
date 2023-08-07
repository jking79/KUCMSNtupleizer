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

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSGenObjectHeader
#define KUCMSGenObjectHeader

//#define GenEBUG true
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
    void PostProcessEvent( ItemManager<float>& geVar ){};

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
    int llpGenChaseP( const reco::GenParticle kid, int depth );
    void kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz );
    std::vector<float> kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  );
    std::vector<float> getGenPartMatch( const reco::SuperCluster* scptr, std::vector<reco::GenParticle> fgenparts, float pt );
    std::vector<float> getGenPartMatch( const reco::SuperCluster* scptr, float pt );
    std::vector<float> getGenJetInfo( float jetEta, float jetPhi, float jetPt );

    private:

    std::vector<reco::GenParticle> fgenparts;
    std::vector<int> fgenpartllp;
    std::vector<reco::GenJet> fgenjets;
    std::vector<reco::GenParticle> fgjquarks;
    std::vector<int> fgjquarksllp;

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
    Branches.makeBranch("genPdgId","Gen_pdgId",VUINT);
    Branches.makeBranch("genWgt","Evt_genWgt",FLOAT);
    Branches.makeBranch("genSusId","Gen_susId",VINT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSGen::InitObject( TTree* fOutTree )

void KUCMSGenObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // ECAL RECHITS

    iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);
    iEvent.getByToken(gent0Token_, gent0_);
    iEvent.getByToken(genxyz0Token_, genxyz0_);
    iEvent.getByToken(pileupInfosToken_, pileupInfos_);
    iEvent.getByToken(genParticlesToken_, genParticles_);
    iEvent.getByToken(genJetsToken_, genJets_);

    if( GenDEBUG ) std::cout << "Collecting Gens" << std::endl;

    fgenparts.clear();
    fgenpartllp.clear();
    fgenjets.clear();
    fgjquarks.clear();
    fgjquarksllp.clear();
    int nGenPart(0);

    if( GenDEBUG ) std::cout << "Collecting Gen Particles" << std::endl;
    for( const auto &genPart : *genParticles_ ){

        if( genPart.energy() < cfPrm("minGenE") ) continue;
        auto genPartID = llpGenChaseP(genPart,0);
        nGenPart++;
        if( genPart.status() == 1 ){
            fgenparts.push_back(genPart);
            fgenpartllp.push_back(genPartID);
        }//<<>>if( genPart.isLastCopy() )
		else if( genPart.pdgId() < 9 ){
			fgjquarks.push_back(genPart);
            fgjquarksllp.push_back(genPartID);
		}//<<>>else if( genpart.pdgId() < 9 && genPartID < 22 )

    }//<<>>for(const auto& genPart : *genParticles_ )

    if( GenDEBUG ) std::cout << "Collecting GenJets" << std::endl;
    for(const auto &genjet : *genJets_ ){ fgenjets.push_back(genjet); }

}//<<>>void KUCMSGen::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSGenObject::ProcessEvent( ItemManager<float>& geVar ){

    if( GenDEBUG ) std::cout << "Processing Gens" << std::endl;

    Branches.clearBranches();

    if( GenDEBUG ) std::cout << " - enetering Gen loop" << std::endl;

    for (const auto & genpart : fgenparts ){

        const float genPt = genpart.pt();
        const float genEnergy = genpart.energy();
        const float genPhi = genpart.phi();
        const float genEta = genpart.eta();
        const uInt genPdgId = genpart.pdgId();
        const float genPx = genpart.px();
        const float genPy = genpart.py();
        const float genPz = genpart.pz();
		const int genSusId = llpGenChaseP( genpart, 0 );
		if( GenDEBUG ) std::cout << "GenPart : genSusId = " << genSusId << std::endl;;

        Branches.fillBranch("genPt",genPt);
        Branches.fillBranch("genEnergy",genEnergy);
        Branches.fillBranch("genPhi",genPhi);
        Branches.fillBranch("genEta",genEta);
        Branches.fillBranch("genPdgId",genPdgId);
        Branches.fillBranch("genPx",genPx);
        Branches.fillBranch("genPy",genPy);
        Branches.fillBranch("genPz",genPz);
        Branches.fillBranch("genSusId",genSusId);

    }//<<>> for (const auto genpart : fgenparts )

    float wgt = genEvtInfo_->weight();
    Branches.fillBranch("genWgt",wgt);

}//<<>>void KUCMSGen::ProcessEvent()

void KUCMSGenObject::EndJobs(){}

std::vector<float> KUCMSGenObject::getGenPartMatch( const reco::SuperCluster* scptr, float pt ){

    // GenParticle Info   -------------------------------------------------------------------
    if( GenDEBUG ) std::cout << "Getting phoGenParton Information" << std::endl;
    // set defaults for no match here
    vector<float> results;
    float minDr(0.3);
    float minSDr(0.2);
    float minDp(0.3);
    float minSDp(0.2);
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
        const auto gpt = genPart.pt();
        auto cphoEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
        auto cphoPhi = std::atan2(rhY-gnY,rhX-gnX);
        auto dr = std::sqrt(reco::deltaR2(genPart.eta(), genPart.phi(), cphoEta, cphoPhi ));
        auto dp = std::abs(gpt-pt)/gpt;
        if( dr < minDr && dp <= minDp ){
            minDr = dr;
            minDp = dp;
            matchedIdx = index;
			llpId = fgenpartllp[index];
        } else if( dr < minSDr && dp <= minSDp ){
            minSDr = dr;
            minSDp = dp;
            matchedSIdx = index;
            llpSId = fgenpartllp[index];
        }//<<>>if( dr < minDr && dp < minDp )
        index++;

    }//<<>>for(const auto& genPart : fgenparts  ) 

    float matDr(-1.0), matDp(-1.0), matSDr(-1.0), matSDp(-1.0);
    if( matchedIdx >=0 ){ matDr = minDr; matDp = minDp; matSDr = minSDr; matSDp = minDp; }

    results.push_back( matchedIdx );
    results.push_back( matDr );
    results.push_back( matDp );
    results.push_back( matchedSIdx );
    results.push_back( matSDr );
    results.push_back( matSDp );
    results.push_back( llpId );
    results.push_back( llpSId );
    return results;

}//<<>>getGenPartMatch( reco::SuperClusterCollection *scptr, std::vector<reco::GenParticle> fgenparts )


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

int KUCMSGenObject::llpGenChaseP( const reco::GenParticle kid, int depth ){

    auto kidPdgID = std::abs(kid.pdgId());
    bool wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? true : false;
    bool gmFlag = ( kidPdgID == 22 ) ? true : false;
    bool qFlag = ( kidPdgID < 7 ) ? true : false;
    bool pFlag = ( kidPdgID == 2212 ) ? true : false;

    if( pFlag ) return 29;
    int nMoms = kid.numberOfMothers();
    if( ! nMoms ) return 39;
    int genPartID = 40;
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
        bool gli = ( momPdgID == 1000021 ) ? true : false;
		bool nti = ni2 || ni3 || ni4;
        bool chi = ci1 || ci2;
		bool sqk = sq1 || sq2;

        if( gmFlag && ni1 ){ if( depth == 0 ) return 22; else return 32; }
		else if( gmFlag && nti ){ if( depth == 0 ) return 25; else return 35; }
		else if( wzFlag && ( ni1 || nti ) ) return 23;
        else if( wzFlag && chi ) return 24;
        else if( qFlag && gli ){ if( depth == 0 ) return 21; else return 31; }
		else if( qFlag && sqk ){ if( depth == 0 ) return 20; else return 30; }

        auto mGenPartID = llpGenChaseP(kid.mother(gmit),depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return genPartID;

}//<<>> void llpChase( Candidate* kid )

int KUCMSGenObject::llpGenChaseP( const reco::Candidate* kid, int depth ){

    auto kidPdgID = std::abs(kid->pdgId());
    bool wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? true : false;
    bool gmFlag = ( kidPdgID == 22 ) ? true : false;
    bool qFlag = ( kidPdgID < 7 ) ? true : false;
    bool pFlag = ( kidPdgID == 2212 ) ? true : false;

    if( pFlag ) return 29;
    int nMoms = kid->numberOfMothers();
    if( ! nMoms ) return 39;
    int genPartID = 40;
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
        bool gli = ( momPdgID == 1000021 ) ? true : false;
        bool nti = ni1 || ni2 || ni3 || ni4;
        bool chi = ci1 || ci2;
        bool sqk = sq1 || sq2;

        if( gmFlag && ni1 ){ if( depth == 0 ) return 22; else return 32; }
        else if( gmFlag && nti ){ if( depth == 0 ) return 25; else return 35; }
        else if( wzFlag && ( ni1 || nti ) ) return 23;
        else if( wzFlag && chi ) return 24;
        else if( qFlag && gli ){ if( depth == 0 ) return 21; else return 31; }
        else if( qFlag && sqk ){ if( depth == 0 ) return 20; else return 30; }

        auto mGenPartID = llpGenChaseP(kid->mother(gmit),depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
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

std::vector<float> KUCMSGenObject::getGenJetInfo( float jetEta, float jetPhi, float jetPt ){

    float genImpactAngle(-10.0);
    float genTime(-99.9);
    float genPt(-10.0);
    float genEta(-10.0);
    float genPhi(-10.0);
    float genEnergy(-10.0);
    float genEMFrac(-1.0);
    float genDrMatch(-1.0);
    float genDptMatch(-1.0);
    float genTimeLLP(-99.9);
    float genTOF(-99.9);
    float genJetLLP(-1.0);

    float goodDr(5.0);
    float goodDpt(5.0);
    int matchedIdx(-1);
    int index(0);
    for(const auto &genJet : fgenjets ){

        auto gjeta = genJet.eta();
        auto gjphi = genJet.phi();
        auto gjpt = genJet.pt();
        auto jtgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jetEta, jetPhi ));
        auto jtgjdp = std::abs(jetPt-gjpt)/jetPt;
        if( jtgjdr < goodDr && jtgjdp <= goodDpt ){
            goodDr = jtgjdr;
            goodDpt = jtgjdp;
            matchedIdx = index;
        }//<<>>if( jtgjdr <= goodDr )
        index++;

    }//<<>>for(const auto& genJet : fgenjets ) 

    float gDr(5.0);
    float gDpt(5.0);
    int mIdx(-1);
    int idx(0);
    for(const auto &genJet : fgjquarks ){

        auto gjeta = genJet.eta();
        auto gjphi = genJet.phi();
        auto gjpt = genJet.pt();
        auto jgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jetEta, jetPhi ));
        auto jgjdp = std::abs(jetPt-gjpt)/jetPt;
        if( jgjdr < gDr && jgjdp <= gDpt ){
            gDr = jgjdr;
            gDpt = jgjdp;
            mIdx = idx;
        }//<<>>if( jtgjdr <= goodDr )
        idx++;

    }//<<>>for(const auto& genJet : fgenjets ) 

    if( mIdx >= 0 ){ genJetLLP = fgjquarksllp[mIdx]; }
    else{ genJetLLP = 35; }

    if( matchedIdx >= 0 ){

        auto genJet = fgenjets[matchedIdx];
        if( GenDEBUG ) std::cout << " --- Jet-GenJet dR match : " << goodDr << std::endl;

        //auto nSources = genJet.numberOfSourceCandidatePtrs();
        auto kids = genJet.daughterPtrVector();
        auto theta = 2*std::atan(std::exp(-1*genJet.eta()));
        auto cx = 120*std::sin(genJet.phi());
        auto cy = 120*std::cos(genJet.phi());
        auto cz = 120/std::tan(theta);
        auto tofcor = hypo( cx, cy, cz )/SOL;
        //if( GenDEBUG ) kidChase( kids, vtxX, vtxY, vtxZ );
        auto genKidInfo = kidTOFChain( kids, cx, cy, cz );
        if( GenDEBUG ) std::cout << " - genJet GenTime noTOF : " << genKidInfo[0] << " rhPos: " << cx;
        if( GenDEBUG ) std::cout << "," << cy << "," << cz << std::endl;
        genEta = genJet.eta();
        genPhi = genJet.phi();
        if( genKidInfo[0] > 25.0 ) genTime = -28.0;
        else if( genKidInfo[0] > -25.0 ) genTime = genKidInfo[0]-tofcor;
        else genTime = -27.0;
        genImpactAngle = genKidInfo[1];
        if( GenDEBUG ) std::cout << " - genJet GenTime : " << genTime << " Angle: " << genImpactAngle << std::endl;
        genPt = genJet.pt();
        genEnergy = genJet.energy();
        genEMFrac = (genJet.chargedEmEnergy() + genJet.neutralEmEnergy())/genEnergy;
        genDrMatch = goodDr; //std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), genJet.eta(), genJet.phi()));
        genDptMatch = goodDpt;
        //genTimeVar = genKidInfo[2];
        //genNextBX = genKidInfo[3];
        genTimeLLP = genKidInfo[4];
        //genLLPPurity = genKidInfo[5];
        //genNKids = genKidInfo[6];
        genTOF = tofcor;
        if( GenDEBUG ) std::cout << " -- Energy : " << genEnergy << " Pt : " << genPt << " EMfrac : " << genEMFrac << std::endl;

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
    results.push_back(genDptMatch);
    results.push_back(genTimeLLP);
    results.push_back(genTOF);

    results.push_back(goodDr);
    results.push_back(goodDpt);	

    results.push_back(genJetLLP);
    results.push_back(gDr);
    results.push_back(gDpt);

	return results;

}//<<>>std::vector<float> KUCMSGenObject::getGenJetInfo( float jetEta, float jetPhi, float jetPt )

#endif
