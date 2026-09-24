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

#include <cmath>
#include <limits>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

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
#include "KUCMSDisplacedVertex_mini.hh"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSGenObjectMiniHeader
#define KUCMSGenObjectMiniHeader

//#define GenDEBUG true
#define GenDEBUG false

using namespace edm; 

class KUCMSGenObjectMini : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSGenObjectMini( const edm::ParameterSet& iConfig );
    ~KUCMSGenObjectMini(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadGenEvtInfoTokens( edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken ){ genEvtInfoToken_ = genEvtInfoToken;};
    void LoadGenT0Tokens( edm::EDGetTokenT<float> gent0Token ){ gent0Token_ = gent0Token; };
    void LoadGenXYZ0Tokens( edm::EDGetTokenT<Point3D> genxyz0Token ){ genxyz0Token_ = genxyz0Token; };
    void LoadGenPUInfosTokens( edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfosToken ){ pileupInfosToken_ = pileupInfosToken; };
    void LoadgpTokens( edm::EDGetTokenT<std::vector<reco::GenParticle> > genPartPToken ){ genParticlesPToken_ = genPartPToken; };
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
    //void motherChase( const reco::Candidate* kid, string & depth );

    void motherChase( const reco::Candidate* kid, string depth );
    void motherChase( const reco::GenParticle & kid, string  depth );
    int typeChase( const reco::Candidate* kid );
    int typeChase( const reco::GenParticle & kid );
	int getEvtType( int genPdgId );

    int llpGenChase( const reco::Candidate* kid );
    int llpGenChase( const reco::GenParticle & kid );
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
    //GlobalPoint GenVertex() {return GlobalPoint(genxyz0_->x(), genxyz0_->y(), genxyz0_->z()); };

    // new exclusive gen matching  v3fPoint
    std::vector<int> getGenMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& energies );
    std::vector<int> getGenPhoMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& energies );
    std::vector<int> getGenEleMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& energies, std::vector<reco::TransientTrack>& trakcs );
    std::vector<int> getGenMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& energies, int select, std::vector<reco::TransientTrack>& trakcs );
	std::vector<int> getGenMuonMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& reco_e, std::vector<reco::TransientTrack>& trakcs );
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

    std::vector<reco::GenParticle> fpgenparts;
    std::vector<int> fpgpLlp;
    std::vector<int> fpgpMomIdx;
    std::vector<int> fpgpGMomIdx;
    std::vector<int> fpgpKidIdx;

    std::vector<reco::GenParticle> tpgenparts;
    std::vector<int> tpgpLlp;
    std::vector<bool> tpgpGetMom;
    std::vector<bool> tpgpKeep;

	std::vector<int> sourceToFinal;
	std::vector<int> finalToSource;

    std::vector<reco::GenJet> fgenjets;
    std::vector<int> fgenjetllp;

    std::vector<int> nueEvntId;
    std::vector<int> chrEvntId;

    std::vector<float> pmatdr;
    std::vector<float> pmatde;

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

    const edm::InputTag genParticlesPTag;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesPToken_;
    edm::Handle<std::vector<reco::GenParticle>> genParticlesP_;

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

	float pvx, pvy, pvz;

	int sourceMotherIndex( int sourceIndex, unsigned int motherNumber = 0 ) const;
	int normalizedMotherSourceIndex( int sourceIndex ) const;
	void retainSourceAncestry( int sourceIndex, std::vector<bool>& keep, std::unordered_set<int>& activePath, int depth = 0 ) const;
	int llpGenChasePImpl( const reco::Candidate* kid, int childPdgId,
		std::unordered_set<const reco::Candidate*>& activePath, int depth ) const;


	public:

		// gen particle info accessors
    int getGenLlpId( int genIndex ) const { 
			if( genIndex < 0 || static_cast<std::size_t>(genIndex) >= fpgpLlp.size() ) return -1;
			return fpgpLlp[genIndex];
		};//<<>>int getGenLlpId( uInt genIndex ) const
    int getGenMomIdx( int genIndex ) const {
		if( genIndex < 0 || static_cast<std::size_t>(genIndex) >= fpgpMomIdx.size() ) return -1;
		return fpgpMomIdx[genIndex];
    };//<<>>int getGenLlpId( uInt genIndex ) const 
    int getGenPdgId( int genIndex ) const {
		if( genIndex < 0 || static_cast<std::size_t>(genIndex) >= fpgenparts.size() ) return 0;
		return fpgenparts[genIndex].pdgId();
    };//<<>>int getGenLlpId( uInt genIndex ) const 

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

};//<<>>class KUCMSg : public KUCMSObjectBase

KUCMSGenObjectMini::KUCMSGenObjectMini( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : false );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfFlag.set( "motherChase", iConfig.existsAs<bool>("doGenMotherChase") ? iConfig.getParameter<bool>("doGenMotherChase") : false );
    cfPrm.set( "minGenE", iConfig.existsAs<double>("minGenE") ? iConfig.getParameter<double>("minGenE") : 0.0 );

}//<<>>KUCMSg::KUCMSGen( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSGenObjectMini::InitObject( TTree* fOutTree ){

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
    Branches.makeBranch("genMomVx","Gen_MomVx",VFLOAT);
    Branches.makeBranch("genMomVy","Gen_MomVy",VFLOAT);
    Branches.makeBranch("genMomVz","Gen_MomVz",VFLOAT);
    Branches.makeBranch("genMomPx","Gen_MomPx",VFLOAT);
    Branches.makeBranch("genMomPy","Gen_MomPy",VFLOAT);
    Branches.makeBranch("genMomPz","Gen_MomPz",VFLOAT);
    Branches.makeBranch("genMomEnergy","Gen_MomEnergy",VFLOAT);
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
    Branches.makeBranch("isGG","Evt_isGG", BOOL);
    Branches.makeBranch("isZZ", "Evt_isZZ", BOOL);
    Branches.makeBranch("isGZ","Evt_isGZ", BOOL);
    Branches.makeBranch("fsType","Evt_fsType", INT);
    Branches.makeBranch("nXs","Evt_nXs", INT);

    Branches.makeBranch("Xa_MomDisplacment","Xa_Displacment",FLOAT);
    Branches.makeBranch("Xa_PdgId","Xa_pdgId",UINT);
    Branches.makeBranch("Xa_Vx","Xa_vx",FLOAT);
    Branches.makeBranch("Xa_Vy","Xa_vy",FLOAT);
    Branches.makeBranch("Xa_Vz","Xa_vz",FLOAT);
    Branches.makeBranch("Xa_Pt","Xa_pt",FLOAT);
    Branches.makeBranch("Xa_P","Xa_p",FLOAT);
    Branches.makeBranch("Xa_beta","Xa_beta",FLOAT);
    Branches.makeBranch("Xa_ctau","Xa_ctau",FLOAT);
    Branches.makeBranch("Xa_Mass","Xa_mass",FLOAT);
    Branches.makeBranch("Xa_Phi","Xa_phi",FLOAT);
    Branches.makeBranch("Xa_Eta","Xa_eta",FLOAT);
    Branches.makeBranch("Xa_Energy","Xa_energy",FLOAT);

    Branches.makeBranch("Xb_MomDisplacment","Xb_Displacment",FLOAT);
    Branches.makeBranch("Xb_PdgId","Xb_pdgId",UINT);
    Branches.makeBranch("Xb_Vx","Xb_vx",FLOAT);
    Branches.makeBranch("Xb_Vy","Xb_vy",FLOAT);
    Branches.makeBranch("Xb_Vz","Xb_vz",FLOAT);
    Branches.makeBranch("Xb_Pt","Xb_pt",FLOAT);
    Branches.makeBranch("Xb_P","Xb_p",FLOAT);
    Branches.makeBranch("Xb_beta","Xb_beta",FLOAT);
    Branches.makeBranch("Xb_ctau","Xb_ctau",FLOAT);
    Branches.makeBranch("Xb_Mass","Xb_mass",FLOAT);
    Branches.makeBranch("Xb_Phi","Xb_phi",FLOAT);
    Branches.makeBranch("Xb_Eta","Xb_eta",FLOAT);
    Branches.makeBranch("Xb_Energy","Xb_energy",FLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSGen::InitObject( TTree* fOutTree )

int KUCMSGenObjectMini::sourceMotherIndex( int sourceIndex, unsigned int motherNumber ) const {

    if( !genParticlesP_.isValid() ) return -1;
    if( sourceIndex < 0 || static_cast<std::size_t>(sourceIndex) >= genParticlesP_->size() ) return -1;

    const auto & particle = (*genParticlesP_)[sourceIndex];
    if( motherNumber >= particle.numberOfMothers() ) return -1;

    const auto motherRef = particle.motherRef(motherNumber);
    if( motherRef.isNull() ) return -1;
    if( motherRef.id() != genParticlesP_.id() ) return -1;

    const auto motherKey = motherRef.key();
    if( motherKey >= genParticlesP_->size() ) return -1;
    return static_cast<int>(motherKey);

}//<<>>int KUCMSGenObjectMini::sourceMotherIndex

int KUCMSGenObjectMini::normalizedMotherSourceIndex( int sourceIndex ) const {

    if( !genParticlesP_.isValid() ) return -1;
    if( sourceIndex < 0 || static_cast<std::size_t>(sourceIndex) >= genParticlesP_->size() ) return -1;

    std::unordered_set<int> visited;
    int currentIndex = sourceIndex;
    for( int depth = 0; depth < 100; depth++ ){
        if( !visited.insert(currentIndex).second ) return -1;

        const int motherIndex = sourceMotherIndex(currentIndex, 0);
        if( motherIndex < 0 ) return -1;

        const auto & current = (*genParticlesP_)[currentIndex];
        const auto & mother = (*genParticlesP_)[motherIndex];
        if( mother.pdgId() != current.pdgId() ) return motherIndex;
        currentIndex = motherIndex;
    }

    return -1;

}//<<>>int KUCMSGenObjectMini::normalizedMotherSourceIndex

void KUCMSGenObjectMini::retainSourceAncestry( int sourceIndex, std::vector<bool>& keep,
        std::unordered_set<int>& activePath, int depth ) const {

    if( !genParticlesP_.isValid() ) return;
    if( sourceIndex < 0 || static_cast<std::size_t>(sourceIndex) >= genParticlesP_->size() ) return;
    if( depth > 100 || !activePath.insert(sourceIndex).second ) return;

    const int motherIndex = normalizedMotherSourceIndex(sourceIndex);
    if( motherIndex >= 0 ){
        keep[motherIndex] = true;
        retainSourceAncestry(motherIndex, keep, activePath, depth+1);
    }

    activePath.erase(sourceIndex);

}//<<>>void KUCMSGenObjectMini::retainSourceAncestry

void KUCMSGenObjectMini::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);
    iEvent.getByToken(genParticlesPToken_, genParticlesP_);
    iEvent.getByToken(genJetsToken_, genJets_);

    fpgenparts.clear();
    fpgpLlp.clear();
    fpgpMomIdx.clear();
    fpgpGMomIdx.clear();
    fpgpKidIdx.clear();
    tpgenparts.clear();
    tpgpLlp.clear();
    tpgpGetMom.clear();
    tpgpKeep.clear();
    sourceToFinal.clear();
    finalToSource.clear();
    pmatdr.clear();
    pmatde.clear();
    fgenjets.clear();
    fgenjetllp.clear();
    nueEvntId.clear();
    chrEvntId.clear();

    const int nSource = static_cast<int>(genParticlesP_->size());
    tpgenparts.reserve(nSource);
    tpgpLlp.reserve(nSource);
    tpgpGetMom.assign(nSource, false);
    tpgpKeep.assign(nSource, false);

    for( int sourceIndex = 0; sourceIndex < nSource; sourceIndex++ ){
        const auto & genPart = (*genParticlesP_)[sourceIndex];
        const int absPdgId = std::abs(genPart.pdgId());
        const int genPartSusId = llpGenChaseP(genPart, 0);

        tpgenparts.push_back(genPart);
        tpgpLlp.push_back(genPartSusId);

        if( genPart.pt() < cfPrm("minGenE") ) continue;

        if( cfFlag("motherChase") ) motherChase(genPart, "");

        const bool goodPdgId = absPdgId < 100 || absPdgId == 1000022 || absPdgId == 1000023;
        if( genPart.status() == 1 && goodPdgId ){
            tpgpKeep[sourceIndex] = true;
            tpgpGetMom[sourceIndex] = true;
        } else {
            const bool isTree = genPart.status() < 40;
            const bool isQuark = absPdgId < 9;
            const bool isLepton = absPdgId > 10 && absPdgId < 19;
            const bool isPhoton = absPdgId == 22;
            const bool isWZ = absPdgId == 23 || absPdgId == 24;
            const bool isSusy = absPdgId > 1000000 && absPdgId < 3000000;
            const bool isNeutralino1 = absPdgId == 1000022;
            const bool isNeutralino2 = absPdgId == 1000023;
            const bool isChargino1 = absPdgId == 1000024;

            if( isQuark && isTree &&
                    ( genPartSusId == 41 || genPartSusId == 42 || genPartSusId == 43 ||
                      genPartSusId == 33 || genPartSusId == 34 ) ){
                tpgpKeep[sourceIndex] = true;
            }
            if( isWZ && isTree && ( genPartSusId == 23 || genPartSusId == 24 ) ){
                tpgpKeep[sourceIndex] = true;
                tpgpGetMom[sourceIndex] = true;
            }
            if( ( isSusy || isLepton || isPhoton ) && isTree ){
                tpgpKeep[sourceIndex] = true;
                tpgpGetMom[sourceIndex] = true;
            }

            if( isNeutralino1 || isNeutralino2 || isChargino1 ){
                const int motherIndex = sourceMotherIndex(sourceIndex, 0);
                const bool lastInChain = motherIndex < 0 ||
                    (*genParticlesP_)[motherIndex].pdgId() != genPart.pdgId();
                if( lastInChain ){
                    if( isNeutralino1 || isNeutralino2 ) nueEvntId.push_back(typeChase(genPart));
                    else chrEvntId.push_back(typeChase(genPart));
                }
            }
        }
    }

    for( int sourceIndex = 0; sourceIndex < nSource; sourceIndex++ ){
        if( !tpgpKeep[sourceIndex] || !tpgpGetMom[sourceIndex] ) continue;
        std::unordered_set<int> activePath;
        retainSourceAncestry(sourceIndex, tpgpKeep, activePath);
    }

    sourceToFinal.assign(nSource, -1);
    for( int sourceIndex = 0; sourceIndex < nSource; sourceIndex++ ){
        if( !tpgpKeep[sourceIndex] ) continue;

        const int finalIndex = static_cast<int>(fpgenparts.size());
        sourceToFinal[sourceIndex] = finalIndex;
        finalToSource.push_back(sourceIndex);
        fpgenparts.push_back(tpgenparts[sourceIndex]);
        fpgpLlp.push_back(tpgpLlp[sourceIndex]);
        fpgpMomIdx.push_back(-1);
        fpgpGMomIdx.push_back(-1);
        fpgpKidIdx.push_back(-1);
        pmatdr.push_back(-1.f);
        pmatde.push_back(1000.f);
    }

    for( int finalIndex = 0; finalIndex < static_cast<int>(finalToSource.size()); finalIndex++ ){
        const int sourceIndex = finalToSource[finalIndex];
        const int motherSourceIndex = normalizedMotherSourceIndex(sourceIndex);
        if( motherSourceIndex < 0 ) continue;

        const int motherFinalIndex = sourceToFinal[motherSourceIndex];
        if( motherFinalIndex < 0 ) continue;
        fpgpMomIdx[finalIndex] = motherFinalIndex;
        if( fpgpKidIdx[motherFinalIndex] == -1 ) fpgpKidIdx[motherFinalIndex] = finalIndex;

        const int grandmotherSourceIndex = normalizedMotherSourceIndex(motherSourceIndex);
        if( grandmotherSourceIndex >= 0 ){
            fpgpGMomIdx[finalIndex] = sourceToFinal[grandmotherSourceIndex];
        }
    }

    for( const auto & genJet : *genJets_ ){
        fgenjets.push_back(genJet);
        fgenjetllp.push_back(llpGenChaseP(genJet, 0));
    }



}//<<>>void KUCMSGen::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSGenObjectMini::PostProcessEvent( ItemManager<float>& geVar ){

    const int nGenPart = static_cast<int>(fpgenparts.size());
    if( static_cast<int>(pmatdr.size()) != nGenPart || static_cast<int>(pmatde.size()) != nGenPart ){
        std::cout << " !!!! -- gen match vector size mismatch !!!! " << std::endl;
        return;
    }

    for( int idx = 0; idx < nGenPart; idx++ ){
        if( pmatdr[idx] >= 0.f && pmatde[idx] < 1000.f ){
            Branches.fillBranch("genMatchDr",pmatdr[idx]);
            Branches.fillBranch("genMatchDe",pmatde[idx]);
        } else {
            Branches.fillBranch("genMatchDr",-1.f);
            Branches.fillBranch("genMatchDe",-10.f);
        }
    }

}//<<>>void KUCMSGenObjectMini::PostProcessEvent( ItemManager<float>& geVar )

void KUCMSGenObjectMini::ProcessEvent( ItemManager<float>& geVar ){

    if( GenDEBUG ) std::cout << "Processing Gens" << std::endl;

    Branches.clearBranches();

	bool isGG = false;
	bool isZZ = false;
	bool isGZ = false;
	int nXG = 0;
	int nXZ = 0;
	bool hasLWZX = false;
    bool hasLWZQ = false;
	int nLWZX = 0;
	int nXs = 0;

    int nGenParts = 0;
    const int nFgenParts = static_cast<int>(fpgenparts.size());

    for (const auto & genpart : fpgenparts ){

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
        const int genSusId = fpgpLlp[nGenParts];
        const int genMomIdx = fpgpMomIdx[nGenParts];
        const int genGMomIdx = fpgpGMomIdx[nGenParts];
        const int genCharge = genpart.charge();
        float genVx = genpart.vx();
        float genVy = genpart.vy();
        float genVz = genpart.vz();
        const float genMass = genpart.mass();

        if( genSusId == 22 or genSusId == 25 ) nXG++;
        if( genSusId == 23 or genSusId == 24 ) nXZ++;

        const int partMomIdx = genMomIdx;
        const int masterMomIndex = genMomIdx;
        const bool validMom = partMomIdx >= 0 && partMomIdx < nFgenParts;
        float momVx = validMom ? fpgenparts[partMomIdx].vx() : -999;
        float momVy = validMom ? fpgenparts[partMomIdx].vy() : -999;
        float momVz = validMom ? fpgenparts[partMomIdx].vz() : -999;
  
        if( validMom && genVx == 0 && momVx != 0 ){
            genVx = momVx; genVy = momVy; genVz = momVz;
            momVx = geVar("vtxX"); momVy = geVar("vtxY"); momVz = geVar("vtxZ");
        }//<<>>if( genVx == 0 && momVx != 0 )
        if( genVx == 0 ){ genVx = geVar("vtxX"); genVy = geVar("vtxY"); genVz = geVar("vtxZ"); }
        if( validMom && momVx == 0 ){ momVx = geVar("vtxX"); momVy = geVar("vtxY"); momVz = geVar("vtxZ"); }

        const float displacment = validMom ? hypo( genVx-momVx, genVy-momVy, genVz-momVz ) : -10;
        const float momPx = validMom ? fpgenparts[partMomIdx].px() : -999;
        const float momPy = validMom ? fpgenparts[partMomIdx].py() : -999;
        const float momPz = validMom ? fpgenparts[partMomIdx].pz() : -999;
        //const float momP = hypo( momPx, momPy, momPz );
        const float momE = validMom ? fpgenparts[partMomIdx].energy() : -999;
        //const float momBeta = momP/momE;

        int genMomPdgId = validMom ? fpgenparts[partMomIdx].pdgId() : 0;
        genMomPdgId = std::abs( genMomPdgId );

        const int partGMomIdx = genGMomIdx;
        const bool validGMom = partGMomIdx >= 0 && partGMomIdx < nFgenParts;
        int genGMomPdgId = validGMom ? fpgenparts[partGMomIdx].pdgId() : 0;
        genGMomPdgId = std::abs( genGMomPdgId );

        bool isLWZX = ( genPdgId == 11 || genPdgId == 13 ) && ( genMomPdgId == 23 || genMomPdgId == 24 ) && genGMomPdgId == 1000023;
        bool isLWZQ = ( genPdgId == 11 || genPdgId == 13 ) && ( genMomPdgId == 23 || genMomPdgId == 24 ) && genGMomPdgId < 10;
        if( isLWZX ){ hasLWZX = true; nLWZX++; }
        if( isLWZQ ) hasLWZQ = true;

        const int genKidIdx = fpgpKidIdx[nGenParts];
        const int partKidIdx = genKidIdx;
        const bool validKid = partKidIdx >= 0 && partKidIdx < nFgenParts;
        int kidPdgId = validKid ? fpgenparts[partKidIdx].pdgId() : 0;
        kidPdgId = std::abs( kidPdgId );

        //std::cout << "Part : " << nGenParts << " pdg " << genPdgId << " st " << genStatus; 
		//std::cout << " v " << genVx << ", " << genVy << ", " << genVz;
		//std::cout << " 4v e " << genEnergy << " eta " << genEta << " phi " << genPhi << " pt " << genPt << " mass " << genMass;
        //std::cout << " mom " << genMomIdx << " gmom "  << genGMomIdx << " kid " << genKidIdx << std::endl;

		bool isX2 = ( genPdgId == 1000023 );
        //bool isX = ( genPdgId > 1000021 ) and ( genPdgId < 1000038 );
        //bool fsGZ = ( kidPdgId > 0 ) ? ( ( kidPdgId == 22 ) or ( kidPdgId == 23 ) ) : false;
        //bool fsLSP = ( kidPdgId > 0 ) ? ( kidPdgId == 1000022 ) : false;
        //if( isX and ( fsGZ or fsLSP ) ){
        if( isX2 ){
            nXs++;
            float kVx = -999; float kVy = -999; float kVz = -999;
            if( validKid ){
                kVx = fpgenparts[partKidIdx].vx();
                kVy = fpgenparts[partKidIdx].vy();
                kVz = fpgenparts[partKidIdx].vz();
                if( kVx == 0 ){
                    kVx = genVx; kVy = genVy; kVz = genVz;
                    genVx = geVar("vtxX"); genVy = geVar("vtxY"); genVz = geVar("vtxZ");
                }//<<>>if( kVx == 0 )
            }//<<>>if( genKidIdx > -1 )

        	//std::cout << "Part : " << nGenParts << " pdg " << genPdgId << " st " << genStatus; 
        	//std::cout << " v " << genVx << ", " << genVy << ", " << genVz;
			//std::cout << " k " << kVx << ", " << kVy << ", " << kVz;
        	//std::cout << " 4v e " << genEnergy << " eta " << genEta << " phi " << genPhi << " pt " << genPt << " mass " << genMass;
        	//std::cout << " mom " << genMomIdx << " gmom "  << genGMomIdx << " kid " << genKidIdx;
			//std::cout << std::endl;

            const float dis = validKid ? hypo( genVx-kVx, genVy-kVy, genVz-kVz ) : -10;
            float xp = hypo( genPx, genPy, genPz );
            float beta = genEnergy != 0.f ? xp/genEnergy : -1.f;
            float gbeta = genMass != 0.f ? xp/genMass : -1.f;
            float ct = gbeta > 0.f ? dis/gbeta : -1.f;
            if( nXs == 1 ){
                Branches.fillBranch("Xa_MomDisplacment",dis);
                Branches.fillBranch("Xa_beta",beta);
                Branches.fillBranch("Xa_ctau",ct);
                Branches.fillBranch("Xa_PdgId",genPdgId);
                Branches.fillBranch("Xa_Vx",genVx);
                Branches.fillBranch("Xa_Vy",genVy);
                Branches.fillBranch("Xa_Vz",genVz);
                Branches.fillBranch("Xa_Pt",genPt);
                Branches.fillBranch("Xa_P",xp);
                Branches.fillBranch("Xa_Mass",genMass);
                Branches.fillBranch("Xa_Phi",genPhi);
                Branches.fillBranch("Xa_Eta",genEta);
                Branches.fillBranch("Xa_Energy",genEnergy);
            }//<<>>if( nXs == 1 )
            if( nXs == 2 ){
                Branches.fillBranch("Xb_MomDisplacment",dis);
                Branches.fillBranch("Xb_PdgId",genPdgId);
                Branches.fillBranch("Xb_Vx",genVx);
                Branches.fillBranch("Xb_Vy",genVy);
                Branches.fillBranch("Xb_Vz",genVz);
                Branches.fillBranch("Xb_Pt",genPt);
                Branches.fillBranch("Xb_P",xp);
                Branches.fillBranch("Xb_beta",beta);
                Branches.fillBranch("Xb_ctau",ct);
                Branches.fillBranch("Xb_Mass",genMass);
                Branches.fillBranch("Xb_Phi",genPhi);
                Branches.fillBranch("Xb_Eta",genEta);
                Branches.fillBranch("Xb_Energy",genEnergy);
            }//<<>>if( nXs == 2 )
        }//<<>>if( genPdgId > 1000021 and genPdgId < 1000038 )

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
        Branches.fillBranch("genMomIdx",masterMomIndex);
        Branches.fillBranch("genSusId",genSusId);
        Branches.fillBranch("genStatus",genStatus);
        Branches.fillBranch("genCharge",genCharge);
        Branches.fillBranch("genMass",genMass);
        Branches.fillBranch("genMomDisplacment",displacment);
        Branches.fillBranch("genMomPdgId",genMomPdgId);
        Branches.fillBranch("genMomVx",momVx);
        Branches.fillBranch("genMomVy",momVy);
        Branches.fillBranch("genMomVz",momVz);
        Branches.fillBranch("genMomPx",momPx);
        Branches.fillBranch("genMomPy",momPy);
        Branches.fillBranch("genMomPz",momPz);
        Branches.fillBranch("genMomEnergy",momE);
        Branches.fillBranch("genGMomPdgId",genGMomPdgId);
        Branches.fillBranch("genLWZX",isLWZX);
        Branches.fillBranch("genLWXQ",isLWZQ);

        nGenParts++;
    }//<<>> for (const auto & genpart : fpgenparts )


    if( GenDEBUG ) std::cout << " - Gen Event type Endjobs " << nXs << std::endl;
    Branches.fillBranch("genNtotal", unsigned( fpgenparts.size() ));

	int fsType = 0;
	if( nXG == 2 ){ isGG = true; fsType = 1; }
	if( nXZ == 2 ){ isZZ = true; fsType = 2; }
	if( nXG == 1 and nXZ == 1 ){ isGZ = true; fsType = 3; }

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
    Branches.fillBranch("isGG",isGG);
    Branches.fillBranch("isZZ", isZZ);
    Branches.fillBranch("isGZ",isGZ);
    Branches.fillBranch("fsType",fsType);
	Branches.fillBranch("nXs",nXs);

	geVar.fill("genWgt",wgt);
    if( GenDEBUG ) std::cout << "GenPart : Done " << std::endl;

}//<<>>void KUCMSGen::ProcessEvent()

///////  -----------------------------   End of ProcessEvent post and main ------------------------------------------

void KUCMSGenObjectMini::EndJobs(){}

///////  -----------------------------   helper functions --------------------------------------------------------------

std::vector<int> KUCMSGenObjectMini::getGenPhoMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& reco_e ){

	// gen matching to supcer clusters froom reeco phtons to gen photons
	std::vector<reco::TransientTrack> null;
	return getGenMatch( sc, reco_e, 22, null );

}//<<>>std::vector<int> KUCMSGenObjectMini::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObjectMini::getGenEleMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& reco_e, std::vector<reco::TransientTrack>& trakcs ){

    // gen matching to supcer clusters froom reco electrons to gen electrons
    return getGenMatch( sc, reco_e, 11, trakcs );

}//<<>>std::vector<int> KUCMSGenObjectMini::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObjectMini::getGenMuonMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& reco_e, std::vector<reco::TransientTrack>& trakcs ){

    // gen matching to supcer clusters froom reco electrons to gen electrons
    return getGenMatch( sc, reco_e, 13, trakcs );

}//<<>>std::vector<int> KUCMSGenObjectMini::getGenPhoMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObjectMini::getGenMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& reco_e ){

	// matchoing any particle :: matching valid for nutral particles only
	std::vector<reco::TransientTrack> null;
	return getGenMatch( sc, reco_e, 0, null );

}//<<>>std::vector<int> KUCMSGenObjectMini::getGenMatch( const std::vector<reco::SuperCluster*> sc, std::vectro<float> reco_e )

std::vector<int> KUCMSGenObjectMini::getGenMatch( const std::vector<v3fPoint>& sc, const std::vector<float>& en, int select, std::vector<reco::TransientTrack>& tracks ){

    bool thisDEBUG = GenDEBUG;
	//bool thisDEBUG = true;

    // GenParticle Info   -------------------------------------------------------------------
    if( thisDEBUG ) std::cout << "Getting GenParton Match Information for type " << select << std::endl;
	// code caclutes eta & phi of the momentum vector requird for a nuetral particle 
	// produced at a gen particles production vertix to cause the given supercluster and compares
	// this to the actual gen particles mometum vectors eta & phi to find best, exclusive, match
	// between a reco particle with that supercluster and a gen particle

    const float drthres = 0.4f;
    const float demax = 1.0f;
    const int nSCPtrs = static_cast<int>(sc.size());
    std::vector<int> recoPartIndx(nSCPtrs, -1);

    if( en.size() != sc.size() ){
        std::cout << "KUCMSGenObjectMini::getGenMatch: sc and energy vector sizes differ ("
                  << sc.size() << " != " << en.size() << "). No matches returned." << std::endl;
        return recoPartIndx;
    }
    if( ( select == 11 || select == 13 ) && tracks.size() != sc.size() ){
        std::cout << "KUCMSGenObjectMini::getGenMatch: reco and track vector sizes differ ("
                  << sc.size() << " != " << tracks.size() << "). No matches returned." << std::endl;
        return recoPartIndx;
    }

    std::vector<int> genPartIndx;
    const int nGenPart = static_cast<int>(fpgenparts.size());
    for( int idx = 0; idx < nGenPart; idx++ ){
        if( fpgenparts[idx].status() != 1 ) continue;
        if( select != 0 && std::abs(fpgenparts[idx].pdgId()) != select ) continue;
        genPartIndx.push_back(idx);
        if( static_cast<std::size_t>(idx) < pmatdr.size() && static_cast<std::size_t>(idx) < pmatde.size() ){
            pmatdr[idx] = drthres;
            pmatde[idx] = 1000.f;
        }
    }

    const int nGPIs = static_cast<int>(genPartIndx.size());
    if( thisDEBUG ) std::cout << " - Finding Matches nSC " << nSCPtrs << " nGP " << nGPIs << std::endl;
    if( nSCPtrs == 0 || nGPIs == 0 ) return recoPartIndx;

    // Build the bipartite graph using only pairs which pass both thresholds.
    // A min-cost maximum flow then maximizes the number of matches first and,
    // among assignments of that size, minimizes the total delta-R.
    std::vector<std::vector<float>> matchDR(nSCPtrs, std::vector<float>(nGPIs, drthres));
    std::vector<std::vector<float>> matchDE(nSCPtrs, std::vector<float>(nGPIs, 1000.f));
    std::vector<std::vector<bool>> validMatch(nSCPtrs, std::vector<bool>(nGPIs, false));

    for( int it = 0; it < nSCPtrs; it++ ){
        const auto rhX = sc[it].x();
        const auto rhY = sc[it].y();
        const auto rhZ = sc[it].z();
        for( int pgidx = 0; pgidx < nGPIs; pgidx++ ){
            const int partidx = genPartIndx[pgidx];
            const auto & genPart = fpgenparts[partidx];
            const auto gnEta = genPart.eta();
            const auto gnPhi = genPart.phi();
            const auto gnX = genPart.vx();
            const auto gnY = genPart.vy();
            const auto gnZ = genPart.vz();
            const double transverseDistance = std::hypot(rhX-gnX, rhY-gnY);
            if( transverseDistance == 0.0 ) continue;

            double cmEta = std::asinh((rhZ-gnZ)/transverseDistance);
            double cmPhi = std::atan2(rhY-gnY, rhX-gnX);

            if( select == 11 || select == 13 ){
                SteppingHelixPropagator propagator( tracks[it].field(), anyDirection );
                const GlobalPoint genLocation( gnX, gnY, gnZ );
                const FreeTrajectoryState propLocation = propagator.propagateWithPath(tracks[it].initialFreeState(), genLocation).first;
                const GlobalPoint pca( propLocation.position() );
                const reco::Track track(tracks[it].track());
                const double denominator = 87.78*std::sin(track.lambda());
                if( denominator == 0.0 ) continue;
                const double tmin = (pca.z() - track.vz())*track.qoverp()/denominator;
                cmPhi -= tmin;
            }

            const float dr = std::sqrt(reco::deltaR2(gnEta, gnPhi, cmEta, cmPhi));
            const float genenergy = genPart.energy();
            if( genenergy <= 0.f ) continue;
            const float de = std::abs(en[it] - genenergy)/genenergy;
            if( !std::isfinite(dr) || !std::isfinite(de) ) continue;
            if( dr < drthres && de < demax ){
                validMatch[it][pgidx] = true;
                matchDR[it][pgidx] = dr;
                matchDE[it][pgidx] = de;
            }
        }
    }

    struct MatchEdge {
        int to;
        int reverse;
        int capacity;
        double cost;
    };

    const int source = 0;
    const int recoBase = 1;
    const int genBase = recoBase + nSCPtrs;
    const int sink = genBase + nGPIs;
    const int nNodes = sink + 1;
    std::vector<std::vector<MatchEdge>> graph(nNodes);
    auto addEdge = [&graph]( int from, int to, int capacity, double cost ){
        const int forwardReverse = static_cast<int>(graph[to].size());
        const int backwardReverse = static_cast<int>(graph[from].size());
        graph[from].push_back({to, forwardReverse, capacity, cost});
        graph[to].push_back({from, backwardReverse, 0, -cost});
    };

    for( int recoIdx = 0; recoIdx < nSCPtrs; recoIdx++ ) addEdge(source, recoBase+recoIdx, 1, 0.0);
    for( int genIdx = 0; genIdx < nGPIs; genIdx++ ) addEdge(genBase+genIdx, sink, 1, 0.0);
    for( int recoIdx = 0; recoIdx < nSCPtrs; recoIdx++ ){
        for( int genIdx = 0; genIdx < nGPIs; genIdx++ ){
            if( validMatch[recoIdx][genIdx] ) addEdge(recoBase+recoIdx, genBase+genIdx, 1, matchDR[recoIdx][genIdx]);
        }
    }

    const double infinity = std::numeric_limits<double>::infinity();
    while( true ){
        std::vector<double> distance(nNodes, infinity);
        std::vector<int> previousNode(nNodes, -1);
        std::vector<int> previousEdge(nNodes, -1);
        std::vector<bool> inQueue(nNodes, false);
        std::queue<int> nodes;
        distance[source] = 0.0;
        nodes.push(source);
        inQueue[source] = true;

        while( !nodes.empty() ){
            const int node = nodes.front();
            nodes.pop();
            inQueue[node] = false;
            for( int edgeIdx = 0; edgeIdx < static_cast<int>(graph[node].size()); edgeIdx++ ){
                const MatchEdge & edge = graph[node][edgeIdx];
                if( edge.capacity == 0 ) continue;
                const double newDistance = distance[node] + edge.cost;
                if( newDistance >= distance[edge.to] ) continue;
                distance[edge.to] = newDistance;
                previousNode[edge.to] = node;
                previousEdge[edge.to] = edgeIdx;
                if( !inQueue[edge.to] ){
                    nodes.push(edge.to);
                    inQueue[edge.to] = true;
                }
            }
        }

        if( previousNode[sink] == -1 ) break;
        for( int node = sink; node != source; node = previousNode[node] ){
            MatchEdge & edge = graph[previousNode[node]][previousEdge[node]];
            edge.capacity--;
            graph[node][edge.reverse].capacity++;
        }
    }

    for( int recoIdx = 0; recoIdx < nSCPtrs; recoIdx++ ){
        const int recoNode = recoBase + recoIdx;
        for( const MatchEdge & edge : graph[recoNode] ){
            if( edge.to < genBase || edge.to >= sink || edge.capacity != 0 ) continue;
            const int filteredGenIdx = edge.to - genBase;
            const int genIndex = genPartIndx[filteredGenIdx];
            if( genIndex < 0 || genIndex >= nGenPart ) continue;
            recoPartIndx[recoIdx] = genIndex;
            if( static_cast<std::size_t>(genIndex) < pmatdr.size() && static_cast<std::size_t>(genIndex) < pmatde.size() ){
                pmatdr[genIndex] = matchDR[recoIdx][filteredGenIdx];
                pmatde[genIndex] = matchDE[recoIdx][filteredGenIdx];
            }
            break;
        }
    }

	if( thisDEBUG ){ std::cout << " -- IdLIST: #" << nSCPtrs << " "; for( auto it : recoPartIndx ) std::cout << it << " "; std::cout << std::endl; }
	if( thisDEBUG ) std::cout << " - Done with Gen Matching" << std::endl;
    return recoPartIndx;

}//<<>>getGenPhoMatch( reco::SuperClusterCollection *scptr, std::vector<reco::GenParticle> fgenparts )

std::vector<float> KUCMSGenObjectMini::getGenPartMatch( const reco::SuperCluster* scptr, float energy ){

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
    for( int index = 0; index < static_cast<int>(fpgenparts.size()); index++ ){

        const auto & genPart = fpgenparts[index];

		if( genPart.status() != 1 ) continue;

        const auto rhX = scptr->x();
        const auto rhY = scptr->y();
        const auto rhZ = scptr->z();
        auto gnX = genPart.vx();
        auto gnY = genPart.vy();
        auto gnZ = genPart.vz();
        const auto ge = genPart.energy();
		//if( gnX == 0 ){ gnX = pvx; gnY = pvy; gnZ = pvz; }
		//if( GenDEBUG ) std::cout << " GenPart: " << gnX << " " << gnY << " " << gnZ << std::endl;
		//if( ge < 1.0 ) continue;
        const auto transverseDistance = hypo(rhX-gnX,rhY-gnY);
        if( transverseDistance == 0.0 || ge <= 0.0 ) continue;
        auto cphoEta = std::asinh((rhZ-gnZ)/transverseDistance);
        auto cphoPhi = std::atan2(rhY-gnY,rhX-gnX);
        auto dr = std::sqrt(reco::deltaR2(genPart.eta(), genPart.phi(), cphoEta, cphoPhi ));
        auto re = energy/ge;
        if( dr < minDr ){
            minSDr = minDr;
            minSRe = minRe;
            matchedSIdx = matchedIdx;
            llpSId = llpId;
            minDr = dr;
            minRe = re;
            matchedIdx = index;
            llpId = fpgpLlp[index];
        } else if( dr < minSDr ){
            minSDr = dr;
            minSRe = re;
            matchedSIdx = index;
            llpSId = fpgpLlp[index];
        }//<<>>if( dr < minDr && dp < minDp )

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

std::pair<int,std::pair<float,float>> KUCMSGenObjectMini::getGenLeptonMatch( const reco::RecoCandidate & parton, int type ){

	// basic dr match for best match to a single gen particle   - thrshold 0.8

	float eta = parton.eta();
    float phi = parton.phi();
	float en = parton.energy();
	int parPdgId = std::abs(parton.pdgId());
	if( parPdgId == 0 ) parPdgId = std::abs(type);
	//std::cout << " -- type/pdgid: " << type << " ?= " << parPdgId << std::endl;
    float minRe(100.0);
    float minDr(0.8);
	int matchedIdx(-9);
	//std::cout << " -- Getting fgenparts size " << std::endl;

    int nGenPart = fpgenparts.size();
    //std::cout << " - Starting fgenparts loop for : " << nGenPart << std::endl;
    for( int idx = 0; idx < nGenPart; idx++ ){

        //std::cout << " -- in loop with :: " << std::endl;
        if( fpgenparts[idx].status() != 1 ) continue;

		int pdgId = std::abs( fpgenparts[idx].pdgId() );
		if( pdgId != parPdgId ) continue;
		//std::cout << " -- pdgID: " << pdgId << std::endl;		

		int genMomIdx = fpgpMomIdx[idx];
        const bool validMom = genMomIdx >= 0 && genMomIdx < nGenPart;
        int genMomPdgId = validMom ? std::abs(fpgenparts[genMomIdx].pdgId()) : 0;
		if( genMomPdgId != 23 && genMomPdgId != 24 && genMomPdgId != 1000023  ) continue;
		//std::cout << " -- momPdgID: " << genMomPdgId << std::endl;

		float geta = fpgenparts[idx].eta();
		float gphi = fpgenparts[idx].phi();
		float gen = fpgenparts[idx].energy();
        if( gen <= 0.f ) continue;
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

}//<<>> int KUCMSGenObjectMini::getGenPartMatch( const T &lepton )

int KUCMSGenObjectMini::getGenSigPhoXMother( uInt genIndex ){

    return getGenSigPhoXMother(genIndex, 0);

}//<<>>int KUCMSGenObjectMini::getGenSigPhoInfo( uInt genIndex )

int KUCMSGenObjectMini::getGenSigPhoXMother( uInt genIndex, int loopcnt ){

    if( loopcnt > 100 ) return -1;
    if( genIndex >= fpgenparts.size() || genIndex >= fpgpMomIdx.size() ) return -1;

    const int motherIndex = fpgpMomIdx[genIndex];
    if( motherIndex < 0 || motherIndex >= static_cast<int>(fpgenparts.size()) ) return -1;
    if( std::abs(fpgenparts[motherIndex].pdgId()) > 1000000 ) return motherIndex;
    return getGenSigPhoXMother(static_cast<uInt>(motherIndex), loopcnt+1);

}//<<>>int KUCMSGenObjectMini::getGenSigPhoInfo( uInt genIndex )

int KUCMSGenObjectMini::getGenSigEleXMother( uInt genIndex ){

	return getGenSigPhoXMother( genIndex );

}//<<>>std::vector<float> KUCMSGenObjectMini::getGenSigEleInfo( uInt genIndex )


std::vector<float> KUCMSGenObjectMini::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  ){
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

}//>>>>vector<float> KUCMSGenObjectMini::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  )

void KUCMSGenObjectMini::kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz ){

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

int KUCMSGenObjectMini::getEvtType( int genPdgId ){

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

}//<<>>int KUCMSGenObjectMini::getEvtType( int genPdgId )

const reco::Candidate* KUCMSGenObjectMini::llpGenSearch( const reco::Candidate* kid, int target ){

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

int KUCMSGenObjectMini::llpGenChaseP( const reco::GenJet & kid, int depth ){


    int genPartID = 100;
    int nMoms = kid.numberOfMothers();
    for( int gmit(0); gmit < nMoms; gmit++ ){
        auto mGenPartID = llpGenChaseP(kid.mother(gmit),0);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;
    }//<<>>for( int gmit(0); gmit < nMoms; gmit++ )
    return genPartID;

}//<<>>int KUCMSGenObjectMini::llpGenChaseP( const reco::GenJet & kid, int depth )

int KUCMSGenObjectMini::llpGenChaseP( const reco::GenParticle & kid, int depth ){

    int genPartSusID = 100;
    int nMoms = kid.numberOfMothers();
    int genPartId = std::abs(kid.pdgId());
    //if( GenDEBUG ) std::cout << "Genpart motherCase for : " << genPartId << std::endl;
    for( int gmit(0); gmit < nMoms; gmit++ ){
        auto mGenPartSusID = llpGenChaseP(kid.mother(gmit),genPartId);
        if( mGenPartSusID < genPartSusID ) genPartSusID = mGenPartSusID;
    }//<<>>for( int gmit(0); gmit < nMoms; gmit++ )
    return genPartSusID;

}//<<>>int KUCMSGenObjectMini::llpGenChaseP( const reco::GenJet & kid, int depth )

int KUCMSGenObjectMini::llpGenChaseP( const reco::Candidate* kid, int childPdgId ){

    std::unordered_set<const reco::Candidate*> activePath;
    return llpGenChasePImpl(kid, childPdgId, activePath, 0);

}//<<>>int KUCMSGenObjectMini::llpGenChaseP

int KUCMSGenObjectMini::llpGenChasePImpl( const reco::Candidate* kid, int childPdgId,
        std::unordered_set<const reco::Candidate*>& activePath, int depth ) const {

    if( kid == nullptr || depth > 100 ) return 100;
    if( !activePath.insert(kid).second ) return 100;
    int nMoms = kid->numberOfMothers();
	if( nMoms == 1 && kid->mother(0) != nullptr && kid->pdgId() == kid->mother(0)->pdgId() ){
        return llpGenChasePImpl(kid->mother(0), childPdgId, activePath, depth+1);
    }

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

    if( nMoms == 0 ) return 99;

    int genPartID = 100;
    for( int gmit(0); gmit < nMoms; gmit++ ){

        auto genmom = kid->mother(gmit);
        if( genmom == nullptr ) continue;
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

        auto branchPath = activePath;
        auto mGenPartID = llpGenChasePImpl(genmom, 0, branchPath, depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;

    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    //if( GenDEBUG ) std::cout << " --- Found llpgen chase canadate: " << kidPdgID << " susid: " << genPartID << std::endl;
    return genPartID;

}//<<>>int KUCMSGenObjectMini::llpGenChasePImpl

void KUCMSGenObjectMini::motherChase( const reco::Candidate* kid, string  depth ){

    //if( ( kid->numberOfMothers() == 1 ) && ( kid->pdgId() == (kid->mother(0))->pdgId() ) ){ motherChase( kid->mother(0), depth ); }
    if( false );
    else{
	    if( true ){ //status()
		//if( kid->status() < 40 ){
    		std::cout <<  depth  << " CkidID : " << kid->pdgId() << " pt : " << kid->pt() << " status : " << kid->status();
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

void KUCMSGenObjectMini::motherChase( const reco::GenParticle & kid, string  depth ){
 
	//if( ( kid.numberOfMothers() == 1 ) && ( kid.pdgId() == (kid.mother(0))->pdgId() ) ){ motherChase( kid.mother(0), depth ); }
	if( false );
	else{
	    if( true ){ //status()
	    //if( kid.status() < 40 ){
    		std::cout <<  depth  << " GkidID : " << kid.pdgId() << " pt : " << kid.pt() << " status : " << kid.status();
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

int KUCMSGenObjectMini::typeChase( const reco::Candidate* kid ){

    std::unordered_set<const reco::Candidate*> visited;
    const reco::Candidate* current = kid;
    for( int depth = 0; current != nullptr && depth <= 100; depth++ ){
        if( !visited.insert(current).second ) return 0;
        const auto nMothers = current->numberOfMothers();
        if( nMothers == 2 ) return std::abs(current->pdgId());
        if( nMothers != 1 ) return 0;
        current = current->mother(0);
    }
    return 0;

}//<<>> void MotherChase( Candidate* kid, string depth  )

int KUCMSGenObjectMini::typeChase( const reco::GenParticle & kid ){

	if( kid.numberOfMothers() == 2 ) return std::abs(kid.pdgId());
    if( kid.numberOfMothers() == 1 ) return typeChase(kid.mother(0));
    return 0;

}//<<>> void MotherChase( Candidate* kid, string depth  )

int KUCMSGenObjectMini::llpGenChase( const reco::Candidate* kid ){

    auto match = ( std::abs(kid->pdgId()) > 1000020 ) && ( std::abs(kid->pdgId()) < 1000040 );
    if( match ) return kid->pdgId();
    for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
        if( kid->mother(gmit)->numberOfMothers() > 0 ){ auto r = llpGenChase(kid->mother(gmit)); if(r != 0) return r; }
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return 0;

}//<<>> void llpChase( Candidate* kid )

int KUCMSGenObjectMini::llpGenChase( const reco::GenParticle & kid ){

    auto match = ( std::abs(kid.pdgId()) > 1000020 ) && ( std::abs(kid.pdgId()) < 1000040 );
    auto hard = kid.isLastCopy();
    if( match && hard ) return kid.pdgId();
    for( long unsigned int gmit(0); gmit < kid.numberOfMothers(); gmit++ ){
        //if( kid->pdgId() == 6000113 ) return true;
        if( kid.mother(gmit)->numberOfMothers() > 0 ){ auto r = llpGenChase(kid.mother(gmit)); if(r != 0) return r; }
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return 0;

}//<<>> void llpChase( Candidate* kid )

string KUCMSGenObjectMini::bigKidChase( std::vector<reco::CandidatePtr> kids, float vx ){

    std::string result("");
    float maxpt(0.0);
    const reco::Candidate* bigKid(nullptr);
    for( auto &kid : kids ){
        if( kid->pt() > maxpt ){
            maxpt = kid->pt();
            bigKid = kid.get();
        }//<<>>if( kid->pt() > maxpt )
    }//<<>>for( auto kid : kids )
    if( bigKid == nullptr ) return " -- BigKid : none";
    result += std::string(" -- BigKid : pdgID : ") + std::to_string(bigKid->pdgId()) +
        std::string(" pt : ") + std::to_string(bigKid->pt());
    if( bigKid->numberOfMothers() > 0 && bigKid->mother(0) != nullptr ){
        result += std::string(" vtx : ") + std::to_string(bigKid->mother(0)->vx()) +
            std::string(" / ") + std::to_string(vx);
    }
    return result;

}//<<>>string bigKidChase( std::vector<reco::CandidatePtr> kids, float vx )

std::vector<float> KUCMSGenObjectMini::getGenJetInfo( float jetEta, float jetPhi, float jetEnergy ){

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
    for(const auto &genPart : fpgenparts ){

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
		genPartLlpId = float(fpgpLlp[mGenPartIdx]);
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

}//<<>>std::vector<float> KUCMSGenObjectMini::getGenJetInfo( float jetEta, float jetPhi, float jetPt )

// Gen Electron methods

/*
std::vector<int>KUCMSGenObjectMini:: MomIDs(const reco::GenParticle &genElectron) const {

  auto mother = genElectron.mother();


  std::vector<int> motherIDs;

  //if( genElectron.numberOfMothers() < 1 ) return motherIDs; 
  auto mother = genElectron.mother(0);
  while(mother->pt() > 0) {
    const int motherID = mother->pdgId();

    mother = mother->mother();

    if(motherID == mother->pdgId())
      continue;

    motherIDs.push_back(motherID);
  }

  return motherIDs;

}

bool KUCMSGenObjectMini::IsMotherZ(const reco::GenParticle &genElectron) const {
  LepMomType momType = ClassifyGenElectron(genElectron);
  return (momType == kZ);
}

bool KUCMSGenObjectMini::isSignalGenElectron(const reco::GenParticle &genElectron) const {

  LepMomType momType = ClassifyGenElectron(genElectron);
  return (momType == kZ || momType == kSusy);

}

LepMomType KUCMSGenObjectMini::ClassifyGenElectron(const reco::GenParticle &genElectron) const {

  std::vector<int> motherIDs(MomIDs(genElectron));

  LepMomType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignLeptonMomType(id);

    if(momType != kUnmatched)
      break;
  }
  return momType;
}

void KUCMSGenObjectMini::PrintMother(const LepMomType &momType) const {

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

void KUCMSGenObjectMini::GenElectronContent() const {

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

LepMomType KUCMSGenObjectMini::AssignLeptonMomType(const int motherID) const {
  
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
GenClassifiedElectrons<T> KUCMSGenObjectMini::GetGenClassifiedElectrons(const std::vector<T> &candidateTracks) const {
  
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
