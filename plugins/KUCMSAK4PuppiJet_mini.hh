// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//
// KUCMS AK4 PUPPI Jet Object
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

// Includes for all object types
#include <algorithm>
#include <cmath>
#include <memory>

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
//#include "DataFormats/AK4Jet/interface/AK4Jet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

// JECS
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JERs
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit_mini.hh"
#include "KUCMSGenObjects_mini.hh"
#include "KUCMSGenJet_mini.hh"
#include "KUCMSElectron_mini.hh"
#include "KUCMSPhoton_mini.hh"

#ifndef KUCMSAK4PuppiJetObjectMiniHeader
#define KUCMSAK4PuppiJetObjectMiniHeader

//#define AK4PuppiJetDEBUG true
#define AK4PuppiJetDEBUG false

using namespace edm;

// Run 3 AK4 PUPPI PF jet ID, CMS-DP-2024/028:
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsDP24028
// Energy fractions are defined using the uncorrected jet energy (as in
// pat::Jet's fraction accessors); multiplicities are the PF multiplicities.
// |eta| <= 2.6: NHF < .99, NEF < .90, Nconst > 1, CHF > .01, Nch > 0.
// 2.6 < |eta| <= 2.7: NHF < .90, NEF < .99.
// 2.7 < |eta| <= 3.0: NHF < .99.
// 3.0 < |eta| < 5.0: NEF < .40, Nneutral >= 2.
// Tight lepton veto adds MUF < .80 and CEF < .80 for |eta| <= 2.7.
// The eta=3.0 boundary follows the CMSSW selector's region convention.
// Return -1 for unsupported eta or invalid PF inputs, not a failed decision.
static void run3PuppiJetId( double eta, double nhf, double nef, double chf,
                            double cef, double muf, int nch, int nneutral,
                            int& tight, int& tightLepVeto ){
  tight = tightLepVeto = -1;
  const double absEta = std::abs(eta);
  if( !std::isfinite(absEta) || absEta >= 5.0 || !std::isfinite(nhf) ||
      !std::isfinite(nef) || !std::isfinite(chf) || !std::isfinite(cef) ||
      !std::isfinite(muf) || nch < 0 || nneutral < 0 ) return;
  bool passes = false;
  if( absEta <= 2.6 ) passes = nhf < .99 && nef < .90 && nch + nneutral > 1 && chf > .01 && nch > 0;
  else if( absEta <= 2.7 ) passes = nhf < .90 && nef < .99;
  else if( absEta <= 3.0 ) passes = nhf < .99;
  else passes = nef < .40 && nneutral >= 2;
  tight = int(passes);
  tightLepVeto = int(passes && (absEta > 2.7 || (muf < .80 && cef < .80)));
}

/** Collects AK4 jets and derives their object and generator associations. */
class KUCMSAK4PuppiJetObjectMini : public KUCMSObjectBase {

  public:

  // use base class constructor
  KUCMSAK4PuppiJetObjectMini( const edm::ParameterSet& iConfig );
  ~KUCMSAK4PuppiJetObjectMini(){};

  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadAK4JetTokens( edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken ){ jetsToken_ = jetsToken; inputEnabled_ = true; };
  void LoadPfcandTokens( edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcandToken ){ pfcandToken_ = pfcandToken; };
  // sets up branches, do preloop jobs
  void InitObject( TTree* fOutTree );
  // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
  // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed
  void LoadRecHitObject( KUCMSEcalRecHitObjectMini* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed
  void LoadGenObject( KUCMSGenObjectMini* genObjs_ ){ genObjs = genObjs_; };
  void LoadGenJetObject( KUCMSGenJetObjectMini* genJetObj_ ){ genJetObj = genJetObj_; };
  void LoadElectronObject( KUCMSElectronObjectMini* electronObj_ ){ electronObj = electronObj_; };
  void LoadPhotonObject( KUCMSPhotonObjectMini* photonObj_ ){ photonObj = photonObj_; };

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

  private:

  //std::vector<examnple> fexamples;
  std::vector<pat::Jet> fjets;
  std::vector<int> jetIds;
  std::unique_ptr<PFJetIDSelectionFunctor> tightJetId;
  std::unique_ptr<PFJetIDSelectionFunctor> tightLepVetoJetId;
  bool run3JetId_ = false;

  //edm::EDGetTokenT<example> exampleToken_;
  //edm::Handle<example> example_;
  // Jets ----------------------------------------------------------

  //const edm::InputTag jetsTag;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::Handle<std::vector<pat::Jet>> jets_;
  bool inputEnabled_ = false;
  bool inputAvailable_ = false;

  // PF Candidates View
  //const edm::InputTag pfcandTag;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcandToken_;
  edm::Handle<edm::View<pat::PackedCandidate>> pfcands_;

  // Jet corrections  ----------------------------------------------------------

  // JEC & JER  -> From miniAOD code not converted
  // JECs
  edm::ESHandle<JetCorrectorParametersCollection> jetCorrH;
  // JERs
  JME::JetResolution jetRes;
  JME::JetResolutionScaleFactor jetRes_sf;

  // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
  // exampleObject* otherObjectPtr;

  KUCMSEcalRecHitObjectMini* rhObj;
  KUCMSGenObjectMini* genObjs;
  KUCMSGenJetObjectMini* genJetObj;
  KUCMSElectronObjectMini* electronObj;
  KUCMSPhotonObjectMini* photonObj;

};//<<>>class KUCMSAK4Jet : public KUCMSObjectBase

KUCMSAK4PuppiJetObjectMini::KUCMSAK4PuppiJetObjectMini( const edm::ParameterSet& iConfig ){
  // ---- end constructor initilizations  --------------------------

  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : false );
  const std::string jetRunEra = iConfig.existsAs<std::string>("RunEra") ? iConfig.getParameter<std::string>("RunEra") : "";
  PFJetIDSelectionFunctor::Version_t version = PFJetIDSelectionFunctor::N_VERSIONS;
  if( jetRunEra == "Run2_2016_Pre" || jetRunEra == "Run2_2016_Post" ) version = PFJetIDSelectionFunctor::RUN2UL16PUPPI;
  else if( jetRunEra == "Run2_2017" || jetRunEra == "Run2_2018" ) version = PFJetIDSelectionFunctor::RUN2ULPUPPI;
  // 2022 inputs are assumed reprocessed data or MC, never early prompt data.
  // For 2024-2026, use the 2023 convention as an analysis assumption.
  else if( jetRunEra == "Run3_2022" || jetRunEra == "Run3_2023" || jetRunEra == "Run3_2024" || jetRunEra == "Run3_2025" || jetRunEra == "Run3_2026" ){
    run3JetId_ = true;
  }
  if( version != PFJetIDSelectionFunctor::N_VERSIONS ){
    tightJetId.reset(new PFJetIDSelectionFunctor(version, PFJetIDSelectionFunctor::TIGHT));
    tightLepVetoJetId.reset(new PFJetIDSelectionFunctor(version, PFJetIDSelectionFunctor::TIGHTLEPVETO));
  }
  cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
  cfPrm.set( "jetPTmin", iConfig.existsAs<double>("jetPTmin") ? iConfig.getParameter<double>("jetPTmin") : 10.0 );
  cfPrm.set( "jetEtaMax", iConfig.existsAs<double>("jetEtaMax") ? iConfig.getParameter<double>("jetEtaMax") : 3.5 );
  cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.44 );

}//<<>>KUCMSAK4Jet::KUCMSAK4Jet( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSAK4PuppiJetObjectMini::InitObject( TTree* fOutTree ){

  Branches.makeBranch("Available","PuppiJet_available",BOOL,"PUPPI jet input collection is available");
  Branches.makeBranch("E","PuppiJet_energy",VFLOAT);
  Branches.makeBranch("M","PuppiJet_mass",VFLOAT);
  Branches.makeBranch("Pt","PuppiJet_pt",VFLOAT);
  Branches.makeBranch("Eta","PuppiJet_eta",VFLOAT);
  Branches.makeBranch("Phi","PuppiJet_phi",VFLOAT);
  Branches.makeBranch("Area","PuppiJet_area",VFLOAT,"jet catchment area, for JECs");
  Branches.makeBranch("LegacyCoverage","PuppiJet_passLegacyCoverage",VBOOL,"jet passes the pre-expansion pt and eta selection");
  Branches.makeBranch("RawFactor","PuppiJet_rawFactor",VFLOAT,"1 - uncorrected JEC factor; -999 if unavailable");
  Branches.makeBranch("MuonSubtrRawPt","PuppiJet_muonSubtrRawPt",VFLOAT,"muon-subtracted raw jet pt; -999 if unavailable");
  Branches.makeBranch("MuonSubtrRawEta","PuppiJet_muonSubtrRawEta",VFLOAT,"muon-subtracted raw jet eta; -999 if unavailable");
  Branches.makeBranch("MuonSubtrRawPhi","PuppiJet_muonSubtrRawPhi",VFLOAT,"muon-subtracted raw jet phi; -999 if unavailable");
  Branches.makeBranch("MuonSubtrFactor","PuppiJet_muonSubtrFactor",VFLOAT,"1 - muon-subtracted raw pt / raw pt; -999 if unavailable");
  Branches.makeBranch("MuonSubtrDeltaEta","PuppiJet_muonSubtrDeltaEta",VFLOAT,"muon-subtracted raw eta minus jet eta; -999 if unavailable");
  Branches.makeBranch("MuonSubtrDeltaPhi","PuppiJet_muonSubtrDeltaPhi",VFLOAT,"muon-subtracted raw phi minus jet phi; -999 if unavailable");
  Branches.makeBranch("IdLoose","PuppiJet_idLoose",VINT,"PUPPI PF jet ID: 1 pass, 0 fail, -1 unavailable");
  Branches.makeBranch("IdTight","PuppiJet_idTight",VINT,"era-dependent PUPPI PF jet ID: 1 pass, 0 fail, -1 unavailable");
  Branches.makeBranch("IdTightLepVeto","PuppiJet_idTightLepVeto",VINT,"era-dependent PUPPI PF jet ID: 1 pass, 0 fail, -1 unavailable");
  Branches.makeBranch("PuIdLoose","PuppiJet_puIdLoose",VINT,"PUPPI PU ID decision disabled: -1 unavailable");
  Branches.makeBranch("PuIdMedium","PuppiJet_puIdMedium",VINT,"PUPPI PU ID decision disabled: -1 unavailable");
  Branches.makeBranch("PuIdTight","PuppiJet_puIdTight",VINT,"PUPPI PU ID decision disabled: -1 unavailable");
  Branches.makeBranch("PuIdDisc","PuppiJet_puIdDisc",VFLOAT,"embedded PUPPI PU ID discriminant, training depends on input; -999 if unavailable");
  Branches.makeBranch("PuIdBitmap","PuppiJet_puIdBitmap",VINT,"embedded PUPPI PU ID raw bitmap, interpretation depends on input; -1 if unavailable");

  Branches.makeBranch("nKids","PuppiJet_nConstituents",VINT,"Number of particles in the jet");
  Branches.makeBranch("NHF","PuppiJet_neHEF",VFLOAT,"neutral Hadron Energy Fraction");
  Branches.makeBranch("NEMF","PuppiJet_neEmEF",VFLOAT,"neutral Electromagnetic Energy Fraction");
  Branches.makeBranch("CHF","PuppiJet_chHEF",VFLOAT,"charged Hadron Energy Fraction");
  Branches.makeBranch("CEMF","PuppiJet_chEmEF",VFLOAT,"charged Electromagnetic Energy Fraction");
  Branches.makeBranch("MUF","PuppiJet_muEF",VFLOAT,"muon Energy Fraction");
  Branches.makeBranch("NHM","PuppiJet_neHM",VFLOAT,"neutral hadron multiplicity");
  Branches.makeBranch("CHM","PuppiJet_chHM",VFLOAT,"charged hadron multiplicity");

  Branches.makeBranch("Parts","PuppiJet_egIndxs",VVUINT);
  Branches.makeBranch("DrRhIds","PuppiJet_drRhIds",VVUINT);

  if( cfFlag("hasGenInfo") ){

    Branches.makeBranch("GenJetIdx","PuppiJet_genJetIdx",VINT,"PAT GenJet index; -1 no match, -9 match not saved");
    Branches.makeBranch("GenImpactAngle","PuppiJet_genImpactAngle",VFLOAT);
    Branches.makeBranch("GenTime","PuppiJet_genTime",VFLOAT);
    Branches.makeBranch("GenPt","PuppiJet_genPt",VFLOAT);
    Branches.makeBranch("GenEta","PuppiJet_genEta",VFLOAT);
    Branches.makeBranch("GenPhi","PuppiJet_genPhi",VFLOAT);
    Branches.makeBranch("GenEnergy","PuppiJet_genEnergy",VFLOAT);
    Branches.makeBranch("GenDrMatch","PuppiJet_genDrMatch",VFLOAT);
    Branches.makeBranch("GenDptMatch","PuppiJet_genDptMatch",VFLOAT);
    Branches.makeBranch("GenTimeLLP","PuppiJet_genTimeLLP",VFLOAT);
    Branches.makeBranch("GenTOF","PuppiJet_genTOF",VFLOAT);

    Branches.makeBranch("GenQrkLlpDr","PuppiJet_genQrkLlpDr",VFLOAT);
    Branches.makeBranch("GenQrkLlpRp","PuppiJet_genQrkLlpRp",VFLOAT);
    Branches.makeBranch("GenQrkLlpId","PuppiJet_genQrkLlpId",VINT);

    Branches.makeBranch("GejJetDr","PuppiJet_genJetLlpDr",VFLOAT);
    Branches.makeBranch("GenJetRp","PuppiJet_genJetLlpRp",VFLOAT);
    Branches.makeBranch("GenJetLlpId","PuppiJet_genJetLlpId",VINT);

  }//<<>>if( cfFlag("hasGenInfo") )

  Branches.attachBranches(fOutTree);

}//<<>>void KUCMSAK4Jet::InitObject( TTree* fOutTree )

void KUCMSAK4PuppiJetObjectMini::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  fjets.clear();
  jetIds.clear();
  inputAvailable_ = false;
  if( !inputEnabled_ ) return;
  iEvent.getByToken(jetsToken_, jets_);
  if( !jets_.isValid() )
    throw cms::Exception("ProductNotFound") << "Configured puppiJets product is unavailable; use notAvailable when it is intentionally absent.";
  inputAvailable_ = true;
  iEvent.getByToken(pfcandToken_, pfcands_);

  if( AK4PuppiJetDEBUG ) std::cout << "Collecting AK4 PUPPI jets" << std::endl;
  for(const auto &jet : *jets_ ){ // Preserve all input jets in their original order.
    fjets.emplace_back(jet);
    auto jetID = 0;
    jetIds.push_back(jetID);

  }//<<>>for(const auto& jet : *jets_ )
  //std::sort(fjets.begin(),fjets.end(),sortByPt);


}//<<>>void KUCMSAK4Jet::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSAK4PuppiJetObjectMini::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSAK4PuppiJetObjectMini::ProcessEvent( ItemManager<float>& geVar ){

  if( AK4PuppiJetDEBUG ) std::cout << "Processing AK4Jets" << std::endl;

  Branches.clearBranches();
  Branches.fillBranch("Available",inputAvailable_);

  if( AK4PuppiJetDEBUG ) std::cout << " - enetering AK4Jet loop" << std::endl;

  vector<uInt> noIdGrp(1,0);
  // **** const reco::TrackRefVector& associatedTracks() const;
  // **** const reco::PixelClusterTagInfo* tagInfoPixelCluster(const std::string& label = "") const; // from b tagging info methods
  // **** reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }

  auto nJets = fjets.size();
  float jetHt = 0;
  if(AK4PuppiJetDEBUG ) std::cout << "Starting Jet Loop for " << nJets << " jets " << std::endl;
  for ( uInt ijet(0); ijet < nJets; ijet++ ){

    const auto & jet = fjets[ijet];
    //const auto jetid = fjetsID[ijet];
    if(AK4PuppiJetDEBUG ) std::cout << " --- Proccesssing : " << jet  << std::endl;

    jetHt += jet.pt();

    const float jetE = jet.energy();
    const float jetM = jet.mass();
    const float jetPt = jet.pt();
    const float jetEta = jet.eta();
    const float jetPhi = jet.phi();
    const float jetArea = jet.jetArea();
    const bool legacyCoverage = !(jet.pt() < cfPrm("jetPTmin")) &&
        !(cfFlag("onlyEB") && std::abs(jet.eta()) > cfPrm("ebMaxEta")) &&
        !(std::abs(jet.eta()) > cfPrm("jetEtaMax"));

    const float jetNHF = jet.neutralHadronEnergyFraction();
    const float jetNEMF = jet.neutralEmEnergyFraction();
    const float jetCHF = jet.chargedHadronEnergyFraction();
    const float jetCEMF = jet.chargedEmEnergyFraction();
    const float jetMUF = jet.muonEnergyFraction();
    const float jetNHM = jet.neutralMultiplicity();
    const float jetCHM = jet.chargedMultiplicity();

    Branches.fillBranch("E",jetE);
    Branches.fillBranch("M",jetM);
    Branches.fillBranch("Pt",jetPt);
    Branches.fillBranch("Eta",jetEta);
    Branches.fillBranch("Phi",jetPhi);
    Branches.fillBranch("Area",jetArea);
    Branches.fillBranch("LegacyCoverage",legacyCoverage);

    constexpr float unavailable = -999.f;
    float rawFactor = unavailable;
    float muonSubtrRawPt = unavailable, muonSubtrRawEta = unavailable, muonSubtrRawPhi = unavailable;
    float muonSubtrFactor = unavailable, muonSubtrDeltaEta = unavailable, muonSubtrDeltaPhi = unavailable;
    // Match PATJetCleanerForType1MET: start from the uncorrected jet and
    // subtract the p4 of global or standalone muon constituents.
    if( jet.jecSetsAvailable() ){
      try {
        const auto levels = jet.availableJECLevels();
        if( std::find(levels.begin(), levels.end(), "Uncorrected") != levels.end() ){
          const float rawScale = jet.jecFactor("Uncorrected");
          if( std::isfinite(rawScale) && rawScale > 0.f ){
            rawFactor = 1.f - rawScale;
            auto rawP4 = jet.p4() * rawScale;
            auto cleanedP4 = rawP4;
            bool constituentsAvailable = true;
            for( const auto& cand : jet.daughterPtrVector() ){
              if( cand.isNull() || !cand.isAvailable() ){ constituentsAvailable = false; break; }
              const auto* particle = cand.get();
              const auto* pfCandidate = dynamic_cast<const reco::PFCandidate*>(particle);
              const reco::Candidate* muon = particle;
              if( pfCandidate ){
                const auto& muonRef = pfCandidate->muonRef();
                if( muonRef.isNonnull() && !muonRef.isAvailable() ){ constituentsAvailable = false; break; }
                muon = muonRef.isNonnull() ? muonRef.get() : nullptr;
              }
              if( muon && (muon->isGlobalMuon() || muon->isStandAloneMuon()) ) cleanedP4 -= particle->p4();
            }
            const double rawPt = rawP4.pt();
            const double cleanedPt = cleanedP4.pt();
            if( constituentsAvailable && std::isfinite(rawPt) && rawPt > 0. &&
                std::isfinite(cleanedPt) && cleanedPt > 0. &&
                std::isfinite(cleanedP4.eta()) && std::isfinite(cleanedP4.phi()) ){
              muonSubtrRawPt = cleanedPt;
              muonSubtrRawEta = cleanedP4.eta();
              muonSubtrRawPhi = cleanedP4.phi();
              muonSubtrFactor = 1. - cleanedPt/rawPt;
              muonSubtrDeltaEta = cleanedP4.eta() - jetEta;
              muonSubtrDeltaPhi = cleanedP4.phi() - jetPhi;
            }
          }
        }
      } catch( const cms::Exception& ) {
        // Retain the documented sentinels when correction information is invalid.
      }
    }
    Branches.fillBranch("RawFactor",rawFactor);
    Branches.fillBranch("MuonSubtrRawPt",muonSubtrRawPt);
    Branches.fillBranch("MuonSubtrRawEta",muonSubtrRawEta);
    Branches.fillBranch("MuonSubtrRawPhi",muonSubtrRawPhi);
    Branches.fillBranch("MuonSubtrFactor",muonSubtrFactor);
    Branches.fillBranch("MuonSubtrDeltaEta",muonSubtrDeltaEta);
    Branches.fillBranch("MuonSubtrDeltaPhi",muonSubtrDeltaPhi);

    int idTight = -1, idTightLepVeto = -1;
    if( jet.isPFJet() && run3JetId_ ){
      run3PuppiJetId(jetEta, jetNHF, jetNEMF, jetCHF, jetCEMF, jetMUF,
                     jet.chargedMultiplicity(), jet.neutralMultiplicity(), idTight, idTightLepVeto);
    } else if( jet.isPFJet() && tightJetId && tightLepVetoJetId ){
      pat::strbitset tightBits = tightJetId->getBitTemplate();
      pat::strbitset vetoBits = tightLepVetoJetId->getBitTemplate();
      idTight = int((*tightJetId)(jet, tightBits));
      idTightLepVeto = int((*tightLepVetoJetId)(jet, vetoBits));
    }
    Branches.fillBranch("IdLoose",-1);
    Branches.fillBranch("IdTight",idTight);
    Branches.fillBranch("IdTightLepVeto",idTightLepVeto);
    // Preserve the embedded payload. Its training is assumed, not verified.
    float puDisc = -999.f;
    int puBitmap = -1;
    if( jet.hasUserFloat("pileupJetIdPuppi:fullDiscriminant") ){
      const float value = jet.userFloat("pileupJetIdPuppi:fullDiscriminant");
      if( std::isfinite(value) ) puDisc = value;
    }
    if( jet.hasUserInt("pileupJetIdPuppi:fullId") ) puBitmap = jet.userInt("pileupJetIdPuppi:fullId");
    // Keep an inclusive schema while PUPPI PU-ID decisions are disabled.
    Branches.fillBranch("PuIdLoose",-1);
    Branches.fillBranch("PuIdMedium",-1);
    Branches.fillBranch("PuIdTight",-1);
    Branches.fillBranch("PuIdDisc",puDisc);
    Branches.fillBranch("PuIdBitmap",puBitmap);

    Branches.fillBranch("NHF",jetNHF);
    Branches.fillBranch("NEMF",jetNEMF);
    Branches.fillBranch("CHF",jetCHF);
    Branches.fillBranch("CEMF",jetCEMF);
    Branches.fillBranch("MUF",jetMUF);
    Branches.fillBranch("NHM",jetNHM);
    Branches.fillBranch("CHM",jetCHM);

    // dR matched RH group  -----------------------------------------------
    // --------------------------------------------------------------------

    if(AK4PuppiJetDEBUG ) std::cout << "Getting jet dR rechit group" << std::endl;
    auto deltaRmin = 0.8;
    auto minRhE = 0.5;
    const auto jetDrRhGroup = rhObj->getRHGroup( jetEta, jetPhi, deltaRmin, minRhE );
    auto jetDrRhIdsGroup = rhObj->getRhGrpIDs( jetDrRhGroup );
    Branches.fillBranch("DrRhIds",jetDrRhIdsGroup);
    if( legacyCoverage ) rhObj->setRecHitUsed( jetDrRhIdsGroup );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // !!!!!!!!  create list of associated photons and electrons with jet via daughter particles of jet !!!!!!!!!!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const int nKids = jet.nConstituents();
    Branches.fillBranch("nKids",nKids);

    std::vector<uInt> jetParts;
    for( const auto &kid : jet.daughterPtrVector() ){

      // Only use packed-candidate keys when the daughter refers to this collection.
      if( kid.isNull() || kid.id() != pfcands_.id() || kid.key() >= pfcands_->size() ) continue;
      auto kidcand = (pfcands_->ptrAt(kid.key())).get();
      auto kideta = kidcand->eta();
      auto kidphi = kidcand->phi();
      auto idx = photonObj->getIndex( kideta, kidphi );
      if( idx >= 0 ){ jetParts.push_back(idx); }
      else {
        idx = electronObj->getIndex( kideta, kidphi );
        if( idx >= 0 ) jetParts.push_back(1000+idx);
      }//<<>>else

    }//<<>.for( const auto kid : jet.daughterPtrVector() )

    Branches.fillBranch("Parts",jetParts);

    // GenJet Info for MC  -------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------

    if(AK4PuppiJetDEBUG ) std::cout << "Getting jetGenParton Information" << std::endl;


    if( cfFlag("hasGenInfo") ){

      const int genJetIdx = genJetObj->savedIndex(jet);
      auto genJetInfo = genObjs->getGenJetInfo( jetEta, jetPhi, jetE, genJetIdx );
      Branches.fillBranch("GenJetIdx",genJetIdx);
      // load event level vectors for this jet with gen info

      Branches.fillBranch("GenImpactAngle",genJetInfo[0]);
      Branches.fillBranch("GenTime",genJetInfo[1]);
      Branches.fillBranch("GenPt",genJetInfo[2]);
      Branches.fillBranch("GenEta",genJetInfo[3]);
      Branches.fillBranch("GenPhi",genJetInfo[4]);
      Branches.fillBranch("GenEnergy",genJetInfo[5]);
      //Branches.fillBranch("genEMFrac",genJetInfo[6]);
      Branches.fillBranch("GenDrMatch",genJetInfo[7]);
      Branches.fillBranch("GenDptMatch",genJetInfo[8]);
      Branches.fillBranch("GenTimeLLP",genJetInfo[9]);
      Branches.fillBranch("GenTOF",genJetInfo[10]);
      Branches.fillBranch("GejJetDr",genJetInfo[11]);
      Branches.fillBranch("GenJetRp",genJetInfo[12]);
      Branches.fillBranch("GenJetLlpId",int(genJetInfo[13]));
      Branches.fillBranch("GenQrkLlpDr",genJetInfo[14]);
      Branches.fillBranch("GenQrkLlpRp",genJetInfo[15]);
      Branches.fillBranch("GenQrkLlpId",int(genJetInfo[16]));

    }//<<>>if( hasGenInfo )
    if(AK4PuppiJetDEBUG ) std::cout << "Next Jet .......................... " << std::endl;


  }//<<>>for ( uInt ijet = 0; ijet < nJets; ijet++ )
  // ** end of jets   ***************************************************************************************************

}//<<>>void KUCMSAK4Jet::ProcessEvent()

void KUCMSAK4PuppiJetObjectMini::EndJobs(){}

#endif
