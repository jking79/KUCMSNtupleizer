#ifndef KUCMSNtupleizer_KUCMSNtupleizer_GenLeptonInfo_h
#define KUCMSNtupleizer_KUCMSNtupleizer_GenLeptonInfo_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenTools.h"

//enum LepMomType {kW, kZ, kTau, kConversion, kLight, kHeavy, kSusy, kOther, kUnmatched};

//===============================================================================================//                                                                               
//                                   class: GenLeptonInfo                                        //                                                                               
//-----------------------------------------------------------------------------------------------//                                                                               
// Class for storing results of a gen-matched lepton including its mother.                       //                                                                               
//-----------------------------------------------------------------------------------------------//                                                                               
//===============================================================================================//
class GenLeptonInfo {

 public:

  GenLeptonInfo()
    : index_(-1), genLepton_(reco::GenParticle()), type_(kUnmatched), deltaR_(-1.) {}

  GenLeptonInfo(const int index, const reco::GenParticle &genLepton, const LepMomType &type, const double deltaR = -1)
    : index_(index), genLepton_(genLepton), type_(type), deltaR_(deltaR) {}

  int GetIndex() const { return index_; }
  reco::GenParticle GetGenLepton() const { return genLepton_; }
  LepMomType GetLepMomType() const { return type_; }
  double GetDeltaR() const { return deltaR_; }
  bool isValid() const { return (index_ > 0 && deltaR_ > 0); }
  int GetPdgID() const {
    if(isValid())
      return int(genLepton_.pdgId());
    else
      return -999;}

 private:

  int index_;
  reco::GenParticle genLepton_;
  LepMomType type_;
  double deltaR_;
};

#endif
