#ifndef KUCMSNtupleizer_KUCMSNtupleizer_GenTools_h
#define KUCMSNtupleizer_KUCMSNtupleizer_GenTools_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

enum LepMomType {kW, kZ, kTau, kConversion, kLight, kHeavy, kSusy, kOther, kUnmatched};

inline std::vector<int> MomIDs(const reco::GenParticle &genLepton) {

  auto mother = genLepton.mother();

  if(!mother)
    return std::vector<int>();

  std::vector<int> motherIDs;
  while(mother->pt() > 0) {
    const int motherID = mother->pdgId();

    if(auto *mom = mother->mother())
      mother = mom;
    else break;
    
    if(motherID == mother->pdgId())
      continue;

    motherIDs.push_back(motherID);
  }

  return motherIDs;
}

inline LepMomType AssignGenLeptonMomType(const int motherID) {

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

inline int GenParticleMomPdgID(const reco::GenParticle &genLepton) {

  int pdgID = -1;
  //std::cout << "here" << std::endl;
  std::vector<int> motherIDs(MomIDs(genLepton));

  LepMomType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignGenLeptonMomType(id);

    if(momType != kUnmatched) {
      pdgID = id;
      break;
     }
  }
  return pdgID;
}

inline LepMomType ClassifyGenLeptonMomType(const reco::GenParticle &genLepton) {
  int momID(GenParticleMomPdgID(genLepton));
  return AssignGenLeptonMomType(momID);
}

inline bool isSignalGenElectron(const reco::GenParticle &gen) {
  if(abs(gen.pdgId()) != 11 || gen.status() != 1) return false;
  LepMomType momType = ClassifyGenLeptonMomType(gen);
  return (momType == kZ || momType == kSusy);
}

inline bool isSignalGenMuon(const reco::GenParticle &gen) {
  if(abs(gen.pdgId()) != 13 || gen.status() != 1) return false;
  LepMomType momType = ClassifyGenLeptonMomType(gen);
  return (momType == kZ || momType == kSusy);
}

inline bool isSignalGenJet(const reco::GenParticle &gen) {
  if(abs(gen.pdgId()) > 6 || gen.status() != 23) return false;
  LepMomType momType = ClassifyGenLeptonMomType(gen);
  return (momType == kZ || momType == kSusy);
}

#endif
