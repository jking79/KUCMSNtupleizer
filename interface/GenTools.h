#ifndef KUCMSNtupleizer_KUCMSNtupleizer_GenTools_h
#define KUCMSNtupleizer_KUCMSNtupleizer_GenTools_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

enum LepMomType {kW, kZ, kTau, kConversion, kLight, kHeavy, kSusy, kOther, kUnmatched};

inline std::vector<int> MomIDs(const reco::GenParticle &genLepton) {

  auto mother = genLepton.mother();

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

inline LepMomType ClassifyGenLeptonMomType(const reco::GenParticle &genLepton) {

  std::vector<int> motherIDs(MomIDs(genLepton));

  LepMomType momType = kUnmatched;
  for(auto const& id : motherIDs) {
    momType = AssignGenLeptonMomType(id);

    if(momType != kUnmatched)
      break;
  }
  return momType;
}

inline bool isSignalGenElectron(const reco::GenParticle &genLepton) {
  LepMomType momType = ClassifyGenLeptonMomType(genLepton);
  return (momType == kZ || momType == kSusy);
}

#endif
