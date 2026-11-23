// MiniAOD AK4 generator jets and PAT jet association lookup.
#ifndef KUCMSGenJetObjectMiniHeader
#define KUCMSGenJetObjectMiniHeader

#include <vector>

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "KUCMSObjectBase.hh"

class KUCMSGenJetObjectMini : public KUCMSObjectBase {
public:
  explicit KUCMSGenJetObjectMini(const edm::ParameterSet&) {}

  void LoadGenJetsTokens(edm::EDGetTokenT<std::vector<reco::GenJet>> token) { genJetsToken_ = token; }
  void InitObject(TTree* tree) override;
  void LoadEvent(const edm::Event& event, const edm::EventSetup&, ItemManager<float>&) override;
  void ProcessEvent(ItemManager<float>&) override;

  // -1: no PAT association; -9: associated jet cannot be mapped to this collection.
  int savedIndex(const pat::Jet& jet) const;

private:
  edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
  edm::Handle<std::vector<reco::GenJet>> genJets_;
};

void KUCMSGenJetObjectMini::InitObject(TTree* tree) {
  Branches.makeBranch("Pt", "GenJet_pt", VFLOAT);
  Branches.makeBranch("Eta", "GenJet_eta", VFLOAT);
  Branches.makeBranch("Phi", "GenJet_phi", VFLOAT);
  Branches.makeBranch("Energy", "GenJet_energy", VFLOAT);
  Branches.makeBranch("Mass", "GenJet_mass", VFLOAT);
  Branches.attachBranches(tree);
}

void KUCMSGenJetObjectMini::LoadEvent(const edm::Event& event, const edm::EventSetup&, ItemManager<float>&) {
  event.getByToken(genJetsToken_, genJets_);
  if (!genJets_.isValid()) {
    throw cms::Exception("ProductNotFound") << "Configured genjets collection is unavailable with hasGenInfo enabled.";
  }
}

void KUCMSGenJetObjectMini::ProcessEvent(ItemManager<float>&) {
  Branches.clearBranches();
  for (const auto& jet : *genJets_) {
    Branches.fillBranch("Pt", float(jet.pt()));
    Branches.fillBranch("Eta", float(jet.eta()));
    Branches.fillBranch("Phi", float(jet.phi()));
    Branches.fillBranch("Energy", float(jet.energy()));
    Branches.fillBranch("Mass", float(jet.mass()));
  }
}

int KUCMSGenJetObjectMini::savedIndex(const pat::Jet& jet) const {
  const auto& association = jet.genJetFwdRef();
  if (association.isNull()) {
    // Older PAT jets can hold an embedded or legacy reference instead.
    try {
      return jet.genJet() ? -9 : -1;
    } catch (const cms::Exception&) {
      return -9;
    }
  }

  const auto& forward = association.ref();
  const auto& back = association.backRef();
  if (forward.isNonnull() && forward.id() == genJets_.id() && forward.key() < genJets_->size())
    return static_cast<int>(forward.key());
  if (back.isNonnull() && back.id() == genJets_.id() && back.key() < genJets_->size())
    return static_cast<int>(back.key());
  return -9;
}

#endif
