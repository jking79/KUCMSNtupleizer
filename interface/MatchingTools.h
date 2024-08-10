#ifndef KUCMSNtupleizer_KUCMSNtupleizer_MatchingTools_h
#define KUCMSNtupleizer_KUCMSNtupleizer_MatchingTools_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

inline reco::GenParticleCollection CleanGenParticles(const reco::GenParticleCollection &genParticles, std::map<int, int> &indexMap) {

  reco::GenParticleCollection statusOne;
  for(size_t i = 0; i < genParticles.size(); i++) {

    if(genParticles[i].status() == 1 && genParticles[i].charge() != 0) {
      indexMap.emplace(statusOne.size(), i);
      statusOne.emplace_back(genParticles[i]);
    }
  }

  return statusOne;
}

template <typename T>
inline std::vector<T> RemoveDataAtIndices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete) {
  if(indicesToDelete.empty())
    return data;
  
  std::vector<T> ret;
  ret.reserve(data.size() - indicesToDelete.size());

  std::sort(indicesToDelete.begin(), indicesToDelete.end());

  typename std::vector<T>::const_iterator itBlockBegin = data.begin();
  for(std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++ it) {
    typename std::vector<T>::const_iterator itBlockEnd = data.begin() + *it;
    if(itBlockBegin != itBlockEnd) {
      std::copy(itBlockBegin, itBlockEnd, std::back_inserter(ret));
    }
    itBlockBegin = itBlockEnd + 1;
  }

  if(itBlockBegin != data.end()) {
    std::copy(itBlockBegin, data.end(), std::back_inserter(ret));
  }

  return ret;
}

#endif
