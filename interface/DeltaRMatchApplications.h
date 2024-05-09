#ifndef KUCMSNtupleizer_KUCMSNtupleizer_DeltaRMatchApplications_h
#define KUCMSNtupleizer_KUCMSNtupleizer_DeltaRMatchApplications_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenTools.h"

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                               Templated Function: RemoveOverlap                               //
//-----------------------------------------------------------------------------------------------//
// Uses DeltaRMatchHungarian to remove the overlap between two sets of objects. The list that it //
// returns is the reduced version of the set given to the toReduce argument of the function. The //
// objects remove are the ones that have a DR below the specified deltaRcut.                     //
//===============================================================================================//
template <typename T, typename S>
  inline std::vector<T> RemoveOverlap(const std::vector<T> &toReduce, const std::vector<S> &toMatch, const double &deltaRcut) {

  DeltaRMatchHungarian<T, S> assigner(toReduce, toMatch);

  std::vector<size_t> matchedIndexes;

  for(auto const &pair : assigner.GetPairedObjects()) {
    if(pair.GetDeltaR() < deltaRcut) {
      matchedIndexes.push_back(size_t(pair.GetIndexA()));
    }
  }

  return RemoveDataAtIndices<T>(toReduce, matchedIndexes);
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                                 Templated Function: GetOverlap                                //
//-----------------------------------------------------------------------------------------------//
// In contrast to the RemoveOverlap function, this one returns the pairs of matched objects that //
// that are below the specified deltaRcut. The return type is PairedObjects class.               //
//===============================================================================================//
template <typename T, typename S>
  inline std::vector<PairedObjects<T, S>> GetOverlap(const std::vector<T> &setA, const std::vector<S> &setB, const double &deltaRcut) {

  DeltaRMatchHungarian<T, S> assigner(setA, setB);

  return assigner.GetPairedObjects();
}

//===============================================================================================//
//-----------------------------------------------------------------------------------------------//
//                             Templated Function: GetSignalElectrons                            //
//-----------------------------------------------------------------------------------------------//
// Matches a set of reco::Electron objects to a set of reco::GenParticle and saves only the ones //
// that correspond to electrons decaying from a Z boson or directly from a neutralino (GMSB).    //
//===============================================================================================//
inline reco::ElectronCollection GetSignalElectrons(const reco::ElectronCollection &electrons, const reco::GenParticleCollection &genParticles) {

  reco::ElectronCollection signalElectrons;

  DeltaRGenMatchHungarian<reco::Electron> genAssigner(electrons, genParticles);

  for(const auto &pair : genAssigner.GetPairedObjects()) 
    if(isSignalGenElectron(pair.GetObjectB()))
       signalElectrons.emplace_back(pair.GetObjectA());

  return signalElectrons;  
}

template <typename T>
inline std::vector<T> GetSignalMatches(const std::vector<T> &toMatch, const reco::GenParticleCollection &genParticles) {

  std::vector<T> signalMatches;
  DeltaRGenMatchHungarian<T> genAssigner(toMatch, genParticles);

  for(const auto &pair : genAssigner.GetPairedObjects()) {
    const reco::GenParticle genParticle(pair.GetObjectB());
    if(abs(pair.GetObjectB().pdgId()) == 11) {
      if(isSignalGenElectron(pair.GetObjectB())) {
	signalMatches.emplace_back(pair.GetObjectA());
	//std::cout << "Found a signal match in the function with deltaR: "<< pair.GetDeltaR() << std::endl;
	//std::cout << "gen pt: " << pair.GetObjectB().pt() << std::endl;
	//std::cout << "track pt: "<< pair.GetObjectA().pt() <<std::endl;
	//std::cout << "track quality: "<< pair.GetObjectA().qualityMask() <<std::endl;
      }
    }
  }
  return signalMatches;
} 

#endif
