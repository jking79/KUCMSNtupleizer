#ifndef KUCMSNtupleizer_KUCMSNtupleizer_ArbitrationInfo_h
#define KUCMSNtupleizer_KUCMSNtupleizer_ArbitrationInfo_h

#include <boost/math/distributions/chi_squared.hpp>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

class ArbitrationInfo {

 public:
  ArbitrationInfo(const TransientVertex &vGold, const TransientVertex &vBronze);
  virtual ~ArbitrationInfo() = default;

  double goldDecision() const {return goldDecision_;}
  double bronzeDecision() const {return bronzeDecision_;}
  double goldPValue() const {return goldPVal_;}
  double bronzePValue() const {return bronzePVal_;}
  double cosThetaGold() const {return cosThetaGold_;}
  double cosThetaBronze() const {return cosThetaBronze_;}
  
  void setTrackCosTheta(const reco::Vertex &primaryVertex);
  
  
 private:
  TransientVertex vGold_, vBronze_;
  reco::Track commonTrack_;
  
  double goldDecision_, bronzeDecision_, goldPVal_, bronzePVal_;
  double cosThetaGold_, cosThetaBronze_;

};

class MergeInfo {

public:
  MergeInfo(const TransientVertex &vGold, const TransientVertex &vBronze, bool hasSignal = false);

  virtual ~MergeInfo() = default;
  
  double distance() const { return distance_; }
  double error() const { return error_; }
  double distanceSignificance() const { return distanceSignificance_; }
  double vertex1NormChi2() const { return vertex1NormChi2_; }
  double vertex2NormChi2() const { return vertex2NormChi2_; }
  double mergedNormChi2() const { return mergedNormChi2_; }

  bool isValidMerge() const { return isValidMerge_; }
  bool hasSignal() const { return hasSignal_; }

  
private:
  bool hasSignal_, isValidMerge_;
  double distance_, error_, distanceSignificance_;
  double vertex1NormChi2_, vertex2NormChi2_, mergedNormChi2_;
};

/*
class InitialSeedInfo {

public:
  InitialSeedInfo(const TransientVertex &vertex, const int nSignal = 0) :
    isBronze_(nSignal == 1), isGold_(nSignal == 2) {}

  bool isBronze() const { return isBronze_; }
  bool isGold() const { return isGold_; }

  double cosTheta(const reco::Vertex &primaryVertex) const;
  double bronzeSignalCosTheta(const reco::Vertex &primaryVertex) const;
  
  
private:
  bool isBronze_, isGold_;
}
*/
#endif
