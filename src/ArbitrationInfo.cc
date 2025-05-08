#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/ArbitrationInfo.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

ArbitrationInfo::ArbitrationInfo(const TransientVertex &vGold, const TransientVertex &vBronze) :
  vGold_(vGold), vBronze_(vBronze), commonTrack_(VertexHelper::GetCommonTrack(vGold, vBronze).track()) {
  
  const double chi2Gold(vGold.totalChiSquared()), chi2Bronze(vBronze.totalChiSquared());
  const double ndofGold(vGold.degreesOfFreedom()), ndofBronze(vBronze.degreesOfFreedom());
  const double normChi2Gold(vGold.normalisedChiSquared()), normChi2Bronze(vBronze.normalisedChiSquared());
  
  goldDecision_ = normChi2Gold;
  bronzeDecision_ = normChi2Bronze;

  auto calculateChiSquaredPValue = [](double chiSquaredValue, int degreesOfFreedom) {
				     // Define the chi-squared distribution with the specified degrees of freedom
				     boost::math::chi_squared chiSqDist(degreesOfFreedom);

				     // Compute the p-value as 1 - CDF(chiSquaredValue)
				     return 1.0 - boost::math::cdf(chiSqDist, chiSquaredValue);
				   };

  
  goldPVal_ = calculateChiSquaredPValue(chi2Gold, ndofGold);
  bronzePVal_ = calculateChiSquaredPValue(chi2Bronze, ndofBronze);

}

void ArbitrationInfo::setTrackCosTheta(const reco::Vertex &primaryVertex) {
    
  cosThetaGold_ = TrackHelper::CalculateCosTheta(primaryVertex, vGold_, commonTrack_);
  cosThetaBronze_ = TrackHelper::CalculateCosTheta(primaryVertex, vBronze_, commonTrack_);
  
}

MergeInfo::MergeInfo(const TransientVertex &vGold, const TransientVertex &vBronze, bool hasSignal) :
  hasSignal_(hasSignal),
  isValidMerge_(VertexHelper::Merge(vGold, vBronze).isValid()),
  distance_(VertexHelper::GetDistance(vGold, vBronze)),
  error_(VertexHelper::GetDistanceError(vGold, vBronze)),
  distanceSignificance_(VertexHelper::GetDistanceSignificance(vGold, vBronze)),
  vertex1NormChi2_(vGold.normalisedChiSquared()),
  vertex2NormChi2_(vBronze.normalisedChiSquared()),
  mergedNormChi2_(isValidMerge_? VertexHelper::Merge(vGold, vBronze).normalisedChiSquared() : -999.) {}
