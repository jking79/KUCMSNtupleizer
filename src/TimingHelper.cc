#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TimingHelper.h"

double TimingHelper::PathLength(const reco::TransientTrack &track, const GlobalPoint &vertexPosition, const GlobalPoint &positionAtECAL) {

  //Calculate path length with helical extrapolation
  SteppingHelixPropagator propagator(track.field());

  //Calculate path length using the stepping helix propagator
  double pathLength = propagator.propagateWithPath(track.initialFreeState(), positionAtECAL).second;

  return pathLength;
}

double TimingHelper::Beta(const reco::Track &track, const double mass) {
  return track.p() / sqrt(track.p()*track.p() + mass*mass);
}

double TimingHelper::Time(const reco::Track &track, const double mass, const double pathLength) {
  return pathLength/( Beta(track, mass) * 30);
}

reco::Track TimingHelper::TimeStampTrack(const reco::Track &track, const double pathLength) { 

  // Track parameters for creating new track instance with path length info                                                                                                       
  const double chi2 = track.chi2();
  const int charge = track.charge();
  const double ndof = track.ndof();
  const math::XYZPoint trackRefPoint = track.referencePoint();
  const math::XYZVector trackMomentum = track.momentum();

  // Declare different mass hypotheses (pion, kaon, proton and electron) and calculate their tofs
  const double mass_pi(0.13957018); const double tof_pi(Time(track, mass_pi, pathLength));
  const	double mass_k(0.493677); const double tof_k(Time(track, mass_k, pathLength));
  const double mass_p(0.9382720813); const double tof_p(Time(track, mass_p, pathLength));
  const	double mass_e(0.000511); const double tof_e(Time(track, mass_e, pathLength));
  
  // Create track and save the tofs of the four different mass hypotheses
  const reco::Track updatedTrack(chi2, ndof, trackRefPoint, trackMomentum, charge, track.covariance(),
                                 //reco::TrackBase::undefAlgorithm, reco::TrackBase::undefQuality,
				 track.algo(), reco::TrackBase::undefQuality,
                                 tof_pi, tof_k, tof_p, tof_e);

  return updatedTrack;
}
