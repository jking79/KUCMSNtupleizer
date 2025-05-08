#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DisplacedAVR.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoVertex/VertexTools/interface/DummyVertexSmoother.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"

#include <algorithm>

using namespace std;

TransientVertex DisplacedAVR::cleanUp(const TransientVertex& old) const {
  vector<reco::TransientTrack> trks = old.originalTracks();
  vector<reco::TransientTrack> newtrks;
  TransientVertex::TransientTrackToFloatMap mp;
  static const float minweight = 1.e-8;  // discard all tracks with lower weight
  //std::cout << "Track Cleaning" << std::endl;
  for (vector<reco::TransientTrack>::const_iterator i = trks.begin(); i != trks.end(); ++i) {
    //std::cout << "Track " << TrackHelper::FindTrackIndex(*i, trks) << " weight: " << old.trackWeight(*i) << std::endl;
    if (old.trackWeight(*i) > minweight) {
      newtrks.push_back(*i);
      mp[*i] = old.trackWeight(*i);
    }
  }

  TransientVertex ret;

  if (old.hasPrior()) {
    VertexState priorstate(old.priorPosition(), old.priorError());
    ret = TransientVertex(priorstate, old.vertexState(), newtrks, old.totalChiSquared(), old.degreesOfFreedom());
  } else {
    ret = TransientVertex(old.vertexState(), newtrks, old.totalChiSquared(), old.degreesOfFreedom());
  }
  ret.weightMap(mp);  // set weight map

  if (old.hasRefittedTracks()) {
    // we have refitted tracks -- copy them!
    vector<reco::TransientTrack> newrfs;
    vector<reco::TransientTrack> oldrfs = old.refittedTracks();
    vector<reco::TransientTrack>::const_iterator origtrkiter = trks.begin();
    for (vector<reco::TransientTrack>::const_iterator i = oldrfs.begin(); i != oldrfs.end(); ++i) {
      if (old.trackWeight(*origtrkiter) > minweight) {
        newrfs.push_back(*i);
      }
      origtrkiter++;
    }
    if (!newrfs.empty())
      ret.refittedTracks(newrfs);  // copy refitted tracks
  }

  if (ret.refittedTracks().size() > ret.originalTracks().size()) {
    edm::LogError("AdaptiveVertexReconstructor") << "More refitted tracks than original tracks!";
  }

  return ret;
}

void DisplacedAVR::erase(const TransientVertex& newvtx,
                                        set<reco::TransientTrack>& remainingtrks,
                                        float w) const {
  /*
   * Erase tracks that are in newvtx from remainingtrks 
   * But erase only if trackweight > w
   */
  const vector<reco::TransientTrack>& origtrks = newvtx.originalTracks();
  for (vector<reco::TransientTrack>::const_iterator i = origtrks.begin(); i != origtrks.end(); ++i) {
    double weight = newvtx.trackWeight(*i);
    if (weight > w) {
      remainingtrks.erase(*i);
    };
  };
}

DisplacedAVR::DisplacedAVR(float primcut, float seccut, float min_weight, bool smoothing)
  : thePrimaryFitter(nullptr), theSecondaryFitter(nullptr), theMinWeight(min_weight), theWeightThreshold(0.001) {
  setupFitters(primcut, 256., 0.25, seccut, 256., 0.25, smoothing);
}

DisplacedAVR::~DisplacedAVR() {
  if (thePrimaryFitter)
    delete thePrimaryFitter;
  if (theSecondaryFitter)
    delete theSecondaryFitter;
}

void DisplacedAVR::setupFitters(
    float primcut, float primT, float primr, float seccut, float secT, float secr, bool smoothing) {
  VertexSmoother<5>* smoother;
  if (smoothing) {
    smoother = new KalmanVertexSmoother();
  } else {
    smoother = new DummyVertexSmoother<5>();
  }

  if (thePrimaryFitter)
    delete thePrimaryFitter;
  if (theSecondaryFitter)
    delete theSecondaryFitter;

  /*
  edm::LogError ("AdaptiveVertexReconstructor" )
    << "Tini and r are hardcoded now!";
    */
  thePrimaryFitter = new AdaptiveVertexFitter(GeometricAnnealing(primcut, primT, primr),
                                              DefaultLinearizationPointFinder(),
                                              KalmanVertexUpdator<5>(),
                                              KalmanVertexTrackCompatibilityEstimator<5>(),
                                              *smoother);
  thePrimaryFitter->setWeightThreshold(theWeightThreshold);
  // if the primary fails, sth is wrong, so here we set a threshold on the weight.
  theSecondaryFitter = new AdaptiveVertexFitter(GeometricAnnealing(seccut, secT, secr),
                                                DefaultLinearizationPointFinder(),
                                                KalmanVertexUpdator<5>(),
                                                KalmanVertexTrackCompatibilityEstimator<5>(),
                                                *smoother);
  theSecondaryFitter->setWeightThreshold(0.);
  // need to set it or else we have
  // unwanted exceptions to deal with.
  // cleanup can come later!
  delete smoother;
}

DisplacedAVR::DisplacedAVR(const edm::ParameterSet& m)
    : thePrimaryFitter(nullptr), theSecondaryFitter(nullptr), theMinWeight(0.5), theWeightThreshold(0.001) {
  float primcut = 2.0;
  float seccut = 6.0;
  bool smoothing = false;
  // float primT = 4096.;
  // float primr = 0.125;
  float primT = 256.;
  float primr = 0.25;
  float secT = 256.;
  float secr = 0.25;

  try {
    primcut = m.getParameter<double>("primcut");
    primT = m.getParameter<double>("primT");
    primr = m.getParameter<double>("primr");
    seccut = m.getParameter<double>("seccut");
    secT = m.getParameter<double>("secT");
    secr = m.getParameter<double>("secr");
    theMinWeight = m.getParameter<double>("minweight");
    theWeightThreshold = m.getParameter<double>("weightthreshold");
    smoothing = m.getParameter<bool>("smoothing");
  } catch (edm::Exception& e) {
    edm::LogError("AdaptiveVertexReconstructor") << e.what();
  }

  setupFitters(primcut, primT, primr, seccut, secT, secr, smoothing);
}

vector<TransientVertex> DisplacedAVR::vertices(const vector<reco::TransientTrack>& t,
                                                              const reco::BeamSpot& s) const {
  return vertices(vector<reco::TransientTrack>(), t, s, false, true);
}

vector<TransientVertex> DisplacedAVR::vertices(const vector<reco::TransientTrack>& primaries,
                                                              const vector<reco::TransientTrack>& tracks,
                                                              const reco::BeamSpot& s) const {
  return vertices(primaries, tracks, s, true, true);
}

vector<TransientVertex> DisplacedAVR::vertices(const vector<reco::TransientTrack>& tracks) const {
  return vertices(vector<reco::TransientTrack>(), tracks, reco::BeamSpot(), false, false);
}

vector<TransientVertex> DisplacedAVR::vertices(const vector<reco::TransientTrack>& primaries,
					       const vector<reco::TransientTrack>& tracks,
					       const reco::BeamSpot& s,
					       bool has_primaries,
					       bool usespot) const {

  vector<TransientVertex> ret;
  set<reco::TransientTrack> remainingtrks;
  vector<reco::TransientTrack> rejectTracks;
  
  //std::cout << "\nWe start with " << tracks.size() << " tracks" << std::endl;
  
  copy(tracks.begin(), tracks.end(), inserter(remainingtrks, remainingtrks.begin()));

  unsigned int n_tracks = remainingtrks.size();

  // cout << "[AdaptiveVertexReconstructor] DEBUG ::vertices!!" << endl;
  try {
    while (remainingtrks.size() > 1) {
  
      //std::cout << "  There are " << remainingtrks.size() << " remaining tracks: " << std::endl; 
      //for(const auto &track : remainingtrks)
      //std::cout << TrackHelper::FindTrackIndex(track, tracks) << ", ";
      //std::cout << std::endl;
  
      updateRejectTracks(rejectTracks, tracks, remainingtrks);
      //std::cout << "Rejected tracks: " << rejectTracks.size() << std::endl;
      
      const AdaptiveVertexFitter* fitter = theSecondaryFitter;
      if (ret.empty()) {
        fitter = thePrimaryFitter;
      };
      vector<reco::TransientTrack> fittrks;
      fittrks.reserve(remainingtrks.size());

      copy(remainingtrks.begin(), remainingtrks.end(), back_inserter(fittrks));

      TransientVertex tmpvtx;
      if ((ret.empty()) && has_primaries) {
        // add the primaries to the fitted tracks.
        copy(primaries.begin(), primaries.end(), back_inserter(fittrks));
      }
      if ((ret.empty()) && usespot) {
        tmpvtx = fitter->vertex(fittrks, s);
      } else {
        tmpvtx = fitter->vertex(fittrks);
      }
      TransientVertex newvtx = cleanUp(tmpvtx);
      ret.push_back(newvtx);
      erase(newvtx, remainingtrks, theMinWeight);
      if (n_tracks == remainingtrks.size()) {
        if (usespot) {
          // try once more without beamspot constraint!
          usespot = false;
          LogDebug("AdaptiveVertexReconstructor") << "no tracks in vertex. trying again without beamspot constraint!";
          continue;
        }
        LogDebug("AdaptiveVertexReconstructor")
            << "all tracks (" << n_tracks << ") would be recycled for next fit. Trying with low threshold!";
        erase(newvtx, remainingtrks, 1.e-5);
        if (n_tracks == remainingtrks.size()) {
          LogDebug("AdaptiveVertexReconstructor") << "low threshold didnt help! "
                                                  << "Discontinue procedure!";
          break;
        }
      };

      // cout << "[AdaptiveVertexReconstructor] erased" << endl;
      n_tracks = remainingtrks.size();
    };
  } catch (VertexException& v) {
    // Will catch all (not enough significant tracks exceptions.
    // in this case, the iteration can safely terminate.
  };

  return cleanUpVertices(ret);
}

vector<TransientVertex> DisplacedAVR::cleanUpVertices(const vector<TransientVertex>& old) const {
  vector<TransientVertex> ret;
  for (vector<TransientVertex>::const_iterator i = old.begin(); i != old.end(); ++i) {
    if (!(i->hasTrackWeight())) {  // if we dont have track weights, we take the vtx
      ret.push_back(*i);
      std::cout << "This happens!" << std::endl;
      continue;
    }


    // maybe this should be replaced with asking for the ndf ...
    // e.g. if ( ndf > - 1. )
    int nsig = 0;  // number of significant tracks.
    //if(i->degreesOfFreedom() > -1 && (i->originalTracks().size() > 1 || i->refittedTracks().size() > 1))
    //ret.push_back(*i);
    
    TransientVertex::TransientTrackToFloatMap wm = i->weightMap();
    for (TransientVertex::TransientTrackToFloatMap::const_iterator w = wm.begin(); w != wm.end(); ++w) {
      if (w->second > theWeightThreshold)
        nsig++;
    }
    if (nsig > 1)
      ret.push_back(*i);
    
  }
  
  //for (vector<TransientVertex>::const_iterator i = old.begin(); i != old.end(); ++i)
  //if(i->originalTracks().size() > 2 || i->refittedTracks().size() > 2)
  /*
  std::cout << "There were " << old.size() <<  " vertices before cleaning up and " << ret.size() << " after." << std::endl;

  std::cout << "Old vertex info: " << std::endl;
  for(const auto &vertex : old) {
    if(vertex.isValid()) std::cout << "\tThe vertex is valid" << std::endl;
    std::cout << "\tndof: " << vertex.degreesOfFreedom() << std::endl;
    std::cout << "\tchi2: " << vertex.totalChiSquared() << std::endl;
     std::cout << "\tnormChi2: " << vertex.normalisedChiSquared() << std::endl;
    std::cout << "\tThis vertex has " << vertex.originalTracks().size() << " tracks and " << vertex.refittedTracks().size() << " refitted tracks." << std::endl;
    if(vertex.originalTracks().size() > 0) {
      TransientVertex::TransientTrackToFloatMap wm = vertex.weightMap();
      for(TransientVertex::TransientTrackToFloatMap::const_iterator w = wm.begin(); w != wm.end(); ++w)
	std::cout << "\t\tTrack weight: " << w->second << std::endl;
    }
  }
  
  std::cout << "End of AVR" << std::endl;
  */
  return ret;
}

// Function to update rejectTrackIndices based on removed tracks
void DisplacedAVR::updateRejectTracks(
    std::vector<reco::TransientTrack>& rejectTracks,   // Removed tracks
    const std::vector<reco::TransientTrack>& originalTracks, // Original set of tracks
    const std::set<reco::TransientTrack>& remainingTracks    // Remaining set of tracks
) const {
  // Iterate through the original track list and check if each track exists in remainingTracks
  for(const auto& track : remainingTracks) {
    if(remainingTracks.find(track) == remainingTracks.end()) {
      rejectTracks.push_back(track); // Add it to the rejectTracks vector
    }
  }
}
