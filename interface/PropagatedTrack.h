#ifndef KUCMSNtupleizer_KUCMSNtupleizer_PropagatedTrack_h
#define KUCMSNtupleizer_KUCMSNtupleizer_PropagatedTrack_h

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"

//===============================================================================================//
//                                      class: PropagatedTrack                                   //
//-----------------------------------------------------------------------------------------------//
// Class for storing a track that was succesfully propagated to the ECAL with its index and its  //
// associated TrackDetInfo.                                                                      //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <class T>
class PropagatedTrack {

 public:

  PropagatedTrack(const unsigned int index, const T &track, const TrackDetMatchInfo &detInfo)
    : index_(index), track_(track), detInfo_(detInfo), isoInfo_(IsolationInfo()) {}

  virtual ~PropagatedTrack() = default;

  unsigned int GetIndex() const { return index_; }
  T GetTrack() const { return track_; }
  TrackDetMatchInfo GetDetInfo() const { return detInfo_; }
  IsolationInfo GetIsoInfo() const { return isoInfo_; }
  void SetIsoInfo(const IsolationInfo &isoInfo) { isoInfo_ = isoInfo; }

 private:

  const unsigned int index_;
  const T track_;
  const TrackDetMatchInfo detInfo_;
  IsolationInfo isoInfo_;
};

// Define alias for a vector of ProgagatedTrack objects
template <typename T> using PropagatedTracks = std::vector<PropagatedTrack<T> >;

#endif
