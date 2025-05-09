#ifndef KUCMSNtupleizer_KUCMSNtupleizer_TTBuilderWrapper_h
#define KUCMSNtupleizer_KUCMSNtupleizer_TTBuilderWrapper_h

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

class TTBuilderWrapper {

 public:

 TTBuilderWrapper() : ttBuilder_(nullptr) {}
  
 TTBuilderWrapper(const TransientTrackBuilder* ttBuilder) :
  ttBuilder_(ttBuilder) {}

  reco::TransientTrack build(const reco::Track &track) const { return ttBuilder_->build(track); }

  const TransientTrackBuilder* getBuilder() const {return ttBuilder_;}
  
 private:
  const TransientTrackBuilder* ttBuilder_;
};

#endif
