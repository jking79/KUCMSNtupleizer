#ifndef ElectronTools_hh
#define ElectronTools_hh

class TrackInfo {

 public:

  TrackInfo(const double &pt, const double &eta, const double &phi)
    : pt_(pt), eta_(eta), phi_(phi) {}

  // Getters
  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}

  bool isEmpty() {
    return pt_ < 0 && eta_ < 0 && phi_ < 0;
  }

  bool operator == (const TrackInfo& other) const {
    return pt_ == other.pt_ && eta_ == other.eta_ && phi_ == other.phi_;
  }

 private:

  double pt_, eta_, phi_;
};

template <typename T>
double GetDXY(T &object) {
  const double objectDx = object.vx();
  const double objectDy = object.vy();
  const double objectDxy = sqrt(objectDx*objectDx + objectDy*objectDy);
  return objectDxy;
}

#endif
