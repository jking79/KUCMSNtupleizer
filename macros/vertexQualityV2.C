//c++ includes
#include <vector>

// ROOT includes
#include <TH1D.h>
#include <TH2D.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TLegend.h>
#include <Math/Vector4D.h>
#include <TEfficiency.h>

// local includes
//#include "DeltaRMatch.h"

// nTuple class
#include "vtxtree.h"
//#include "../interface/TrackTools.h"

enum particleID {kPhoton, kPion, kProtons, kNeutrons, kKaon, kKlong, kMuon, kElectrons, kSusy, kNeutrino, kUnmatched};
particleID ClassifyParticle(const int pdgID);
template <typename T>
void FillIDHists(std::vector<TH1D>& histID, const int pdgID, const T quantity);
std::string GetParticleTypeString(particleID id);

class DeltaRMatchContainer {

 public:

  DeltaRMatchContainer(const double &pt, const double &eta, const double &phi)
    : pt_(pt), eta_(eta), phi_(phi) {}

  // Getters
  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}

 private:

  double pt_, eta_, phi_;
};

typedef ROOT::Math::PtEtaPhiMVector LorentzVec;

int GenMatchIndex(const LorentzVec &trackVec, const std::vector<LorentzVec> &genVecs, float &deltaR);
float GetDXY(const float x, const float y);
//std::vector<LorentzVec> GetLorentzVecCollection(const llpgtree &tree);
//std::vector<LorentzVec> KeepOnlySS(const int recoCharge, const llpgtree &tree);
TCanvas* plotHistogram(TH1D& hist1D, const std::string& xlabel, const std::string& ylabel = "Events");
TCanvas* plotHistogram(TH2D& hist2D, const std::string& xlabel, const std::string& ylabel, const std::string& zlabel = "Events");
TCanvas* plotHistograms(TH1D &signal, TH1D &background, const std::string& name, const std::string& xlabel, const std::string& ylabel = "Events");
TCanvas* plotHistograms(std::vector<TH1D> &histograms, const std::vector<std::string> &labels, const std::string name, std::string xlabel);
TCanvas* plotEfficiency(std::vector<TEfficiency*> &efficiencies, std::vector<string> &labels, string name, string xlabel);

void vertexQualityV2() {

  using namespace std;

  std::vector<std::string> particleTypeNames = {
    "Photon",
    "Pion",
    "Proton",
    "Neutron",
    "Kaon",
    "K-long",
    "Muon",
    "Electron",
    "SUSY-Particle",
    "Neutrino",
    "Unmatched"
  };

  // Histograms
  int nbins = 100;
  double xmin = 0.;
  double xmax = 10;

  // Gen Particle
  TH1D genParticle_sigPt("genParticle_sigPt", "genParticle_sigPt", 100, 0., 500);
  TH1D genParticle_sigDeltaR("genParticle_sigDeltaR", "genParticle_sigDeltaR", 100, 0., 1.);
  TH1D genParticle_dxy("genParticle_dxy", "genParticle_dxy", 100, 0., 10.);

  TH1D track_ptSigMatch("track_ptSigMatch", "track_ptSigMatch", 100, 0., 500);
  TH1D track_ptSigMatchPV("track_ptSigMatchPV", "track_ptSigMatchPV", 100, 0., 500);
  TH1D track_ptBkgMatch("track_ptBkgMatch", "track_ptBkgMatch", 100, 0., 500);
  TH1D track_ptBkgMatchPV("track_ptBkgMatchPV", "track_ptBkgMatchPV", 100, 0., 500);

  TH1D track_etaSigMatch("track_etaSigMatch", "track_ptSigMatch", 100, -4., 4.);
  TH1D track_etaSigMatchPV("track_etaSigMatchPV", "track_etaSigMatchPV", 100, -4., 4.);
  TH1D track_etaBkgMatch("track_etaBkgMatch", "track_etaBkgMatch", 100, -4., 4.);
  TH1D track_etaBkgMatchPV("track_etaBkgMatchPV", "track_etaBkgMatchPV", 100, -4., 4.);

  TH1D track_dxySigMatch("track_dxySigMatch", "track_dxySigMatch", 70, 0., 70.);
  TH1D track_dxySigMatchPV("track_dxySigMatchPV", "track_dxySigMatchPV", 70, 0., 70.);
  TH1D track_dxyBkgMatch("track_dxyBkgMatch", "track_dxyBkgMatch", 70, 0., 70.);
  TH1D track_dxyBkgMatchPV("track_dxyBkgMatchPV", "track_dxyBkgMatchPV", 70, 0., 70.);

  const double deltaRange(0.1);
  TH1D track_pvDeltaZSigMatch("track_pvDeltaZSigMatch", "track_pvDeltaZSigMatch", 100, 0., deltaRange);
  TH1D track_pvDeltaZSigMatchPV("track_pvDeltaZSigMatchPV", "track_pvDeltaZSigMatchPV", 100, 0., deltaRange);
  TH1D track_pvDeltaZBkgMatch("track_pvDeltaZBkgMatch", "track_pvDeltaZBkgMatch", 100, 0., deltaRange);
  TH1D track_pvDeltaZBkgMatchPV("track_pvDeltaZBkgMatchPV", "track_pvDeltaZBkgMatchPV", 100, 0., deltaRange);

  TH1D track_pvDelta3DSigMatch("track_pvDelta3DSigMatch", "track_pvDelta3DSigMatch", 100, 0., deltaRange);
  TH1D track_pvDelta3DSigMatchPV("track_pvDelta3DSigMatchPV", "track_pvDelta3DSigMatchPV", 100, 0., deltaRange);
  TH1D track_pvDelta3DBkgMatch("track_pvDelta3DBkgMatch", "track_pvDelta3DBkgMatch", 100, 0., deltaRange);
  TH1D track_pvDelta3DBkgMatchPV("track_pvDelta3DBkgMatchPV", "track_pvDelta3DBkgMatchPV", 100, 0., deltaRange);
  
  TH1D track_normChi2SigMatch("track_normChi2SigMatch", "track_normChi2SigMatch", 100, 0., 4.);
  TH1D track_normChi2SigMatchPV("track_normChi2SigMatchPV", "track_normChi2SigMatchPV", 100, 0., 4.);
  TH1D track_normChi2BkgMatch("track_normChi2BkgMatch", "track_normChi2BkgMatch", 100, 0., 4.);
  TH1D track_normChi2BkgMatchPV("track_normChi2BkgMatchPV", "track_normChi2BkgMatchPV", 100, 0., 4.);

  //TH1D track_ptSigMatchBelow("track_ptSigMatchBelow", "track_ptSigMatchBelow", 100, 0., 500);
  //TH1D track_ptSigMatchAbove("track_ptSigMatchAbove", "track_ptSigMatchAbove", 100, 0., 500);
  TH1D track_ptBkgMatchBelow("track_ptBkgMatchBelow", "track_ptBkgMatchBelow", 100, 0., 500);
  TH1D track_ptBkgMatchAbove("track_ptBkgMatchAbove", "track_ptBkgMatchAbove", 100, 0., 500);

  TH1D track_isECALsigOutPV("track_isECALsigOutPV", "track_isECALsigOutPV", 3, -1, 2);
  TH1D track_isECALbkgOutPV("track_isECALbkgOutPV", "track_isECALbkgOutPV", 3, -1, 2);
  TH1D track_isECALsigInPV("track_isECALsigInPV", "track_isECALsigInPV", 3, -1, 2);
  TH1D track_isECALbkgInPV("track_isECALbkgInPV", "track_isECALbkgInPV", 3, -1, 2);

  // Tracks not in PV that make it to the ECAL 
  TH1D track_ptSig("track_ptSig", "track_ptSig", 100, 0., 500.);
  TH1D track_ptBkg("track_ptBkg", "track_ptBkg", 100, 0., 500.);

  TH1D track_etaSig("track_etaSig", "track_etaSig", 100, -4., 4.);
  TH1D track_etaBkg("track_etaBkg", "track_etaBkg", 100, -4., 4.);

  TH1D track_phiSig("track_phiSig", "track_phiSig", 100, -4., 4.);
  TH1D track_phiBkg("track_phiBkg", "track_phiBkg", 100, -4., 4.);

  TH1D track_chargeSig("track_chargeSig", "track_chargeSig", 3, -1, 2);
  TH1D track_chargeBkg("track_chargeBkg", "track_chargeBkg", 3, -1, 2);

  TH1D track_normChi2Sig("track_normChi2Sig", "track_normChi2Sig", 100, 0., 4.);
  TH1D track_normChi2Bkg("track_normChi2Bkg", "track_normChi2Bkg", 100, 0., 4.);

  TH1D track_dxySig("track_dxySig", "track_dxySig", 120, 0., 60.);
  TH1D track_dxyBkg("track_dxyBkg", "track_dxyBkg", 120, 0., 60.);

  TH1D track_ip2DSig("track_ip2DSig", "track_ip2DSig", 100, 0., 20.);
  TH1D track_ip2DBkg("track_ip2DBkg", "track_ip2DBkg", 100, 0., 20.);

  TH1D track_sip2DSig("track_sip2DSig", "track_sip2DSig", 100, -20., 20.);
  TH1D track_sip2DBkg("track_sip2DBkg", "track_sip2DBkg", 100, -20., 20.);

  TH1D track_ip3DSig("track_ip3DSig", "track_ip3DSig", 100, 0., 20.);
  TH1D track_ip3DBkg("track_ip3DBkg", "track_ip3DBkg", 100, 0., 20.);
  
  TH1D track_sip3DSig("track_sip3DSig", "track_sip3DSig", 100, -20., 20.);
  TH1D track_sip3DBkg("track_sip3DBkg", "track_sip3DBkg", 100, -20., 20.);
  
  TH2D genSignal_DxyVsTrackWeight("genSignal_DxyVsTrackWeight", "genSignal_DxyVsTrackWeight", 100, -1., 1., 100, 0., 1.);

  TH2D genBkg_DxyVsTrackWeight("genBkg_DxyVsTrackWeight", "genBkg_DxyVsTrackWeight", 100, -1., 1., 100, 0., 10.);
  TH2D genParticle_PdgIdVsDxy("genParticle_PdgIdVsDxy", "genParticle_PdgIdVsDxy", 100, 0., 20, 350, 0., 350.);
  TH2D genParticle_PtVsDxy("genParticle_PtVsDxy", "genParticle_PtVsDxy", 150, 0., 15, 300, 0., 300.);
  TH2D genParticle_EtaVsDxy("genParticle_EtaVsDxy", "genParticle_EtaVsDxy", 150, 0., 15, 100, -4., 4.);
  TH2D genParticle_normChi2VsDxy("genParticle_normChi2VsDxy", "genParticle_normChi2VsDxy", 150, 0., 15, 100, 0.,2.);

  TH2D genParticle_PtVsTrackWeight("genParticle_PtVsTrackWeight", "genParticle_PtVsTrackWeight", 100, -1., 1., 300, 0., 300.);
  TH2D genParticle_EtaVsTrackWeight("genParticle_EtaVsTrackWeight", "genParticle_EtaVsTrackWeight", 100, -1., 1., 100, -4., 4.);
  TH2D genParticle_normChi2VsTrackWeight("genParticle_normChi2VsTrackWeight", "genParticle_normChi2VsTrackWeight", 100, -1., 1., 100, 0.,4.);

  TH2D track_pvDelta3DVsTrackWeightSig("track_pvDelta3DVsTrackWeightSig", "track_pvDelta3DVsTrackWeightSig", 100, -1., 1., 100, 0.,10.);
  TH2D track_pvDelta3DzoomVsTrackWeightSig("track_pvDelta3DzoomVsTrackWeightSig", "track_pvDelta3DzoomVsTrackWeightSig", 100, -1., 1., 100, 0., 0.1);
  TH2D track_pvDelta3DVsTrackWeightBkg("track_pvDelta3DVsTrackWeightBkg", "track_pvDelta3DVsTrackWeightBkg", 100, -1., 1., 100, 0.,10.);
  TH2D track_pvDelta3DzoomVsTrackWeightBkg("track_pvDelta3DzoomVsTrackWeightBkg", "track_pvDelta3DzoomVsTrackWeightBkg", 100, -1., 1., 100, 0., 0.1);
  
  TH2D electron_DxyVsTrackWeight("electron_DxyVsTrackWeight", "electron_DxyVsTrackWeight", 100, -1., 1., 100, 0., 1.);
  TH2D electron_PtVsDxy("electron_PtVsDxy", "electron_PtVsDxy", 150, 0., 15, 300, 0., 300.);
  TH2D electron_EtaVsDxy("electron_EtaVsDxy", "electron_EtaVsDxy", 150, 0., 15, 100, -4., 4.);
  TH2D electron_normChi2VsDxy("electron_normChi2VsDxy", "electron_normChi2VsDxy", 150, 0., 15, 100, 0., 4.);

  TH2D electron_PtVsTrackWeight("electron_PtVsTrackWeight", "electron_PtVsTrackWeight", 100, -1., 1., 300, 0., 300.);
  TH2D electron_EtaVsTrackWeight("electron_EtaVsTrackWeight", "electron_EtaVsTrackWeight", 100, -1., 1., 100, -4., 4.);
  TH2D electron_normChi2VsTrackWeight("electron_normChi2VsTrackWeight", "electron_normChi2VsTrackWeight", 100, -1., 1., 100, 0., 4.);
  
  // Vertex
  TH1D vertex_nTotal("vertex_nTotal", "vertex_nTotal", 20, 0., 20.);
  TH1D vertex_nTracks("vertex_nTracks", "vertex_nTracks", 100, 0., 100.);
  
  TH2D vertex_nVtxVsNtrks("vertex_nVtxVsNtrks", "vertex_nVtxVsNtrks", 110, 0., 110., 20, 0., 20.);
  //TH2D vertex_nVtxVsNtrksInVtx("vertex_nVtxVsNtrksInVtx", "vertex_nVtxVsNtrksInVtx", 100, 0., 100., 20, 0., 20.);
  TH2D vertex_nTracksVsNsig("vertex_nTracksVsNsig", "vertex_nTracksVsNsig", 5, 0., 5., 100, 0., 100.);
  TH2D vertex_nTracksVsSigWeight("vertex_nTracksVsSigWeight", "vertex_nTracksVsSigWeight", 101, 0., 1.01, 80, 0., 80.);
  TH2D vertex_TrackWeightVsSigWeight("vertex_TrackWeightVsSigWeight", "vertex_TrackWeightVsSigWeight", 101, 0., 2.01, 80, 0., 80.);
  TH2D vertex_TrackRatioVsSigWeight("vertex_TrackRatioVsSigWeight", "vertex_TrackRatioVsSigWeight", 101, 0., 2.01, 80, 0., 80.);
  
  // Tracks
  TH1D tracks_nTotal("tracks_nTotal", "tracks_nTotal", 110, 0., 110.);
  
  
  TH2D testHist("testHist", "testHist", 60, 0., 120., 40, 0., 200.);


  vector<string> labelsAll = {"signal out PV", "signal in PV", "background out PV", "background in PV"};
  vector<string> labelsOutPV = {"above", "below"};
  vector<TH1D> trackPtOverlay;
  vector<TH1D> trackEtaOverlay;
  vector<TH1D> trackDxyOverlay;
  vector<TH1D> trackNormChi2Overlay;
  vector<TH1D> trackPVdeltaZOverlay;
  vector<TH1D> trackPVdelta3DOverlay;
  vector<TH1D> trackPtOverlayOutPV;
  
  std::string filename = "root/GMSB_L-150TeV_Ctau-0_1cm_Fall17.root";                                                                                          
  //std::string filename = "root/GMSB_L-150TeV_Ctau-200cm_Fall17.root";
  //std::string filename = "gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root";
  //std::string filename = "/eos/uscms/store/group/lpcsusylep/anazario/LLP/electronMatching/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/electronMatching_DYJetsToLL_M-50_AODSIM_RunIISummer20UL18RECO-106X_upgrade2018/240710_144043/0000/output_*.root";

  TChain* tch = new TChain("tree/llpgtree");
  tch->Add(filename.c_str());

  llpgtree handle = llpgtree(tch);

  // Loop over ntuple entries
  const int nEntries = tch->GetEntries();
  //const int nEntries = 500000;
  //std::cout << nEntries << std::endl;

  for(int e = 0; e < nEntries; e++){

    if(e%1000 == 0) {
      fprintf(stdout, "\r  Processed events: %5d of %5d ", e, nEntries);
      fflush(stdout);
    }

    handle.fChain->GetEntry(e);

    const int nTracks(handle.Track_nTotal);
    const int nGen(handle.GenParticle_nMatches);
    const int nSig(handle.GenParticle_nSignal);
    const int nVtx(handle.Vertex_nTotal);
    const double bsx(handle.Beamspot_x), bsy(handle.Beamspot_y), bsz(handle.Beamspot_z);
    const double beamspotDxy(sqrt(handle.Beamspot_x*handle.Beamspot_x + handle.Beamspot_y*handle.Beamspot_y));

    vertex_nTotal.Fill(nVtx);
    tracks_nTotal.Fill(nTracks);
    vertex_nVtxVsNtrks.Fill(nTracks, nVtx);
    
    std::vector<int> signalIndex;
    
    for(int g = 0; g < nGen; g++) {

      const bool isSignal(handle.GenParticle_isSignal->at(g));
      const int trackIndex(handle.GenParticle_matchedTrackIndex->at(g));
      const int trackCharge(handle.Track_charge->at(trackIndex));
      const double trackPt(handle.Track_pt->at(trackIndex));
      const double trackEta(handle.Track_eta->at(trackIndex));
      const double trackPhi(handle.Track_phi->at(trackIndex));
      const double trackNormChi2(fabs(handle.Track_normalizedChi2->at(trackIndex)));
      const double trackDxy(handle.Track_dxy->at(trackIndex));
      const double trackX(handle.Track_vx->at(trackIndex));
      const double trackY(handle.Track_vy->at(trackIndex));
      const double trackZ(handle.Track_vz->at(trackIndex));
      const double trackIp2D(handle.Track_ip2D->at(trackIndex));
      const double trackSip2D(handle.Track_sip2D->at(trackIndex));
      const double trackIp3D(handle.Track_ip3D->at(trackIndex));
      const double trackSip3D(handle.Track_sip3D->at(trackIndex));
      const double genDxy(handle.GenParticle_dxy->at(g));
      const int pdgId(handle.GenParticle_pdgId->at(g));
      const double dxyFromBS(fabs(beamspotDxy-genDxy));
      const bool isECAL(handle.Track_isECAL->at(trackIndex));

      //if(trackNormChi2 >= 1.5 || trackNormChi2 < 0.5 || trackDxy < 0.5 || fabs(trackSip3D) <= 8 || trackPt < 20) continue;
      //if(trackNormChi2 >= 2 || trackNormChi2 < 0.5 || fabs(trackEta) > 2.4 || fabs(trackSip3D) <= 8 || fabs(trackSip2D) <= 5 || trackPt < 20 || trackDxy < 0.5) continue; 
      
      genParticle_dxy.Fill(genDxy);
      genParticle_PdgIdVsDxy.Fill(genDxy, abs(pdgId));

      if(isSignal) {
	genParticle_sigPt.Fill(handle.GenParticle_pt->at(g));
	genParticle_sigDeltaR.Fill(handle.GenParticle_deltaR->at(g));
	
	signalIndex.push_back(trackIndex);

	if(handle.GenParticle_deltaR->at(g) >= 0.1) continue;
	
	int idx(0);
	double trackWeight(-1.);
	double pvDeltaZ(999.);
	double pvDelta3D(999.);
	for(const auto &index : *(handle.PVCollection_trackIndex)) {

	  if(index == trackIndex && trackWeight < handle.PVCollection_signalWeight->at(idx))
	    trackWeight = handle.PVCollection_signalWeight->at(idx);

	  // Z-difference between track and current PV
	  const double pvZ(handle.PVCollection_z->at(idx));
	  const double currentDeltaZ = fabs(trackZ - pvZ);

	  if(currentDeltaZ < pvDeltaZ)
	    pvDeltaZ = currentDeltaZ;

	  // 3D-difference between track and current PV
	  const double pvX(handle.PVCollection_x->at(idx));
	  const double pvY(handle.PVCollection_y->at(idx));
          const double currentDelta3D = sqrt((trackX - pvX)*(trackX - pvX) + (trackY - pvY)*(trackY - pvY) + (trackZ - pvZ)*(trackZ - pvZ));

	  if(currentDelta3D < pvDelta3D)
            pvDelta3D = currentDelta3D;
	  
	  idx++;
	}
	genSignal_DxyVsTrackWeight.Fill(trackWeight, dxyFromBS);
	track_pvDelta3DVsTrackWeightSig.Fill(trackWeight, pvDelta3D);
	track_pvDelta3DzoomVsTrackWeightSig.Fill(trackWeight, pvDelta3D);
	
	int madeItToECAL = -1;
	if(isECAL) madeItToECAL = 1;
	
	if(trackWeight < 0.) {
	  track_ptSigMatch.Fill(trackPt);
	  track_etaSigMatch.Fill(trackEta);
	  track_dxySigMatch.Fill(trackDxy);
	  track_normChi2SigMatch.Fill(trackNormChi2);
	  track_pvDeltaZSigMatch.Fill(pvDeltaZ);
	  track_pvDelta3DSigMatch.Fill(pvDelta3D);
	  track_isECALsigOutPV.Fill(madeItToECAL);
	  if(madeItToECAL > 0) {
	    track_ptSig.Fill(trackPt);
	    track_etaSig.Fill(trackEta);
	    track_phiSig.Fill(trackPhi);
	    track_chargeSig.Fill(trackCharge);
	    track_normChi2Sig.Fill(trackNormChi2);
	    track_dxySig.Fill(trackDxy);
	    track_ip2DSig.Fill(trackIp2D);
	    track_sip2DSig.Fill(trackSip2D);
	    track_ip3DSig.Fill(trackIp3D);
            track_sip3DSig.Fill(trackSip3D);
	  }
	}
	else {
	  track_ptSigMatchPV.Fill(trackPt);
	  track_etaSigMatchPV.Fill(trackEta);
	  track_dxySigMatchPV.Fill(trackDxy);
          track_normChi2SigMatchPV.Fill(trackNormChi2);
	  track_pvDeltaZSigMatchPV.Fill(pvDeltaZ);
	  track_pvDelta3DSigMatchPV.Fill(pvDelta3D);
	  track_isECALsigInPV.Fill(madeItToECAL);
	}
      }

      else {

	if(handle.GenParticle_deltaR->at(g) >= 0.1) continue;

        int idx(0);
        double trackWeight(-1.);
	double pvDeltaZ(999.);
	double pvDelta3D(999.);
        for(const auto &index : *(handle.PVCollection_trackIndex)) {

          if(index == trackIndex && trackWeight < handle.PVCollection_signalWeight->at(idx))
            trackWeight = handle.PVCollection_signalWeight->at(idx);

	  const double pvZ(handle.PVCollection_z->at(idx));
          const double currentDeltaZ = fabs(trackZ - pvZ);

          if(currentDeltaZ < pvDeltaZ)
            pvDeltaZ = currentDeltaZ;

	  const double pvX(handle.PVCollection_x->at(idx));
	  const double pvY(handle.PVCollection_y->at(idx));
	  const double currentDelta3D = sqrt((trackX - pvX)*(trackX - pvX) + (trackY - pvY)*(trackY - pvY) + (trackZ - pvZ)*(trackZ - pvZ));
	  
	  if(currentDelta3D < pvDelta3D)
	    pvDelta3D = currentDelta3D;
	  
          idx++;
        }

	if(pvDeltaZ < 0.3) {
	  track_ptBkgMatchBelow.Fill(trackPt);
	}
	else {
	  track_ptBkgMatchAbove.Fill(trackPt);
	}

	int madeItToECAL = -1;
        if(isECAL) madeItToECAL = 1;

	if(trackWeight < 0.) {
          track_ptBkgMatch.Fill(trackPt);
          track_etaBkgMatch.Fill(trackEta);
	  track_dxyBkgMatch.Fill(trackDxy);
          track_normChi2BkgMatch.Fill(trackNormChi2);
	  track_pvDeltaZBkgMatch.Fill(pvDeltaZ);
	  track_pvDelta3DBkgMatch.Fill(pvDelta3D);
	  track_isECALbkgOutPV.Fill(madeItToECAL);
	  if(madeItToECAL > 0) {
            track_ptBkg.Fill(trackPt);
            track_etaBkg.Fill(trackEta);
            track_phiBkg.Fill(trackPhi);
	    track_chargeBkg.Fill(trackCharge);
            track_normChi2Bkg.Fill(trackNormChi2);
	    track_dxyBkg.Fill(trackDxy);
	    track_ip2DBkg.Fill(trackIp2D);
            track_sip2DBkg.Fill(trackSip2D);
	    track_ip3DBkg.Fill(trackIp3D);
            track_sip3DBkg.Fill(trackSip3D);
          }
	}
        else {
          track_ptBkgMatchPV.Fill(trackPt);
          track_etaBkgMatchPV.Fill(trackEta);
	  track_dxyBkgMatchPV.Fill(trackDxy);
          track_normChi2BkgMatchPV.Fill(trackNormChi2);
	  track_pvDeltaZBkgMatchPV.Fill(pvDeltaZ);
	  track_pvDelta3DBkgMatchPV.Fill(pvDelta3D);
	  track_isECALbkgInPV.Fill(madeItToECAL);
	}

	genBkg_DxyVsTrackWeight.Fill(trackWeight, dxyFromBS);
	genParticle_PtVsDxy.Fill(dxyFromBS, trackPt);
	genParticle_PtVsTrackWeight.Fill(trackWeight, trackPt);
	genParticle_EtaVsDxy.Fill(dxyFromBS, trackEta);
	genParticle_EtaVsTrackWeight.Fill(trackWeight, trackEta);
	genParticle_normChi2VsDxy.Fill(dxyFromBS, trackNormChi2);
	genParticle_normChi2VsTrackWeight.Fill(trackWeight, trackNormChi2);
	track_pvDelta3DVsTrackWeightBkg.Fill(trackWeight, pvDelta3D);
	track_pvDelta3DzoomVsTrackWeightBkg.Fill(trackWeight, pvDelta3D);

	if(abs(pdgId) == 11) {
	  electron_DxyVsTrackWeight.Fill(trackWeight, dxyFromBS);
	  electron_PtVsDxy.Fill(dxyFromBS, trackPt);
	  electron_EtaVsDxy.Fill(dxyFromBS, trackEta);
	  electron_normChi2VsDxy.Fill(dxyFromBS, trackNormChi2);
	  electron_PtVsTrackWeight.Fill(trackWeight, trackPt);
	  electron_EtaVsTrackWeight.Fill(trackWeight, trackEta);
	  electron_normChi2VsTrackWeight.Fill(trackWeight, trackNormChi2);
	}
      }
    }

    trackPtOverlay = {track_ptSigMatch, track_ptSigMatchPV, track_ptBkgMatch, track_ptBkgMatchPV};
    trackEtaOverlay = {track_etaSigMatch, track_etaSigMatchPV, track_etaBkgMatch, track_etaBkgMatchPV};
    trackDxyOverlay = {track_dxySigMatch, track_dxySigMatchPV, track_dxyBkgMatch, track_dxyBkgMatchPV};
    trackNormChi2Overlay = {track_normChi2SigMatch, track_normChi2SigMatchPV, track_normChi2BkgMatch, track_normChi2BkgMatchPV};
    trackPVdeltaZOverlay = {track_pvDeltaZSigMatch, track_pvDeltaZSigMatchPV, track_pvDeltaZBkgMatch, track_pvDeltaZBkgMatchPV};
    trackPVdelta3DOverlay = {track_pvDelta3DSigMatch, track_pvDelta3DSigMatchPV, track_pvDelta3DBkgMatch, track_pvDelta3DBkgMatchPV};
    trackPtOverlayOutPV = {track_ptBkgMatchAbove, track_ptBkgMatchBelow};
    
    for(int v = 0; v < nVtx; v++) {

      const int nTracksInVertex(handle.Vertex_nTracks->at(v));
      const bool hasSignal(handle.Vertex_hasSignal->at(v));
      const double vertexDxy(handle.Vertex_dxy->at(v));

      if(hasSignal) {
	const int nSig(handle.Vertex_signalCount->at(v));
	vertex_nTracksVsNsig.Fill(nSig, nTracksInVertex);
	//if(nSig == 2 && nTracksInVertex == 2)
	//cout << "found a signal vertex candidate with dxy = " << vertexDxy << endl;	
      }
      
      vertex_nTracks.Fill(nTracksInVertex);

      bool isSignal(false);
      double totalWeight(0.), signalWeight(0.);
      for(int t = 0; t < nTracksInVertex; t++) {
	const double trackWeight(handle.Vertex_trackWeight->at(t));
	totalWeight += trackWeight;
	for(const auto &idx : signalIndex) {
	  if(idx == t && handle.Vertex_vertexIndex->at(t) == v) {
	    isSignal = true;
	    signalWeight += trackWeight;
	    vertex_nTracksVsSigWeight.Fill(trackWeight, nTracksInVertex);
	    break;
	  }
	}
      }
      vertex_TrackWeightVsSigWeight.Fill(signalWeight, totalWeight);
      vertex_TrackRatioVsSigWeight.Fill(signalWeight, totalWeight/nTracksInVertex);
    }
  }// End event loop
  cout << endl;

  TFile outfile("vertex_histograms.root", "RECREATE");
  outfile.cd();
  plotHistogram(genParticle_sigPt, "gen pT")->Write();
  plotHistogram(genParticle_sigDeltaR, "delta R")->Write();
  plotHistogram(genParticle_dxy, "gen dxy [cm]")->Write();
  plotHistogram(vertex_nTotal, "nVertex")->Write();
  plotHistogram(vertex_nTracks, "nTracks")->Write();
  plotHistogram(vertex_nVtxVsNtrks, "nTracks", "nVertex")->Write();
  plotHistogram(vertex_nTracksVsNsig, "nSig", "nTracks in vertex")->Write();
  plotHistogram(vertex_nTracksVsSigWeight, "Signal Track Weight in Vertex", "nTracks in Vertex")->Write();
  plotHistogram(vertex_TrackWeightVsSigWeight, "Signal Track Weight in Vertex", "Total Track Weight in Vertex")->Write();
  plotHistogram(vertex_TrackRatioVsSigWeight, "Signal Track Weight in Vertex", "Total Track Weight/nTracks in Vertex")->Write();
  plotHistogram(genSignal_DxyVsTrackWeight, "weight", "dxy [cm]")->Write();
  plotHistogram(genBkg_DxyVsTrackWeight, "weight", "dxy [cm]")->Write();
  plotHistogram(genParticle_PdgIdVsDxy, "dxy [cm]", "pdg ID")->Write();
  plotHistogram(genParticle_PtVsDxy, "dxy [cm]", "pt")->Write();
  plotHistogram(genParticle_EtaVsDxy, "dxy [cm]", "eta")->Write();
  plotHistogram(genParticle_normChi2VsDxy, "dxy [cm]", "norm Chi2")->Write();
  plotHistogram(genParticle_PtVsTrackWeight, "weight", "pt")->Write();
  plotHistogram(genParticle_EtaVsTrackWeight, "weight", "eta")->Write();
  plotHistogram(genParticle_normChi2VsTrackWeight, "weight", "norm Chi2")->Write();
  plotHistogram(track_pvDelta3DVsTrackWeightSig, "weight", "pvDelta3D")->Write();
  plotHistogram(track_pvDelta3DzoomVsTrackWeightSig, "weight", "pvDelta3D")->Write();
  plotHistogram(track_pvDelta3DVsTrackWeightBkg, "weight", "pvDelta3D")->Write();
  plotHistogram(track_pvDelta3DzoomVsTrackWeightBkg, "weight", "pvDelta3D")->Write();
  //plotHistogram(electron_PtVsTrackWeight, "weight", "pt")->Write();
  //plotHistogram(electron_EtaVsTrackWeight, "weight", "eta")->Write();
  //plotHistogram(electron_normChi2VsTrackWeight, "weight", "norm Chi2")->Write();
  plotHistogram(electron_DxyVsTrackWeight, "weight", "dxy [cm]")->Write();
  plotHistogram(electron_PtVsDxy, "dxy [cm]", "pt")->Write();
  plotHistogram(electron_EtaVsDxy, "dxy [cm]", "eta")->Write();
  plotHistogram(electron_normChi2VsDxy, "dxy [cm]", "norm Chi2")->Write();
  plotHistogram(tracks_nTotal, "nTracks")->Write();
  plotHistograms(trackPtOverlay, labelsAll, "trackPtOverlay", "pT")->Write();
  plotHistograms(trackEtaOverlay, labelsAll, "trackEtaOverlay", "#eta")->Write();
  plotHistograms(trackDxyOverlay, labelsAll, "trackDxyOverlay", "dxy")->Write();
  plotHistograms(trackNormChi2Overlay, labelsAll, "trackNormChi2Overlay", "normChi2")->Write();
  plotHistograms(trackPVdeltaZOverlay, labelsAll, "trackPVdeltaZOverlay", "pvDeltaZ")->Write();
  plotHistograms(trackPVdelta3DOverlay, labelsAll, "trackPVdelta3DOverlay", "pvDelta3D")->Write();
  plotHistograms(track_isECALsigOutPV, track_isECALbkgOutPV, "track_isECALoutPV", "isECAL")->Write();
  plotHistograms(track_isECALsigInPV, track_isECALbkgInPV, "track_isECALInPV", "isECAL")->Write();
  plotHistograms(track_ptSig, track_ptBkg, "track_pt", "pt")->Write();
  plotHistograms(track_etaSig, track_etaBkg, "track_eta", "eta")->Write();
  plotHistograms(track_phiSig, track_phiBkg, "track_phi", "phi")->Write();
  plotHistograms(track_chargeSig, track_chargeBkg, "track_charge", "charge")->Write();
  plotHistograms(track_normChi2Sig, track_normChi2Bkg, "track_normChi2", "normChi2")->Write();
  plotHistograms(track_dxySig, track_dxyBkg, "track_dxy", "dxy")->Write();
  plotHistograms(track_ip2DSig, track_ip2DBkg, "track_ip2D", "ip2D")->Write();
  plotHistograms(track_sip2DSig, track_sip2DBkg, "track_sip2D", "sip2D")->Write();
  plotHistograms(track_ip3DSig, track_ip3DBkg, "track_ip3D", "ip3D")->Write();
  plotHistograms(track_sip3DSig, track_sip3DBkg, "track_sip3D", "sip3D")->Write();
  //plotHistograms(trackPtOverlayOutPV, labelsOutPV, "trackPtOverlayOutPVbkg", "pT")->Write();
  //plotHistogram(testHist, "pT", "M")->Write();

  outfile.Close();
}

float GetDXY(const float x, const float y) {
  return sqrt(x*x + y*y);
}

int GenMatchIndex(const LorentzVec &trackVec, const std::vector<LorentzVec> &genVecs, float &deltaR) {

  int matchIndex = -1;
  float minDeltaR = 10;//0.05;
  
  const int nGenElectrons = genVecs.size();
  for(int e = 0; e < nGenElectrons; e++) {
    const LorentzVec genEle = genVecs[e];
    deltaR = ROOT::Math::VectorUtil::DeltaR(trackVec, genEle);

    if(deltaR < minDeltaR) {
      minDeltaR = deltaR;
      matchIndex = e;
    }
  }

  deltaR = minDeltaR;
  
  return matchIndex;
}

TCanvas* plotHistograms(TH1D &signal, TH1D &background, const std::string& name, const std::string& xlabel, const std::string& ylabel) {
  // Create a canvas
  TCanvas* canvas = new TCanvas((name+"_canvas").c_str(), (name+"_canvas").c_str(), 800, 600);

  // Find the maximum value of the histograms
  double maxSignal = signal.GetMaximum();
  double maxBackground = background.GetMaximum();
  double maxY = std::max(maxSignal, maxBackground);
  
  // Set the range in y accordingly
  signal.SetMaximum(1.5 * maxY);
  signal.SetMinimum(0.9);
  
  // Draw signal histogram
  signal.SetLineColor(kBlue);
  signal.Draw("hist");
  
  // Draw background histogram
  background.SetLineColor(kRed);
  background.SetLineWidth(2);
  background.Draw("histsame");
  
  // Create a legend
  TLegend legend(0.7, 0.7, 0.9, 0.9);
  legend.AddEntry(&signal, "Signal", "l");
  legend.AddEntry(&background, "Background", "l");
  //legend.Draw();
  
  // Set labels and titles
  signal.GetXaxis()->SetTitle(xlabel.c_str());
  signal.GetYaxis()->SetTitle(ylabel.c_str());
  signal.SetTitle("");
  
  // Remove the title and the stat box
  signal.SetStats(0);
  
  // Add grids in x and y
  signal.SetMarkerStyle(20);
  signal.SetMarkerSize(0.5);
  signal.SetMarkerColor(kBlue);
  signal.SetLineColor(kBlue);
  signal.SetLineWidth(2);
  
  background.SetMarkerStyle(20);
  background.SetMarkerSize(0.5);
  background.SetMarkerColor(kRed);
  background.SetLineColor(kRed);
  
  signal.Draw("");
  //signal.Draw("hist same axis");
  background.Draw("histsame");
  legend.Draw("same");
  
  // Save the plot as a .pdf
  std::string pdfName = name + ".pdf";
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy();
  //canvas->SaveAs(pdfName.c_str());

  return canvas;
}

TCanvas* plotHistogram(TH1D& hist1D, const std::string& xlabel, const std::string& ylabel) {
    // Extract the name from the TH1D
    std::string name = hist1D.GetName();

    // Create a canvas
    TCanvas* canvas = new TCanvas((name+"_canvas").c_str(), (name+"_canvas").c_str(), 800, 600);

    // Draw the 1D histogram
    hist1D.SetStats(0);
    hist1D.Draw();

    // Set labels and titles
    hist1D.GetXaxis()->SetTitle(xlabel.c_str());
    hist1D.GetYaxis()->SetTitle(ylabel.c_str());
    hist1D.SetTitle("");

    // Save the plot as a .pdf
    std::string pdfName = name + ".pdf";
    canvas->SetGridx();
    canvas->SetGridy();
    //canvas->SaveAs(pdfName.c_str());

    return canvas;
}


TCanvas* plotHistogram(TH2D& hist2D, const std::string& xlabel, const std::string& ylabel, const std::string& zlabel = "Events") {
    // Extract the name from the TH2D
    std::string name = hist2D.GetName();

    // Create a canvas
    TCanvas* canvas = new TCanvas((name+"_canvas").c_str(), (name+"_canvas").c_str(), 800, 600);

    // Draw the 2D histogram
    hist2D.SetStats(0);
    hist2D.Draw("colz");

    // Set labels and titles
    hist2D.GetXaxis()->SetTitle(xlabel.c_str());
    hist2D.GetYaxis()->SetTitle(ylabel.c_str());
    hist2D.GetZaxis()->SetTitle(zlabel.c_str());
    hist2D.SetTitle("");

    // Save the plot as a .pdf
    std::string pdfName = name + ".pdf";
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->SetLogz();
    //canvas->SaveAs(pdfName.c_str());

    return canvas;
}

TCanvas* plotEfficiency(std::vector<TEfficiency*> &efficiencies, std::vector<string> &labels, string name, string xlabel) {

  vector<int> markers = {20,21,22,23,29,33,34,43,49};
  std::vector<int> colors = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-2,kBlack,kAzure+2,kGray};

  TMultiGraph* mg = new TMultiGraph();
  TCanvas *canvas = new TCanvas(("can_"+name).c_str(), ("can_"+name).c_str(), 800, 600);
  TLegend* legend = new TLegend(0.75, 0.18, 0.89, 0.32);
  canvas->SetLeftMargin(0.15);
  canvas->SetGrid();

  for(int i = 0; i < efficiencies.size(); i++) {
    efficiencies[i]->SetMarkerStyle(markers[i]);
    efficiencies[i]->SetStatisticOption(TEfficiency::kFNormal);
    efficiencies[i]->SetMarkerColor(colors[i]);
    efficiencies[i]->SetLineColor(colors[i]);
    mg->Add(efficiencies[i]->CreateGraph());
  }

  //efficiencies[0]->Draw();
  mg->Draw("ap");
  mg->SetMinimum(0.); 
  mg->SetMaximum(1.); 
  mg->GetXaxis()->SetTitle(xlabel.c_str());
  mg->GetXaxis()->SetTitleOffset(1.25);
  mg->GetYaxis()->SetTitle("Efficiency");
  legend->AddEntry(efficiencies[0], labels[0].c_str(), "lep");

  for(int i = 1; i < efficiencies.size(); i++) {
    //efficiencies[i]->Draw("same");
    legend->AddEntry(efficiencies[i], labels[i].c_str(), "lep");
  }

  legend->Draw("same");

  //gPad->Update(); 
  //efficiencies[0]->GetPaintedGraph()->GetHistogram()->GetYaxis()->SetRangeUser(0.,1.);
  //efficiencies[0]->GetPaintedGraph()->SetMaximum(1);
  gPad->Update(); 

  return canvas;
}

particleID ClassifyParticle(const int pdgID) {
  switch (pdgID) {
  case 22:
    return particleID::kPhoton;
  case 211:
    return particleID::kPion;
  case 2212:
    return particleID::kProtons;
  case 2112:
    return particleID::kNeutrons;
  case 321:
    return particleID::kKaon;
  case 130:
    return particleID::kKlong;
  case 11:
    return particleID::kElectrons;
  case 13:
    return particleID::kMuon;
  case 1000039:
    return particleID::kSusy;
  case 12:
  case 14:
  case 16:
    return particleID::kNeutrino;
  default:
    return particleID::kUnmatched;
  }
}

template <typename T>
void FillIDHists(std::vector<TH1D>& histID, const int pdgID, const T quantity) {
  particleID id = ClassifyParticle(pdgID);
  switch (id) {
  case particleID::kPhoton:
    histID[0].Fill(quantity);
    break;
  case particleID::kPion:
    histID[1].Fill(quantity);
    break;
  case particleID::kProtons:
    histID[2].Fill(quantity);
    break;
  case particleID::kNeutrons:
    histID[3].Fill(quantity);
    break;
  case particleID::kKaon:
    histID[4].Fill(quantity);
    break;
  case particleID::kKlong:
    histID[5].Fill(quantity);
    break;
  case particleID::kMuon:
    histID[6].Fill(quantity);
    break;
  case particleID::kElectrons:
    histID[7].Fill(quantity);
    break;
  case particleID::kSusy:
    histID[8].Fill(quantity);
    break;
  case particleID::kNeutrino:
    histID[9].Fill(quantity);
    break;
  case particleID::kUnmatched:
    histID[10].Fill(quantity);
    break;
  }
}

std::string GetParticleTypeString(particleID id) {
  switch (id) {
  case particleID::kPhoton:
    return "Photon";
  case particleID::kPion:
    return "Pion";
  case particleID::kProtons:
    return "Proton";
  case particleID::kNeutrons:
    return "Neutron";
  case particleID::kKaon:
    return "Kaon";
  case particleID::kKlong:
    return "K-long";
  case particleID::kMuon:
    return "Muon";
  case particleID::kElectrons:
    return "Electron";
  case particleID::kSusy:
    return "SUSY Particle";
  case particleID::kNeutrino:
    return "Neutrino";
  default:
    return "Unmatched";
  }
}

TCanvas* plotHistograms(std::vector<TH1D> &histograms, const std::vector<std::string> &labels, const std::string name, std::string xlabel) {

  histograms[0].SetStats(0);
  histograms[0].SetTitle("");
  // Define the color vector
  std::vector<int> colors = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-2,kBlack,kAzure+2,kGray};

  // Create the canvas
  TCanvas* canvas = new TCanvas(("can_"+name).c_str(), "Histogram Plot", 800, 600);

  // Create the legend
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

  histograms[0].GetXaxis()->SetTitle(xlabel.c_str());

  // Loop through the histograms and plot them
  for (size_t i = 0; i < histograms.size(); i++) {
    histograms[i].SetLineWidth(3);
    histograms[i].SetLineColor(colors[i]);
    histograms[i].Draw(i == 0 ? "" : "same");
    legend->AddEntry(&histograms[i], labels[i].c_str(), "l");
  }

  // Find the maximum y-value
  double maxY = 0;
  for (const auto& h : histograms) {
    double thisMax = h.GetMaximum();
    if (thisMax > maxY) {
      maxY = thisMax;
    }
  }

  // Set the y-axis range and add a grid
  canvas->cd();
  canvas->SetGrid();
  canvas->SetLogy();
  double yMax = maxY * 1.2; // 20% padding
  histograms[0].GetYaxis()->SetRangeUser(0.01, yMax); 

  // Draw the legend
  legend->Draw();

  return canvas;
}
