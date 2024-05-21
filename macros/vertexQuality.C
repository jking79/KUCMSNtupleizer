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
TCanvas* plotHistograms(std::vector<TH1D> &histograms, const std::vector<std::string> &labels, std::string xlabel);
TCanvas* plotEfficiency(std::vector<TEfficiency*> &efficiencies, std::vector<string> &labels, string name, string xlabel);

void vertexQuality() {

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

  // General Tracks
  TH1D generalTracks_ip2D("generalTracks_ip2D", "generalTracks_ip2D", 120, -60, 60);
  TH1D generalTracksSig_ip2D("generalTracksSig_ip2D", "generalTracksSig_ip2D", 120, -60, 60); 
  TH1D generalTracks_sip2D("generalTracks_sip2D", "generalTracks_sip2D", 120, -60, 60);
  TH1D generalTracksSig_sip2D("generalTracksSig_sip2D", "generalTracksSig_sip2D", 120, -60, 60);

  TH2D generalTracks_ip2DVsSip2D("generalTracks_ip2DVsSip2D", "generalTracks_ip2DVsSip2D",  120, -60, 60,  120, -60, 60);
  TH2D generalTracksSig_ip2DVsSip2D("generalTracksSig_ip2DVsSip2D", "generalTracksSig_ip2DVsSip2D",  120, -60, 60,  120, -0.2, 0.2);

  // General Vertices
  TH1D generalVertices_ecalness("generalVertices_ecalness", "generalVertices_ecalness", 100, 0, 1);
  TH1D generalVerticesSig_ecalness("generalVerticesSig_ecalness", "generalVerticesSig_ecalness", 100, 0, 1);
  TH1D generalVertices_trackOverlapPV("generalVertices_trackOverlapPV", "generalVertices_trackOverlapPV", 100, 0, 1);
  TH1D generalVerticesSig_nTracks("generalVerticesSig_nTracks", "generalVerticesSig_nTracks", 100, 0, 100);
  TH1D generalVerticesSig_dxyPVnoSig("generalVerticesSig_dxyPVnoSig", "generalVerticesSig_dxyPVnoSig", 100, 0, 1);
  TH1D generalVerticesSig_dxyPVSig("generalVerticesSig_dxyPVSig", "generalVerticesSig_dxyPVSig", 100, 0, 1);

  TH2D generalVertices_nTracksVsEcalness("generalVertices_nTracksVsEcalness", "generalVertices_nTracksVsEcalness", 100, 0, 1, 100, 0, 200);
  TH2D generalVerticesSig_nTracksVsEcalness("generalVerticesSig_nTracksVsEcalness", "generalVerticesSig_nTracksVsEcalness", 100, 0, 1, 100, 0, 200);
  TH2D generalVerticesSig_nTracksVsDxy("generalVerticesSig_nTracksVsDxy", "generalVerticesSig_nTracksVsDxy", 70, 0, 70, 100, 0, 200);
  TH2D generalVertices_signalCountPVvsAVR("generalVertices_signalCountPVvsAVR", "generalVertices_signalCountPVvsAVR", 20, 0, 5, 20, 0 ,5);

  // Signal Vertices
  TH1D signalVertices_nTotal("signalVertices_nTotal", "signalVertices_nTotal", 3, 2, 5);
  TH1D signalVertices_nTracks("signalVertices_nTracks", "signalVertices_nTracks", 3, 2, 5);
  TH1D signalVertices_dxy("signalVertices_dxy", "signalVertices_dxy", 280, 0 , 70);
  TH1D signalVertices_genPtRelDiff("signalVertices_genPtRelDiff", "signalVertices_genPtRelDiff", 100, 0., 10.);
  TH1D signalVertices_trackPtBadMatches("signalVertices_trackPtBadMatches", "signalVertices_trackPtBadMatches", 100, 0., 200.);
  TH1D signalVertices_trackPtGoodMatches("signalVertices_trackPtGoodMatches", "signalVertices_trackPtGoodMatches", 100, 0., 200.);

  TH2D signalVertices_nTracksVsDxy("signalVertices_nTracksVsDxy", "signalVertices_nTracksVsDxy", 350, 0, 70, 3, 2, 5);
  TH2D signalVertices_nTracksVsDxyZoom("signalVertices_nTracksVsDxyZoom", "signalVertices_nTracksVsDxyZoom", 100, 0, 0.1, 3, 2, 5);

  TH2D signal_ptDiffVsDeltaR("signal_ptDiffVsDeltaR", "signal_ptDiffVsDeltaR", 60, 0, 0.1, 200, -20, 100);

  TH2D signal_genDxyVsSVdxy("signal_genDxyVsSVdxy", "signal_genDxyVsSVdxy", 100, 0, 100, 100, 0, 100);
  TH2D signal_dxyDiffVsNormChi2("signal_dxyDiffVsNormChi2", "signal_dxyDiffVsNormChi2", 100, 0, 5, 100, 0, 30);
  TH2D signal_dxyDiffVsGenPtRelDiff("signal_dxyDiffVsGenPtRelDiff", "signal_dxyDiffVsGenPtRelDiff", 100, -10., 10., 100, 0., 30.);
  TH2D signal_dxyDiffVsWeightRatio("signal_dxyDiffVsWeightRatio", "signal_dxyDiffVsWeightRatio", 100, 0., 1.,  100, 0., 30.);
  
  TH2D testHist("testHist", "testHist", 60, 0., 120., 40, 0., 200.);

  //std::string filename = "root/GMSB_L-150TeV_Ctau-0_1cm_Fall17.root";                                                                                          
  std::string filename = "root/GMSB_L-150TeV_Ctau-200cm_Fall17.root";
  //std::string filename = "root/GMSB_L-150TeV_Ctau-0_1cm_Fall17_GeneralVertices.root";
  //std::string filename = "root/GMSB_L-150TeV_Ctau-0_1cm_Fall17_GeneralVerticesRefitTracks.root"; 
  //std::string filename = "root/GMSB_L-150TeV_Ctau-200cm_Fall17_GeneralVertices.root";
  //std::string filename = "root/GMSB_L-150TeV_Ctau-200cm_Fall17_GeneralVerticesRefitTracks.root";
  //std::string filename = "gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root";

  TChain* tch = new TChain("tree/llpgtree");
  tch->Add(filename.c_str());

  llpgtree handle = llpgtree(tch);

  // Loop over ntuple entries
  const int nEntries = tch->GetEntries();
  
  //std::cout << nEntries << std::endl;

  for(int e = 0; e < nEntries; e++){

    if(e%1000 == 0) {
      fprintf(stdout, "\r  Processed events: %5d of %5d ", e, nEntries);
      fflush(stdout);
    }

    handle.fChain->GetEntry(e);

    const int nTracks(handle.Track_nTotal);
    const int nGen(handle.GenParticle_nMatches);
    const int nVtx(handle.Vertex_nTotal);
    const double beamspotDxy(sqrt(handle.Beamspot_x*handle.Beamspot_x + handle.Beamspot_y*handle.Beamspot_y));
    //std::cout << nTracks << std::endl;
    
    std::vector<int> signalIndex;

    for(int g = 0; g < nGen; g++) {

      const bool isSignal(handle.GenParticle_isSignal->at(g));
      const int trackIndex(handle.GenParticle_matchedTrackIndex->at(g));
      const double ip2D(handle.Track_ip2D->at(trackIndex));
      const double sip2D(handle.Track_sip2D->at(trackIndex));

      if(isSignal) {
	signalIndex.push_back(trackIndex);
	generalTracksSig_ip2D.Fill(ip2D);
	generalTracksSig_sip2D.Fill(sip2D);
	generalTracksSig_ip2DVsSip2D.Fill(sip2D, ip2D);
      }
      else {
	generalTracks_ip2D.Fill(handle.Track_ip2D->at(trackIndex));
	generalTracks_sip2D.Fill(handle.Track_ip2D->at(trackIndex));
	generalTracks_ip2DVsSip2D.Fill(sip2D, ip2D);
      }
    }

    generalVertices_trackOverlapPV.Fill(handle.Vertex_trackRatioWithPV->at(0));
    //generalVertices_signalCountPVvsAVR.Fill(handle.Vertex_signalCount->at(0), handle.PV_signalCount);

    double pvWeightedNtracks(0);
    for(int t = 0; t < handle.PV_trackIndex->size(); t++) {
      const int trackIndex(handle.PV_trackIndex->at(t));
      
      for(const auto &idx : signalIndex)
          if(trackIndex == idx)
	    pvWeightedNtracks += handle.PV_trackWeight->at(t);
    }

    double avrPVWeightedNtracks(0);
    for(int v = 0; v < nVtx; v++) {
      
      bool hasSignal(false);

      int nTrk = 0;//(handle.Vertex_nTracks->at(v));
      
      for(int t = 0; t < handle.Vertex_nTracks->size(); t++) {
	const int vertexIndex(handle.Vertex_vertexIndex->at(t));
	const int trackIndex(handle.Vertex_trackIndex->at(t));
	const double weight(handle.Vertex_trackWeight->at(t));
	nTrk += weight;
	
	for(const auto &idx : signalIndex) {
	  if(vertexIndex == v && trackIndex == idx)
	    hasSignal = true;

	  if(vertexIndex == 0 && trackIndex == idx)
	    avrPVWeightedNtracks += weight;
	}
      }

      const double vertexDxy(fabs(handle.Vertex_dxy->at(v)-beamspotDxy));
      const double ecalness(handle.Vertex_ecalness->at(v));

      if(hasSignal) {
	generalVerticesSig_ecalness.Fill(ecalness);
	generalVerticesSig_nTracksVsEcalness.Fill(ecalness, nTrk);
	generalVerticesSig_nTracksVsDxy.Fill(vertexDxy, nTrk);
	generalVerticesSig_nTracks.Fill(handle.Vertex_signalCount->at(v));
	if(!handle.PV_hasSignal && v != 0)
	  generalVerticesSig_dxyPVnoSig.Fill(vertexDxy);
	else if(handle.PV_hasSignal && v != 0)
	  generalVerticesSig_dxyPVSig.Fill(vertexDxy);
      }
      else {
	generalVertices_ecalness.Fill(ecalness);
	generalVertices_nTracksVsEcalness.Fill(ecalness, nTrk);
      }
    }
    if(avrPVWeightedNtracks > 0 || pvWeightedNtracks > 0)
      generalVertices_signalCountPVvsAVR.Fill(avrPVWeightedNtracks, pvWeightedNtracks);
    
    
    signalVertices_nTotal.Fill(handle.SignalSV_nTotal);
    for(int v = 0; v < handle.SignalSV_nTotal; v++) {
      
      const int nTracks(handle.SignalSV_nTracks->at(v));
      const double dxy(fabs(handle.SignalSV_dxy->at(v)-beamspotDxy));
      
      signalVertices_nTracks.Fill(nTracks);
      signalVertices_dxy.Fill(dxy);
      signalVertices_nTracksVsDxy.Fill(dxy, nTracks);
      signalVertices_nTracksVsDxyZoom.Fill(dxy, nTracks);

      double totalWeight(0.);
      double genDxy(0.);
      for(int t = 0; t < handle.SignalSV_trackIndex->size(); t++) {
	const int vertexIndex(handle.SignalSV_vertexIndex->at(t));
	const int trackIndex(handle.SignalSV_trackIndex->at(t));
	const int genIndex(handle.SignalSV_genIndex->at(t));
	//const double genDxy(handle.GenParticle_dxy->at(genIndex));
	if(!(genDxy > 0))
	  genDxy = handle.GenParticle_dxy->at(genIndex);

	const double trackPt(handle.Track_pt->at(trackIndex));
	const double genPt(handle.GenParticle_pt->at(genIndex));
	const double genPtRelDiff((trackPt-genPt)/genPt);

	if(genPtRelDiff < -0.5)
	  signalVertices_trackPtBadMatches.Fill(trackPt);
	else
	  signalVertices_trackPtGoodMatches.Fill(trackPt);
	
	if(vertexIndex == v) {
	  totalWeight += handle.SignalSV_trackWeight->at(t);
	  signal_genDxyVsSVdxy.Fill(dxy, genDxy);
	  signalVertices_genPtRelDiff.Fill(genPtRelDiff);
	  signal_dxyDiffVsNormChi2.Fill(handle.SignalSV_normalizedChi2->at(v), fabs(genDxy-dxy));
	  signal_dxyDiffVsGenPtRelDiff.Fill(genPtRelDiff, fabs(genDxy-dxy));
	}
      }
      signal_dxyDiffVsWeightRatio.Fill(totalWeight/nTracks, fabs(genDxy-dxy));
    }

    /*
    for(int t = 0; t < nTracks; t++) {

      generalTracks_1p2D.Fill(handle.Track_ip2D->at(t));
    }
    */
  }// End event loop
  cout << endl;

  TFile outfile("generalVertex_histograms.root", "RECREATE");
  outfile.cd();

  plotHistograms(generalTracksSig_ip2D, generalTracks_ip2D, "generalTracks_ip2D", "ip2D")->Write();
  plotHistograms(generalTracksSig_sip2D, generalTracks_sip2D, "generalTracks_sip2D", "sip2D")->Write();
  plotHistograms(generalVerticesSig_dxyPVnoSig, generalVerticesSig_dxyPVSig, "generalVerticesSig_dxyPV", "dxy")->Write();
  plotHistogram(generalTracksSig_ip2DVsSip2D, "sip2D", "ip2D")->Write();
  plotHistogram(generalTracks_ip2DVsSip2D, "sip2D", "ip2D")->Write();
  plotHistograms(generalVerticesSig_ecalness, generalVertices_ecalness, "generalVertices_ecalness", "ecalness")->Write();
  plotHistogram(generalVerticesSig_nTracksVsEcalness, "ecalness", "nTracks")->Write();
  plotHistogram(generalVertices_nTracksVsEcalness, "ecalness", "nTracks")->Write();
  plotHistogram(generalVerticesSig_nTracksVsDxy, "dxy", "nTracks")->Write();
  plotHistogram(generalVerticesSig_nTracks, "nTracks")->Write();
  plotHistogram(generalVertices_trackOverlapPV, "overlap")->Write();
  plotHistogram(generalVertices_signalCountPVvsAVR, "nSignal AVR", "nSignal PV")->Write();
  plotHistogram(generalVerticesSig_dxyPVnoSig, "dxy")->Write();

  plotHistograms(signalVertices_trackPtGoodMatches, signalVertices_trackPtBadMatches, "signalVertices_trackPt", "pt")->Write();
  
  plotHistogram(signalVertices_nTotal, "nVertices")->Write();
  plotHistogram(signalVertices_nTracks, "nTracks")->Write();
  plotHistogram(signalVertices_dxy, "dxy")->Write();
  plotHistogram(signalVertices_genPtRelDiff, "(genPt - trackPt) / genPt")->Write();
  plotHistogram(signalVertices_nTracksVsDxy, "dxy", "nTracks")->Write();
  plotHistogram(signalVertices_nTracksVsDxyZoom, "dxy", "nTracks")->Write();
  plotHistogram(signal_genDxyVsSVdxy, "vertex dxy", "gen dxy")->Write();
  plotHistogram(signal_dxyDiffVsNormChi2, "norm Chi2", "|gen dxy - vertex dxy|")->Write();
  plotHistogram(signal_dxyDiffVsGenPtRelDiff, "(genPt - trackPt) / genPt", "|gen dxy - vertex dxy|")->Write();
  plotHistogram(signal_dxyDiffVsWeightRatio, "totalWeight/nTracks", "dxy")->Write();
  
  //plotHistogram(signal_ptDiffVsDeltaR, "deltaR", "p_T^{Gsf}-p_T^{General}")->Write();
 
  //plotEfficiency(idPtEfficiencies, idLabels, "UnmatchedComposition", "pT")->Write();

  // Vertex Quality
  //plotHistograms(vertexPosXSig, vertexPosXBkg, "vertexPosX", "x [cm]")->Write();

  // Test
  //plotHistograms(test1DSig, test1DBkg, "test1D", test1D_xlabel)->Write();
  //plotHistograms( , , "", "", "");
  /*
  //1D plots
  plotHistogram(genElectron_pt, "pt")->Write();
  plotHistogram(genElectron_eta, "eta")->Write();
  plotHistogram(generalTracksGenMatch_deltaR, "deltaR")->Write();
  plotHistogram(gsfTracksGenMatch_deltaR, "deltaR")->Write();
  plotHistograms(generalTracksGoodGenMatch_pt, generalTracksBadGenMatch_pt, "generalTracksGenMatch_pt", "pt")->Write();
  plotHistograms(generalTracksGoodGenMatch_genPt, generalTracksBadGenMatch_genPt, "generalTracksGenMatch_genPt", "genPt")->Write();
  plotHistograms(generalTracksGoodGenMatch_eta, generalTracksBadGenMatch_eta, "generalTracksGenMatch_eta", "eta")->Write();
  plotHistograms(generalTracksGoodGenMatch_genEta, generalTracksBadGenMatch_genEta, "generalTracksGenMatch_genEta", "genEta")->Write();
  plotHistograms(gsfTracksGoodGenMatch_pt, gsfTracksBadGenMatch_pt, "gsfTracksGenMatch_pt", "pt")->Write();
  plotHistograms(gsfTracksGoodGenMatch_genPt, gsfTracksBadGenMatch_genPt, "gsfTracksGenMatch_genPt", "genPt")->Write();
  plotHistograms(gsfTracksGoodGenMatch_eta, gsfTracksBadGenMatch_eta, "gsfTracksGenMatch_eta", "eta")->Write();
  plotHistograms(gsfTracksGoodGenMatch_genEta, gsfTracksBadGenMatch_genEta, "gsfTracksGenMatch_genEta", "genEta")->Write();
  plotEfficiency(matchEfficiencies, labels, "trackEfficiency")->Write();
  
  //2D plots
  plotHistogram(generalTracksGoodGenMatch_ptVsEta, "eta", "pt")->Write();
  plotHistogram(generalTracksBadGenMatch_ptVsEta, "eta", "pt")->Write();
  plotHistogram(generalTracksBadGenMatch_ptVsDeltaR, "deltaR", "pt")->Write();
  plotHistogram(generalTracksBadGenMatch_etaVsDeltaR, "deltaR", "eta")->Write();
  plotHistogram(gsfTracksGoodGenMatch_ptVsEta, "eta", "pt")->Write();
  plotHistogram(gsfTracksBadGenMatch_ptVsEta, "eta", "pt")->Write();
  plotHistogram(gsfTracksBadGenMatch_ptVsDeltaR, "deltaR", "pt")->Write();
  plotHistogram(gsfTracksBadGenMatch_etaVsDeltaR, "deltaR", "eta")->Write();
  plotHistogram(generalTracksGoodGenMatch_ptVsgenPt, "genPt", "pt")->Write();
  plotHistogram(generalTracksGoodGenMatch_etaVsgenEta, "genEta", "eta")->Write();
  plotHistogram(generalTracksBadGenMatch_ptVsgenPt, "genPt", "pt")->Write();
  plotHistogram(generalTracksBadGenMatch_etaVsgenEta, "genEta", "eta")->Write();
  plotHistogram(gsfTracksGoodGenMatch_ptVsgenPt, "genPt", "pt")->Write();
  plotHistogram(gsfTracksGoodGenMatch_etaVsgenEta, "genEta", "eta")->Write();
  plotHistogram(gsfTracksBadGenMatch_ptVsgenPt, "genPt", "pt")->Write();
  plotHistogram(gsfTracksBadGenMatch_etaVsgenEta, "genEta", "eta")->Write();
  //plotHistogram(testHist, "pT", "M")->Write();
  */
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

TCanvas* plotHistograms(std::vector<TH1D> &histograms, const std::vector<std::string> &labels, std::string xlabel) {

  histograms[0].SetStats(0);
  histograms[0].SetTitle("");
  // Define the color vector
  std::vector<int> colors = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-2,kBlack,kAzure+2,kGray};

  // Create the canvas
  TCanvas* canvas = new TCanvas(("can_"+xlabel).c_str(), "Histogram Plot", 800, 600);

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
