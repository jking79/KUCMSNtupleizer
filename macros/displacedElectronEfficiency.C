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
#include "llpgtree.h"
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
TCanvas* plotHistograms(TH1D &signal, TH1D &background, const std::string& name, const std::string& xlabel, const std::string& ylabel = "Events");
TCanvas* plotHistogram(TH1D& hist1D, const std::string& xlabel, const std::string& ylabel = "Events");
TCanvas* plotHistogram(TH2D& hist2D, const std::string& xlabel, const std::string& ylabel, const std::string& zlabel = "Events");
TCanvas* plotEfficiency(std::vector<TEfficiency*> &efficiencies, std::vector<string> &labels, string name, string xlabel);

void displacedElectronEfficiency() {

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
  // track gen electron matching
  //TH1D deltaRGenMatchedGeneralTracksSC("deltaRGenMatchedGeneralTracksSC", "deltaRGenMatchedGeneralTracksSC", 100, 0, 5);
  //TH2D ptPercentDiffGsfTracksVsDeltaR("ptPercentDiffGsfTracksVsDeltaR", "ptPercentDiffGsfTracksVsDeltaR", 100, 0, 0.1, 100, -1, 3);
  const int nbinsEff(60);
  const float maxPtEff(150);
  TH1D genElectron_pt("genElectron_pt", "genElectron_pt", 75, 0., 150);
  TH1D unmatched_deltaR("unmatched_deltaR", "unmatched_deltaR", 100, 0., 1);

  std::vector<TH1D> unmatchedID_deltaR(11);
  for(int i = 0; i < 11; i++) {
    unmatchedID_deltaR[i] = TH1D(("unmatchedID"+particleTypeNames[i]+"_deltaR").c_str(), ("unmatchedID"+particleTypeNames[i]+"_deltaR").c_str(), 100, 0., 0.1);
  }

  TH1D unmatchedIDpion_pt("unmatchedIDpion_pt", "unmatchedIDpion_pt", nbinsEff, 0., maxPtEff);
  TH1D unmatchedIDproton_pt("unmatchedIDproton_pt", "unmatchedIDproton_pt", nbinsEff, 0., maxPtEff);
  TH1D unmatchedIDkaon_pt("unmatchedIDkaon_pt", "unmatchedIDkaon_pt", nbinsEff, 0., maxPtEff);
  TH1D unmatchedIDmuon_pt("unmatchedIDmuon_pt", "unmatchedIDmuon_pt", nbinsEff, 0., maxPtEff);
  TH1D unmatchedIDelectron_pt("unmatchedIDelectron_pt", "unmatchedIDelectron_pt", nbinsEff, 0., maxPtEff);
  TH1D unmatchedID_pt("unmatchedID_pt", "unmatchedID_pt", nbinsEff, 0., maxPtEff);
  
  // Pt efficiencies
  TH1D allElectron_pt("allElectron_pt","allElectron_pt", nbinsEff, 0., maxPtEff);
  TH1D sigElectron_pt("sigElectron_pt","sigElectron_pt", nbinsEff, 0., maxPtEff);
  TH1D tauElectron_pt("tauElectron_pt","tauElectron_pt", nbinsEff, 0., maxPtEff);
  TH1D lightElectron_pt("lightElectron_pt","lightElectron_pt", nbinsEff, 0., maxPtEff);
  TH1D heavyElectron_pt("heavyElectron_pt","heavyElectron_pt", nbinsEff, 0., maxPtEff);
  TH1D conversionElectron_pt("conversionElectron_pt","conversionElectron_pt", nbinsEff, 0., maxPtEff);
  TH1D unmatchedElectron_pt("unmatchedElectron_pt","unmatchedElectron_pt", nbinsEff, 0., maxPtEff);

  const int nbinsEnergyEff(80);
  const float maxEnergyEff(200);
  // Eta efficiencies
  TH1D allElectron_E("allElectron_E","allElectron_E", nbinsEnergyEff, 0., maxPtEff);
  TH1D sigElectron_E("sigElectron_E","sigElectron_E", nbinsEnergyEff, 0., maxPtEff);
  TH1D tauElectron_E("tauElectron_E","tauElectron_E", nbinsEnergyEff, 0., maxPtEff);
  TH1D lightElectron_E("lightElectron_E","lightElectron_E", nbinsEnergyEff, 0., maxPtEff);
  TH1D heavyElectron_E("heavyElectron_E","heavyElectron_E", nbinsEnergyEff, 0., maxPtEff);
  TH1D conversionElectron_E("conversionElectron_E","conversionElectron_E", nbinsEnergyEff, 0., maxPtEff);
  TH1D unmatchedElectron_E("unmatchedElectron_E","unmatchedElectron_E", nbinsEnergyEff, 0., maxPtEff);  

  TH2D gsfTracksBadGenMatch_etaVsDeltaR("gsfTracksBadGenMatch_etaVsDeltaR", "gsfTracksBadGenMatch_etaVsDeltaR", 100, 0., 1., 100, -4., 4);

  // Efficiency histograms
  const int nBins = 50;
  const double xmin = 0.;
  const double xmax = 100;
  
  // numerator
  TH1D genMatchedTrack_dxy("genMatchedTrack_dxy", "genMatchedTrack_dxy", nBins, xmin, xmax); 
  TH1D genMatchedGeneralTrack_dxy("genMatchedGeneralTrack_dxy", "genMatchedGeneralTrack_dxy", nBins, xmin, xmax);
  TH1D genMatchedGsfTrack_dxy("genMatchedGsfTrack_dxy", "genMatchedGsfTrack_dxy", nBins, xmin, xmax);

  // denominator
  TH1D acceptedGenElectron_dxy("acceptedGenElectron_dxy", "acceptedGenElectron_dxy", nBins, xmin, xmax);
  
  // Test
  //TH1D test1DSig("test1DSig", "test1DSig", nBins, x0, xf);
  //TH1D test1DBkg("test1DBkg", "test1DBkg", nBins, x0, xf);

  TH2D testHist("testHist", "testHist", 60, 0., 120., 40, 0., 200.);
  //TH1D vertexXSignal("", "", 100, 0, 10);
  //TH1D ("", "", 100, -10, 10);
  
  //std::string filename = "test_ntuple.root";
  //std::string filename = "test2.root";

  //std::string filename = "test.root";
  //std::string filename = "root/GMSB_L-150TeV_Ctau*_Fall17.root";
  std::string filename = "root/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-*_Fall17.root";
  //std::string filename = "gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root";//"root/GMSB_L-150TeV_Ctau*_Fall17.root";
  TChain* tch = new TChain("tree/llpgtree");
  tch->Add(filename.c_str());

  llpgtree handle = llpgtree(tch);

  int totalGenAcceptanceAndPt(0), totalGen(0), totalGenInAcceptance(0), totalDYToEE(0), totalGenPtCut(0);
  int passDeltaRGeneralTrack(0), passDeltaRGsfTrack(0), totalGeneralTracks(0), totalGsfTracks(0);
  
  // Loop over ntuple entries
  const int nEntries = tch->GetEntries();
  for(int i = 0; i < nEntries; i++){

    if(i%1000 == 0) {
      fprintf(stdout, "\r  Processed events: %5d of %5d ", i, nEntries);
      fflush(stdout);
    }

    //cout << "\nThis is event " << i << endl;

    handle.fChain->GetEntry(i);

    const int nGenElectrons(handle.DisplacedElectron_nGenSig);

    if(nGenElectrons < 2)
      continue;

    for(int g = 0; g < nGenElectrons; g++) {
      const int genIndex(handle.DisplacedElectron_genSigIndex->at(g));
      const float genVx(handle.Gen_vx->at(genIndex));
      const float genVy(handle.Gen_vy->at(genIndex));
      const float genDxy(GetDXY(genVx, genVy));
      acceptedGenElectron_dxy.Fill(genDxy);
    }

    const int nElectrons(handle.DisplacedElectron_nTotal);

    for(int e = 0; e < nElectrons; e++) {

      const bool isSignal(handle.DisplacedElectron_isSignal->at(e));
      const bool isTau(handle.DisplacedElectron_isTau->at(e));
      const bool isLight(handle.DisplacedElectron_isLight->at(e));
      const bool isHeavy(handle.DisplacedElectron_isHeavy->at(e));
      const bool isConversion(handle.DisplacedElectron_isConversion->at(e));
      const bool isUnmatched(handle.DisplacedElectron_isUnmatched->at(e));
      const bool isGsfTrack(handle.DisplacedElectron_isGsfTrack->at(e));
      const bool isGeneralTrack(handle.DisplacedElectron_isGeneralTrack->at(e));

      const int genIndex(handle.DisplacedElectron_indexGen->at(e));
      const float pt(handle.DisplacedElectron_pt->at(e));
      const float energy(handle.DisplacedElectron_energy->at(e));
      const float deltaR(handle.DisplacedElectron_dRGen->at(e));

      allElectron_pt.Fill(pt);
      allElectron_E.Fill(energy);

      if(genIndex > 0) {
	const int genPdgID(handle.Gen_pdgId->at(genIndex));
	const float genPt(handle.Gen_pt->at(genIndex));
	const float genVx(handle.Gen_vx->at(genIndex));
	const float genVy(handle.Gen_vy->at(genIndex));
	const float genDxy(GetDXY(genVx, genVy));

	if(isSignal) {
	  sigElectron_pt.Fill(pt);
	  sigElectron_E.Fill(energy);
	  genMatchedTrack_dxy.Fill(genDxy);
	  if(isGsfTrack)
	    genMatchedGsfTrack_dxy.Fill(genDxy);
	  if(isGeneralTrack)
            genMatchedGeneralTrack_dxy.Fill(genDxy);

	  const int trackIndex(handle.DisplacedElectron_indexTrack->at(e));
	  const int trackGenIndex(handle.ECALTrack_genIndex->at(trackIndex));
	  const int pdgID(handle.Gen_pdgId->at(trackGenIndex));
	  const float trackGenDeltaR(handle.ECALTrack_genDeltaR->at(trackIndex));
	  const float trackPt(handle.ECALTrack_pt->at(trackIndex));
	  /*
	  if(genIndex != trackGenIndex) {
	    cout << "event " << e << endl;
	    cout << "\nelectron genIndex: " << genIndex << ", genIndex from track: " << trackGenIndex << endl; 
	    cout << "electron pdgID: " << genPdgID << ", track pdgID: " << pdgID << endl;
	    cout << "electron DR = " << deltaR << ", track DR = " << trackGenDeltaR << endl;
	    cout << "electron pT = " << pt << ", track pt = " << trackPt << endl;
	    cout << "electron gen pt = " << handle.Gen_pt->at(genIndex) << ", track gen pt = " << handle.Gen_pt->at(trackGenIndex) << endl;
	  }
	  */
	}
	if(isTau) {
          tauElectron_pt.Fill(pt);
	  tauElectron_E.Fill(energy);
	}
	if(isLight) {
          lightElectron_pt.Fill(pt);
	  lightElectron_E.Fill(energy);
	}
	if(isHeavy) {
          heavyElectron_pt.Fill(pt);
	  heavyElectron_E.Fill(energy);
	}
	if(isConversion) {
          conversionElectron_pt.Fill(pt);
	  conversionElectron_E.Fill(energy);
	}
      }
      if(isUnmatched) {
	unmatchedElectron_pt.Fill(pt);
	unmatchedElectron_E.Fill(energy);
      
	const int trackIndex(handle.DisplacedElectron_indexTrack->at(e));
	const int trackGenIndex(handle.ECALTrack_genIndex->at(trackIndex));
	const float trackGenDeltaR(handle.ECALTrack_genDeltaR->at(trackIndex));
	const int pdgID(handle.Gen_pdgId->at(trackGenIndex));

	particleID id(ClassifyParticle(pdgID));
	//if(ClassifyParticle(pdgID) == kUnmatched)
	//cout << "unmatched track pdgID = " << pdgID << ", match deltaR = " << trackGenDeltaR << endl;

	if(trackGenDeltaR < 0.005) {
	  if(id == kPion)
	    unmatchedIDpion_pt.Fill(pt);
	  else if(id == kProtons)
	    unmatchedIDproton_pt.Fill(pt);
	  else if(id == kKaon)
	    unmatchedIDkaon_pt.Fill(pt);
	  else if(id == kElectrons)
	    unmatchedIDelectron_pt.Fill(pt);
	  else if(id == kMuon)
	    unmatchedIDmuon_pt.Fill(pt);
	  else 
	    unmatchedID_pt.Fill(pt);
	}
	else unmatchedID_pt.Fill(pt);
	    
	FillIDHists<float>(unmatchedID_deltaR, pdgID, trackGenDeltaR);
	unmatched_deltaR.Fill(trackGenDeltaR);
      }
    }
    
  }// End event loop

  // unmatched pt composition
  vector<string> idLabels = {"electrons", "muons", "pions", "kaon", "protons", "unmatched"};
  TEfficiency *unmatchedIDelectronEff_pt(new TEfficiency(unmatchedIDelectron_pt, unmatchedElectron_pt));
  TEfficiency *unmatchedIDmuonEff_pt(new TEfficiency(unmatchedIDmuon_pt, unmatchedElectron_pt));
  TEfficiency *unmatchedIDpionEff_pt(new TEfficiency(unmatchedIDpion_pt, unmatchedElectron_pt));
  TEfficiency *unmatchedIDkaonEff_pt(new TEfficiency(unmatchedIDkaon_pt, unmatchedElectron_pt));
  TEfficiency *unmatchedIDprotonEff_pt(new TEfficiency(unmatchedIDproton_pt, unmatchedElectron_pt));
  TEfficiency *unmatchedIDEff_pt(new TEfficiency(unmatchedID_pt, unmatchedElectron_pt));
  vector<TEfficiency*> idPtEfficiencies;
  idPtEfficiencies.emplace_back(unmatchedIDelectronEff_pt);
  idPtEfficiencies.emplace_back(unmatchedIDmuonEff_pt);
  idPtEfficiencies.emplace_back(unmatchedIDpionEff_pt);
  idPtEfficiencies.emplace_back(unmatchedIDkaonEff_pt);
  idPtEfficiencies.emplace_back(unmatchedIDprotonEff_pt);
  idPtEfficiencies.emplace_back(unmatchedIDEff_pt);

  // Pt Efficiency
  TEfficiency *signalPtEff(new TEfficiency(sigElectron_pt, allElectron_pt));
  TEfficiency *tauPtEff(new TEfficiency(tauElectron_pt, allElectron_pt));
  TEfficiency *lightPtEff(new TEfficiency(lightElectron_pt, allElectron_pt));
  TEfficiency *heavyPtEff(new TEfficiency(heavyElectron_pt, allElectron_pt));
  TEfficiency *conversionPtEff(new TEfficiency(conversionElectron_pt, allElectron_pt));
  TEfficiency *unmatchedPtEff(new TEfficiency(unmatchedElectron_pt, allElectron_pt));

  vector<string> labels = {"signal", "#tau", "light", "heavy", "conversion", "unmatched"}; 
  vector<TEfficiency*> ptEfficiencies;
  ptEfficiencies.emplace_back(signalPtEff);
  ptEfficiencies.emplace_back(tauPtEff);
  ptEfficiencies.emplace_back(lightPtEff);
  ptEfficiencies.emplace_back(heavyPtEff);
  ptEfficiencies.emplace_back(conversionPtEff);
  ptEfficiencies.emplace_back(unmatchedPtEff);
  
  // Energy Efficiency
  TEfficiency *signalEnergyEff(new TEfficiency(sigElectron_E, allElectron_E));
  TEfficiency *tauEnergyEff(new TEfficiency(tauElectron_E, allElectron_E));
  TEfficiency *lightEnergyEff(new TEfficiency(lightElectron_E, allElectron_E));
  TEfficiency *heavyEnergyEff(new TEfficiency(heavyElectron_E, allElectron_E));
  TEfficiency *conversionEnergyEff(new TEfficiency(conversionElectron_E, allElectron_E));
  TEfficiency *unmatchedEnergyEff(new TEfficiency(unmatchedElectron_E, allElectron_E));

  vector<TEfficiency*> EnergyEfficiencies;
  EnergyEfficiencies.emplace_back(signalEnergyEff);
  EnergyEfficiencies.emplace_back(tauEnergyEff);
  EnergyEfficiencies.emplace_back(lightEnergyEff);
  EnergyEfficiencies.emplace_back(heavyEnergyEff);
  EnergyEfficiencies.emplace_back(conversionEnergyEff);
  EnergyEfficiencies.emplace_back(unmatchedEnergyEff);

  // Tranverse Impact Parameter Efficiency
  std::vector<string> genLabels = {"Displaced Electron", "Gsf Track", "General Track"};
  TEfficiency *genMatchedTrackEfficiency = new TEfficiency(genMatchedTrack_dxy, acceptedGenElectron_dxy);
  TEfficiency *genMatchedGsfTrackEfficiency = new TEfficiency(genMatchedGsfTrack_dxy, acceptedGenElectron_dxy);
  TEfficiency *genMatchedGeneralTrackEfficiency = new TEfficiency(genMatchedGeneralTrack_dxy, acceptedGenElectron_dxy);
  vector<TEfficiency*> matchEfficiencies;
  matchEfficiencies.emplace_back(genMatchedTrackEfficiency);
  matchEfficiencies.emplace_back(genMatchedGsfTrackEfficiency);
  matchEfficiencies.emplace_back(genMatchedGeneralTrackEfficiency);

  TFile outfile("matching_histograms.root", "RECREATE");
  outfile.cd();

  plotEfficiency(idPtEfficiencies, idLabels, "UnmatchedComposition", "pT")->Write();
  plotEfficiency(ptEfficiencies, labels, "ptEfficiency", "pT")->Write();
  plotEfficiency(EnergyEfficiencies, labels, "EnergyEfficiency", "E")->Write();
  plotEfficiency(matchEfficiencies, genLabels, "dxyEfficiency", "dxy [cm]")->Write();

  plotHistogram(unmatched_deltaR, "deltaR")->Write();
  for(int i = 0; i < 9; i++)
    plotHistogram(unmatchedID_deltaR[i], "deltaR")->Write();

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
