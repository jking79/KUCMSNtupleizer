#ifndef KUCMSNtupleizer_KUCMSNtupleizer_Plot_h
#define KUCMSNtupleizer_KUCMSNtupleizer_Plot_h
//c++ includes
#include <vector>
#include <iostream>

// ROOT includes
#include <TH1D.h>
#include <TH2D.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TLegend.h>
#include <Math/Vector4D.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>

class Plot {

 public:
  
  static TCanvas* plotHistogram(TH1D& hist1D, const std::string& xlabel, const std::string& ylabel = "Events");
  static TCanvas* plotHistogram(TH2D& hist2D, const std::string& xlabel, const std::string& ylabel, const std::string& zlabel = "Events");
  static TCanvas* plotHistograms(TH1D &signal, TH1D &background, const std::string& name, const std::string& xlabel, const std::string& ylabel = "Events");
  static TCanvas* plotHistograms(std::vector<TH1D> &histograms, const std::vector<std::string> &labels, const std::string name, std::string xlabel);
  static TCanvas* plotEfficiency(TEfficiency *efficiency, std::string name, std::string xlabel);
  static TCanvas* plotEfficiency(std::vector<TEfficiency*> &efficiencies, std::vector<std::string> &labels, std::string name, std::string xlabel);

};

#endif
