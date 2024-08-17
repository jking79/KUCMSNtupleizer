#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/Plot.h"

TCanvas* Plot::plotHistogram(TH1D& hist1D, const std::string& xlabel, const std::string& ylabel) {
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
  canvas->SetGrid();
  //canvas->SetGridy();
   
  return canvas;
}

TCanvas* Plot::plotHistogram(TH2D& hist2D, const std::string& xlabel, const std::string& ylabel, const std::string& zlabel) {

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

TCanvas* Plot::plotHistograms(TH1D &signal, TH1D &background, const std::string& name, const std::string& xlabel, const std::string& ylabel) {

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
  background.Draw("histsame");
  legend.Draw("same");

  // Save the plot as a .pdf
  std::string pdfName = name + ".pdf";
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy();

  return canvas;
}

TCanvas* Plot::plotHistograms(std::vector<TH1D> &histograms, const std::vector<std::string> &labels, const std::string name, std::string xlabel) {

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

TCanvas* Plot::plotEfficiency(TEfficiency *efficiency, std::string name, std::string xlabel) {

  TCanvas *canvas = new TCanvas(("can_"+name).c_str(), ("can_"+name).c_str(), 800, 600);
  canvas->SetLeftMargin(0.15);
  canvas->SetGrid();

  TMultiGraph* mg = new TMultiGraph();
  efficiency->SetMarkerStyle(20);
  efficiency->SetStatisticOption(TEfficiency::kFNormal);
  mg->Add(efficiency->CreateGraph());
  mg->Draw("ap");
  mg->SetMinimum(0.);
  mg->SetMaximum(1.);
  mg->GetXaxis()->SetTitle(xlabel.c_str());
  mg->GetXaxis()->SetTitleOffset(1.25);
  mg->GetYaxis()->SetTitle("Efficiency");
  gPad->Update();
  return canvas;
  
}

TCanvas* Plot::plotEfficiency(std::vector<TEfficiency*> &efficiencies, std::vector<std::string> &labels, std::string name, std::string xlabel) {

  std::vector<int> markers = {20,21,22,23,29,33,34,43,49};
  std::vector<int> colors = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-2,kBlack,kAzure+2,kGray};

  TMultiGraph* mg = new TMultiGraph();
  TCanvas *canvas = new TCanvas(("can_"+name).c_str(), ("can_"+name).c_str(), 800, 600);
  TLegend* legend = new TLegend(0.75, 0.18, 0.89, 0.32);
  canvas->SetLeftMargin(0.15);
  canvas->SetGrid();

  for(size_t i = 0; i < efficiencies.size(); i++) {
    efficiencies[i]->SetMarkerStyle(markers[i]);
    efficiencies[i]->SetStatisticOption(TEfficiency::kFNormal);
    efficiencies[i]->SetMarkerColor(colors[i]);
    efficiencies[i]->SetLineColor(colors[i]);
    mg->Add(efficiencies[i]->CreateGraph());
  }

  mg->Draw("ap");
  mg->SetMinimum(0.);
  mg->SetMaximum(1.);
  mg->GetXaxis()->SetTitle(xlabel.c_str());
  mg->GetXaxis()->SetTitleOffset(1.25);
  mg->GetYaxis()->SetTitle("Efficiency");
  legend->AddEntry(efficiencies[0], labels[0].c_str(), "lep");

  for(size_t i = 1; i < efficiencies.size(); i++) {
    legend->AddEntry(efficiencies[i], labels[i].c_str(), "lep");
  }

  legend->Draw("same");

  gPad->Update();

  return canvas;
  
}
