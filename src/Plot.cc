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
  CMSmark();
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

  CMSmark();
  // Save the plot as a .pdf                                                                                                                                                    
  std::string pdfName = name + ".pdf";
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogz();
  //canvas->SaveAs(pdfName.c_str());                                                                                                                                            
  
  return canvas;
    
}

TCanvas* Plot::plotHistogramYields(TH2D& hist2D, const std::string& xlabel, const std::string& ylabel, const std::string& zlabel) {
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

    // Add text labels for bin content
    TText* t = new TText();
    t->SetTextAlign(22);
    t->SetTextFont(42);
    t->SetTextSize(0.05);

    Int_t nBinsX = hist2D.GetNbinsX();
    Int_t nBinsY = hist2D.GetNbinsY();

    for (Int_t i = 1; i <= nBinsX; i++) {
        for (Int_t j = 1; j <= nBinsY; j++) {
            Double_t content = hist2D.GetBinContent(i, j);
            if (content > 0) {  // Only print non-zero content
                Double_t x = hist2D.GetXaxis()->GetBinCenter(i);
                Double_t y = hist2D.GetYaxis()->GetBinCenter(j);
                t->DrawText(x, y, Form("%.0f", content));
            }
        }
    }
    CMSmark();
    // Save the plot as a .pdf
    std::string pdfName = name + ".pdf";
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->SetLogz();
    //canvas->SaveAs(pdfName.c_str());

    return canvas;
}

/*
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
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(signal, "Signal", "l");
  legend->AddEntry(background, "Background", "l");

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
  legend->Draw("SAME");
  CMSmark();

  // Save the plot as a .pdf
  std::string pdfName = name + ".pdf";
  canvas->SetGrid();
  canvas->SetLogy();

  return canvas;
}
*/

TCanvas* Plot::plotHistograms(TH1D &signal, TH1D &background, const std::string& name, const std::string& xlabel, const std::string& ylabel) {

   signal.SetTitle("");
   signal.SetStats(0);
   // Create a canvas
   TCanvas* canvas = new TCanvas((name+"_canvas").c_str(), (name+"_canvas").c_str(), 800, 600);
   
   // Find the maximum value of the histograms
   double maxSignal = signal.GetMaximum();
   double maxBackground = background.GetMaximum();
   double maxY = std::max(maxSignal, maxBackground);
    
   // Set the range in y accordingly
   signal.SetMaximum(1.5 * maxY);
   signal.SetMinimum(0.9);
   
   // Style signal histogram
   //signal.SetMarkerStyle(20);
   //signal.SetMarkerSize(0.5);
   //signal.SetMarkerColor(kBlue+2);
   signal.SetLineColor(kBlue+2);
   signal.SetLineWidth(2);
   
   // Style background histogram
   //background.SetMarkerStyle(20);
   //background.SetMarkerSize(0.5);
   //background.SetMarkerColor(kRed);
   background.SetLineColor(kRed+1);
   background.SetLineWidth(2);
   
   // Set labels and titles
   signal.GetXaxis()->SetTitle(xlabel.c_str());
   signal.GetYaxis()->SetTitle(ylabel.c_str());
   
   // Draw histograms
   signal.Draw("hist");
   background.Draw("histsame");
   
   // Create and draw legend
   TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
   legend->AddEntry(&signal, "Signal", "l");
   legend->AddEntry(&background, "Background", "l");
   //legend->Draw();
   
   // Add CMS mark and set canvas properties
   canvas->SetGrid();
   canvas->SetLogy();
   legend->Draw();
   CMSmark();
   
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
    histograms[i].SetLineWidth(2);
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
  CMSmark();
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
  CMSmark();
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

  CMSmark("", 0.16);
  
  return canvas;  
}

TCanvas* Plot::plotGraph(const TH1D& hist1D, const std::string& xlabel, const std::string& ylabel) {
  // Extract the name from the TH1D
  std::string name = hist1D.GetName();

  // Create a canvas
  TCanvas* canvas = new TCanvas((name+"_graph_canvas").c_str(), (name+"_graph_canvas").c_str(), 800, 600);

  // Create a TGraphErrors from the histogram to display uncertainties
  int numPoints = hist1D.GetNbinsX();
  TGraphErrors* graph = new TGraphErrors(numPoints);
  
  // Fill the graph with histogram bin centers, contents, and uncertainties
  for (int i = 1; i <= numPoints; i++) {
    double x = hist1D.GetBinCenter(i);
    double y = hist1D.GetBinContent(i);
    double ex = 0; // No x error
    double ey = hist1D.GetBinError(i);
    graph->SetPoint(i-1, x, y);
    graph->SetPointError(i-1, ex, ey);
  }
  
  // Set graph properties
  graph->SetMarkerStyle(20);  // Circular marker
  graph->SetMarkerSize(1.0);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);
  graph->SetLineWidth(2);
  
  // Set graph title and axis labels
  graph->SetTitle("");
  graph->GetXaxis()->SetTitle(xlabel.c_str());
  graph->GetYaxis()->SetTitle(ylabel.c_str());
  
  // Draw the graph with only markers and error bars (no connecting lines)
  graph->Draw("AP");  // "A" for axes, "P" for points
  
  // Add CMS mark
  CMSmark();
  
  // Set grid
  canvas->SetGrid();
  
  return canvas;
}

TCanvas* Plot::plotRatio(const TH1D& hist1D_num, const TH1D& hist1D_denom, const std::string& name, const std::string& xlabel) {

  // Create a canvas
  TCanvas* canvas = new TCanvas((name+"_ratio_canvas").c_str(), (name+"_ratio_canvas").c_str(), 800, 600);
  
  // Create a new histogram for the ratio
  TH1D* ratio = (TH1D*)hist1D_num.Clone((name+"_ratio").c_str());
  
  // Loop through bins and manually calculate ratio and propagate errors correctly
  for (int i = 1; i <= ratio->GetNbinsX(); i++) {
    // Get bin values and errors
    double num_val = hist1D_num.GetBinContent(i);
    double num_err = hist1D_num.GetBinError(i);
    double denom_val = hist1D_denom.GetBinContent(i);
    double denom_err = hist1D_denom.GetBinError(i);
    
    // Calculate ratio
    double ratio_val = 0;
    double ratio_err = 0;
    
    // Avoid division by zero
    if (denom_val != 0) {
      ratio_val = num_val / denom_val;
      
      ratio_err = ratio_val * sqrt(
	  pow(num_err / num_val, 2) + 
	  pow(denom_err / denom_val, 2)
      );
    }
      // Set the bin content and error in the ratio histogram
    ratio->SetBinContent(i, ratio_val);
    ratio->SetBinError(i, ratio_err);
  }

  // Create a TGraphErrors from the ratio histogram to display uncertainties
  int numPoints = ratio->GetNbinsX();
  TGraphErrors* graph = new TGraphErrors(numPoints);
  
  // Fill the graph with histogram bin centers, contents, and uncertainties
  for (int i = 1; i <= numPoints; i++) {
    double x = ratio->GetBinCenter(i);
    double y = ratio->GetBinContent(i);
    double ex = 0; // No x error
    double ey = ratio->GetBinError(i);
    
    graph->SetPoint(i-1, x, y);
    graph->SetPointError(i-1, ex, ey);
  }
  
  // Set graph properties
  graph->SetMarkerStyle(20); // Circular marker
  graph->SetMarkerSize(1.0);
  graph->SetMarkerColor(kOrange+1);
  graph->SetLineColor(kOrange+1);
  graph->SetLineWidth(2);
  
  // Set graph title and axis labels
  graph->SetTitle("");
  graph->GetXaxis()->SetTitle(xlabel.c_str());
  graph->GetYaxis()->SetTitle("Ratio");
  
  // Draw the graph with only markers and error bars (no connecting lines)
  graph->Draw("AP"); // "A" for axes, "P" for points
  
  // Add CMS mark
  CMSmark();
  
  // Set grid
  canvas->SetGrid();
  
  // Clean up
  delete ratio;
    
  return canvas;
}

void Plot::CMSmark(TString plotTitle, const double xpos) {
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(0.51,0.91,plotTitle);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(xpos, 0.91, "#bf{CMS} #it{Preliminary}");
}
