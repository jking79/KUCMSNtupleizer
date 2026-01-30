#include <iostream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TColorWheel.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TError.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TApplication.h>

#include "../include/ReducedBase_V2.hh"
#include "../include/SampleTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/Leptonic.hh"
#include "../include/FitPlotter.hh"
#include "../include/L_Cand.hh"

#include "RestFrames/RestFrames.hh"

#include "../include/mt2_bisect.hh"
mt2_bisect::mt2 mt2calc;

using namespace std;
using namespace RestFrames;

// declare global vars
string g_PlotTitle;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX = 64.;
string g_Yname;
double g_Ymin;
double g_Ymax;
double g_NY = 64.;
string g_dm_Yname;
double g_dm_Ymin;
double g_dm_Ymax;
double g_dm_YN;
string g_Label;
// root file to store output of plots
string output_root_file = "output_Plot_";
bool SavePDF = true; // whether or not to save pdfs of plots
string folder_name = "";

using namespace RestFrames;

LabRecoFrame LAB("LAB","lab");
DecayRecoFrame CM("CM","CM");
DecayRecoFrame S("S","S");
DecayRecoFrame Pa("Pa","P_{a}");
DecayRecoFrame Pb("Pb","P_{b}");
VisibleRecoFrame Ja("Ja","J_{a}");
VisibleRecoFrame Jb("Jb","J_{b}");
VisibleRecoFrame Va("Va","V_{a}");
VisibleRecoFrame Vb("Vb","V_{b}");
VisibleRecoFrame La("La","L_{a}");
VisibleRecoFrame Lb("Lb","L_{b}");
VisibleRecoFrame ISR("ISR","ISR");
InvisibleRecoFrame Ia("Ia","I_{a}");
InvisibleRecoFrame Ib("Ib","I_{b}");
SelfAssemblingRecoFrame sLa("sLa","sLa");
SelfAssemblingRecoFrame sLb("sLb","sLb");
InvisibleGroup INV("INV","Invisible System");
SetMassInvJigsaw InvM("InvM", "Set inv. system mass");
SetRapidityInvJigsaw InvEta("InvEta", "Set inv. system rapidity");
MinMassesSqInvJigsaw InvSplit("InvSplit", "INV -> I_{a} + I_{b}", 2);
CombinatoricGroup COMB_L("COMB_L","Lepton Jigsaws");
MinMassesSqCombJigsaw CombSplitSq_L("CombSplitSq_L", "Minimize M_{Sa}^{2} + M_{Sb}^{2}",2,2);

void InitRJRtree(){

  LAB.SetChildFrame(CM);
  CM.AddChildFrame(S);
  CM.AddChildFrame(ISR);
  S.AddChildFrame(Pa);
  S.AddChildFrame(Pb);
  Pa.AddChildFrame(Ia);
  Pb.AddChildFrame(Ib);
  Pa.AddChildFrame(sLa);
  Pb.AddChildFrame(sLb);
  sLa.AddChildFrame(La);
  sLb.AddChildFrame(Lb);

  LAB.InitializeTree();

  COMB_L.AddFrame(La);
  COMB_L.SetNElementsForFrame(La, 0);
  COMB_L.AddFrame(Lb);
  COMB_L.SetNElementsForFrame(Lb, 0);
  COMB_L.AddJigsaw(CombSplitSq_L);
  CombSplitSq_L.AddCombFrame(La, 0);
  CombSplitSq_L.AddCombFrame(Lb, 1);
  CombSplitSq_L.AddObjectFrame(La, 0);
  CombSplitSq_L.AddObjectFrame(Lb, 1);

  INV.AddFrame(Ia);
  INV.AddFrame(Ib);
  INV.AddJigsaw(InvM);
  INV.AddJigsaw(InvEta);
  InvEta.AddVisibleFrames(S.GetListVisibleFrames());
  
  INV.AddJigsaw(InvSplit);
  InvSplit.AddVisibleFrames(Pa.GetListVisibleFrames(), 0);
  InvSplit.AddVisibleFrames(Pb.GetListVisibleFrames(), 1);
  InvSplit.AddInvisibleFrame(Ia, 0);
  InvSplit.AddInvisibleFrame(Ib, 1);

  LAB.InitializeAnalysis();
}

void DrawRJRtree(bool treeplot_INV = true){ // treeplot_INV: whether to flip treeplots to dark mode
  TreePlot tree_plot("TreePlot","TreePlot");
  tree_plot.SetTree(LAB);
  tree_plot.Draw("ISR_LEP_tree", "Analysis Tree", treeplot_INV);
  tree_plot.SetTree(INV);
  tree_plot.Draw("ISR_LEP_inv", "Invisible Jigsaws", treeplot_INV);
  tree_plot.SetTree(COMB_L);
  tree_plot.Draw("ISR_LEP_leps", "Lep Jigsaws", treeplot_INV);

  tree_plot.WriteOutput("trees.root");
  std::cout << "Writing trees to trees.root" << endl;
}



    Long64_t start = gSystem->Now();
  RestFrames::SetStyle();
  InitRJRtree();

    DrawRJRtree(false); gApplication->Terminate(0); // end macro if making tree plot

