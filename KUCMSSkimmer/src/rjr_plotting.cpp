#include <iostream>
#include <string>

// ROOT
#include <TSystem.h>
#include <TApplication.h>

// RestFrames
#include "RestFrames/RestFrames.hh"

// If you truly need these later, keep them; otherwise remove to speed compile.
// #include <TFile.h>
// #include <TTree.h>
// ...

using namespace std;
using namespace RestFrames;

// ------------------------------------------------------------
// Global plotting/config vars (kept from your snippet)
// ------------------------------------------------------------
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

string output_root_file = "output_Plot_";
bool SavePDF = true;
string folder_name = "";

// ------------------------------------------------------------
// RestFrames objects (must be declared BEFORE InitRJRtree uses them)
// ------------------------------------------------------------
LabRecoFrame   LAB("LAB","Lab");
DecayRecoFrame CM("CM","CM");
DecayRecoFrame S ("S","S");

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

SetMassInvJigsaw      InvM  ("InvM",   "Set inv. system mass");
SetRapidityInvJigsaw  InvEta("InvEta", "Set inv. system rapidity");
MinMassesSqInvJigsaw  InvSplit("InvSplit", "INV -> I_{a} + I_{b}", 2);

CombinatoricGroup     COMB_L("COMB_L","Lepton Jigsaws");
MinMassesSqCombJigsaw CombSplitSq_L("CombSplitSq_L",
                                    "Minimize M_{Sa}^{2} + M_{Sb}^{2}", 2, 2);

// ------------------------------------------------------------
// Functions
// ------------------------------------------------------------
void InitRJRtree() {

  // Build decay tree
  LAB.SetChildFrame(S);
  //CM.AddChildFrame(S);
  //CM.AddChildFrame(ISR);

  S.AddChildFrame(Pa);
  S.AddChildFrame(Pb);

  Pa.AddChildFrame(Ia);
  Pb.AddChildFrame(Ib);

  Pa.AddChildFrame(Va);
  Pb.AddChildFrame(Vb);

  //sLa.AddChildFrame(La);
  //sLb.AddChildFrame(Lb);

  // Initialize tree
  LAB.InitializeTree();

  // Lepton combinatorics
  COMB_L.AddFrame(Va);
  COMB_L.SetNElementsForFrame(Va, 0);
  COMB_L.AddFrame(Vb);
  COMB_L.SetNElementsForFrame(Vb, 0);

  COMB_L.AddJigsaw(CombSplitSq_L);
  CombSplitSq_L.AddCombFrame(Va, 0);
  CombSplitSq_L.AddCombFrame(Vb, 1);
  CombSplitSq_L.AddObjectFrame(Va, 0);
  CombSplitSq_L.AddObjectFrame(Vb, 1);

  // Invisible system
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

  // Initialize analysis
  LAB.InitializeAnalysis();
}

void DrawRJRtree(bool treeplot_INV = true) {
  TreePlot tree_plot("TreePlot","TreePlot");

  tree_plot.SetTree(LAB);
  tree_plot.Draw("RJR_tree", "Analysis Tree", treeplot_INV);

  tree_plot.SetTree(INV);
  tree_plot.Draw("RJR_inv", "Invisible Jigsaws", treeplot_INV);

  tree_plot.SetTree(COMB_L);
  tree_plot.Draw("RJR_leps", "Lep Jigsaws", treeplot_INV);

  tree_plot.WriteOutput("trees.root");
  std::cout << "Writing trees to trees.root" << std::endl;
}

// ------------------------------------------------------------
// MAIN
// ------------------------------------------------------------
int main(int argc, char** argv) {
  // Needed because you call gApplication->Terminate()
  TApplication app("rjr_treeplot", &argc, argv);

  Long64_t start = gSystem->Now();

  RestFrames::SetStyle();
  InitRJRtree();
  DrawRJRtree(false);

  Long64_t stop = gSystem->Now();
  std::cout << "Done. Elapsed time [ms]: " << (stop - start) << std::endl;

  // End program immediately after making tree plots
  gApplication->Terminate(0);
  return 0; // (won't be reached, but keeps compilers happy)
}

