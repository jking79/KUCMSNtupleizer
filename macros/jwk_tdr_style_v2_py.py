from ROOT import *

can_x_width = 600
can_y_width = 800

def fixOverlay():
        gPad.RedrawAxis()


def setTDRStyle():
        #gStyle = TStyle("tdrStyle","Style for P-TDR")

##// For the canvas:
        gStyle.SetCanvasBorderMode(0)
        gStyle.SetCanvasColor(kWhite)
        gStyle.SetCanvasDefH(can_x_width)  #//Height of canvas
        gStyle.SetCanvasDefW(can_y_width)  #//Width of canvas
        gStyle.SetCanvasDefX(0)    #//POsition on screen
        gStyle.SetCanvasDefY(0)

#// For the Pad:
        gStyle.SetPadBorderMode(0)
        #// gStyle.SetPadBorderSize(Width_t size = 1) 
        gStyle.SetPadColor(kWhite)
        gStyle.SetPadGridX(False)
        gStyle.SetPadGridY(False)
        gStyle.SetGridColor(0)
        gStyle.SetGridStyle(3)
        gStyle.SetGridWidth(1)

#// For the frame:
        gStyle.SetFrameBorderMode(0)
        gStyle.SetFrameBorderSize(1)
        gStyle.SetFrameFillColor(0)
        gStyle.SetFrameFillStyle(0)
        gStyle.SetFrameLineColor(1)
        gStyle.SetFrameLineStyle(1)
        gStyle.SetFrameLineWidth(1)


#// set the paper & margin sizes
        gStyle.SetPaperSize(20,26);
        #gStyle.SetPadTopMargin(0.09);
        gStyle.SetPadTopMargin(0.095);
        gStyle.SetPadRightMargin(0.05);
        gStyle.SetPadBottomMargin(0.12);
        gStyle.SetPadLeftMargin(0.15);

#// For the histo:
        #gStyle.SetHistFillColor(1) 
        #gStyle.SetHistFillStyle(0) 
        gStyle.SetHistLineColor(1)
        gStyle.SetHistLineStyle(0)
        gStyle.SetHistLineWidth(1)
        #gStyle.SetLegoInnerR(Float_t rad = 0.5) 
        #gStyle.SetNumberContours(Int_t number = 20) 
        gStyle.SetMarkerStyle(8)
        gStyle.SetMarkerSize(0.8)
        gStyle.SetEndErrorSize(2)
        #// gStyle.SetErrorMarker(20) 
        gStyle.SetErrorX(0.5)

#//For the fit/function:
        gStyle.SetOptFit(0)
        gStyle.SetFitFormat("5.4g")
        gStyle.SetFuncColor(2)
        gStyle.SetFuncStyle(1)
        gStyle.SetFuncWidth(1)

#//For the date:
        gStyle.SetOptDate(0)
        #// gStyle.SetDateX(Float_t x = 0.01) 
        #// gStyle.SetDateY(Float_t y = 0.01) 

#// For the statistics box:
        gStyle.SetOptFile(0)
        gStyle.SetOptStat(0)  #// To display the mean and RMS:   SetOptStat("mr");
        gStyle.SetStatColor(kWhite)
        gStyle.SetStatFont(42)
        gStyle.SetStatFontSize(0.025)
        gStyle.SetStatTextColor(1)
        gStyle.SetStatFormat("6.4g")
        gStyle.SetStatBorderSize(1)
        gStyle.SetStatH(0.1)
        gStyle.SetStatW(0.15)
        #// gStyle.SetStatStyle(Style_t style = 1001) 
        #gStyle.SetStatX( 0.825 ) 
        #gStyle.SetStatY( 0.85 ) 


#// For the Global title:
        gStyle.SetOptTitle(0)
        gStyle.SetTitleFont(42)
        gStyle.SetTitleColor(1)
        gStyle.SetTitleTextColor(1)
        gStyle.SetTitleFillColor(10)
        gStyle.SetTitleFontSize(0.045)
        #// gStyle.SetTitleH(0)  // Set the height of the title box
        #// gStyle.SetTitleW(0)  // Set the width of the title box
        #gStyle.SetTitleX(0.5)  #// Set the position of the title box
        #gStyle.SetTitleAlign(23) 
        #gStyle.SetTitleY(0.985)  #// Set the position of the title box
        #//  gStyle.SetTitleStyle(1001) 
        #gStyle.SetTitleBorderSize(0) 

#// For the axis titles:

        gStyle.SetTitleColor(1, "XYZ")
        gStyle.SetTitleFont(42, "XYZ")
        gStyle.SetTitleSize(0.04, "XYZ")
        #gStyle.SetTitleXSize(0.06)  #// Another way to set the size?  0.02
        #gStyle.SetTitleYSize(0.06) 
        #gStyle.SetTitleXOffset(0.9)
        gStyle.SetTitleXOffset(1.0)
        #gStyle.SetTitleYOffset(1.1)
        gStyle.SetTitleYOffset(1.4)
        #//  gStyle.SetTitleOffset(1.1, "Y")  // Another way to set the Offset

#// For the axis labels:

        gStyle.SetLabelColor(1, "XYZ")
        gStyle.SetLabelFont(42, "XYZ")
        gStyle.SetLabelOffset(0.007, "XYZ")
        #gStyle.SetLabelSize(0.05, "XYZ")
        gStyle.SetLabelSize(0.02, "XYZ")

#// For the axis:

        gStyle.SetAxisColor(1, "XYZ")
        gStyle.SetStripDecimals(kTRUE)
        gStyle.SetTickLength(0.03, "XYZ")
        gStyle.SetNdivisions(505, "X")
        gStyle.SetPadTickX(1)   #// To get tick marks on the opposite side of the frame
        gStyle.SetPadTickY(1)

#// Change for log plots:
        #gStyle.SetOptLogx(0)
        #gStyle.SetOptLogy(0)
        #gStyle.SetOptLogz(0)

#// Postscript options:
        #//#gStyle.SetPaperSize(20.,20.) 
        #// gStyle.SetLineScalePS(Float_t scale = 3) 
        #// gStyle.SetLineStyleString(Int_t i, const char* text) 
        #// gStyle.SetHeaderPS(const char* header) 
        #// gStyle.SetTitlePS(const char* pstitle) 

        #// gStyle.SetBarOffset(Float_t baroff = 0.5) 
        #// gStyle.SetBarWidth(Float_t barwidth = 0.5) 
        #// gStyle.SetPaintTextFormat(const char* format = "g") 
        #// gStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0) 
        #// gStyle.SetTimeOffset(Double_t toffset) 
        #// gStyle.SetHistMinimumZero(kTRUE) 

        #//gStyle.SetHatchesLineWidth(5) 
        #//gStyle.SetHatchesSpacing(0.05)

        #// color palette
        gStyle.SetPalette(kBird)
        gStyle.SetNumberContours(255)

        #gStyle.cd()
        gROOT.ForceStyle()
