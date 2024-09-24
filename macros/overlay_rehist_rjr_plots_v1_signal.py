#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t, norm ):

    first = True
    paramn = []
    parnerror = []
    paramc = []
    parcerror = []
    params = []
    parserror = []
    thebinmid = []
    thebinerror = []
    f1 = []
    h1 = []
    n = 0

    setTDRStyle()
    gROOT.SetBatch(True)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTopMargin(0.09)

    if layout['logx'] : gStyle.SetOptLogx(1)
    if layout['logy'] : gStyle.SetOptLogy(1)
    c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
    c1.cd()
    #if layout['logx'] : c1.SetLogx()
    #if layout['logy'] : c1.SetLogy()
    c1.SetGridx(1)
    c1.SetGridy(1)
    c1.Update()
    c1.Draw()

    legend = TLegend(l[0],l[1],l[2],l[3]);
#    legend.SetHeader(layout['legtitle'], '')
    legend.SetName('legend')
    gStyle.SetLegendFont(42)
    gStyle.SetLegendTextSize(0.03)
    legend.SetBorderSize(0)
    #legend.SetLineColor(kBlack)

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    for histname, tree, infile, lego  in hist_list :
    
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        orighist = f1[n].Get(hist)
        #htitle = 'hist' + str(n)
        #lenmybins = int(orighist.GetNbinsX())
        scale = orighist.Integral()
        if scale == 0 : scale == 1
        if norm : orighist.Scale(1/scale)
        h1.append(orighist)
        #h1.append(orighist.Rebin(2))
        #norm = 1;
        #if norm == 0 : norm = 1

        #for bn in range(1,lenmybins-1):
        #    #binmid1 = float(orighist1.GetBinCenter(bn))
        #    obinval1 = float(h1[n].GetBinContent(bn)/norm)
        #    #binwidth1 = float(orighist1.GetBinWidth(bn))
        #    #nextbin = float(h1[n].GetBinContent(bn+1)/norm)
        #    h1[n].SetBinContent(bn,obinval1)

        h1[n].UseCurrentStyle()
        #h1[n].SetMarkerStyle(n+25)
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        h1[n].SetLineWidth(2)
        h1[n].SetFillColor(k[n])
        #h1[n].SetFillColorAlpha(k[n],0.2)
        h1[n]. SetFillStyle(3002)
        #msz = 0.8
        #if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        #elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        #else : h1[n].SetMarkerSize(msz)
        #h1[n].SetMarkerSize(msz)
        h1[n].SetTitle(layout['title'])
        h1[n].GetXaxis().CenterTitle(True)
        h1[n].GetXaxis().SetTitle(layout['xtitle'])
        h1[n].GetYaxis().CenterTitle(True)
        h1[n].GetYaxis().SetTitle(layout['ytitle'])
        h1[n].SetMinimum(y[0])
        h1[n].SetMaximum(y[1])
        h1[n].GetXaxis().SetRangeUser(x[0],x[1])
        if( first ) :
             first = False
             h1[n].Draw("hist")
        else :
             h1[n].Draw("hist same")

        legend.AddEntry(h1[n],lego,'lf');
        n += 1

        #End of loop

    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch
    gPad.Modified()

    #lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    lat_cms = '#bf{CMS} #it{WkInPgrs}' + ptitle[0]
    lat_title = ptitle[1]
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.DrawLatex(0.16,0.925,lat_cms);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.045);
    lat.DrawLatex(t[0],t[1],ptitle[2]);

    if layout['logx'] : c1.SetLogx()
    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.SaveAs( outname + '_' + date + '.C' )
    #c1.Show()
    c1.Close()

from overlay_hist_defs_v3 import *

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run2 AODSIM'
#Ic_legtitle = ''
xtitle = ''
#xtitle = '#Delta_{Run}'
#xtitle = 'GeV'
#xtitle = '[ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#xtitle = 'HcalTowerSumEtBcConeDR04 [GeV]'
ytitle = 'a.u.'
htitle = ''
#islogx = True
islogx = False
islogy = True
#islogy = False

#---------------------------------------------------------------

#model = "GluGlu"
#model = "SqkSqk"
#model = "SqkGlu"
#model = "XinoXino"
model = "Bkgrd"

#rjrtype = "ootmet_phomet"
rjrtype = "ootmet_phojet"

sig100 = "rjr_multihist_files/KUCMS_GMSB_L100_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
sig150 = "rjr_multihist_files/KUCMS_GMSB_L150_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
sig200 = "rjr_multihist_files/KUCMS_GMSB_L200_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
sig250 = "rjr_multihist_files/KUCMS_GMSB_L250_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
sig300 = "rjr_multihist_files/KUCMS_GMSB_L300_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
sig350 = "rjr_multihist_files/KUCMS_GMSB_L350_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
sig400 = "rjr_multihist_files/KUCMS_GMSB_L400_"+model+"_v21_genSigPerfect_nSigPho1_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"

#bgQCD = "rjr_multihist_files/KUCMS_QCD_v19_BG_nSigPho0_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
#bgGJets = "rjr_multihist_files/KUCMS_GJets_v19_BG_nSigPho0_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
bgQCD = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v22_"+rjrtype+"_multiHists.root"
bgGJets = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v22_"+rjrtype+"_multiHists.root"

#y = [ 0.0001, 1000000000000.0 ]
#y = [ 0.0000001, 10000000.0 ] # vdiff, QNSum
#y = [ 0.0001, 100000000.0 ]
#y = [ 0.0000001, 10000.0 ] # vsum, asmass, QSum 
#y = [ 0.00000001, 10.0 ] # normizlized to 1 for shape comparision
#y = [ 0.00001, 1.0 ] # normizlized to 1 for shape comparision
y = [ 0.0001, 5 ] # normizlized to 1 for shape comparision

#l = [ 0.7,0.6,0.925,0.875 ] # legend position top right
l = [ 0.6,0.65,0.94,0.875 ] # legend position top right short
#l = [ 0.475,0.6,0.925,0.875 ] # legend position top middle
#l = [ 0.2,0.65,0.425,0.9 ] # legend position top left
t = [0.2,0.85,0.0,0.175,0.225] # titles position

#rhname = "rjrVSum"
rhname = "rjrVDiff"
#rhname = "ASMass"
#rhname = "rjrAX2QSum"
#rhname = "rjrAX2NQSum"
#rhname = "rjrNJetsJa"
#rhname = "rjrNJetsJb"

#x = [ 0.0, 2600 ]
x = [ -1.0, 1.0 ]
#x = [ 0.0, 12000 ]
#x = [ 0.0, 2600 ]
#x = [ 0, 1.0 ]
#x = [ 0.0, 20 ]

fhname = rhname
xtitle = fhname
outname = 'llpa_rjr_' + model + '_' + rjrtype + '_' + fhname
#outname = 'llpa_rjr_op_' + fhname
layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
ptitle=[' 2017', '137 fb', model + ' ' + rjrtype ]
#ptitle=[' 2017','137 fb','']

#leg100 = "M_{\Chi_{1}} = 139 GeV"
#leg200 = "M_{\Chi_{1}} = 285 GeV"
#leg300 = "M_{\Chi_{1}} = 430 GeV"
#leg400 = "M_{\Chi_{1}} = 576 GeV"
leg100 = "M_{glu},M_{\Chi_{1}} = 838, 139  GeV"
leg150 = "M_{glu},M_{\Chi_{1}} = 1207, 212  GeV"
leg200 = "M_{glu},M_{\Chi_{1}} = 1565, 285  GeV"
leg250 = "M_{glu},M_{\Chi_{1}} = 1915, 358  GeV"
leg300 = "M_{glu},M_{\Chi_{1}} = 2260, 430  GeV"
leg350 = "M_{glu},M_{\Chi_{1}} = 2599, 503  GeV"
leg400 = "M_{glu},M_{\Chi_{1}} = 2935, 576  GeV"
#leg100 = "M_{sqk, \Chi_{1}} = 1119, 139 GeV"
#leg200 = "M_{sqk, \Chi_{1}} = 2122, 285 GeV"
#leg300 = "M_{sqk, \Chi_{1}} = 3084, 430 GeV"
#leg400 = "M_{sqk, \Chi_{1}} = 4025, 576 GeV"
#leg100 = "M_{glu, sqk, \Chi_{1}} = 838, 1119, 139 GeV"
#leg200 = "M_{glu, sqk, \Chi_{1}} = 1565, 2122, 285 GeV"
#leg300 = "M_{glu, sqk, \Chi_{1}} = 2260, 3084, 430 GeV"
#leg400 = "M_{glu, sqk, \Chi_{1}} = 2935, 4025, 576 GeV"

inhistlist = [

            #[rhname,'',sig100,leg100],
            #[rhname,'',sig150,leg150],
            #[rhname,'',sig200,leg200],
            #[rhname,'',sig250,leg250],
            ##[rhname,'',sig300,leg300],

            [rhname,'',bgGJets,"GJets"],
            [rhname,'',bgQCD,"QCD"],

            ]

#norm = False
norm = True

dostack(inhistlist, outname, date, layout, ptitle,  y, x, l, t, norm)

#ptitle=[' 2022 IOV5 359421-360089','','#splitline{EBEB}{CC Ave RH Time by Channel}'] #{GT 106X_dataRun2_v28}'
#y = [ 4.5, 0.5 ]
#x = [ 200.0, 700.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#t = [0.2,0.825,0.0,0.175,0.225]
#outname = 'downloads/tr_hl_r3_iov5cali_v7'
#dostack(hl_r3_iov5cali_v7, outname, date, Ic_layout, ptitle,  y, x, l, t)

#    legend = TLegend(0.25,0.20,0.52,0.525); # bottom left
#    legend = TLegend(0.4,0.205,0.6,0.525);   # bottom middle
#    legend = TLegend(0.4,0.60,0.6,0.90);   # top middle
#    legend = TLegend(0.645,0.50,0.825,0.9);   # top mid right
#    legend = TLegend(0.605,0.50,0.945,0.9);   # top right very wide
#    legend = TLegend(0.705,0.50,0.945,0.9);   # top right wide 
#    legend = TLegend(0.745,0.50,0.925,0.9);   # top right
#    legend = TLegend(0.745,0.40,0.925,0.9);   # top right tall
#    legend = TLegend(0.650,0.375,0.925,0.875);   # top mid right wide
#    legend = TLegend(0.62,0.60,0.8,0.9);   # top right
#    legend = TLegend(0.65,0.60,0.9,0.90);   # top right large

#      g_2->GetYaxis()->SetMoreLogLabels();
#      g_2->GetYaxis()->SetNoExponent();


