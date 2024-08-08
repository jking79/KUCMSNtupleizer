#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t ):

    first = True
    f1 = []
    f2 = []
    h1 = []
    n = 0

    setTDRStyle()
    gROOT.SetBatch(True)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
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

    mg = TMultiGraph();

    for histname, infile1, infile2, lego, cutside in hist_list :
    
        tree = ''
        f1.append(TFile.Open(infile1))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        f2.append(TFile.Open(infile2))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname

        orighist1 = f1[n].Get(hist)
        orighist2 = f2[n].Get(hist)
        htitle = 'hist' + str(n)
        lenmybins = int(orighist1.GetNbinsX())
        h1.append(TGraphErrors(lenmybins))
        #h1.append(TGraphErrors(lenmybins))
        norm1 = float(orighist1.Integral())
        norm2 = float(orighist2.Integral())
        if norm1 == 0 : norm1 = 1.0
        if norm2 == 0 : norm2 = 1.0
        if cutside == "no" : 
            norm1 = 1.0
            norm2 = 1.0

        for bn in range(1,lenmybins):
            binmid1 = float(orighist1.GetBinCenter(bn))
            #binmid2 = float(orighist2.GetBinCenter(bn))
            obinval1 = float(orighist1.GetBinContent(bn)/norm1)
            obinval2 = float(orighist2.GetBinContent(bn)/norm2)
            binerr1 = float(orighist1.GetBinError(bn)/norm1)
            binerr2 = float(orighist2.GetBinError(bn)/norm2)
            binwidth1 = float(orighist1.GetBinWidth(bn))
            #binwidth2 = float(orighist2.GetBinWidth(bn))
            rval = 0
            if obinval2 == 0 and obinval1 == 0 : rval = 0.0000000001
            if obinval2 == 0 and obinval1 > 0 : rval = 0.0000000001
            if obinval2 > 0 and obinval1 > 0 : rval = float(obinval1/obinval2)
            h1[n].SetPoint(bn,binmid1,rval)
            #h1[n].SetPoint(bn,binmid1,obinval1)
            #h1[n].SetPoint(bn,binmid1,obinval2)
            #if obinval1 == 0 : obinval1 = 1;
            #if obinval2 == 0 : obinval2 = 1; 
            #widtherr1 = 0 #binwidth1/(2*sqrt(obinval1))
            #widtherr2 = 0 #binwidth2/(2*sqrt(obinval2))
            #if obinval2 == 0 : obinval2 = 1
            #t1 = binerr1/obinval2
            #t2 = obinval1*binerr2/(obinval2*obinval2)
            #rerror = sqrt(t1*t1+t2*t2) 
            #h1[n].SetPointError(bn,widtherr1,rerror)
            #h1[n+1].SetPointError(bn,widtherr2,binerr2)

        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #h1[n+1].UseCurrentStyle()
        #h1[n+1].SetMarkerStyle(n+26)
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        #h1[n+1].SetLineColor(k[n+1])
        #h1[n+1].SetMarkerColor(k[n+1])
        msz = 0.8
        if( n == 1 ) : 
           h1[n].SetMarkerSize(msz+0.3)
           #h1[n+1].SetMarkerSize(msz+0.5)
        elif( n == 2 ) : 
           h1[n].SetMarkerSize(msz+0.5)
           #h1[n+1].SetMarkerSize(msz)
        elif( n == 4 ) :
           h1[n].SetMarkerSize(msz+0.5)
        else : 
           h1[n].SetMarkerSize(msz)
           #h1[n+1].SetMarkerSize(msz)
        #lego0 = lego + ' Sig Eff'
        legend.AddEntry(h1[n],lego,'epl');
        #lego1 = lego + ' Bkg Eff'
        #legend.AddEntry(h1[n+1],lego1,'epl');
        n += 1

        #End of loop
    
    for h in range(0,n):
        mg.Add(h1[h])

    mg.Draw('AP')

    #mg.UseCurrentStyle()
    #mg.SetMarkerStyle(n+25)
    #mg.SetMarkerStyle(6)
    mg.SetTitle(layout['title'])
#+';X axis '+layout['xtitle']+';Y Axis '+layout['ytitle']+';')
    mg.GetXaxis().CenterTitle(True)
    mg.GetXaxis().SetTitle(layout['xtitle'])
    mg.GetYaxis().CenterTitle(True)
    mg.GetYaxis().SetTitle(layout['ytitle'])
    mg.SetMinimum(y[0])
    mg.SetMaximum(y[1])
#   mg.GetXaxis().SetRangeUser(200.0,1100.0)
    mg.GetXaxis().SetRangeUser(x[0],x[1])
#    if layout['logx'] : mg.GetXaxis().SetMoreLogLabels()
#    if layout['logy'] : mg.GetYaxis().SetMoreLogLabels()

    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch
    gPad.Modified()

    #lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    lat_cms = '#bf{CMS} #it{WkInPgrs}' + ptitle[0]
    #lat_title = ptitle[1]+' (13 TeV)'
    lat_title = ptitle[1]
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (N/Eeff)^{2} + 2C^{2}'
    lat.SetTextSize(0.035);
    lat.SetTextFont(42);
    lat.DrawLatex(0.15,0.965,lat_cms);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.035);
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

legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run2 AODSIM'
#Ic_legtitle = ''
#xtitle = 'Signal Eff.'
#xtitle = '#Delta_{Run}'
#xtitle = 'GeV'
#xtitle = '[ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#xtitle = 'HcalTowerSumEtBcConeDR04 [GeV]'
ytitle = 'a.u.'
#ytitle = 'Eff.'
htitle = ''
#islogx = True
islogx = False
#islogy = True
islogy = False

#---------------------------------------------------------------
#hl_mc_fms_loc = [
#     #['hist_name","tree_name",hist_file_location","legend_name"],
#     #["Data_sigma","",mc_full_loc+pcal+lstfr,"Full"],
#     #["Data_sigma","",mc_multi_loc+pcal+lstfr,"Multi"],
#     #["Data_sigma","",mc_single_loc+pcal+lstfr,"Single"],
#]

#y = [ 0.001, 0.25 ]
#x = [ 0.0, 1.2 ]
l = [ 0.7,0.725,0.925,0.9 ] # legend position top right
#l = [ 0.2,0.65,0.425,0.9 ] # legend position top left
#l = [ 0.25,0.20,0.52,0.525 ] # legend position bottom left
t = [0.2,0.825,0.0,0.175,0.225] # titles position

y1 = [ 0.000001, 100000000 ]
y1l = [ 0.00001, 4 ]
#yn = [0.75,1.25 ]
#yn = [0.95,1.05 ]
yn = [0.98,1.02 ]
y = [0.95,1.05 ]
norm = 'no'
norml = 'yes'

rhname1 = [ "" ]
#rhname1 = [ "selPhoHcalTowerSumEtBcConeDR04", "pho30ptHcalTowerSumEtBcConeDR04", "pho100ptHTSumEtBcConeDR04" ]
#rhname1 = [ "rjrX2NGMean" ]
#rhname1a = [ "rjrAX2NGMean" ]
#rhname1 = [ "rjrAX2NAve" ]
#rhname1 = [ "rjrX2NVSum" ]
#rhname1a = [ "rjrAX2NVSum" ]
#x1 = [ 0.0, 1.0 ]

#rhname1 = [ "SMass" ]
#rhname1a = [ "ASMass" ]
#x1 = [ 0, 15000 ]

#rhname1 = [ "X2aMass" ]
#rhname1a = [ "AX2aMass" ]
#rhname1 = [ "X2bMass" ]
#rhname1a = [ "AX2bMass" ]
#rhname1 = [ "rjrX2GMean" ]
#rhname1a = [ "rjrAX2GMean" ]
#rhname1 = [ "rjrX2Ave" ]
#rhname1 = [ "rjrX2VSum" ]
#rhname1a = [ "rjrAX2VSum" ]
#rhname1 = [ "rjrX2QSum" ] ##  AX2Mass
rhname1a = [ "rjrAX2QSum" ]
#rhname1 = [ "rjrVSVis" ]
#rhname1a = [ "rjrVSVis" ]
#rhname1 = [ "rjrVSum" ]
#rhname1a = [ "rjrVSum" ]
#x1 = [ 0, 3600 ]
x1 = [ 0, 400 ]
#x1 = [ 0, 500 ]

#rhname1 = [ "rjrX2aPtS" ]
#rhname1 = [ "rjrX2bPtS" ]
#x1 = [ 0, 10000 ]

#rhname1 = [ "rjrNJetsJa" ]
#rhname1 = [ "rjrNJetsJb" ]
#y1 = [ 0.00001, 1.0 ]
#x1 = [ 0, 15 ]

#rhname1 = [ "rjrX2Diff" ]
#rhname1a = [ "rjrAX2Diff" ]
#rhname1 = [ "rjrVDiff" ]
#rhname1a = [ "rjrVDiff" ]
#x1 = [ -1.0, 1.0 ]

#rhname1 = [ "SCosA" ]
#x1 = [ -1.5, 1.5 ]
#y1 = [ 0.01, 1 ]

#rhname1 = [ "evtMetPt" ]
#rhname1a = [ "evtMetE" ]
#rhname1 = [ "evtRjrMetPt" ]
#rhname1a = [ "evtRawMetPt" ]
#x1 = [ 0, 10000 ]

#rhname1 = [ "rjrDphiMETV" ]
#rhname1 = [ "rjrSdphiBV" ]
#rhname1 = [ "rjrSdphiDA" ]
#x1 = [ -1.5, 1.5 ]
#y1 = [ 0.01, 1 ]

#rhname1 = [ "selPhoMBetaEql" ]
#rhname1a = [ "selPhoMBetaPmt" ]
#x1 = [ 0, 1200.0 ]

#rhname1 = [ "selPhoNMBetaEql" ]
#rhname1a = [ "selPhoNMBetaPmt" ]
#x1 = [ 0, 2 ]

#rjrgjpj = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v21_phojet_alt2_ztc4_multiHists.root"
#rjrgjpm = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v21_phomet_c1a_ztc4_multiHists.root"
#rjrgjpj = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v21_phojet_met0_ztc4_multiHists.root"
#rjrgjpm = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v21_phomet_met0_ztc4_multiHists.root"
#rjrgjpj = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v21a_phojet_ztc4_multiHists.root"
#rjrgjpm = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v21a_phomet_ztc4_multiHists.root"

rjrgjpjr = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v22_rawmet_phojet_ztc4_multiHists.root"
rjrgjpmr = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v22_rawmet_phomet_ztc4_multiHists.root"
rjrgjpjo = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v22_ootmet_phojet_ztc4_multiHists.root"
rjrgjpmo = "rjr_multihist_files/KUCMS_GJets_v19_wt2_RjrSkim_v22_ootmet_phomet_ztc4_multiHists.root"

#rjrg1pj = "rjr_multihist_files/KUCMS_GMSB_L100_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21a_phojet_ztc4_multiHists.root"
#rjrg1pm = "rjr_multihist_files/KUCMS_GMSB_L100_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21a_phomet_ztc4_multiHists.root"
#rjrg3pj = "rjr_multihist_files/KUCMS_GMSB_L350_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21a_phojet_ztc4_multiHists.root"
#rjrg3pm = "rjr_multihist_files/KUCMS_GMSB_L350_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21a_phomet_ztc4_multiHists.root"

rjrg1pjr = "rjr_multihist_files/KUCMS_GMSB_L100_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_rawmet_phojet_ztc4_multiHists.root"
rjrg1pmr = "rjr_multihist_files/KUCMS_GMSB_L100_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_rawmet_phomet_ztc4_multiHists.root"
rjrg3pjr = "rjr_multihist_files/KUCMS_GMSB_L350_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_rawmet_phojet_ztc4_multiHists.root"
rjrg3pmr = "rjr_multihist_files/KUCMS_GMSB_L350_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_rawmet_phomet_ztc4_multiHists.root"
rjrg1pjo = "rjr_multihist_files/KUCMS_GMSB_L100_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_ootmet_phojet_ztc4_multiHists.root"
rjrg1pmo = "rjr_multihist_files/KUCMS_GMSB_L100_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_ootmet_phomet_ztc4_multiHists.root"
rjrg3pjo = "rjr_multihist_files/KUCMS_GMSB_L350_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_ootmet_phojet_ztc4_multiHists.root"
rjrg3pmo = "rjr_multihist_files/KUCMS_GMSB_L350_T30_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v22_ootmet_phomet_ztc4_multiHists.root"

#rjrg1pj = "rjr_multihist_files/KUCMS_GMSB_L100_T11_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg1pm = "rjr_multihist_files/KUCMS_GMSB_L100_T11_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"
#rjrg3pj = "rjr_multihist_files/KUCMS_GMSB_L350_T11_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg3pm = "rjr_multihist_files/KUCMS_GMSB_L350_T11_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"

#rjrg1pj = "rjr_multihist_files/KUCMS_GMSB_L100_T22_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg1pm = "rjr_multihist_files/KUCMS_GMSB_L100_T22_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"
#rjrg3pj = "rjr_multihist_files/KUCMS_GMSB_L350_T22_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg3pm = "rjr_multihist_files/KUCMS_GMSB_L350_T22_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"

#rjrg1pj = "rjr_multihist_files/KUCMS_GMSB_L100_T12_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg1pm = "rjr_multihist_files/KUCMS_GMSB_L100_T12_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"
#rjrg3pj = "rjr_multihist_files/KUCMS_GMSB_L350_T12_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg3pm = "rjr_multihist_files/KUCMS_GMSB_L350_T12_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"

#rjrg1pj = "rjr_multihist_files/KUCMS_GMSB_L100_T50_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg1pm = "rjr_multihist_files/KUCMS_GMSB_L100_T50_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"
#rjrg3pj = "rjr_multihist_files/KUCMS_GMSB_L350_T50_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phojet_ztc4_multiHists.root"
#rjrg3pm = "rjr_multihist_files/KUCMS_GMSB_L350_T50_v19_genSigPerfect_nSigPho1_wt2_RjrSkim_v21_phomet_ztc4_multiHists.root"


#rjrqcpj = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v21_phojet_alt2_ztc4_multiHists.root"
#rjrqcpm = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v21_phomet_c1a_ztc4_multiHists.root"
#rjrqcpj = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v21_phojet_met0_ztc4_multiHists.root"
#rjrqcpm = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v21_phomet_met0_ztc4_multiHists.root"
#rjrqcpj = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v21a_phojet_ztc4_multiHists.root"
#rjrqcpm = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v21a_phomet_ztc4_multiHists.root"

rjrqcpjr = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v22_rawmet_phojet_ztc4_multiHists.root"
rjrqcpmr = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v22_rawmet_phomet_ztc4_multiHists.root"
rjrqcpjo = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v22_ootmet_phojet_ztc4_multiHists.root"
rjrqcpmo = "rjr_multihist_files/KUCMS_QCD_v19_wt2_RjrSkim_v22_ootmet_phomet_ztc4_multiHists.root"


rjrdpjr = "rjr_multihist_files/KUCMS_DEG_v19_wt2_RjrSkim_v22_rawmet_phojet_ztc4_multiHists.root"
rjrdpmr = "rjr_multihist_files/KUCMS_DEG_v19_wt2_RjrSkim_v22_rawmet_phomet_ztc4_multiHists.root"
rjrdpjo = "rjr_multihist_files/KUCMS_DEG_v19_wt2_RjrSkim_v22_ootmet_phojet_ztc4_multiHists.root"
rjrdpmo = "rjr_multihist_files/KUCMS_DEG_v19_wt2_RjrSkim_v22_ootmet_phomet_ztc4_multiHists.root"

boutname = 'llpa_rjr'
title = ' 2017 '
ptnorm = ''
#isonorm = 'hadTowOverEM < 0.02'
isonorm = ''
ptlist = [ '' ]
#ptlist = [ '100','300' ]

for rhname,rhnamea,pt in zip(rhname1,rhname1a,ptlist) :
  ptitle=[ title, ptnorm+pt, isonorm ]
  soutname1 = boutname+pt+'_pm_'
  soutname2 = boutname+pt+'_pj_'

#  ifname = rjrgjetsc02;
#  inhistlist = [[rhnamea,"",ifname,"GJet",norm ],[rhnameb,"",ifname,"GJets2",norm],
#                [rhnamec,"",ifname,"GJetsC",norm],
#                #[rhnamea,"",rjrqcdv1,"QCD",norm],
#                #[rhname,"",rjrqcdv1,"QCD2",norm],[rhname,"",rjrqcdc02,"QCDC",norm]]
#                ]

#  inhistlist = [[rhname,"",rjrgjwt,"GJet",norm ],[rhname,"",rjrqwt,"QCD",norm],
#                [rhname,"",rjrg1wt,"L100",norm],[rhname,"",rjrg3wt,"L300",norm],
#                #[rhname,"",rjrgmsbl100t33p1c02,"L100t33",norm],[rhname,"",rjrgmsbl300t33p1c02,"L300t33",norm]]
#                ]

#  inhistlist = [[rhname,"",rjrgjc2wt,"GJet",norm ],[rhname,"",rjrqc2wt,"QCD",norm],
#                [rhname,"",rjrg1c2wt,"L100",norm],[rhname,"",rjrg3c2wt,"L300",norm],
#                #[rhname,"",rjrgmsbl100t33p1c02,"L100t33",norm],[rhname,"",rjrgmsbl300t33p1c02,"L300t33",norm]]
#                ]

#  inhistlist1 = [[rhname,"",rjrgjpm,"GJet PM",norm ],#[rhname,"",rjrgjpj,"GJets PJ",norm],
#                [rhname,"",rjrqcpm,"QCD PM",norm ],#[rhname,"",rjrqcpj,"QCD PJ",norm],
#                [rhname,"",rjrg1pm,"L100 PM",norm],#[rhname,"",rjrg1pj,"L100 PJ",norm],
#                [rhname,"",rjrg3pm,"L350 PM",norm],#[rhname,"",rjrg3pj,"L350 PJ",norm],
#                [rhname,"",rjrdpm,"DEG PM",norm],
#                ]

#  inhistlist2 = [[rhname,"",rjrgjpj,"GJets PJ",norm],
#                [rhname,"",rjrqcpj,"QCD PJ",norm],
#                [rhname,"",rjrg1pj,"L100 PJ",norm],
#                [rhname,"",rjrg3pj,"L350 PJ",norm],
#                [rhname,"",rjrdpj,"DEG PJ",norm],
#                ]

  inhistlist1a = [[rhnamea,rjrgjpmr,rjrgjpmo,"GJet PM",norm ],#[rhname,"",rjrgjpj,"GJets PJ",norm],
                [rhnamea,rjrqcpmr,rjrqcpmo,"QCD PM",norm ],#[rhname,"",rjrqcpj,"QCD PJ",norm],
                [rhnamea,rjrg1pmr,rjrg1pmo,"L100 PM",norm],#[rhname,"",rjrg1pj,"L100 PJ",norm],
                [rhnamea,rjrg3pmr,rjrg3pmo,"L350 PM",norm],#[rhname,"",rjrg3pj,"L350 PJ",norm],
                [rhnamea,rjrdpmr,rjrdpmo,"DEG PM",norm],
                ]

  inhistlist2a = [[rhnamea,rjrgjpjr,rjrgjpjo,"GJets PJ",norm],
                [rhnamea,rjrqcpjr,rjrqcpjo,"QCD PJ",norm],
                [rhnamea,rjrg1pjr,rjrg1pjo,"L100 PJ",norm],
                [rhnamea,rjrg3pjr,rjrg3pjo,"L350 PJ",norm],
                [rhnamea,rjrdpjr,rjrdpjo,"DEG PJ",norm],
                ]

  date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')
  outname1 = soutname1 + rhname
  outname2 = soutname2 + rhname
  layout = { 'xtitle' : rhname, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  #dostack(inhistlist1, outname1, date, layout, ptitle,  y1, x1, l, t)
  #dostack(inhistlist2, outname2, date, layout, ptitle,  y1, x1, l, t)

  date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')
  outname1 = soutname1 + rhnamea
  outname2 = soutname2 + rhnamea
  layout = { 'xtitle' : rhnamea, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  dostack(inhistlist1a, outname1, date, layout, ptitle,  y, x1, l, t)
  dostack(inhistlist2a, outname2, date, layout, ptitle,  y, x1, l, t)

#  inhistlist1 = [[rhname,"",rjrgjpm,"GJet PM",norml ],#[rhname,"",rjrgjpj,"GJets PJ",norm],
#                [rhname,"",rjrqcpm,"QCD PM",norml ],#[rhname,"",rjrqcpj,"QCD PJ",norm],
#                [rhname,"",rjrg1pm,"L100 PM",norml ],#[rhname,"",rjrg1pj,"L100 PJ",norm],
#                [rhname,"",rjrg3pm,"L350 PM",norml ],#[rhname,"",rjrg3pj,"L350 PJ",norm],
#                [rhname,"",rjrdpm,"DEG PM",norml ],
#                ]

#  inhistlist2 = [[rhname,"",rjrgjpj,"GJets PJ",norml],
#                [rhname,"",rjrqcpj,"QCD PJ",norml],
#                [rhname,"",rjrg1pj,"L100 PJ",norml],
#                [rhname,"",rjrg3pj,"L350 PJ",norml],
#                [rhname,"",rjrdpj,"DEG PJ",norml ],
#                ]

  inhistlist1a = [[rhnamea,rjrgjpmr,rjrgjpmo,"GJet PM",norml ],#[rhname,"",rjrgjpj,"GJets PJ",norm],
                [rhnamea,rjrqcpmr,rjrqcpmo,"QCD PM",norml ],#[rhname,"",rjrqcpj,"QCD PJ",norm],
                [rhnamea,rjrg1pmr,rjrg1pmo,"L100 PM",norml],#[rhname,"",rjrg1pj,"L100 PJ",norm],
                [rhnamea,rjrg3pmr,rjrg3pmo,"L350 PM",norml],#[rhname,"",rjrg3pj,"L350 PJ",norm],
                [rhnamea,rjrdpmr,rjrdpmo,"DEG PM",norml],
                ]

  inhistlist2a = [[rhnamea,rjrgjpjr,rjrgjpjo,"GJets PJ",norml],
                [rhnamea,rjrqcpjr,rjrqcpjo,"QCD PJ",norml],
                [rhnamea,rjrg1pjr,rjrg1pjo,"L100 PJ",norml],
                [rhnamea,rjrg3pjr,rjrg3pjo,"L350 PJ",norml],
                [rhnamea,rjrdpjr,rjrdpjo,"DEG PJ",norml ],
                ]

  date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')
  outname1 = soutname1 + rhname + "_norm"
  outname2 = soutname2 + rhname + "_norm"
  layout = { 'xtitle' : rhname, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  #dostack(inhistlist1, outname1, date, layout, ptitle,  y1l, x1, l, t)
  #dostack(inhistlist2, outname2, date, layout, ptitle,  y1l, x1, l, t)

  date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')
  outname1 = soutname1 + rhnamea + "_norm"
  outname2 = soutname2 + rhnamea + "_norm"
  layout = { 'xtitle' : rhnamea, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  #dostack(inhistlist1a, outname1, date, layout, ptitle,  yn, x1, l, t)
  #dostack(inhistlist2a, outname2, date, layout, ptitle,  yn, x1, l, t)



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


