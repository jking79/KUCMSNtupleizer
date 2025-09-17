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
    f1b = []
    f2b = []
    h1 = []
    n = 0

    setTDRStyle()
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

    mg = TMultiGraph()
    first = True

    for histname, tree, infile1s, infile2b, lego, cutside in hist_list :
    
        print( 'using : ', infile1s, ' for ', histname )
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname

        f1.append(TFile.Open(infile1s))
        f2.append(TFile.Open(infile2b))
        orighist1 = f1[n].Get(hist)
        orighist2 = f2[n].Get(hist)
        htitle = 'hist' + str(n)
        lenmybins = int(orighist1.GetNbinsX())
        h1.append(TGraphErrors(lenmybins))
        totbg = orighist2.Integral()
        for bn in range( 2, lenmybins):
            sigpass = 0
            bkgfail = 0
            rint1 = orighist1.Integral(bn,lenmybins)
            lint1 = orighist1.Integral(1,bn-1)
            rint2 = orighist2.Integral(bn,lenmybins)
            lint2 = orighist2.Integral(1,bn-1)
            cutv = orighist1.GetBinCenter(bn)
            #print( "`Intergrals : ", lint1, lint2 )
            if( cutside == "<" ):
                sigpass = lint1
                bkgpass = lint2
            else :
                sigpass = rint1
                bkgpass = rint2

            non = sigpass+bkgpass
            muoff = bkgpass
            #non = 140
            #muoff = 83.33
            #noff = 100
            if muoff != 0 :
                sigma = 0.1
                #tau = muoff/(muoff*sigma*muoff*sigma)
                #####tau = m/(m*s*m*s) = m/(mmss) = 1/mss 
                tau = 1/(muoff*sigma*sigma)
                #tau = 1./Nbkg/(deltaNbkg*deltaNbkg); 
                #tau = (1/muoff)/(sigma*sigma)
                #tau = 1.2
                invtau = 1.0/(1.0+tau)
                noff = tau*muoff
                pbi = TMath.BetaIncomplete(invtau,non,noff+1)
                zbi = sqrt(2)*TMath.ErfcInverse(2*pbi)
                #zbi = sqrt(2)*TMath.ErfInverse( 1 - 2*pbi )
                #print( " set point : ", bn, non, noff, tau, pbi, zbi)
                #if zbi != 0 :
                h1[n].SetPoint(bn,cutv,zbi)
                print( round(cutv,2), round(zbi,2), round(sigpass,2), round(bkgpass,2) )
        #ms = [26,25,26,25,26,25,26,25,26,25]
        ms = [24,25,26,30,28,24,25,26,30,28,24,25,26,30]
        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(ms[n])
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kMagenta+2,kMagenta+2,kGreen+2,kGreen+2,kBlue+2,kBlue+2,kRed+2,kRed+2,kAzure+4,kAzure+4]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        msz = 0.8
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        legend.AddEntry(h1[n],lego+' '+cutside,'epl');
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
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    lat.DrawLatex(0.16,0.96,lat_cms);
    lat.DrawLatex((0.828-t[2]),0.93,lat_title);
    lat.SetTextSize(0.03);
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

rtitle = 'Run2UL AODSIM'
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
#ytitle = 'a.u.'
ytitle = 'Zbi'
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

#l = [ 0.7,0.7,0.925,0.9 ] # legend position top right
l = [ 0.2,0.60,0.4,0.90 ] # legend position top left
#l = [ 0.25,0.20,0.52,0.525 ] # legend position bottom left
#t = [0.2,0.48,0.23,0.175,0.22] # titles position
t = [0.45,0.825,0.23,0.175,0.22] # titles position

#x = [ 0.0, 1.05 ]
x = [ 0.0, 0.3 ]
#x = [ 0.8, 1.05 ]
#x = [ 0.0, 8000 ]
#y = [ 0.001, 10.0 ]
#y = [ -0.3, 0.25 ]
#y = [ -0.3, 0.5 ]
#y = [ -0.3, 0.75 ]
#y = [ -0.3, 1.0 ]
y = [ -0.3, 2.0 ]
#y = [ -1.0, 5.0 ]
#y = [ -1.0, 15.0 ]
#y = [ -1.0, 25.0 ]
#x = [ 0.0, 0.8 ]
#y = [ 0.001, 10.0 ]
#y = [ 0.0001, 100.0 ]
#x = [ 0.01, 0.4 ]
#y = [ 0.2, 0.8 ]

#title = ' GGG 20-15-5'
title = ' GGG Tot BG'
#title = ' QQG Tot BG'
#isocut = '#splitline{#splitline{nSigPho > 0, phoTime > -0.6}{}}{}'
#isocut = '#splitline{#splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1, Min(R_{xa},R_{xb}) > 0.3}}{}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1, Min(R_{xa},R_{xb}) > 0.3} }{M_{s} > 3400}'
isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1, Min(R_{xa},R_{xb}) > 0.3} }{M_{s} > 2200}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1, Min(R_{xa},R_{xb}) > 0.3} }{M_{s} > 750}'
#isocut = '#splitline{#splitline{M_{s} > 1250}{R_{xa} && R_{xb} != 1}}{nSigPho > 0}'
#isocut = '#splitline{#splitline{R_{s} > 0.3}{R_{xa} && R_{xb} != 1}}{nSigPho > 0}'
#isocut = '#splitline{#splitline{M_{s} > 1250, R_{s} > 0.3}{R_{xa} && R_{xb} != 1}}{nSigPho > 0}'
#isocut = '#splitline{#splitline{M_{s} > 2000}{R_{xa} && R_{xb} != 1}}{nSigPho > 0, Rs > 0.2}'
ptitle=[ title, "", isocut ]
bkgrd1 = 'GJets '
bkgrd2 = 'WJets '
bkgrd3 = 'ZJets '

cut1 = '>'
cut2 = '<'


v37tbg000 = 'KUCMS_BG_llpana_v37_tm06_cv0_TotBG_RjrSkim_Hists.root'
v37tbg001 = 'KUCMS_BG_llpana_v37_tm06_cv1_TotBG_RjrSkim_Hists.root'
v37tbg002 = 'KUCMS_BG_llpana_v37_tm06_cv2_TotBG_RjrSkim_Hists.root'
v37tbg003 = 'KUCMS_BG_llpana_v37_tm06_cv3_TotBG_RjrSkim_Hists.root'
v37tbg004 = 'KUCMS_BG_llpana_v37_tm06_cv4_TotBG_RjrSkim_Hists.root'
v37tbg005 = 'KUCMS_BG_llpana_v37_tm06_cv5_TotBG_RjrSkim_Hists.root'

v37tsig101 = 'KUCMS_llpana_v37_tm06_cv1_PG15dN5dN1dC1_RjrSkim_Hists.root'
v37tsig201 = 'KUCMS_llpana_v37_tm06_cv1_PG25dN15dN1kC1_RjrSkim_Hists.root'

v37tsig301 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN15dN1C1_RjrSkim_Hists.root'
v37tsig401 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN15dN1kC1_RjrSkim_Hists.root'
v37tsig501 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN15dN250C1_RjrSkim_Hists.root'
v37tsig601 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN15dN5dC1_RjrSkim_Hists.root'

v37tsig701 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1950N15dC1_RjrSkim_Hists.root'
v37tsig801 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1950N19dC1_RjrSkim_Hists.root'
v37tsig901 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1950N1C1_RjrSkim_Hists.root'
v37tsig1001 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1950N1kC1_RjrSkim_Hists.root'
v37tsig1101 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1950N250C1_RjrSkim_Hists.root'
v37tsig1201 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1950N5dC1_RjrSkim_Hists.root'

v37tsig1301 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN19dN15dC1_RjrSkim_Hists.root'
v37tsig1401 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN19dN1C1_RjrSkim_Hists.root'
v37tsig1501 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN19dN1kC1_RjrSkim_Hists.root'
v37tsig1601 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN19dN250C1_RjrSkim_Hists.root'
v37tsig1701 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN19dN5dC1_RjrSkim_Hists.root'

v37tsig1801 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1kN1C1_RjrSkim_Hists.root'
v37tsig1901 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1kN250C1_RjrSkim_Hists.root'
v37tsig2001 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN1kN5dC1_RjrSkim_Hists.root'

v37tsig2101 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN5dN1C1_RjrSkim_Hists.root'
v37tsig2201 = 'KUCMS_llpana_v37_tm06_cv1_PG2kN5dN250C1_RjrSkim_Hists.root'

# ------- sqsqG -----------------
v37tqsig100 = 'KUCMS_llpana_v37_tm06_cv0_QG17dN15dN1C1_RjrSkim_Hists.root'
v37tqsig200 = 'KUCMS_llpana_v37_tm06_cv0_QG17dN3dN1C1_RjrSkim_Hists.root'

v37tqsig300 = 'KUCMS_llpana_v37_tm06_cv0_QG1850N1650N1C1_RjrSkim_Hists.root'
v37tqsig400 = 'KUCMS_llpana_v37_tm06_cv0_QG1850N3dN1C1_RjrSkim_Hists.root'

v37tqsig500 = 'KUCMS_llpana_v37_tm06_cv0_QG2kN18dN1C1_RjrSkim_Hists.root'
v37tqsig600 = 'KUCMS_llpana_v37_tm06_cv0_QG2kN3dN1C1_RjrSkim_Hists.root'

v37tqsig700 = 'KUCMS_llpana_v37_tm06_cv0_QG2150N1950N1C1_RjrSkim_Hists.root'
v37tqsig800 = 'KUCMS_llpana_v37_tm06_cv0_QG2150N3dN1C1_RjrSkim_Hists.root'

v37tqsig101 = 'KUCMS_llpana_v37_tm06_cv1_QG17dN15dN1C1_RjrSkim_Hists.root'
v37tqsig201 = 'KUCMS_llpana_v37_tm06_cv1_QG17dN3dN1C1_RjrSkim_Hists.root'

v37tqsig301 = 'KUCMS_llpana_v37_tm06_cv1_QG1850N1650N1C1_RjrSkim_Hists.root'
v37tqsig401 = 'KUCMS_llpana_v37_tm06_cv1_QG1850N3dN1C1_RjrSkim_Hists.root'

v37tqsig501 = 'KUCMS_llpana_v37_tm06_cv1_QG2kN18dN1C1_RjrSkim_Hists.root'
v37tqsig601 = 'KUCMS_llpana_v37_tm06_cv1_QG2kN3dN1C1_RjrSkim_Hists.root'

v37tqsig701 = 'KUCMS_llpana_v37_tm06_cv1_QG2150N1950N1C1_RjrSkim_Hists.root'
v37tqsig801 = 'KUCMS_llpana_v37_tm06_cv1_QG2150N3dN1C1_RjrSkim_Hists.root'

#///////////////////////
v37tsig1904 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN1kN250C1_RjrSkim_Hists.root'
v37tsig504 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN15dN250C1_RjrSkim_Hists.root'
v37tsig1604 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN19dN250C1_RjrSkim_Hists.root'
v37tsig204 = 'KUCMS_llpana_v37_tm06_cv4_PG25dN15dN1kC1_RjrSkim_Hists.root'

v37tsig1302 = 'KUCMS_llpana_v37_tm06_cv2_PG2kN19dN15dC1_RjrSkim_Hists.root'
v37tsig802 = 'KUCMS_llpana_v37_tm06_cv2_PG2kN1950N19dC1_RjrSkim_Hists.root'

v37tsig1305 = 'KUCMS_llpana_v37_tm06_cv5_PG2kN19dN15dC1_RjrSkim_Hists.root'
v37tsig805 = 'KUCMS_llpana_v37_tm06_cv5_PG2kN1950N19dC1_RjrSkim_Hists.root'
#////////////////////////

v37tbg004 = 'KUCMS_BG_llpana_v37_tm06_cv4_TotBG_RjrSkim_Hists.root'
# Mn1 < 1.0 TeV
v37tsig2104 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN15dN1C1_RjrSkim_Hists.root'
v37tsig1904 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN1kN250C1_RjrSkim_Hists.root'
v37tsig604 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN15dN5dC1_RjrSkim_Hists.root'
v37tsig1704 = 'KUCMS_llpana_v37_tm06_cv4_PG2kN19dN5dC1_RjrSkim_Hists.root'
# Mn1 1.0 - 1.5
v37tsig403 = 'KUCMS_llpana_v37_tm06_cv3_PG2kN15dN1kC1_RjrSkim_Hists.root'
v37tsig1503 = 'KUCMS_llpana_v37_tm06_cv3_PG2kN19dN1kC1_RjrSkim_Hists.root'
v37tsig1003 = 'KUCMS_llpana_v37_tm06_cv3_PG2kN1950N1kC1_RjrSkim_Hists.root'
# Mn1 > 1.5
v37tsig1302 = 'KUCMS_llpana_v37_tm06_cv2_PG2kN19dN15dC1_RjrSkim_Hists.root'
v37tsig702 = 'KUCMS_llpana_v37_tm06_cv2_PG2kN1950N15dC1_RjrSkim_Hists.root'
v37tsig802 = 'KUCMS_llpana_v37_tm06_cv2_PG2kN1950N19dC1_RjrSkim_Hists.root'


bg = 'TotBG' 
sig1 = '1.5-0.5-0.001 ' 
sig2 = '2.5-1.5-1.0 ' #

sig3 = '2.0-1.5-0.001 ' 
sig4 = '2.0-1.5-1.0 ' #
sig5 = '2.0-1.5-0.25 ' # 
sig6 = '2.0-1.5-0.5 '

sig7 = '2.0-1.95-1.5 '
sig8 = '2.0-1.95-1.9 ' #
sig9 = '2.0-1.95-0.001 '
sig10 = '2.0-1.95-1.0 '
sig11 = '2.0-1.95-0.25 '
sig12 = '2.0-1.95-0.5 '

sig13 = '2.0-1.9-1.5 ' #
sig14 = '2.0-1.9-0.001 '
sig15 = '2.0-1.9-1.0 '
sig16 = '2.0-1.9-0.25 ' #
sig17 = '2.0-1.9-0.5 '

sig18 = '2.0-1.0-0.001 '
sig19 = '2.0-1.0-0.25 '   #
sig20 = '2.0-1.0-0.5 '

sig21 = '2.0-0.5-0.001 '
sig22 = '2.0-0.5-0.25 '

sigz1 = '1.7-1.5-0.1 '
sigz2 = '1.7-0.3-0.1 '
sigz3 = '1.85-1.65-0.1 '
sigz4 = '1.85-0.3-0.1 '
sigz5 = '2.0-1.8-0.1 '
sigz6 = '2.0-0.3-0.1 '
sigz7 = '2.15-1.95-0.1 '
sigz8 = '2.15-0.3-0.1 '

rhname1 = 'Rs' 
rhname2 = 'R' 
rhname3 = 'Ms'  
rhname4 = 'Mr'   
rhname5 = 'Rx'  
rhname6 = 'Rm'  
rhname7 = 'phoPtsH41'

lint = [ '_' ]
ncomp = 'roc'

rha = rhname7
rhnamea = rha
xtitle = rhnamea

#indir = 'rjr_roc_hist_files/'
indir = 'rjr_zbi_hists/'


outname = 'llpa_met150_multizbi'
layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }

for thing in lint :

    outname = outname+ thing + ncomp

    inhistlist = [
  
        ##[ rhname1, "", base1-sig, base2-bkgrd, rfgmsbroc1-sig, rfgmsbroc2-bkgrd, rhname1, cut1, usebase1 ],
        #[ rhnamea, "", indir+v37tsig1301, indir+v37tbg001, sig13, cut1 ], 
        #[ rhnamea, "", indir+v37tsig701, indir+v37tbg001, sig7, cut1 ],
        #[ rhnamea, "", indir+v37tsig801, indir+v37tbg001, sig8, cut1 ],

        #[ rhnamea, "", indir+v37tsig401, indir+v37tbg001, sig4, cut1 ],
        #[ rhnamea, "", indir+v37tsig1501, indir+v37tbg001, sig15, cut1 ],
        #[ rhnamea, "", indir+v37tsig1001, indir+v37tbg001, sig10, cut1 ],
        ##[ rhnamea, "", indir+v37tsig201, indir+v37tbg001, sig2, cut1 ],

        #[ rhnamea, "", indir+v37tsig2001, indir+v37tbg001, sig20, cut1 ],
        #[ rhnamea, "", indir+v37tsig601, indir+v37tbg001, sig6, cut1 ],
        #[ rhnamea, "", indir+v37tsig1701, indir+v37tbg001, sig17, cut1 ],
        #[ rhnamea, "", indir+v37tsig1201, indir+v37tbg001, sig12, cut1 ],

        #[ rhnamea, "", indir+v37tsig2201, indir+v37tbg001, sig22, cut1 ],
        #[ rhnamea, "", indir+v37tsig1901, indir+v37tbg001, sig19, cut1 ],
        #[ rhnamea, "", indir+v37tsig501, indir+v37tbg001, sig5, cut1 ],  
        #[ rhnamea, "", indir+v37tsig1601, indir+v37tbg001, sig16, cut1 ],
        #[ rhnamea, "", indir+v37tsig1101, indir+v37tbg001, sig11, cut1 ],

        ##[ rhnamea, "", indir+v37tsig101, indir+v37tbg001, sig1, cut1 ],
        #[ rhnamea, "", indir+v37tsig2101, indir+v37tbg001, sig21, cut1 ],
        #[ rhnamea, "", indir+v37tsig1801, indir+v37tbg001, sig18, cut1 ],
        #[ rhnamea, "", indir+v37tsig301, indir+v37tbg001, sig3, cut1 ],
        #[ rhnamea, "", indir+v37tsig1401, indir+v37tbg001, sig14, cut1 ],
        #[ rhnamea, "", indir+v37tsig901, indir+v37tbg001, sig9, cut1 ],

        #[ rhnamea, "", indir+v37tsig1904, indir+v37tbg004, sig19, cut1 ],
        #[ rhnamea, "", indir+v37tsig504, indir+v37tbg004, sig5, cut1 ],
        #[ rhnamea, "", indir+v37tsig1604, indir+v37tbg004, sig16, cut1 ],
        #[ rhnamea, "", indir+v37tsig204, indir+v37tbg004, sig2, cut1 ],

        #[ rhnamea, "", indir+v37tsig1302, indir+v37tbg002, sig13, cut1 ],
        #[ rhnamea, "", indir+v37tsig802, indir+v37tbg002, sig8, cut1 ],

        #  Mn1 > 1.0 TeV
        #[ rhnamea, "", indir+v37tsig2104, indir+v37tbg004, sig21, cut1 ],
        #[ rhnamea, "", indir+v37tsig1904, indir+v37tbg004, sig19, cut1 ],
        #[ rhnamea, "", indir+v37tsig604, indir+v37tbg004, sig6, cut1 ],
        #[ rhnamea, "", indir+v37tsig1704, indir+v37tbg004, sig17, cut1 ],
        #  Mn1 1.0 - 1.5 TeV
        [ rhnamea, "", indir+v37tsig403, indir+v37tbg003, sig4, cut1 ],
        [ rhnamea, "", indir+v37tsig1503, indir+v37tbg003, sig15, cut1 ],
        [ rhnamea, "", indir+v37tsig1003, indir+v37tbg003, sig10, cut1 ],
        # Mn1 > 1.5
        #[ rhnamea, "", indir+v37tsig1301, indir+v37tbg001, sig13, cut1 ],
        #[ rhnamea, "", indir+v37tsig701, indir+v37tbg001, sig7, cut1 ],
        #[ rhnamea, "", indir+v37tsig801, indir+v37tbg001, sig8, cut1 ],

        # sqsq-g
        #[ rhnamea, "", indir+v37tqsig101, indir+v37tbg001, sigz1, cut1 ],
        #[ rhnamea, "", indir+v37tqsig201, indir+v37tbg001, sigz2, cut1 ],
        #[ rhnamea, "", indir+v37tqsig301, indir+v37tbg001, sigz3, cut1 ],
        #[ rhnamea, "", indir+v37tqsig401, indir+v37tbg001, sigz4, cut1 ],
        #[ rhnamea, "", indir+v37tqsig501, indir+v37tbg001, sigz5, cut1 ],
        #[ rhnamea, "", indir+v37tqsig601, indir+v37tbg001, sigz6, cut1 ],
        #[ rhnamea, "", indir+v37tqsig701, indir+v37tbg001, sigz7, cut1 ],
        #[ rhnamea, "", indir+v37tqsig801, indir+v37tbg001, sigz8, cut1 ],

    ]
  
    dostack(inhistlist, outname, date, layout, ptitle,  y, x, l, t)
  
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


