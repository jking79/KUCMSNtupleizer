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

            #if bkgpass == 0 : bkgpass = sigpass*0.001
            non = sigpass+bkgpass
            muoff = bkgpass
            #non = 140
            #muoff = 83.33
            #noff = 100
            #print( "input:",sigpass,bkgpass,non,muoff)
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
    mg.SetTitle(layout['title']) #+';X axis '+layout['xtitle']+';Y Axis '+layout['ytitle']+';')
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

l = [ 0.7,0.6,0.925,0.90 ] # legend position top right
#l = [ 0.2,0.60,0.4,0.90 ] # legend position top left
#l = [ 0.25,0.20,0.52,0.525 ] # legend position bottom left
#t = [0.2,0.48,0.23,0.175,0.22] # titles position
t = [0.45,0.825,0.23,0.175,0.22] # titles position

#x = [ 0.0, 1.05 ]
#x = [ 0.0, 0.3 ]
#x = [ 0.8, 1.05 ]
x = [ 0.0, 8000 ]
#y = [ 0.001, 10.0 ]
#y = [ -0.25, 0.1 ]
#y = [ -0.25, 0.5 ]
#y = [ -0.25, 0.75 ]
#y = [ -1.0, 1.0 ]
#y = [ -1.0, 2.0 ]
#y = [ -1.0, 3.0 ]
#y = [ -1.0, 5.0 ]
#y = [ -1.0, 10.0 ]
y = [ -1.0, 15.0 ]
#y = [ -1.0, 25.0 ]
#x = [ 0.0, 0.8 ]
#y = [ 0.001, 10.0 ]
#y = [ 0.0001, 100.0 ]
#x = [ 0.01, 0.4 ]
#y = [ 0.2, 0.8 ]

#title = ' GGG Tot BG'
#title = ' GGG Tot BG Scale 0.05 M_{g} = 2650'
#title = ' GGG Tot BG Scale 0.1 M_{g} = 2500'
#title = ' GGG Tot BG Scale 0.2 M_{g} = 2345'
title = ' QQG Tot BG'
#isocut = '#splitline{#splitline{nSigPho > 0, phoTime > -0.6}{R_{s} > 0.3}}{}'
#isocut = '#splitline{#splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1}}{}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{R_{s} > 0.34}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{R_{s} > 0.40}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{R_{s} > 0.50}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 3400}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 2200}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 1400}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 750}'
isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 3400, R_{s} > 0.34}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 2200, R_{s} > 0.40}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 1400, R_{s} > 0.50}'
#isocut = '#splitline{ #splitline{nSigPho > 0, phoTime > -0.6}{R_{xa} < 1 && R_{xb} < 1} }{M_{s} > 750, R_{xb} = 1 w/ R_{xa} > 0.7}'
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

scale = 'sc10'
cutset = 'cv13'

Fbg000 = 'KUCMS_BG_llpana_v37_sc10_cv0_TotBG_RjrSkim_Hists.root'
Fbg001 = 'KUCMS_BG_llpana_v37_sc10_cv1_TotBG_RjrSkim_Hists.root'
Fbg002 = 'KUCMS_BG_llpana_v37_sc10_cv2_TotBG_RjrSkim_Hists.root'
Fbg003 = 'KUCMS_BG_llpana_v37_sc10_cv3_TotBG_RjrSkim_Hists.root'
Fbg004 = 'KUCMS_BG_llpana_v37_sc10_cv4_TotBG_RjrSkim_Hists.root'
Fbg005 = 'KUCMS_BG_llpana_v37_sc10_cv5_TotBG_RjrSkim_Hists.root'
Fbg006 = 'KUCMS_BG_llpana_v37_sc10_cv6_TotBG_RjrSkim_Hists.root'
Fbg007 = 'KUCMS_BG_llpana_v37_sc10_cv7_TotBG_RjrSkim_Hists.root'
Fbg008 = 'KUCMS_BG_llpana_v37_sc10_cv8_TotBG_RjrSkim_Hists.root'

Fbg = 'KUCMS_BG_llpana_v37_sc10_'+cutset+'_TotBG_RjrSkim_Hists.root'

Fsig1 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG15dN5dN1dC1_RjrSkim_Hists.root'
Fsig2 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG25dN15dN1kC1_RjrSkim_Hists.root'

Fsig3 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN15dN1C1_RjrSkim_Hists.root'
Fsig304 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN15dN1C1_RjrSkim_Hists.root'
Fsig4 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN15dN1kC1_RjrSkim_Hists.root'
Fsig403 = 'KUCMS_llpana_v37_'+scale+'_cv3_PG2kN15dN1kC1_RjrSkim_Hists.root'
Fsig5 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN15dN250C1_RjrSkim_Hists.root'
Fsig504 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN15dN250C1_RjrSkim_Hists.root'
Fsig6 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN15dN5dC1_RjrSkim_Hists.root'
Fsig604 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN15dN5dC1_RjrSkim_Hists.root'

Fsig7 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1950N15dC1_RjrSkim_Hists.root'
Fsig702 = 'KUCMS_llpana_v37_'+scale+'_cv2_PG2kN1950N15dC1_RjrSkim_Hists.root'
Fsig8 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1950N19dC1_RjrSkim_Hists.root'
Fsig802 = 'KUCMS_llpana_v37_'+scale+'_cv2_PG2kN1950N19dC1_RjrSkim_Hists.root'
Fsig9 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1950N1C1_RjrSkim_Hists.root'
Fsig904 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1950N1C1_RjrSkim_Hists.root'
Fsig10 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1950N1kC1_RjrSkim_Hists.root'
Fsig1003 = 'KUCMS_llpana_v37_'+scale+'_cv3_PG2kN1950N1kC1_RjrSkim_Hists.root'
Fsig11 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1950N250C1_RjrSkim_Hists.root'
Fsig1104 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1950N250C1_RjrSkim_Hists.root'
Fsig12 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1950N5dC1_RjrSkim_Hists.root'
Fsig1204 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1950N5dC1_RjrSkim_Hists.root'

Fsig13 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN19dN15dC1_RjrSkim_Hists.root'
Fsig1302 = 'KUCMS_llpana_v37_'+scale+'_cv2_PG2kN19dN15dC1_RjrSkim_Hists.root'
Fsig14 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN19dN1C1_RjrSkim_Hists.root'
Fsig1404 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN19dN1C1_RjrSkim_Hists.root'
Fsig15 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN19dN1kC1_RjrSkim_Hists.root'
Fsig1503 = 'KUCMS_llpana_v37_'+scale+'_cv3_PG2kN19dN1kC1_RjrSkim_Hists.root'
Fsig16 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN19dN250C1_RjrSkim_Hists.root'
Fsig1604 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN19dN250C1_RjrSkim_Hists.root'
Fsig17 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN19dN5dC1_RjrSkim_Hists.root'
Fsig1704 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN19dN5dC1_RjrSkim_Hists.root'

Fsig18 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1kN1C1_RjrSkim_Hists.root'
Fsig1804 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1kN1C1_RjrSkim_Hists.root'
Fsig19 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1kN250C1_RjrSkim_Hists.root'
Fsig1904 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1kN250C1_RjrSkim_Hists.root'
Fsig20 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN1kN5dC1_RjrSkim_Hists.root'
Fsig2004 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1kN5dC1_RjrSkim_Hists.root'

Fsig21 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN5dN1C1_RjrSkim_Hists.root'
Fsig2104 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN5dN1C1_RjrSkim_Hists.root'
Fsig22 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_PG2kN5dN250C1_RjrSkim_Hists.root'
Fsig2204 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN5dN250C1_RjrSkim_Hists.root'

# ------- sqsqG -----------------
Fqsig1 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG17dN15dN1C1_RjrSkim_Hists.root'
Fqsig2 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG17dN3dN1C1_RjrSkim_Hists.root'

Fqsig3 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG1850N1650N1C1_RjrSkim_Hists.root'
Fqsig4 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG1850N3dN1C1_RjrSkim_Hists.root'

Fqsig5 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG2kN18dN1C1_RjrSkim_Hists.root'
Fqsig6 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG2kN3dN1C1_RjrSkim_Hists.root'

Fqsig7 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG2150N1950N1C1_RjrSkim_Hists.root'
Fqsig8 = 'KUCMS_llpana_v37_'+scale+'_'+cutset+'_QG2150N3dN1C1_RjrSkim_Hists.root'


G1bg = 'KUCMS_BG_llpana_v37_sc10_cv4_TotBG_RjrSkim_Hists.root'
G2bg = 'KUCMS_BG_llpana_v37_sc10_cv3_TotBG_RjrSkim_Hists.root'
G3bg = 'KUCMS_BG_llpana_v37_sc10_cv2_TotBG_RjrSkim_Hists.root'
# Mn1 < 1.0 TeV
G1sig21 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN15dN1C1_RjrSkim_Hists.root'
G1sig19 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN1kN250C1_RjrSkim_Hists.root'
G1sig6 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN15dN5dC1_RjrSkim_Hists.root'
G1sig17 = 'KUCMS_llpana_v37_'+scale+'_cv4_PG2kN19dN5dC1_RjrSkim_Hists.root'
# Mn1 1.0 - 1.5
G2sig4 = 'KUCMS_llpana_v37_'+scale+'_cv3_PG2kN15dN1kC1_RjrSkim_Hists.root'
G2sig15 = 'KUCMS_llpana_v37_'+scale+'_cv3_PG2kN19dN1kC1_RjrSkim_Hists.root'
G2sig10 = 'KUCMS_llpana_v37_'+scale+'_cv3_PG2kN1950N1kC1_RjrSkim_Hists.root'
# Mn1 > 1.5
G3sig13 = 'KUCMS_llpana_v37_'+scale+'_cv2_PG2kN19dN15dC1_RjrSkim_Hists.root'
G3sig7 = 'KUCMS_llpana_v37_'+scale+'_cv2_PG2kN1950N15dC1_RjrSkim_Hists.root'
G3sig8 = 'KUCMS_llpana_v37_'+scale+'_cv2_PG2kN1950N19dC1_RjrSkim_Hists.root'


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
rhname7 = 'sumPhoPtsH21'
rhname8 = 'Rxa'

lint = [ '_' ]
ncomp = 'roc'

rha = rhname7
rhnamea = rha
xtitle = rhnamea

#indir = 'rjr_roc_hist_files/'
#indir = 'rjr_91725_zbi_hists/'
#indir = 'rjr_91825_zbi_hists/'
#indir = 'rjr_91925_zbi_hists/'
#indir = 'rjr_92125_zbi_hists/'
indir = 'rjr_92325_zbi_hists/'

outname = 'llpa_met150_multizbi'
layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }

for thing in lint :

    outname = outname+ thing + ncomp

    inhistlist = [
  
        ##[ rhname1, "", base1-sig, base2-bkgrd, rfgmsbroc1-sig, rfgmsbroc2-bkgrd, rhname1, cut1, usebase1 ],

        # very comp
        #[ rhnamea, "", indir+Fsig13, indir+Fbg, sig13, cut1 ], 
        #[ rhnamea, "", indir+Fsig7, indir+Fbg, sig7, cut1 ],
        #[ rhnamea, "", indir+Fsig8, indir+Fbg, sig8, cut1 ],

        # 1000
        #[ rhnamea, "", indir+Fsig4, indir+Fbg, sig4, cut1 ],
        #[ rhnamea, "", indir+Fsig15, indir+Fbg, sig15, cut1 ],
        #[ rhnamea, "", indir+Fsig10, indir+Fbg, sig10, cut1 ],
        ##[ rhnamea, "", indir+Fsig2, indir+Fbg, sig2, cut1 ],

        # 500
        #[ rhnamea, "", indir+Fsig20, indir+Fbg, sig20, cut1 ],
        #[ rhnamea, "", indir+Fsig6, indir+Fbg, sig6, cut1 ],
        #[ rhnamea, "", indir+Fsig17, indir+Fbg, sig17, cut1 ],
        #[ rhnamea, "", indir+Fsig12, indir+Fbg, sig12, cut1 ],

        # 250
        #[ rhnamea, "", indir+Fsig22, indir+Fbg, sig22, cut1 ],
        #[ rhnamea, "", indir+Fsig19, indir+Fbg, sig19, cut1 ],
        #[ rhnamea, "", indir+Fsig5, indir+Fbg, sig5, cut1 ],  
        #[ rhnamea, "", indir+Fsig16, indir+Fbg, sig16, cut1 ],
        #[ rhnamea, "", indir+Fsig11, indir+Fbg, sig11, cut1 ],

        # 001
        #[ rhnamea, "", indir+Fsig21, indir+Fbg, sig21, cut1 ],
        #[ rhnamea, "", indir+Fsig18, indir+Fbg, sig18, cut1 ],
        #[ rhnamea, "", indir+Fsig3, indir+Fbg, sig3, cut1 ],
        #[ rhnamea, "", indir+Fsig14, indir+Fbg, sig14, cut1 ],
        #[ rhnamea, "", indir+Fsig9, indir+Fbg, sig9, cut1 ],

        # glue mass spread
        #[ rhnamea, "", indir+Fsig19, indir+Fbg, sig19, cut1 ],
        #[ rhnamea, "", indir+Fsig5, indir+Fbg, sig5, cut1 ],
        #[ rhnamea, "", indir+Fsig16, indir+Fbg, sig16, cut1 ],
        #[ rhnamea, "", indir+Fsig2, indir+Fbg, sig2, cut1 ],

        # comp
        #[ rhnamea, "", indir+Fsig13, indir+Fbg, sig13, cut1 ],
        #[ rhnamea, "", indir+Fsig8, indir+Fbg, sig8, cut1 ],

        #  Mn1 < 1.0 TeV
        #[ rhnamea, "", indir+G1sig21, indir+G1bg, sig21, cut1 ],
        #[ rhnamea, "", indir+G1sig19, indir+G1bg, sig19, cut1 ],
        #[ rhnamea, "", indir+G1sig6, indir+G1bg, sig6, cut1 ],
        #[ rhnamea, "", indir+G1sig17, indir+G1bg, sig17, cut1 ],
        #[ rhnamea, "", indir+Fsig21, indir+Fbg, sig21, cut1 ],
        #[ rhnamea, "", indir+Fsig19, indir+Fbg, sig19, cut1 ],
        #[ rhnamea, "", indir+Fsig6, indir+Fbg, sig6, cut1 ],
        #[ rhnamea, "", indir+Fsig17, indir+Fbg, sig17, cut1 ],

        #  Mn1 1.0 - 1.5 TeV
        #[ rhnamea, "", indir+G2sig4, indir+G2bg, sig4, cut1 ],
        #[ rhnamea, "", indir+G2sig15, indir+G2bg, sig15, cut1 ],
        #[ rhnamea, "", indir+G2sig10, indir+G2bg, sig10, cut1 ],
        #[ rhnamea, "", indir+Fsig4, indir+Fbg, sig4, cut1 ],
        #[ rhnamea, "", indir+Fsig15, indir+Fbg, sig15, cut1 ],
        #[ rhnamea, "", indir+Fsig10, indir+Fbg, sig10, cut1 ],

        # Mn1 > 1.5
        #[ rhnamea, "", indir+G3sig13, indir+G3bg, sig13, cut1 ],
        #[ rhnamea, "", indir+G3sig7, indir+G3bg, sig7, cut1 ],
        #[ rhnamea, "", indir+G3sig8, indir+G3bg, sig8, cut1 ],
        #[ rhnamea, "", indir+Fsig13, indir+Fbg, sig13, cut1 ], 
        #[ rhnamea, "", indir+Fsig7, indir+Fbg, sig7, cut1 ],
        #[ rhnamea, "", indir+Fsig8, indir+Fbg, sig8, cut1 ],

        # sqsq-g
        #[ rhnamea, "", indir+Fqsig1, indir+Fbg, sigz1, cut1 ],
        #[ rhnamea, "", indir+Fqsig3, indir+Fbg, sigz3, cut1 ],
        #[ rhnamea, "", indir+Fqsig5, indir+Fbg, sigz5, cut1 ],
        #[ rhnamea, "", indir+Fqsig7, indir+Fbg, sigz7, cut1 ],

        [ rhnamea, "", indir+Fqsig2, indir+Fbg, sigz2, cut1 ],
        [ rhnamea, "", indir+Fqsig4, indir+Fbg, sigz4, cut1 ],
        [ rhnamea, "", indir+Fqsig6, indir+Fbg, sigz6, cut1 ],
        [ rhnamea, "", indir+Fqsig8, indir+Fbg, sigz8, cut1 ],

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


