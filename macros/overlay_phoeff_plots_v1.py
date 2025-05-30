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
    #dofit = True
    dofit = False
    #sxtal = True
    sxtal = False
    paramn = []
    parnerror = []
    paramc = []
    parcerror = []
    params = []
    parserror = []
    thebinmid = []
    thebinerror = []
    f1s = []
    f1b = []
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

    if dofit :
        ns=str(n)
        hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',100,750,3)
        hfit.SetParName(0,'N')
        hfit.SetParameter(0,40.0)
        #hfit.SetParLimits(0,0,50)
        hfit.SetParLimits(0,0.0,100.0)
        #hfit.SetParameter(0,5)
        #hfit.SetParLimits(0,0,10)
        hfit.SetParName(1,'C')
        hfit.SetParameter(1,0.1)
        hfit.SetParLimits(1,0.0,1.0)
        #hfit.SetParLimits(1,0.01,10.0)
        #hfit.SetParLimits(1,0.02,1.0)
        #hfit.SetParameter(1,0.05)
        #hfit.SetParLimits(1,0.001,1.0)
        hfit.SetParName(2,'S')
        hfit.SetParameter(2,5.0)
        hfit.SetParLimits(2,0.0,25.0)

    mg = TMultiGraph();

    for histname, tree, infiles, infileb, pbin, lego  in hist_list :
    
        f1s.append(TFile.Open(infiles))
        f1b.append(TFile.Open(infileb))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname

        orighists = f1s[n].Get(hist)
        orighistb = f1b[n].Get(hist)
        htitle = 'hist' + str(n)
        lenmybins = n+1
        h1.append(TGraphErrors(lenmybins))

        #for bn in range( 1,lenmybins-1):
        obinvals = float(orighists.GetBinContent(pbin+1))
        obinvalb = float(orighistb.GetBinContent(pbin+1))
        #if obinvals > 0 or obinvalb > 0 :
        print( lego, obinvals, obinvalb )
        h1[n].SetPoint(n+1,obinvals,obinvalb)

        h1[n].UseCurrentStyle()
        msize = n + 25
        if n > 23 : msize = n + 2
        if n > 46 : msize = n - 21 
        h1[n].SetMarkerStyle(msize)
        k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kSpring-7,kSpring+3,kAzure+3,kAzure-7,kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kSpring-7,kSpring+3,kAzure+3,kAzure-7,kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kSpring-7,kSpring+3,kAzure+3]
        #k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        if dofit : hfit.SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        msz = 1.2
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        if( first ) :
                if dofit : h1[n].Fit(hfit,'RE')
                #if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                if dofit : h1[n].Fit(hfit,'RE+')
                #if dofit : h1[n].Fit(hfit,'REQ+') 

        if dofit : 
                 paramn.append(str(abs(hfit.GetParameter(0))))
                 paramc.append(str(abs(hfit.GetParameter(1))))
                 params.append(str(abs(hfit.GetParameter(2))))
                 pne = hfit.GetParError(0)
                 pce = hfit.GetParError(1)
                 pse = hfit.GetParError(2)
                 #print('Fit info',paramn[n],pne,paramc[n],pce)
                 #print('Fit info',paramn[n],pne,paramc[n],pce,params[n],pse)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror.append(str(pne))
                 parcerror.append(str(pce))
                 parserror.append(str(pse))
 
        if not "none" in lego :
                 legend.AddEntry(h1[n],lego,'epl');
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
    lat.DrawLatex(0.15,0.9325,lat_cms);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);
    if dofit : 
        lat.SetTextSize(0.03);
        lat.DrawLatex(t[3],t[4]+.075,lat_form);
        for l in range(0,n):
            lat_param =	'#color['+str(k[l])+']{'
            lat_param = lat_param + 'N : '+paramn[l][0:6]+' #pm '+parnerror[l][0:6]+' [ns] '
            lat_param = lat_param + 'S : '+params[l][0:6]+' #pm '+parserror[l][0:6]+' [ns] '
            lat_param = lat_param + 'C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' [ns]}'
            lat.SetTextSize(0.03);
            lat.DrawLatex(t[3],t[4]-l*.035,lat_param);

    
    if layout['logx'] : c1.SetLogx()
    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root"' )
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
xtitle = 'Signal Photons'
#xtitle = '#Delta_{Run}'
#xtitle = 'GeV'
#xtitle = '[ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#xtitle = 'HcalTowerSumEtBcConeDR04 [GeV]'
ytitle = 'Background Photons'
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

#bhJetht0 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso0_Skim_BaseHists.root"
#bhJetht1 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso1_Skim_BaseHists.root"
#bhJetht2 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso2_Skim_BaseHists.root"
#bhJetht3 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso3_Skim_BaseHists.root"
#bhJetht4 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso4_Skim_BaseHists.root"

#bhGMSBs0 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso0_Skim_BaseHists.root"
#bhGMSBs1 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso1_Skim_BaseHists.root"
#bhGMSBs2 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso2_Skim_BaseHists.root"
#bhGMSBs3 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso3_Skim_BaseHists.root"
#bhGMSBs4 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso4_Skim_BaseHists.root"

#bhGMSBb0 = "hist_root_files/KUCMS_GMSB_L100_Met150_notSig_v15_iso0_Skim_BaseHists.root"
#bhGMSBb1 = "hist_root_files/KUCMS_GMSB_L100_Met150_notSig_v15_iso1_Skim_BaseHists.root"
#bhGMSBb2 = "hist_root_files/KUCMS_GMSB_L100_Met150_notSig_v15_iso2_Skim_BaseHists.root"
#bhGMSBb3 = "hist_root_files/KUCMS_GMSB_L100_Met150_notSig_v15_iso3_Skim_BaseHists.root"
#bhGMSBb4 = "hist_root_files/KUCMS_GMSB_L100_Met150_notSig_v15_iso4_Skim_BaseHists.root"

bhGMSBsp200 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_pt100Default2_Skim_BaseHists.root"
bhGMSBsp100 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_pt200Default2_Skim_BaseHists.root"
bhGMSBsp20 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_pt20Default2_Skim_BaseHists.root"
bhGMSBsp50 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_pt50Default2_Skim_BaseHists.root"
bhGMSBspt = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_ptDefault_Skim_BaseHists.root"

bhJethtsp200 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_pt100Default2_Skim_BaseHists.root"
bhJethtsp100 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_pt200Default2_Skim_BaseHists.root"
bhJethtsp20 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_pt20Default2_Skim_BaseHists.root"
bhJethtsp50 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_pt50Default2_Skim_BaseHists.root"
bhJethtspt = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_ptDefault_Skim_BaseHists.root"

bhGMSBspt3 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_pt20Default3_Skim_BaseHists.root"
bhJethtspt3 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_pt20Default3_Skim_BaseHists.root"

bhJetht0 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_005_360_480_245_Skim_BaseHists.root"
bhJetht1 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_50_100_Skim_BaseHists.root"
bhJetht2 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_50_80_Skim_BaseHists.root"
bhJetht2a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_50_110_Skim_BaseHists.root"
bhJetht3 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_50_90_Skim_BaseHists.root"
bhJetht4 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_60_100_Skim_BaseHists.root"
bhJetht5 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_60_80_Skim_BaseHists.root"
bhJetht5a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_60_110_Skim_BaseHists.root"
bhJetht6 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_60_90_Skim_BaseHists.root"
bhJetht7 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_70_100_Skim_BaseHists.root"
bhJetht8 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_70_80_Skim_BaseHists.root"
bhJetht8a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_70_110_Skim_BaseHists.root"
bhJetht9 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_1_70_90_Skim_BaseHists.root"
bhJetht10 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_50_100_Skim_BaseHists.root"
bhJetht11 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_50_80_Skim_BaseHists.root"
bhJetht11a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_50_110_Skim_BaseHists.root"
bhJetht12 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_50_90_Skim_BaseHists.root"
bhJetht13 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_60_100_Skim_BaseHists.root"
bhJetht14 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_60_80_Skim_BaseHists.root"
bhJetht14a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_60_110_Skim_BaseHists.root"
bhJetht15 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_60_90_Skim_BaseHists.root"
bhJetht16 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_70_100_Skim_BaseHists.root"
bhJetht17 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_70_80_Skim_BaseHists.root"
bhJetht17a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_70_110_Skim_BaseHists.root"
bhJetht18 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_2_70_90_Skim_BaseHists.root"
bhJetht19 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_50_100_Skim_BaseHists.root"
bhJetht20 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_50_80_Skim_BaseHists.root"
bhJetht20a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_50_110_Skim_BaseHists.root"
bhJetht21 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_50_90_Skim_BaseHists.root"
bhJetht22 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_60_100_Skim_BaseHists.root"
bhJetht23 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_60_80_Skim_BaseHists.root"
bhJetht23a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_60_110_Skim_BaseHists.root"
bhJetht24 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_60_90_Skim_BaseHists.root"
bhJetht25 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_70_100_Skim_BaseHists.root"
bhJetht26 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_70_80_Skim_BaseHists.root"
bhJetht26a = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_70_110_Skim_BaseHists.root"
bhJetht27 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_3_70_90_Skim_BaseHists.root"

bhJetht1d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_50_100_Skim_BaseHists.root"
bhJetht2d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_50_110_Skim_BaseHists.root"
bhJetht3d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_50_90_Skim_BaseHists.root"
bhJetht4d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_60_100_Skim_BaseHists.root"
bhJetht5d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_60_110_Skim_BaseHists.root"
bhJetht6d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_60_90_Skim_BaseHists.root"
bhJetht7d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_70_100_Skim_BaseHists.root"
bhJetht8d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_70_110_Skim_BaseHists.root"
bhJetht9d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_1_70_90_Skim_BaseHists.root"
bhJetht10d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_50_100_Skim_BaseHists.root"
bhJetht11d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_50_110_Skim_BaseHists.root"
bhJetht12d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_50_90_Skim_BaseHists.root"
bhJetht13d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_60_100_Skim_BaseHists.root"
bhJetht14d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_60_110_Skim_BaseHists.root"
bhJetht15d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_60_90_Skim_BaseHists.root"
bhJetht16d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_70_100_Skim_BaseHists.root"
bhJetht17d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_70_110_Skim_BaseHists.root"
bhJetht18d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_2_70_90_Skim_BaseHists.root"
bhJetht19d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_50_100_Skim_BaseHists.root"
bhJetht20d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_50_110_Skim_BaseHists.root"
bhJetht21d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_50_90_Skim_BaseHists.root"
bhJetht22d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_60_100_Skim_BaseHists.root"
bhJetht23d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_60_110_Skim_BaseHists.root"
bhJetht24d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_60_90_Skim_BaseHists.root"
bhJetht25d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_70_100_Skim_BaseHists.root"
bhJetht26d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_70_110_Skim_BaseHists.root"
bhJetht27d = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v16_iso_3_70_90_Skim_BaseHists.root"

bhJetht28 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_100_110_Skim_BaseHists.root"
bhJetht29 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_100_120_Skim_BaseHists.root"
bhJetht30 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_100_130_Skim_BaseHists.root"
bhJetht31 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_80_110_Skim_BaseHists.root"
bhJetht32 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_80_120_Skim_BaseHists.root"
bhJetht33 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_80_130_Skim_BaseHists.root"
bhJetht34 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_90_110_Skim_BaseHists.root"
bhJetht35 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_90_120_Skim_BaseHists.root"
bhJetht36 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_4_90_130_Skim_BaseHists.root"
bhJetht37 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_100_110_Skim_BaseHists.root"
bhJetht38 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_100_120_Skim_BaseHists.root"
bhJetht39 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_100_130_Skim_BaseHists.root"
bhJetht40 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_80_110_Skim_BaseHists.root"
bhJetht41 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_80_120_Skim_BaseHists.root"
bhJetht42 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_80_130_Skim_BaseHists.root"
bhJetht43 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_90_110_Skim_BaseHists.root"
bhJetht44 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_90_120_Skim_BaseHists.root"
bhJetht45 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_5_90_130_Skim_BaseHists.root"
bhJetht46 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_100_110_Skim_BaseHists.root"
bhJetht47 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_100_120_Skim_BaseHists.root"
bhJetht48 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_100_130_Skim_BaseHists.root"
bhJetht49 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_80_110_Skim_BaseHists.root"
bhJetht50 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_80_120_Skim_BaseHists.root"
bhJetht51 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_80_130_Skim_BaseHists.root"
bhJetht52 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_90_110_Skim_BaseHists.root"
bhJetht53 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_90_120_Skim_BaseHists.root"
bhJetht54 = "hist_root_files/KUCMS_JetHt_18D_Met150_notSig_v15_iso_6_90_130_Skim_BaseHists.root"

bhGMSBs0 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_005_360_480_245_Skim_BaseHists.root"
bhGMSBs1 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_50_100_Skim_BaseHists.root"
bhGMSBs2 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_50_80_Skim_BaseHists.root"
bhGMSBs3 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_50_90_Skim_BaseHists.root"
bhGMSBs4 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_60_100_Skim_BaseHists.root"
bhGMSBs5 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_60_80_Skim_BaseHists.root"
bhGMSBs6 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_60_90_Skim_BaseHists.root"
bhGMSBs7 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_70_100_Skim_BaseHists.root"
bhGMSBs8 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_70_80_Skim_BaseHists.root"
bhGMSBs9 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_1_70_90_Skim_BaseHists.root"
bhGMSBs10 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_50_100_Skim_BaseHists.root"
bhGMSBs11 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_50_80_Skim_BaseHists.root"
bhGMSBs12 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_50_90_Skim_BaseHists.root"
bhGMSBs13 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_60_100_Skim_BaseHists.root"
bhGMSBs14 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_60_80_Skim_BaseHists.root"
bhGMSBs15 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_60_90_Skim_BaseHists.root"
bhGMSBs16 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_70_100_Skim_BaseHists.root"
bhGMSBs17 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_70_80_Skim_BaseHists.root"
bhGMSBs18 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_2_70_90_Skim_BaseHists.root"
bhGMSBs19 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_50_100_Skim_BaseHists.root"
bhGMSBs20 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_50_80_Skim_BaseHists.root"
bhGMSBs21 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_50_90_Skim_BaseHists.root"
bhGMSBs22 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_60_100_Skim_BaseHists.root"
bhGMSBs23 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_60_80_Skim_BaseHists.root"
bhGMSBs24 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_60_90_Skim_BaseHists.root"
bhGMSBs25 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_70_100_Skim_BaseHists.root"
bhGMSBs26 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_70_80_Skim_BaseHists.root"
bhGMSBs27 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_3_70_90_Skim_BaseHists.root"

bhGMSB4s1 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_50_100_Skim_BaseHists.root"
bhGMSB4s2 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_50_110_Skim_BaseHists.root"
bhGMSB4s3 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_50_90_Skim_BaseHists.root"
bhGMSB4s4 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_60_100_Skim_BaseHists.root"
bhGMSB4s5 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_60_110_Skim_BaseHists.root"
bhGMSB4s6 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_60_90_Skim_BaseHists.root"
bhGMSB4s7 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_70_100_Skim_BaseHists.root"
bhGMSB4s8 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_70_110_Skim_BaseHists.root"
bhGMSB4s9 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_1_70_90_Skim_BaseHists.root"
bhGMSB4s10 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_50_100_Skim_BaseHists.root"
bhGMSB4s11 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_50_110_Skim_BaseHists.root"
bhGMSB4s12 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_50_90_Skim_BaseHists.root"
bhGMSB4s13 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_60_100_Skim_BaseHists.root"
bhGMSB4s14 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_60_110_Skim_BaseHists.root"
bhGMSB4s15 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_60_90_Skim_BaseHists.root"
bhGMSB4s16 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_70_100_Skim_BaseHists.root"
bhGMSB4s17 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_70_110_Skim_BaseHists.root"
bhGMSB4s18 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_2_70_90_Skim_BaseHists.root"
bhGMSB4s19 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_50_100_Skim_BaseHists.root"
bhGMSB4s20 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_50_110_Skim_BaseHists.root"
bhGMSB4s21 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_50_90_Skim_BaseHists.root"
bhGMSB4s22 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_60_100_Skim_BaseHists.root"
bhGMSB4s23 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_60_110_Skim_BaseHists.root"
bhGMSB4s24 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_60_90_Skim_BaseHists.root"
bhGMSB4s25 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_70_100_Skim_BaseHists.root"
bhGMSB4s26 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_70_110_Skim_BaseHists.root"
bhGMSB4s27 = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v15_iso_3_70_90_Skim_BaseHists.root"

bhGMSB4s1d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_50_100_Skim_BaseHists.root"
bhGMSB4s2d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_50_110_Skim_BaseHists.root"
bhGMSB4s3d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_50_90_Skim_BaseHists.root"
bhGMSB4s4d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_60_100_Skim_BaseHists.root"
bhGMSB4s5d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_60_110_Skim_BaseHists.root"
bhGMSB4s6d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_60_90_Skim_BaseHists.root"
bhGMSB4s7d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_70_100_Skim_BaseHists.root"
bhGMSB4s8d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_70_110_Skim_BaseHists.root"
bhGMSB4s9d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_1_70_90_Skim_BaseHists.root"
bhGMSB4s10d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_50_100_Skim_BaseHists.root"
bhGMSB4s11d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_50_110_Skim_BaseHists.root"
bhGMSB4s12d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_50_90_Skim_BaseHists.root"
bhGMSB4s13d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_60_100_Skim_BaseHists.root"
bhGMSB4s14d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_60_110_Skim_BaseHists.root"
bhGMSB4s15d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_60_90_Skim_BaseHists.root"
bhGMSB4s16d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_70_100_Skim_BaseHists.root"
bhGMSB4s17d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_70_110_Skim_BaseHists.root"
bhGMSB4s18d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_2_70_90_Skim_BaseHists.root"
bhGMSB4s19d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_50_100_Skim_BaseHists.root"
bhGMSB4s20d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_50_110_Skim_BaseHists.root"
bhGMSB4s21d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_50_90_Skim_BaseHists.root"
bhGMSB4s22d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_60_100_Skim_BaseHists.root"
bhGMSB4s23d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_60_110_Skim_BaseHists.root"
bhGMSB4s24d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_60_90_Skim_BaseHists.root"
bhGMSB4s25d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_70_100_Skim_BaseHists.root"
bhGMSB4s26d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_70_110_Skim_BaseHists.root"
bhGMSB4s27d = "hist_root_files/KUCMS_GMSB_L400_Met150_Signal_v16_iso_3_70_90_Skim_BaseHists.root"

bhGMSBs28 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_100_110_Skim_BaseHists.root"
bhGMSBs29 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_100_120_Skim_BaseHists.root"
bhGMSBs30 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_100_130_Skim_BaseHists.root"
bhGMSBs31 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_80_110_Skim_BaseHists.root"
bhGMSBs32 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_80_120_Skim_BaseHists.root"
bhGMSBs33 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_80_130_Skim_BaseHists.root"
bhGMSBs34 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_90_110_Skim_BaseHists.root"
bhGMSBs35 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_90_120_Skim_BaseHists.root"
bhGMSBs36 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_4_90_130_Skim_BaseHists.root"
bhGMSBs37 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_100_110_Skim_BaseHists.root"
bhGMSBs38 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_100_120_Skim_BaseHists.root"
bhGMSBs39 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_100_130_Skim_BaseHists.root"
bhGMSBs40 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_80_110_Skim_BaseHists.root"
bhGMSBs41 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_80_120_Skim_BaseHists.root"
bhGMSBs42 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_80_130_Skim_BaseHists.root"
bhGMSBs43 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_90_110_Skim_BaseHists.root"
bhGMSBs44 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_90_120_Skim_BaseHists.root"
bhGMSBs45 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_5_90_130_Skim_BaseHists.root"
bhGMSBs46 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_100_110_Skim_BaseHists.root"
bhGMSBs47 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_100_120_Skim_BaseHists.root"
bhGMSBs48 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_100_130_Skim_BaseHists.root"
bhGMSBs49 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_80_110_Skim_BaseHists.root"
bhGMSBs50 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_80_120_Skim_BaseHists.root"
bhGMSBs51 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_80_130_Skim_BaseHists.root"
bhGMSBs52 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_90_110_Skim_BaseHists.root"
bhGMSBs53 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_90_120_Skim_BaseHists.root"
bhGMSBs54 = "hist_root_files/KUCMS_GMSB_L100_Met150_Signal_v15_iso_6_90_130_Skim_BaseHists.root"

#y = [ 0.0001, 0.02 ]
y = [ 0.015, 0.065 ]
x = [ 0.65, 1.0 ]
#x = [ 0.25, 1.0 ]
#x = [ 0.5, 1.0 ]
#x = [ 0.001, 1.0 ]
#x = [ 0.9, 1.0 ]
#l = [ 0.7,0.7,0.925,0.9 ] # legend position top right
#l = [ 0.2,0.55,0.3,0.875 ] # legend position top left
l = [ 0.2,0.15,0.3,0.875 ] # legend position top left
t = [0.4,0.825,0.0,0.175,0.225] # titles position
ptlist= [ ['20',1], ['30',2], ['100',3] ]

outnames = 'llpa_phocuteff_pt'
rhname = 'phoCEffPhoton'
cut1 = "< 3.6"
cut2 = "< 7.0"
cut3 = "< 6.0"
cut4 = "< 5.0"

for pt, ptn in ptlist :

  #cl3 = ', hadTowOverEM < 0.02, ecalRHSumEtConeDR04 < 9.0'
  cl3 = ''
  cl1 = 'L400-JetHT pt > '+pt+' GeV'+cl3
  #cl2 = 'trkSumPtSolidConeDR04 Cuts'
  #cl2 = 'ecalRHSumEtConeDR04 Cuts'
  #cl2 = 'hadTowOverEM Cuts'
  #cl2 = 'Default 0.13x0.011-0.10x0.004'
  cl2 = ''
  cutline = '#splitline{'+cl1+'}{'+cl2+'}'
  #cutline = cl2+" : "+cl1
  ptitle=[' Photon Cut Efficencies ','',cutline]
  layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  
  inhistlist = [
  
      #[ rhname, "", bhGMSBs0, bhGMSBb0, ptn, "GMSB"],
      #[ rhname, "", bhGMSBs1, bhGMSBb1, ptn, cut1],
      #[ rhname, "", bhGMSBs2, bhGMSBb2, ptn, cut2],
      #[ rhname, "", bhGMSBs3, bhGMSBb3, ptn, cut3],
      #[ rhname, "", bhGMSBs4, bhGMSBb4, ptn, cut4],
      #[ rhname, "", bhGMSBs0, bhJetht0, ptn, "Default"],

      [ rhname, "", bhGMSB4s1, bhJetht1, ptn, "0.01-5.0-10.0"],
      #[ rhname, "", bhGMSB4s1d, bhJetht1d, ptn, "0.01-5.0H-10.0"],
      [ rhname, "", bhGMSB4s2, bhJetht2a, ptn, "0.01-5.0-11.0"],
      [ rhname, "", bhGMSB4s3, bhJetht3, ptn, "0.01-5.0-9.0"],
      [ rhname, "", bhGMSB4s4, bhJetht4, ptn, "0.01-6.0-10.0"],
      #[ rhname, "", bhGMSB4s4d, bhJetht4d, ptn, "0.01-6.0H-10.0"],
      [ rhname, "", bhGMSB4s5, bhJetht5a, ptn, "0.01-6.0-11.0"],
      [ rhname, "", bhGMSB4s6, bhJetht6, ptn, "0.01-6.0-9.0"],
      [ rhname, "", bhGMSB4s7, bhJetht7, ptn, "0.01-7.0-10.0"],
      #[ rhname, "", bhGMSB4s7d, bhJetht7d, ptn, "0.01-7.0H-10.0"],
      [ rhname, "", bhGMSB4s8, bhJetht8a, ptn, "0.01-7.0-11.0"],
      [ rhname, "", bhGMSB4s9, bhJetht9, ptn, "0.01-7.0-9.0"],
      [ rhname, "", bhGMSB4s10, bhJetht10, ptn, "0.02-5.0-10.0"],
      #[ rhname, "", bhGMSB4s10d, bhJetht10d, ptn, "0.02-5.0H-10.0"],
      [ rhname, "", bhGMSB4s11, bhJetht11a, ptn, "0.02-5.0-11.0"],
      [ rhname, "", bhGMSB4s12, bhJetht12, ptn, "0.02-5.0-9.0"],
      [ rhname, "", bhGMSB4s13, bhJetht13, ptn, "0.02-6.0-10.0"],
      #[ rhname, "", bhGMSB4s13d, bhJetht13d, ptn, "0.02-6.0H-10.0"],
      [ rhname, "", bhGMSB4s14, bhJetht14a, ptn, "0.02-6.0-11.0"],
      [ rhname, "", bhGMSB4s15, bhJetht15, ptn, "0.02-6.0-9.0"],
      [ rhname, "", bhGMSB4s16, bhJetht16, ptn, "0.02-7.0-10.0"],
      #[ rhname, "", bhGMSB4s16d, bhJetht16d, ptn, "0.02-7.0H-10.0"],
      [ rhname, "", bhGMSB4s17, bhJetht17a, ptn, "0.02-7.0-11.0"],
      [ rhname, "", bhGMSB4s18, bhJetht18, ptn, "0.02-7.0-9.0"],
      [ rhname, "", bhGMSB4s19, bhJetht19, ptn, "0.03-5.0-10.0"],
      #[ rhname, "", bhGMSB4s19d, bhJetht19d, ptn, "0.03-5.0H-10.0"],
      [ rhname, "", bhGMSB4s20, bhJetht20a, ptn, "0.03-5.0-11.0"],
      [ rhname, "", bhGMSB4s21, bhJetht21, ptn, "0.03-5.0-9.0"],
      [ rhname, "", bhGMSB4s22, bhJetht22, ptn, "0.03-6.0-10.0"],
      #[ rhname, "", bhGMSB4s22d, bhJetht22d, ptn, "0.03-6.0H-10.0"],
      [ rhname, "", bhGMSB4s23, bhJetht23a, ptn, "0.03-6.0-11.0"],
      [ rhname, "", bhGMSB4s24, bhJetht24, ptn, "0.03-6.0-9.0"],
      [ rhname, "", bhGMSB4s25, bhJetht25, ptn, "0.03-7.0-10.0"],
      #[ rhname, "", bhGMSB4s25d, bhJetht25d, ptn, "0.03-7.0H-10.0"],
      [ rhname, "", bhGMSB4s26, bhJetht26a, ptn, "0.03-7.0-11.0"],
      [ rhname, "", bhGMSB4s27, bhJetht27, ptn, "0.03-7.0-9.0"],
#      [ rhname, "", bhGMSBs28, bhJetht28, ptn, "0.04-8.0-11.0"],
#      [ rhname, "", bhGMSBs29, bhJetht29, ptn, "0.04-8.0-12.0"],
#      [ rhname, "", bhGMSBs30, bhJetht30, ptn, "0.04-8.0-13.0"],
#      [ rhname, "", bhGMSBs31, bhJetht31, ptn, "0.04-9.0-11.0"],
#      [ rhname, "", bhGMSBs32, bhJetht32, ptn, "0.04-9.0-12.0"],
#      [ rhname, "", bhGMSBs33, bhJetht33, ptn, "0.04-9.0-13.0"],
#      [ rhname, "", bhGMSBs34, bhJetht34, ptn, "0.04-10.0-11.0"],
#      [ rhname, "", bhGMSBs35, bhJetht35, ptn, "0.04-10.0-12.0"],
#      [ rhname, "", bhGMSBs36, bhJetht36, ptn, "0.04-10.0-13.0"],
#      [ rhname, "", bhGMSBs37, bhJetht37, ptn, "0.05-8.0-11.0"],
#      [ rhname, "", bhGMSBs38, bhJetht38, ptn, "0.05-8.0-12.0"],
#      [ rhname, "", bhGMSBs39, bhJetht39, ptn, "0.05-8.0-13.0"],
#      [ rhname, "", bhGMSBs40, bhJetht40, ptn, "0.05-9.0-11.0"],
#      [ rhname, "", bhGMSBs41, bhJetht41, ptn, "0.05-9.0-12.0"],
#      [ rhname, "", bhGMSBs42, bhJetht42, ptn, "0.05-9.0-13.0"],
#      [ rhname, "", bhGMSBs43, bhJetht43, ptn, "0.05-10.0-11.0"],
#      [ rhname, "", bhGMSBs44, bhJetht44, ptn, "0.05-10.0-12.0"],
#      [ rhname, "", bhGMSBs45, bhJetht45, ptn, "0.05-10.0-13.0"],
#      [ rhname, "", bhGMSBs46, bhJetht46, ptn, "0.06-8.0-11.0"],
#      [ rhname, "", bhGMSBs47, bhJetht47, ptn, "0.06-8.0-12.0"],
#      [ rhname, "", bhGMSBs48, bhJetht48, ptn, "0.06-8.0-13.0"],
#      [ rhname, "", bhGMSBs49, bhJetht49, ptn, "0.06-9.0-11.0"],
#      [ rhname, "", bhGMSBs50, bhJetht50, ptn, "0.06-9.0-12.0"],
#      [ rhname, "", bhGMSBs51, bhJetht51, ptn, "0.06-9.0-13.0"],
#      [ rhname, "", bhGMSBs52, bhJetht52, ptn, "0.06-10.0-11.0"],
#      [ rhname, "", bhGMSBs53, bhJetht53, ptn, "0.06-10.0-12.0"],
#      [ rhname, "", bhGMSBs54, bhJetht54, ptn, "0.06-10.0-13.0"],

#      [ rhname, "", bhGMSBsp200, bhJethtsp200, ptn, "Defaut pt = 200"],
#      [ rhname, "", bhGMSBsp100, bhJethtsp100, ptn, "Defaut pt = 100"],
#      [ rhname, "", bhGMSBsp50, bhJethtsp50, ptn, "Defaut pt = 50"],
#      [ rhname, "", bhGMSBsp20, bhJethtsp20, ptn, "Defaut pt = 20"],
#      [ rhname, "", bhGMSBspt, bhJethtspt, ptn, "Defaut pt based"],

  ]

  outname = outnames+pt
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


