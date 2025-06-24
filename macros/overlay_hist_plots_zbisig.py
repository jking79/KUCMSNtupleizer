#  jack w king 3

from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

#from jwk_tdr_style_py import *
from jwk_tdr_style_v2_py import *

xlabel_list = [ 'Evts', 'PSel', ' ', 'LLL', 'LLM', 'LLT', 'LML', 'LMM', 'LMT', 'LTL', 'LTM', 'LTT', ' ', 'MLL', 'MLM', 'MLT', 'MML', 'MMM', 'MMT', 'MTL', 'MTM', 'MTT', ' ', 'TLL', 'TLM', 'TLT', 'TML', 'TMM', 'TMT', 'TTL', 'TTM', 'TTT' ]

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t ):

    first = True
    f1 = []
    h1 = []
    n = 0

    setTDRStyle()
    gROOT.SetBatch(True)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTopMargin(0.09);
    if layout['logx'] : gStyle.SetOptLogx(1)
    if layout['logy'] : gStyle.SetOptLogy(1)
    c1 = TCanvas( 'c1', 'canvas' )
    #c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
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
    gStyle.SetLegendTextSize(0.025)
    legend.SetBorderSize(0)
    #legend.SetLineColor(kBlack)

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    tree = '' 
    histname = 'cutflow' 

    for infile, lego, indir  in hist_list :

        f1.append(TFile.Open(indir+infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
        h1.append(f1[n].Get(hist))
        ##m = [ 24,25,26,32,27,28,46,24,25,26,32,27,28,46,24,25,26,32,27,28,46,24,25,26,32,27,28,46,24,25,26,32,27,28,46]
        #m = [ 20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34]
        m = [ 24,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34]
        ##m = [ 24,25,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34]
        #m = [ 24,25,26,32,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34,20,21,22,23,47,33,34]
        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(m[n])
        h1[n].SetTitle(layout['title'])
        h1[n].GetXaxis().CenterTitle(True)
        h1[n].GetXaxis().SetTitle(layout['xtitle'])
        h1[n].GetYaxis().CenterTitle(True)
        h1[n].GetYaxis().SetTitle(layout['ytitle'])
        #k = [kMagenta+2,kGreen+2,kOrange+1,kBlue+2,kRed+2,kAzure+4,kViolet+5,kOrange+5,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kMagenta+4,kGreen+4,kYellow+3,kBlue+4,kRed+4,kAzure+6,kBlack,kViolet+3,kOrange+3,kGreen+5,kRed+3,kBlue+3,kGreen+4,kAzure+7,kMagenta+4,kGreen+4]
        k = [kBlack,kMagenta+2,kGreen+2,kOrange+1,kBlue+2,kRed+2,kAzure+4,kViolet+5,kOrange+5,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kMagenta+4,kGreen+4,kYellow+3,kBlue+4,kRed+4,kAzure+6,kBlack,kViolet+3,kOrange+3,kGreen+5,kRed+3,kBlue+3,kGreen+4,kAzure+7,kMagenta+4,kGreen+4]
        h1[n].SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        h1[n].SetMarkerSize(1.0)
        if m[n] == 33 : h1[n].SetMarkerSize(1.2)
        #if n == 0 : h1[n].SetMarkerSize(1.2)
        h1[n].SetMinimum(y[1])
        h1[n].SetMaximum(y[0])
        h1[n].GetXaxis().SetRangeUser(x[0],x[1])
        #if layout['logx'] : h1[n].GetXaxis().SetMoreLogLabels()
        #if layout['logy'] : h1[n].GetYaxis().SetMoreLogLabels()

        h1[n].GetYaxis().SetLabelSize(.025);
        if first :
            cfbin = 1;
            h1[n].GetXaxis().SetLabelSize(.025);
            for label in xlabel_list :
                h1[n].GetXaxis().SetBinLabel(cfbin,label);
                cfbin = cfbin + 1

        if first :
            h1[n].Draw('HIST P')
            first = False
        else :
            h1[n].Draw('HIST P SAME')
        #c1.Update()
        #c1.Draw()

        legend.AddEntry(h1[n],lego,'pl');
        c1.Update()
        c1.Draw()
        n += 1

        #End of loop

    legend.Draw('same')
    #c1.SetGridx(1)
    #c1.SetGridy(1)
    c1.cd()
    c1.Update()
    c1.Draw()

    lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    #lat_cms = '#bf{CMS} #it{Work in Progress}' + ptitle[0]
    #lat_title = 'Run2018D 3206730-320824' #   7 fb^{-1} (#sqrt{s} = 13 TeV)'
    #lat_title = 'Run2018D 1Tier miniAOD'
    lat_title = ptitle[1]+' (13 TeV)'
    #lat_form = '#sigma^{2}_{i}=(#frac{N}{A_{eff}/#sigma_{n}})^{2}+2C^{2}'
    #lat_form = '#sigma^{2}_{i}=(N/ETeff)^{2}+2C^{2}'
    lat_form = ''
    #lat.SetTextSize(0.045);
    #lat.SetTextFont(132);
    #lat.DrawLatex(0.48+max(0.,0.32-lat.GetXsize()),0.947,lat_string);
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    #lat.DrawLatex(0.12,0.9325,lat_cms);
    lat.DrawLatex(0.15,0.9325,lat_cms);
    #lat.DrawLatex(0.2,0.9325,lat_cms);
    #lat.SetTextSize(0.04);
    #lat.SetTextFont(42);
    #lat.DrawLatex(0.58,0.9325,lat_title);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);

    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.Show()
    c1.Close()

cutset = 'cf_2-275-35_2-275-35_0-15-3/'
#cutset = 'cf_1-2-3_2-275-35_1-25-40/'
#cutset = 'cf_2-275-35_2-275-35_1-25-40/'

indir = cutset + 'cf_1pho_1jet/'
pjcut = 'mj1_rp1'
pjleg = 'SelPhoton 1'

indir1 = cutset + 'cf_2pho_1jet/'
pjcut1 = 'mj1_rp2'
pjleg1 = 'SelPhoton 2'

vN = 'v24_'

N1 = '1' ###lep1+ && had1+
lgnd1 = 'lep1+ && had1+'
N2 = '2' ##lep1 && had1
lgnd2 = 'lep1 && had1'
N3 = '3' ##lep0 && had2
lgnd3 = 'lep0 && had2'
N4 = '4' ##lep1+ || had1+
lgnd4 = 'lep1+ || had1+'
N5 = '5' ##lep1+
lgnd5 = 'lep1+'
N6 = '6' ##had1+
lgnd6 = 'had1+'
N7 = '7' ##lep1
lgnd7 = 'Leptonic SV 1'
N8 = '8' ##had1
lgnd8 = 'Hadronic SV 1'
N9 = '9' ##lep2+
lgnd9 = 'Leptonic SV 2+'
N10 = '10' ##had2+
lgnd10 = 'Hadronic SV 2+'

lgnd = lgnd7
N = N7

lgnda = lgnd8
Na = N8

#mpnt = ''
mpnt = 'G2N4N2C3'
#mpnt = 'G2N195N19C1'
mleg = ''
#mleg = '4-2 ct 0.1 '
#mleg = '195-19 ct 0.1 '
#zleg = 'ZJetsToNuNu '
zleg = ''

sfle = 'RjrSkim_v24_cfcatHists.root'

cf_hist_list = [

    #[histname, tree, infile, lego],

    #['KUCMS_Sgg10_llpana_v24_mj1_rp1_mGl-2000_mN2-500_mN1-250_'+sfle,'1Pho1Jet Sig',cutset+'cf_1pho_1jet/'],#
    #['KUCMS_Sgg10_llpana_v24_mj2_rp1_mGl-2000_mN2-500_mN1-250_'+sfle,'1Pho2Jet Sig',cutset+'cf_1pho_2jet/'],#
    #['KUCMS_Sgg10_llpana_v24_mj1_rp2_mGl-2000_mN2-500_mN1-250_'+sfle,'2Pho1Jet Sig',cutset+'cf_2pho_1jet/'],#
    #['KUCMS_Sgg10_llpana_v24_mj2_rp2_mGl-2000_mN2-500_mN1-250_'+sfle,'2Pho2Jet Sig',cutset+'cf_2pho_2jet/'],#
    #['KUCMS_BG_llpana_v24_mj0_rp1_TotBG_'+sfle,'1Pho0Jet BG',cutset+'cf_1pho_0jet/'],
    #['KUCMS_BG_llpana_v24_mj1_rp1_TotBG_'+sfle,'1Pho1Jet BG',cutset+'cf_1pho_1jet/'],
    #['KUCMS_BG_llpana_v24_mj2_rp1_TotBG_'+sfle,'1Pho2Jet BG',cutset+'cf_1pho_2jet/'],
    #['KUCMS_BG_llpana_v24_mj3_rp1_TotBG_'+sfle,'1Pho3Jet BG',cutset+'cf_1pho_3jet/'],
    #['KUCMS_BG_llpana_v24_mj1_rp2_TotBG_'+sfle,'2Pho1Jet BG',cutset+'cf_2pho_1jet/'],
    #['KUCMS_BG_llpana_v24_mj2_rp2_TotBG_'+sfle,'2Pho2Jet BG',cutset+'cf_2pho_2jet/'],
    #['KUCMS_BG_llpana_v24_mj3_rp2_TotBG_'+sfle,'2Pho3Jet BG',cutset+'cf_2pho_3jet/'],

    #['KUCMS_BG_llpana_v24_'+pjcut+'_TotBG_'+sfle,'TotBackground',indir],
    #['KUCMS_BG_llpana_v23_'+pjcut+'_TTTGJets_'+sfle,'TTGJets',indir],
    #['KUCMS_BG_llpana_v23_'+pjcut+'_WZNuTGJets_'+sfle,'WZTGJets',indir],
    #['KUCMS_BG_llpana_v23_'+pjcut+'_ttWZ_'+sfle,'ttWZ',indir],

    #['KUCMS_BG_llpana_v24_'+pjcut+'_GJets_'+sfle,'GJets',indir],
    #['KUCMS_BG_llpana_v24_'+pjcut+'_WJetsLNU_'+sfle,'WJetsLNU',indir],
    #['KUCMS_BG_llpana_v24_'+pjcut+'_ZJetsNuNu_'+sfle,'ZJetsNuNu',indir],
    #['KUCMS_BG_llpana_v24_'+pjcut+'_TTTGJets_'+sfle,'TG/TTGJets',indir],
    #['KUCMS_BG_llpana_v24_'+pjcut+'_ttWZ_'+sfle,'tt/WZ/Jets',indir],


    #['KUCMS_BG_llpana_v24_'+pjcut+'_GJets_'+sfle,'GJets',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_WJetsLNU_'+sfle,'WJetsLNU',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_ZJetsNuNu_'+sfle,'ZJetsNuNu',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_TGJets_'+sfle,'TGJets',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_TTGJets_'+sfle,'TTGJets',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_TTJets_'+sfle,'TTJets',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_WZ_'+sfle,'WZ',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_ttWJets_'+sfle,'ttWJets',indir],
    #['KUCMS_BG_llpana_'+vN+pjcut+'_ttZJets_'+sfle,'ttZJets',indir],

    #['KUCMS_BG_llpana_'+vN+pjcut+'_ZJetsNuNu_'+sfle,pjleg,indir],
    #['KUCMS_BG_llpana_'+vN+pjcut1+'_ZJetsNuNu_'+sfle,pjleg1,indir1],

    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-1500_mN2-500_mN1-100_'+sfle,'1500-500-100',indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-500_mN1-1_'+sfle,'2000, 500, 100',indir],
    #['KUCMS_Sgg10_llpana_'+vN+pjcut+'_mGl-2000_mN2-500_mN1-250_'+sfle,'2000-500-250',indir],#
    #['KUCMS_Sgg10_llpana_'+vN+pjcut+'_mGl-2000_mN2-500_mN1-250_'+sfle,pjleg,indir],#
    #['KUCMS_Sgg10_llpana_'+vN+pjcut1+'_mGl-2000_mN2-500_mN1-250_'+sfle,pjleg1,indir1],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1000_mN1-1_'+sfle,'2000-1000-1',indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1000_mN1-250_'+sfle,'2000-1000-250',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1000_mN1-500_'+sfle,'2000, 1000, 500',indir],
    #['KUCMS_Sgg10_llpana_'+vN+pjcut+'_mGl-2000_mN2-1500_mN1-1000_'+sfle,'2000-1500-1000',indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1500_mN1-1_'+sfle,'2000-1500-1',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1500_mN1-250_'+sfle,'2000-1500-250',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1500_mN1-500_'+sfle,'2000-1500-500',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1900_mN1-1000_'+sfle,'2000, 1900, 1000',indir],
    #['KUCMS_Sgg10_llpana_'+vN+pjcut+'_mGl-2000_mN2-1900_mN1-1500_'+sfle,'2000-1900-1500',indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1900_mN1-1_'+sfle,'2000-1900-1',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1900_mN1-250_'+sfle,'2000-1900-250'+pjleg,indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1900_mN1-500_'+sfle,'2000-1900-500',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1950_mN1-1000_'+sfle,'2000-1950-1000',indir],
    #['KUCMS_Sgg10_llpana_'+vN+pjcut+'_mGl-2000_mN2-1950_mN1-1500_'+sfle,'2000,1950, 1500',indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1950_mN1-1900_'+sfle,'2000-1950-1900',indir],#
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1950_mN1-250_'+sfle,'2000-1950-250',indir],
    #['KUCMS_Sgg10_llpana_v24_'+pjcut+'_mGl-2000_mN2-1950_mN1-500_'+sfle,'2000-1950-500',indir],
    #['KUCMS_Sgg10_llpana_'+vN+pjcut+'_mGl-2500_mN2-1500_mN1-1000_'+sfle,'2500, 1500, 1000',indir],#

    #['zcutflow/KUCMS_ZBG_llpana_v24_lep1_had3_SVZJetsNuNu_'+sfle,'lep1+ && had1+',""],
    #['zcutflow/KUCMS_ZBG_llpana_v24_lep2_had3_SVZJetsNuNu_'+sfle,'lep1 && had1',""],
    #['zcutflow/KUCMS_ZBG_llpana_v24_lep3_had3_SVZJetsNuNu_'+sfle,'lep0 && had2',""],

    ['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had7_SVJetsNuNu_'+sfle,'ZJetsToNuNu',""],
    ['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had7_G2N19N2C3_'+sfle,'2000, 1900, 200',""],
    ['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had7_G2N19N35C3_'+sfle,'2000, 1900, 350',""],
    ['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had7_G2N4N2C3_'+sfle,'2000, 400, 200',""],
    ['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had7_G2N4N35C3_'+sfle,'2000, 400, 350',""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N195N19C1_'+sfle,'195-19 ct 0.1 '+lgnd,""],

    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N19N2C001_'+sfle,'19-2 ct 0.001 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N19N2C1_'+sfle,'19-2 ct 0.1 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N19N2C3_'+sfle,'19-2 ct 0.3 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N19N35C001_'+sfle,'19-35 ct 0.001 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N19N35C1_'+sfle,'19-35 ct 0.1 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N19N35C3_'+sfle,'19-35 ct 0.3 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N4N2C001_'+sfle,'4-2 ct 0.001 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N4N2C1_'+sfle,'4-2 ct 0.1 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N4N2C3_'+sfle,'4-2 ct 0.3 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N4N35C001_'+sfle,'4-35 ct 0.001 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N4N35C1_'+sfle,'4-35 ct 0.1 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N4N35C3_'+sfle,'4-35 ct 0.3 '+lgnd,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N+'_G2N195N19C1_'+sfle,'195-19 ct 0.1 '+lgnd,""],

    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N19N2C001_'+sfle,'19-2 ct 0.001 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N19N2C1_'+sfle,'19-2 ct 0.1 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N19N2C3_'+sfle,'19-2 ct 0.3 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N19N35C001_'+sfle,'19-35 ct 0.001 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N19N35C1_'+sfle,'19-35 ct 0.1 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N19N35C3_'+sfle,'19-35 ct 0.3 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N4N2C001_'+sfle,'4-2 ct 0.001 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N4N2C1_'+sfle,'4-2 ct 0.1 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N4N2C3_'+sfle,'4-2 ct 0.3 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N4N35C001_'+sfle,'4-35 ct 0.001 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N4N35C1_'+sfle,'4-35 ct 0.1 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N4N35C3_'+sfle,'4-35 ct 0.3 '+lgnda,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+Na+'_G2N195N19C1_'+sfle,'195-19 ct 0.1 '+lgnda,""],

    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N1+'_'+mpnt+'_'+sfle,mleg+lgnd1,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N2+'_'+mpnt+'_'+sfle,mleg+lgnd2,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N3+'_'+mpnt+'_'+sfle,mleg+lgnd3,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N4+'_'+mpnt+'_'+sfle,mleg+lgnd4,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N5+'_'+mpnt+'_'+sfle,mleg+lgnd5,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N6+'_'+mpnt+'_'+sfle,mleg+lgnd6,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N7+'_'+mpnt+'_'+sfle,mleg+lgnd7,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N8+'_'+mpnt+'_'+sfle,mleg+lgnd8,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N9+'_'+mpnt+'_'+sfle,mleg+lgnd9,""],
    #['zcutflow/KUCMS_SggZ_llpana_v24_lep3_had'+N10+'_'+mpnt+'_'+sfle,mleg+lgnd10,""],

    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had1_SVJetsNuNu_'+sfle,zleg+lgnd1,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had2_SVJetsNuNu_'+sfle,zleg+lgnd2,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had3_SVJetsNuNu_'+sfle,zleg+lgnd3,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had4_SVJetsNuNu_'+sfle,mleg+lgnd4,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had5_SVJetsNuNu_'+sfle,mleg+lgnd5,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had6_SVJetsNuNu_'+sfle,zleg+lgnd6,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had7_SVJetsNuNu_'+sfle,zleg+lgnd7,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had8_SVJetsNuNu_'+sfle,zleg+lgnd8,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had9_SVJetsNuNu_'+sfle,zleg+lgnd9,""],
    #['zcutflow/KUCMS_SVBG_llpana_v24_lep3_had10_SVJetsNuNu_'+sfle,zleg+lgnd10,""],

]

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

legtitle = ''

rtitle = 'Run2 (13TeV)'

xtitle = ''
ytitle = 'Events'
htitle = ''
#islogx = True
islogx = False
islogy = True
#islogy = False

layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }

#pjcut = 'GlGlZ Sig'
#pjcut = 'ct0.1 lep1+ had1+'
#pjcut = 'Gl 2000 N1 1900 N2 200'
#pjcut = '1 SelPho 1+ Jets/Side ImprvPhoJetIso'
#pjcut = '2 SelPho 1+ Jets/Side'
#pjcut = 'SV ZJetsToNuNu'
#pjcut = 'TotBG & mGl-2000_mN2-500_mN1-250'
#pjcut = 'ZJetsNuNu'
pjcut = 'Leptonic SV 1'
#mcuts = 'Mr:  1000, 2000, 3000'
tcuts =  'Cut  Loose  Medium  Tight '
mcuts =  'Mr   2000    2750       3500  '
rcuts =  'R     0.200   0.275      0.350 ' 
rvcuts = 'Rv   0.00     0.15        0.30  '

cutline = '#splitline{'+pjcut+'}{#splitline{'+tcuts+'}{#splitline{'+mcuts+'}{#splitline{'+rcuts+'}{'+rvcuts+'}}}}'
#ptitle=[' DYJetsToLL','','#splitline{#splitline{ECAL Barrel}{pCalo energy > 1.5 GeV}}{4 Crystal}']
#ptitle=[' DYJetsToLL','','#splitline{ECAL Barrel}{4 Crystal}']
#ptitle=[' ','','#splitline{DiJet dTime Width vs}{Z->ee dTime Width}']
ptitle=[' 400 fb^{-1}','',cutline]
y = [ 40000000000, 0.001 ]
#y = [ 10000000, 0 ]
#y = [ 10000, 0 ]
x = [ 1, 36 ]
l = [ 0.6,0.525,0.925,0.9 ]
#l = [ 0.7,0.725,0.925,0.9 ] # legend position top right
#l = [ 0.2,0.65,0.425,0.9 ] # legend position top left
#l = [ 0.25,0.20,0.52,0.525 ] # legend position bottom left
#l = [ 0.7,0.66,0.9,0.9]
t = [0.2,0.825,0.0,0.175,0.225] # titles position
#t = [0.2,0.825,0.0,0.225,0.6]
#t = [0.225,0.825,0.0]
#t = [0.225,0.775,0.0]

outname = 'tr_cfspec_v1_hist'
dostack(cf_hist_list, outname, date, layout, ptitle,  y, x, l, t)
