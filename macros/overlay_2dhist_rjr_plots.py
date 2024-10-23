#  jack w king 3

from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, z, y, x, l, t ):

    first = True
    f1 = []
    h1 = []
    n = 0

    setTDRStyle()
    gROOT.SetBatch(True)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTopMargin(0.09)

    gStyle.SetPadRightMargin(0.15)
    #gStyle.SetTitleYOffset(0.8)
    #gStyle.SetTitleXOffset(0.8)

    #if layout['logx'] : gStyle.SetOptLogx(1)
    #if layout['logy'] : gStyle.SetOptLogy(1)
    c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
    #c1 = TCanvas( 'c1', 'canvas' , 200, 10, 1000, 1000 )
    c1.cd()
    if layout['logx'] : c1.SetLogx()
    if layout['logy'] : c1.SetLogy()
    if layout['logz'] : c1.SetLogz()
    c1.SetGridx(1)
    c1.SetGridy(1)
    c1.Update()
    c1.Draw()

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    for histname, tree, infile, lego  in hist_list :
    
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
        h1.append(f1[n].Get(hist))
        print( '____________________________________________________________________' )
        print( 'Infile Name : ', infile ) 
        print( 'Address : ', f1[n] )
        print( 'Hist Name : ', hist )
        print( 'Address : ', h1[n] )
        print( '____________________________________________________________________' )

        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #h1[n].SetMarkerStyle(6)
        h1[n].SetTitle(layout['title'])
        h1[n].GetXaxis().CenterTitle(True)
        h1[n].GetXaxis().SetTitle(layout['xtitle'])
        h1[n].GetYaxis().CenterTitle(True)
        h1[n].GetYaxis().SetTitle(layout['ytitle'])
        k = [kBlack]
        h1[n].SetLineColor(k[n])
        #if dofit : hfit.SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        #msz = 0.2
        msz = 0.4
        #msz = 0.8
        #msz = 1.2
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        #h1[n].SetMinimum(z[0])
        #h1[n].SetMaximum(z[1])
        h1[n].GetXaxis().SetRangeUser(x[0],x[1])
        h1[n].GetYaxis().SetRangeUser(y[0],y[1])
        if layout['logx'] : h1[n].GetXaxis().SetMoreLogLabels()
        if layout['logy'] : h1[n].GetYaxis().SetMoreLogLabels()
        if( first ) :
                h1[n].Draw('colz')
                #if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                h1[n].Draw('epSAME')
                #if dofit : h1[n].Fit(hfit,'REQ+') 

        c1.Update()
        c1.Draw()
        n += 1

        #End of loop

    c1.cd()
    c1.Update()
    c1.Draw()

    lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    #lat_cms = '#bf{CMS} #it{Work in Progress}' + ptitle[0]
    #lat_title = 'Run2018D 3206730-320824' #   7 fb^{-1} (#sqrt{s} = 13 TeV)'
    #lat_title = 'Run2018D 1Tier miniAOD'
    lat_title = ptitle[1]+' (13 TeV)'
    #lat_form = '#sigma^{2}_{i}=(#frac{N}{A_{eff}/#sigma_{n}})^{2}+2C^{2}'
    lat_form = '#sigma^{2}_{i}=(N/ETeff)^{2}+2C^{2}'
    #lat.SetTextSize(0.045);
    #lat.SetTextFont(132);
    #lat.DrawLatex(0.48+max(0.,0.32-lat.GetXsize()),0.947,lat_string);
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    #lat.DrawLatex(0.12,0.9325,lat_cms);
    lat.DrawLatex(0.16,0.9325,lat_cms);
    #lat.DrawLatex(0.2,0.9325,lat_cms);
    #lat.SetTextSize(0.04);
    #lat.SetTextFont(42);
    #lat.DrawLatex(0.58,0.9325,lat_title);
    lat.DrawLatex((0.685-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);

    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.Show()
    c1.Close()

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

#--------------------------------------------------------------------------------


#islogx = True
islogx = False
#islogy = True
islogy = False
islogz = True
#islogz = False


#for model in ["Bkgrd"] :
for model in ["GluGlu"] :
#for model in ["GluGlu","SqkSqk"] :
#for model in ["GluGlu","SqkSqk","SqkGlu","XinoXino"] :
#for model in ["GluGlu","SqkSqk","SqkGlu","XinoXino","SleptSlept"] :

    #rjrtype = "ootmet_phomet"
    rjrtype = "ootmet_phojet"

    #nsigpho = "genSigPerfect_nSigPho0"
    #nsigphob = "nSigPho0"
    #nsigpho = "genSigPerfect_nSigPho1"
    #nsigphob = "nSigPho1"
    #nsigpho = "genSigPerfect_nSigPho2"
    #nsigphob = "nSigPho2"
    nsigpho = "gsp_nsp1_nvsvdnjabcut"
    nsigphob = "nsp1_nvsvdnjabcut"

    sig100 = "rjr_multihist_files/KUCMS_GMSB_L100_"+model+"_v23_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    sig150 = "rjr_multihist_files/KUCMS_GMSB_L150_"+model+"_v23_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    sig200 = "rjr_multihist_files/KUCMS_GMSB_L200_"+model+"_v23_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    sig250 = "rjr_multihist_files/KUCMS_GMSB_L250_"+model+"_v23_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    #sig300 = "rjr_multihist_files/KUCMS_GMSB_L300_"+model+"_v21_genSigPerfect_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    #sig350 = "rjr_multihist_files/KUCMS_GMSB_L350_"+model+"_v21_genSigPerfect_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    #sig400 = "rjr_multihist_files/KUCMS_GMSB_L400_"+model+"_v21_genSigPerfect_"+nsigpho+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"

    #bgQCD = "rjr_multihist_files/KUCMS_QCD_v19_BG_nSigPho0_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    #bgGJets = "rjr_multihist_files/KUCMS_GJets_v19_BG_nSigPho0_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    bgQCD = "rjr_multihist_files/KUCMS_QCD_BkGrd_v23_"+nsigphob+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"
    bgGJets = "rjr_multihist_files/KUCMS_GJets_BkGrd_v23_"+nsigphob+"_wt2_RjrSkim_v24_"+rjrtype+"_multiHists.root"

    z = [0.001, 10]

    l = [ 0.6,0.65,0.94,0.875 ] # legend position top right short
    t = [0.55,0.85,0.0,0.175,0.225] # titles position

    for rhname in ["ASMvAX2NQSum","ASMvAX2QSum","ASMvVMass","ASMvVDiff","AX2NQSumvVMass","AX2NQSumvVDiff","AX2NQSumvAX2QSum","AX2QSumvVMass","AX2QSumvVDiff","VMassvVDiff","NJetsJavNJetsJb"] :

        x = []
        y = []
        #at = ["rjrVSum","rjrVDiff","ASMass","rjrAX2QSum","rjrAX2NQSum","rjrNJetsJa","rjrNJetsJb"]
        if "ASMvAX2NQSum" in rhname : 
            x = [ 0.0, 12000 ]
            y = [ 0.0, 1.0 ]
            xtitle = "ASMass [GeV]"
            ytitle = "AX2NQSum"
        if "ASMvAX2QSum" in rhname :
            x = [ 0.0, 12000 ]
            y = [ 0.0, 2600 ]
            xtitle = "ASMass [GeV]"
            ytitle = "AX2QSum [GeV]"
        if "ASMvVMass" in rhname :
            x = [ 0.0, 12000 ]
            y = [ 0.0, 2600 ]
            xtitle = "ASMass [GeV]"
            ytitle = "NVSum"
        if "ASMvVDiff" in rhname :
            x = [ 0.0, 12000 ]
            y = [ -1.0, 1.0 ]
            xtitle = "ASMass [GeV]"
            ytitle = "VDiff"
        if "AX2NQSumvVMass" in rhname :
            x = [ 0, 1.0 ]
            y = [ 0.0, 2600 ]
            xtitle = "AX2NQSum"
            ytitle = "NVSum"
        if "AX2NQSumvVDiff" in rhname :
            x = [ 0, 1.0 ]
            y = [ -1.0, 1.0 ]
            xtitle = "AX2NQSum"
            ytitle = "VDiff"
        if "AX2NQSumvAX2QSum" in rhname :
            x = [ 0, 1.0 ]
            y = [ 0, 2600 ]
            xtitle = "AX2NQSum"
            ytitle = "AX2QSum [GeV]"
        if "AX2QSumvVMass" in rhname :
            x = [ 0, 2600 ]
            y = [ 0, 2600 ]
            xtitle = "AX2QSum [GeV]"
            ytitle = "NVSum"
        if "AX2QSumvVDiff" in rhname :
            x = [ 0, 2600 ]
            y = [ -1.0, 1.0 ]
            xtitle = "AX2QSum [GeV]"
            ytitle = "VDiff"
        if "VMassvVDiff" in rhname :
            x = [ 0, 2600 ]
            y = [ -1.0, 1.0 ]
            xtitle = "NVSum"
            ytitle = "VDiff"
        if "NJetsJavNJetsJb" in rhname :
            x = [ 0, 20 ]
            y = [ 0, 20 ]
            xtitle = "NJetsJa"
            ytitle = "NJetsJb"

        leg = ''
        leg100 = 'L100' #"M_{glu},M_{\Chi_{1}} = 838, 139  GeV"
        leg150 = 'L150' #"M_{glu},M_{\Chi_{1}} = 1207, 212  GeV"
        leg200 = 'L200' #"M_{glu},M_{\Chi_{1}} = 1565, 285  GeV"
        leg250 = 'L250' #"M_{glu},M_{\Chi_{1}} = 1915, 358  GeV"
        leg300 = 'L300' #"M_{glu},M_{\Chi_{1}} = 2260, 430  GeV"
        leg350 = 'L350' #"M_{glu},M_{\Chi_{1}} = 2599, 503  GeV"
        leg400 = 'L400' #"M_{glu},M_{\Chi_{1}} = 2935, 576  GeV"
        #leg100 = "M_{sqk, \Chi_{1}} = 1119, 139 GeV"
        #leg200 = "M_{sqk, \Chi_{1}} = 2122, 285 GeV"
        #leg300 = "M_{sqk, \Chi_{1}} = 3084, 430 GeV"
        #leg400 = "M_{sqk, \Chi_{1}} = 4025, 576 GeV"
        #leg100 = "M_{glu, sqk, \Chi_{1}} = 838, 1119, 139 GeV"
        #leg200 = "M_{glu, sqk, \Chi_{1}} = 1565, 2122, 285 GeV"
        #leg300 = "M_{glu, sqk, \Chi_{1}} = 2260, 3084, 430 GeV"
        #leg400 = "M_{glu, sqk, \Chi_{1}} = 2935, 4025, 576 GeV"
        legGJets = 'L400'

        legtitle = ""
        sig = sig200
        leg = leg200
        if "Bkgrd" in model : 
            sig = bgGJets
            leg = 'GJets'

        fhname = rhname
        htitle = fhname
        #xtitle = fhname
        #ytitle = fhname
        outname = 'llpa_rjr_' + model + '_' + leg + '_' + rjrtype + '_' + nsigpho + '_' + fhname
        #outname = 'llpa_rjr_op_' + fhname
        layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'logz' : islogz,'legtitle' : legtitle }
        ptitle=[' 2017', '137 fb', model + ' ' + leg + ' ' + rjrtype ]
        #ptitle=[' 2017','137 fb','']

        inhistlist = [
            [rhname,'',sig,leg],
            #[rhname,'',sig100,leg100],
            #[rhname,'',sig150,leg150],
            #[rhname,'',sig200,leg200],
            #[rhname,'',sig250,leg250],
            ##[rhname,'',sig300,leg300],
        ]

        dostack(inhistlist, outname, date, layout, ptitle, z, y, x, l, t )
        #dostack(hl_other_plots, outname, date, Ic_layout, ptitle, z, y, x, l, t)


