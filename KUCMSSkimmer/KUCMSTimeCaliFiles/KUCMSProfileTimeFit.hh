////////////////////////////////////////////////////////////////////// 
// -*- C++ -*- 
// 
// 
// Original Author:  Jack W King III 
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT 
// 
////////////////////////////////////////////////////////////////////// 
 
#include "KUCMSRootHelperBaseClass.hh" 
#include "KUCMSEcalDetIDFunctions.hh" 
#include "KUCMSTimeCaliStructs.hh" 
 
#include <TRandom.h> 
#include "TChain.h" 
#include "TGraphErrors.h" 
#include "TMultiGraph.h" 
 
#ifndef KUCMSProfileTimeFit 
#define KUCMSProfileTimeFit

class ProfileTimeFit { 

    public: 
     
    ProfileTimeFit(){ profileHist = NULL; form = NULL; fit = NULL; }; 
    ProfileTimeFit( TH1F* hist ) : profileHist(hist) { form = NULL; fit = NULL; }; 
    //~ProfileTimeFit(){ if( profileHist ) delete profileHist; if( form ) delete form; if( fit ) delete fit; };   
 
    // internal data members 
    TH1F* profileHist; 
    TFormula* form; 
    TF1* fit; 
    TF1* fitHigh; 
    TF1* fitLow; 
    TimeFitResult results; 
 
    // helper functions for making fits to variables 
    //bool isEmpty(){ std::cout << " -- ?  is empty " << profileHist << std::endl; return ( profileHist ) ? profileHist->GetEntries() == 0 : true; }; 
    void DoFit(); 
    TimeFitResult GetFitResult(){ return results; }; 
    void deleteHists(){ if( profileHist ) delete profileHist; if( form ) delete form; if( fit ) delete fit; }; 
 
};//<<>>class ProfileTimeFit 

inline void ProfileTimeFit::DoFit(){ 

    std::cout << " -- ProfileTimeFit DoFit " << std::endl; 
    //if( isEmpty() ){ std::cout << " --- No such profile hist !!!!!!! " << std::endl; return; } 
    //std::cout << " -- yes profile exists " << std::endl; 

    float Sqrt2PI   = 2.506628274631001f; 

    // set tmp init vals 
    auto hsum = profileHist->Integral(); 
    if( hsum == 0 ){ 
        std::cout << " -- !!!!!! Intergral is Zero !!!!!!!! " << std::endl; 
        results.mu = 0; 
        results.emu = 0; 
        results.chi2ndf = 0; 
        results.chi2prob = 0; 
        results.sigma = 10; 
        results.sigmaHigh = 0; 
        results.sigmaLow = 0; 
        results.esigma = 10; 
        results.occ = 0; 
        results.rms = 0; 
        results.std = 0; 
        return; 
    }//<<>>if( hsum == 0 ) 

    float mu = profileHist->GetMean(); 
    float sigma = profileHist->GetStdDev(); 
    float rms = profileHist->GetRMS(); 
    //std::cout << " -- getting intgral, mean, std 5" << std::endl; 
    //float norm = hsum/(sigma*Sqrt2PI); 
    float norm = profileHist->GetBinContent(profileHist->GetMaximumBin()); 
    std::cout << " -- Params: sum " << hsum << " mu " << mu << " s " << sigma << " rms " << rms << " norm " << norm << std::endl; 

    std::cout << " -- setting paramters " << std::endl; 
    // range vars   ? names ? 
    //float sigrange = 2.0; 
    float lowerBound = mu - 2*sigma; 
    float upperBound = mu + 2*sigma; 
    float lowerHBound = mu - 4*sigma; 
    float upperHBound = mu + 4*sigma; 
    float lowerLBound = mu - 1*sigma; 
    float upperLBound = mu + 1*sigma; 

    // names for fits and formulas 
    std::string histname = profileHist->GetName(); 
    std::string formname = histname+"_formula"; 
    std::string fitname  = histname+"_fit"; 
    std::string fitnameH  = histname+"_fitHigh"; 
    std::string fitnameL  = histname+"_fitLow"; 

    //if( lowerBound < -1.0 ) lowerBound = -0.5; 
    //if( upperBound > 1.0 ) upperBound = 0.5; 
    //if( abs(mu) > 0.03 ) mu = -0.02; 
    //if( sigma > 1  ) sigma = 0.3; 

    std::cout << " -- setting form and fit " << std::endl; 
    form = new TFormula(formname.c_str(),"[0]*exp(-0.5*((x-[1])/[2])**2)"); 
    fit  = new TF1(fitname.c_str(),form->GetName(),lowerBound,upperBound); 
    fitHigh  = new TF1(fitnameH.c_str(),form->GetName(),lowerHBound,upperHBound); 
    fitLow  = new TF1(fitnameL.c_str(),form->GetName(),lowerLBound,upperLBound); 
 
    fit->SetParName(0,"N");       
    fit->SetParameter(0,norm);  
    fit->SetParLimits(0,norm/10,norm*10); 
    fit->SetParName(1,"#mu");     
    fit->SetParameter(1,mu);  
    //fit->SetParLimits(1,-0.03,0.03); 
    fit->SetParLimits(1,lowerBound,upperBound); 
    fit->SetParName(2,"#sigma");  
    fit->SetParameter(2,sigma);  
    //fit->SetParLimits(2,0,1); 
    fit->SetParLimits(2,sigma/10,sigma*10); 
 
    fitHigh->SetParName(0,"N"); 
    fitHigh->SetParameter(0,norm); 
    fitHigh->SetParLimits(0,norm/10,norm*10); 
    fitHigh->SetParName(1,"#mu"); 
    fitHigh->SetParameter(1,mu); 
    //fitHigh->SetParLimits(1,-0.03,0.03); 
    fitHigh->SetParLimits(1,lowerHBound,upperHBound); 
    fitHigh->SetParName(2,"#sigma"); 
    fitHigh->SetParameter(2,sigma); 
    //fitHigh->SetParLimits(2,0,1); 
    fitHigh->SetParLimits(2,sigma/10,sigma*10); 

    fitLow->SetParName(0,"N"); 
    fitLow->SetParameter(0,norm); 
    fitLow->SetParLimits(0,norm/10,norm*10); 
    fitLow->SetParName(1,"#mu"); 
    fitLow->SetParameter(1,mu); 
    //fitLow->SetParLimits(1,-0.03,0.03); 
    fitLow->SetParLimits(1,lowerLBound,upperLBound); 
    fitLow->SetParName(2,"#sigma"); 
    fitLow->SetParameter(2,sigma); 
    //fitLow->SetParLimits(2,0,1); 
    fitLow->SetParLimits(2,sigma/10,sigma*10); 

    std::cout << " - Fit > Mu: " << mu << " ( " << lowerBound << " to " << upperBound << " )"; 
    std::cout << " Norm: " << norm << "( " << norm/10 << " to " << norm*10 << " )"; 
    std::cout  << " Sigma: " << sigma << "( " << sigma/10 << " to " << sigma*10 << " )" << std::endl; 

    std::cout << " -- fiting " << std::endl; 
    profileHist->Fit(fit->GetName(),"RBQO"); 
    profileHist->Fit(fitHigh->GetName(),"RBQO"); 
    profileHist->Fit(fitLow->GetName(),"RBQO"); 

    std::cout << " -- setting results " << std::endl; 
    results.mu = fit->GetParameter(1); 
    results.emu = fit->GetParError (1); 
    results.chi2ndf = fit->GetChisquare(); 
    results.chi2prob = fit->GetProb(); 
    results.sigma = fit->GetParameter(2); 
    results.sigmaHigh = fitHigh->GetParameter(2); 
    results.sigmaLow = fitLow->GetParameter(2); 
    results.esigma = fit->GetParError (2); 
    results.occ = hsum; 
    results.rms = rms; 
    results.std = sigma; 

    std::cout << " - Fit Sigma: " << results.sigma << " +/- " << results.esigma; 
    std::cout << " ( low " << results.sigmaLow << " high " << results.sigmaHigh << " ) " << std::endl; 

}//<<>>void ProfileTimeFit::DoFit() 

#endif 
