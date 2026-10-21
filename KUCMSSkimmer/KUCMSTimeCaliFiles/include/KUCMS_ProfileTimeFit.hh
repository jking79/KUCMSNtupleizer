//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMS_RootHelperBaseClass.hh"
#include "KUCMS_EcalDetIDFunctions.hh"
#include "KUCMS_TimeCaliStructs.hh"

#include <TRandom.h>
#include "TChain.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <algorithm>
#include <cmath>

#ifndef KUCMS_ProfileTimeFit
#define KUCMS_ProfileTimeFit

// Fits timing profiles and retains the ROOT fit objects needed by the
// resolution workflow.  Ownership is kept with the helper so callers can
// reuse a profile while producing diagnostic plots and fit parameters.

class kucms_ProfileTimeFit {

  public:

  kucms_ProfileTimeFit(){ profileHist = NULL; form = NULL; fit = NULL; fitHigh = NULL; fitLow = NULL; };
  kucms_ProfileTimeFit( TH1F* hist ) : profileHist(hist) { form = NULL; fit = NULL; fitHigh = NULL; fitLow = NULL; };
  //~ProfileTimeFit(){ if( profileHist ) delete profileHist; if( form ) delete form; if( fit ) delete fit; };

  // internal data members
  TH1F* profileHist;
  TFormula* form;
  TF1* fit;
  TF1* fitHigh;
  TF1* fitLow;
  kucms_TimeFitResult results;

  // helper functions for making fits to variables
  //bool isEmpty(){ std::cout << " -- ?  is empty " << profileHist << std::endl; return ( profileHist ) ? profileHist->GetEntries() == 0 : true; };
  //void DoFit();
  kucms_TimeFitResult GetFitResult(){ return results; };
  void deleteHists(){
    if( fitLow ) delete fitLow;
    if( fitHigh ) delete fitHigh;
    if( fit ) delete fit;
    if( form ) delete form;
    if( profileHist ) delete profileHist;
    profileHist = NULL;
    form = NULL;
    fit = NULL;
    fitHigh = NULL;
    fitLow = NULL;
  };

  //inline void kucms_ProfileTimeFit::DoFit(){
  bool DoFit(){

    std::cout << " -- ProfileTimeFit DoFit " << std::endl;
    //if( isEmpty() ){ std::cout << " --- No such profile hist !!!!!!! " << std::endl; return; }
    //std::cout << " -- yes profile exists " << std::endl;

    // set tmp init vals
    if( not profileHist ) return false;
    auto hsum = profileHist->Integral();
    if( hsum <= 0 ){
      std::cout << " -- !!!!!! Intergral is Zero !!!!!!!! " << std::endl;
      return false;
    }//<<>>if( hsum == 0 )

    float mu = profileHist->GetMean();
    float sigma = profileHist->GetStdDev();
    float rms = profileHist->GetRMS();
    //std::cout << " -- getting intgral, mean, std 5" << std::endl;
    //float norm = hsum/(sigma*Sqrt2PI);
    float norm = profileHist->GetBinContent(profileHist->GetMaximumBin());
    std::cout << " -- Params: sum " << hsum << " mu " << mu << " s " << sigma << " rms " << rms << " norm " << norm << std::endl;
    results.occ = hsum;
    results.rms = rms;
    results.std = sigma;
    if( not std::isfinite(mu) || not std::isfinite(sigma) || sigma <= 0 ||
    not std::isfinite(norm) || norm <= 0 ){
      std::cout << " -- Invalid profile moments; skipping Gaussian fits." << std::endl;
      return false;
    }//<<>>if( profile cannot define a Gaussian fit )

    std::cout << " -- setting paramters " << std::endl;
    // range vars   ? names ?
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
    // N keeps the temporary functions out of the histogram ownership list;
    // deleteHists() owns and deletes all three fit objects explicitly.
    TFitResultPtr centralResult = profileHist->Fit(fit,"RBQONS");
    TFitResultPtr highResult = profileHist->Fit(fitHigh,"RBQONS");
    TFitResultPtr lowResult = profileHist->Fit(fitLow,"RBQONS");

    std::cout << " -- setting results " << std::endl;
    results.mu = fit->GetParameter(1);
    results.emu = fit->GetParError (1);
    results.chi2 = fit->GetChisquare();
    results.ndf = fit->GetNDF();
    results.chi2ndf = results.ndf > 0 ? results.chi2/results.ndf : 0;
    results.chi2prob = fit->GetProb();
    results.sigma = fit->GetParameter(2);
    results.sigmaHigh = fitHigh->GetParameter(2);
    results.sigmaLow = fitLow->GetParameter(2);
    results.esigma = fit->GetParError (2);
    results.fitStatus = static_cast<int>(centralResult);
    results.fitStatusHigh = static_cast<int>(highResult);
    results.fitStatusLow = static_cast<int>(lowResult);
    results.covStatus = centralResult.Get() ? centralResult->CovMatrixStatus() : -1;
    results.covStatusHigh = highResult.Get() ? highResult->CovMatrixStatus() : -1;
    results.covStatusLow = lowResult.Get() ? lowResult->CovMatrixStatus() : -1;

    auto parameterAtLimit = []( TF1* function, int parameter ){
      double lower = 0;
      double upper = 0;
      function->GetParLimits(parameter,lower,upper);
      if( not std::isfinite(lower) || not std::isfinite(upper) || upper <= lower ) return false;
      const double value = function->GetParameter(parameter);
      const double tolerance = 1e-6*std::max(1.0,std::abs(upper-lower));
      return value <= lower+tolerance || value >= upper-tolerance;
    };
    results.sigmaAtLimit = parameterAtLimit(fit,2);
    results.sigmaHighAtLimit = parameterAtLimit(fitHigh,2);
    results.sigmaLowAtLimit = parameterAtLimit(fitLow,2);

    results.valid = results.fitStatus == 0 && results.covStatus >= 2 && results.ndf > 0 &&
    std::isfinite(results.mu) && std::isfinite(results.emu) &&
    std::isfinite(results.sigma) && std::isfinite(results.esigma) &&
    results.sigma > 0 && results.esigma > 0 && not results.sigmaAtLimit;
    results.validHigh = results.fitStatusHigh == 0 && results.covStatusHigh >= 2 && fitHigh->GetNDF() > 0 &&
    std::isfinite(results.sigmaHigh) && std::isfinite(fitHigh->GetParError(2)) &&
    results.sigmaHigh > 0 && fitHigh->GetParError(2) > 0 && not results.sigmaHighAtLimit;
    results.validLow = results.fitStatusLow == 0 && results.covStatusLow >= 2 && fitLow->GetNDF() > 0 &&
    std::isfinite(results.sigmaLow) && std::isfinite(fitLow->GetParError(2)) &&
    results.sigmaLow > 0 && fitLow->GetParError(2) > 0 && not results.sigmaLowAtLimit;

    if( results.valid ){
      if( results.validLow ){
        results.modelError = std::abs(results.sigmaLow-results.sigma);
        results.modelVariations++;
      }//<<>>if( narrow fit is valid )
      if( results.validHigh ){
        results.modelError = std::max(results.modelError,std::abs(results.sigmaHigh-results.sigma));
        results.modelVariations++;
      }//<<>>if( wide fit is valid )
    }//<<>>if( central fit is valid )
    results.modelValid = results.modelVariations > 0 && std::isfinite(results.modelError);

    std::cout << " - Fit Sigma: " << results.sigma << " +/- " << results.esigma;
    std::cout << " ( low " << results.sigmaLow << " high " << results.sigmaHigh << " ) " << std::endl;
    return results.valid;

  }//<<>>bool ProfileTimeFit::DoFit()

};//<<>>class ProfileTimeFit

#endif
