// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

#include "KUCMS_HelperBaseClass.hh"

// ROOT/

#include "TROOT.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"

#include "TF1.h"
#include "TFormula.h"
#include "TGraph.h"

#include "TStyle.h"
#include "TString.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "TMath.h"
#include "Math/PositionVector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TMathBase.h"

#ifndef KUCMS_RootHelperBaseClass_new
#define KUCMS_RootHelperBaseClass_new

class KUCMS_RootHelperBaseClass : public KUCMS_HelperBaseClass {

	private:

	bool rootHelperUsed;


	public:

	KUCMS_RootHelperBaseClass(){ rootHelperUsed = true; };
	//~KUCMSRootHelperBaseClass(){};

//----------------------------------------------------------------------
// Fill a TH1 histogram while clamping values to the plotted bin range.
//
// If the input value is below the center of the first X bin, the first
// bin center is filled instead.  If the input value is above the center
// of the last X bin, the last bin center is filled instead.  Otherwise,
// the input value itself is filled.
//
// This prevents entries from going into the underflow or overflow bins
// and keeps all filled values inside the visible histogram range.
//
// Overloads are provided for TH1F and TH1D histograms.
//----------------------------------------------------------------------

	void fillTH1( float val, TH1F* hist ){

		if( !hist ){ std::cerr << "ERROR: null histogram pointer passed to fillTH1F" << std::endl; return;}

		auto nBins = hist->GetNbinsX();
		auto low = hist->GetBinCenter(1);
		auto high = hist->GetBinCenter(nBins);
		if( val < low ) hist->Fill( low );
		else if ( val > high ) hist->Fill( high );
		else hist->Fill( val );

	}//<<>>void fillTH1( float val, TH1F* hist )

	void fillTH1( float val, TH1D* hist ){

        if( !hist ){ std::cerr << "ERROR: null histogram pointer passed to fillTH1D" << std::endl; return;}

		auto nBins = hist->GetNbinsX();
		auto low = hist->GetBinCenter(1);
		auto high = hist->GetBinCenter(nBins);
		if( val < low ) hist->Fill( low );
		else if ( val > high ) hist->Fill( high );
		else hist->Fill( val );

	}//<<>>void fillTH1( float val, TH1D* hist )

//----------------------------------------------------------------------
// Fill a TH1D ratio histogram from numerator and denominator histograms.
//
// For each bin, the ratio numi/denom is computed and stored in result.
// The bin error is computed according to the selected RatioErrorType:
//
//   kBinomialRatio
//      Efficiency-style binomial uncertainty.  Use when the numerator is
//      a true subset of the denominator, such as Npass/Ntotal.
//
//   kIndependentRatio
//      Standard independent-error propagation.  Use when numerator and
//      denominator histograms are statistically independent samples.
//
//   kLegacyBinomial
//      Legacy binomial-style expression used in earlier versions of this
//      analysis code.  This is equivalent to the standard binomial form
//      only when the numerator error is sqrt(Npass).
//
// By default, the function uses kBinomialRatio.  Bins with denominator
// content less than or equal to the required threshold are assigned zero
// content and zero error.
//
// The input histograms are assumed to have compatible X-axis binning.
//----------------------------------------------------------------------

    enum RatioErrorType {
        kBinomialRatio    = 0,
        kIndependentRatio = 1,
        kLegacyBinomial   = 2
    };//<<>>enum RatioErrorType

    void fillRatioHist( TH1D* numi, TH1D* denom, TH1D* result, RatioErrorType errType = kBinomialRatio ){

        if( !numi || !denom || !result ){
            std::cerr << "ERROR: null histogram pointer passed to fillRatioHist" << std::endl;
            return;
        }//<<>>if( !numi || !denom || !result )

        const int nbins = numi->GetNbinsX();

        if( denom->GetNbinsX() != nbins || result->GetNbinsX() != nbins ){
            std::cerr << "ERROR: incompatible binning in fillRatioHist" << std::endl;
            return;
        }//<<>>if( denom->GetNbinsX() != nbins || result->GetNbinsX() != nbins )

        for( int ibin = 1; ibin <= nbins; ibin++ ){

            const double nc   = numi->GetBinContent(ibin);
            const double ncer = numi->GetBinError  (ibin);
            const double dc   = denom->GetBinContent(ibin);
            const double dcer = denom->GetBinError  (ibin);

            double ratio = 0.0;
            double rerr  = 0.0;

            if( dc > 20.0 ){

                ratio = nc/dc;

                if( errType == kBinomialRatio ){

                    // Efficiency-style binomial uncertainty.
                    // Use when numerator is a true subset of denominator:
                    //
                    //     ratio = Npass / Ntotal
                    //     err   = sqrt( ratio*(1-ratio) / Ntotal )
                    //
                    if( ratio >= 0.0 && ratio <= 1.0 ){
                        rerr = std::sqrt( ratio*(1.0-ratio)/dc );
                    } else {
                        std::cerr << "WARNING: binomial ratio outside [0,1] in fillRatioHist"
                                  << " bin=" << ibin
                                  << " num=" << nc
                                  << " denom=" << dc
                                  << " ratio=" << ratio
                                  << std::endl;
                    }//<<>>if( ratio >= 0.0 && ratio <= 1.0 )

                } else if( errType == kIndependentRatio ){

                    // Generic independent-histogram ratio uncertainty.
                    // Use when numerator and denominator are statistically independent:
                    //
                    //     ratio = N / D
                    //     err^2 = (nerr/D)^2 + (N*derr/D^2)^2
                    //
                    rerr = std::sqrt( sq2(ncer/dc) + sq2(nc*dcer/sq2(dc)) );

                } else if( errType == kLegacyBinomial ){

                    // Legacy binomial-style form.
                    // This is equivalent to the binomial expression when:
                    //
                    //     ncer = sqrt(nc)
                    //
                    // Protect the sqrt against small negative values from floating point
                    // precision or from input errors that do not satisfy the assumptions.
                    //
                    const double rerr2 = sq2(ncer/dc) - sq2(ratio)/dc;

                    if( rerr2 >= 0.0 ){
                        rerr = std::sqrt(rerr2);
                    } else {
                        rerr = 0.0;
                        std::cerr << "WARNING: negative legacy ratio error squared in fillRatioHist"
                                  << " bin=" << ibin
                                  << " num=" << nc
                                  << " denom=" << dc
                                  << " ratio=" << ratio
                                  << " rerr2=" << rerr2
                                  << std::endl;
                    }//<<>>if( rerr2 >= 0.0 )

                } else {

                    std::cerr << "WARNING: unknown RatioErrorType in fillRatioHist"
                              << " errType=" << errType
                              << " using binomial error"
                              << std::endl;

                    if( ratio >= 0.0 && ratio <= 1.0 ){
                        rerr = std::sqrt( ratio*(1.0-ratio)/dc );
                    }//<<>>if( ratio >= 0.0 && ratio <= 1.0 )

                }//<<>>if( errType == kBinomialRatio )

            }//<<>>if( dc > 20.0 )

            result->SetBinContent( ibin, ratio );
            result->SetBinError  ( ibin, rerr  );

        }//<<>>for( int ibin = 1; ibin <= nbins; ibin++ )

    }//<<>>void fillRatioHist( TH1D* numi, TH1D* denom, TH1D* result, RatioErrorType errType = kBinomialRatio )


//----------------------------------------------------------------------
// Normalize each X-bin slice of a TH2D histogram to unit area.
//
// For each X bin, the integral over all Y bins is computed.  The content
// and error of every Y bin in that X slice are then divided by this
// integral, so that each vertical Y distribution at fixed X is normalized
// independently.
//
// Underflow and overflow bins are not included in the normalization,
// since the integral is computed over Y bins 1 through nYbins only.
//
// If the integral for a given X-bin slice is zero, that slice is skipped
// and no bin contents or errors are modified for that X bin.
//----------------------------------------------------------------------

	void normTH2D(TH2D* hist){


        if( !hist ){ std::cerr << "ERROR: null histogram pointer passed to normTH2D" << std::endl; return;}
		std::cout << "Normalizing " << " hist: " << hist->GetName() << std::endl;

		const auto nXbins = hist->GetNbinsX();
		const auto nYbins = hist->GetNbinsY();

		for (auto ibinX = 1; ibinX <= nXbins; ibinX++){

			const auto norm = hist->Integral(ibinX,ibinX,1,nYbins);
			if( norm == 0.0 ) continue;
			for (auto ibinY = 1; ibinY <= nYbins; ibinY++){

				// get content/error
				auto content = hist->GetBinContent(ibinX,ibinY);
				auto error   = hist->GetBinError  (ibinX,ibinY);
				// set new contents
				content /= norm;
				error /= norm;
				hist->SetBinContent(ibinX,ibinY,content);
				hist->SetBinError  (ibinX,ibinY,error);

			}//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
		}//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

	}//<<>>void NormTH2D(TH2D* hist){

	void normTH2F(TH2F* hist){

        if( !hist ){ std::cerr << "ERROR: null histogram pointer passed to normTH2F" << std::endl; return;}
		std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;

		const auto nXbins = hist->GetNbinsX();
		const auto nYbins = hist->GetNbinsY();

		for (auto ibinX = 1; ibinX <= nXbins; ibinX++){

			const auto norm = hist->Integral(ibinX,ibinX,1,nYbins);
			if( norm == 0.0 ) continue;
			for (auto ibinY = 1; ibinY <= nYbins; ibinY++){

				// get content/error
				auto content = hist->GetBinContent(ibinX,ibinY);
				auto error   = hist->GetBinError  (ibinX,ibinY);
				// set new contents
				content /= norm;
				error /= norm;
				hist->SetBinContent(ibinX,ibinY,content);
				hist->SetBinError  (ibinX,ibinY,error);

			}//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
		}//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

	}//<<>>void NormTH2F(TH2F* hist){

//----------------------------------------------------------------------
// Normalize a TH1D histogram to unit area.
//
// The bin contents and bin errors are divided by the histogram integral,
// so that the total integral of the histogram becomes 1.0.  Underflow
// and overflow bins are not included in the normalization, since
// TH1::Integral() is called with the default bin range.
//
// If the histogram integral is zero, no bin contents or errors are
// modified.
//----------------------------------------------------------------------

	void normTH1D(TH1D* hist){

        if( !hist ){ std::cerr << "ERROR: null histogram pointer passed to normTH1D" << std::endl; return;}
		std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;

		const auto nBins = hist->GetNbinsX();
		const auto norm = hist->Integral();
        if( norm == 0.0 ) return;
		for (auto ibinX = 1; ibinX <= nBins; ibinX++){

			// get content/error
			auto content = hist->GetBinContent(ibinX);
			auto error   = hist->GetBinError(ibinX);
			// set new contents
			content /= norm;
			error /= norm;
			hist->SetBinContent(ibinX,content);
			hist->SetBinError  (ibinX,error);

		}//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)

	}//<<>>void NormTH1D(TH1D* hist)

//----------------------------------------------------------------------
// Build a fitted profile of a TH2D by fitting each Y projection.
//
// For each X bin of the input TH2D, a Y-axis projection is made and a
// Gaussian fit is performed in a restricted window around the projection
// mean.  The fit range is defined as:
//
//     mean +/- range * stddev
//
// where range is a configurable scale factor with default value 0.2.
//
// If the projection has nonzero width and sufficient population, the
// fitted Gaussian mean and its uncertainty are stored in the profile
// histogram.  The fit chi2/ndf is stored in fithist for diagnostic use.
//
// A fit result is accepted only when the fit has positive NDF and the
// fitted-mean uncertainty is smaller than the projection standard
// deviation.  Failed or poorly constrained fits leave the corresponding
// output bins unchanged.
//
// The input and output histograms are assumed to have compatible X-axis
// binning.
//----------------------------------------------------------------------

    void profileTH2D( TH2D* nhist, TH1D* prof, TH1D* fithist, float range = 0.2 ){

        if( !nhist || !prof || !fithist ){
            std::cerr << "ERROR: null histogram pointer passed to profileTH2D" << std::endl;
            return;
        }//<<>>if( !nhist || !prof || !fithist )

        std::cout << "Profile " << " hist: " << nhist->GetName() << std::endl;

        const int nXBins = nhist->GetNbinsX();

        if( prof->GetNbinsX() != nXBins || fithist->GetNbinsX() != nXBins ){
            std::cerr << "ERROR: incompatible X binning in profileTH2D" << std::endl;
            return;
        }//<<>>if( prof->GetNbinsX() != nXBins || fithist->GetNbinsX() != nXBins )

        for( int ibinX = 1; ibinX <= nXBins; ibinX++ ){

            const std::string pname = "temp_py_" + std::to_string(ibinX);

            TH1D* phist = nhist->ProjectionY( pname.c_str(), ibinX, ibinX );
            if( !phist ) continue;

            phist->SetDirectory(nullptr);

            const double entries = phist->GetEntries();
            const double mean    = phist->GetMean();
            const double stdv    = phist->GetStdDev();
            const double norm    = phist->GetBinContent( phist->GetMaximumBin() );

            if( entries <= 0.0 || stdv <= 0.0 || norm <= 1.0 ){
                delete phist;
                continue;
            }//<<>>if( entries <= 0.0 || stdv <= 0.0 || norm <= 1.0 )

            const double low  = mean - range*stdv;
            const double high = mean + range*stdv;

            TF1* tmp_fit = new TF1( "tmp_fit", "gaus", low, high );

            tmp_fit->SetParameter( 0, norm );
            tmp_fit->SetParameter( 1, mean );
            tmp_fit->SetParameter( 2, stdv );

            phist->Fit( tmp_fit, "RBQ0" );

            const double fmean = tmp_fit->GetParameter(1);
            const double error = tmp_fit->GetParError(1);
            const int    fNdf  = tmp_fit->GetNDF();
            const double fChi2 = tmp_fit->GetChisquare();

            if( fNdf > 0 && error < stdv ){

                const double fChi2Ndf = fChi2/fNdf;

                fithist->SetBinContent( ibinX, fChi2Ndf );
                fithist->SetBinError  ( ibinX, 0.0 );

                prof->SetBinContent( ibinX, fmean );
                prof->SetBinError  ( ibinX, error );

            }//<<>>if( fNdf > 0 && error < stdv )

            delete tmp_fit;
            delete phist;

        }//<<>>for( int ibinX = 1; ibinX <= nXBins; ibinX++ )

    }//<<>>void profileTH2D( TH2D* nhist, TH1D* prof, TH1D* fithist, float range = 0.2 )

//----------------------------------------------------------------------
// Scale a TH2F histogram by its X and/or Y bin widths.
//
// This function is used to convert between bin-content and bin-density
// conventions for histograms with variable-width binning.  If varBinsX is
// true, the X-bin width is included in the scale factor.  If varBinsY is
// true, the Y-bin width is included in the scale factor.
//
// When isUp is false, each bin content and error are multiplied by the
// selected bin-width factor:
//
//     scale = binWidthX * binWidthY
//
// using only the axes requested by varBinsX and varBinsY.
//
// When isUp is true, each bin content and error are divided by the
// selected bin-width factor:
//
//     scale = 1.0 / ( binWidthX * binWidthY )
//
// This allows the same function to apply or remove bin-width scaling.
// Underflow and overflow bins are not modified.
//
// If the histogram pointer is null, the function prints an error and
// returns without modifying anything.
//----------------------------------------------------------------------

	void scaleHist(TH2F* hist, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY){

        if( !hist ){ std::cerr << "ERROR: null histogram pointer passed to scaleHist" << std::endl; return;}
		std::cout << "Scaling " << (isUp?"up":"down") << " hist: " << hist->GetName() << std::endl;

		const auto nXbins = hist->GetNbinsX();
		const auto nYbins = hist->GetNbinsY();
		for (auto ibinX = 1; ibinX <= nXbins; ibinX++){
			const auto binwidthX = hist->GetXaxis()->GetBinWidth(ibinX);
			for (auto ibinY = 1; ibinY <= nYbins; ibinY++){

				const auto binwidthY = hist->GetYaxis()->GetBinWidth(ibinY);
				// get multiplier/divisor
				auto multiplier = 1.f;
				if( varBinsX ) multiplier *= binwidthX;
				if( varBinsY ) multiplier *= binwidthY;
				auto scale = ( not isUp ) ? multiplier : 1/multiplier;
				// get content/error
				auto content = hist->GetBinContent(ibinX,ibinY)*scale;
				auto error = hist->GetBinError(ibinX,ibinY)*scale;
				// set new contents
				hist->SetBinContent(ibinX,ibinY,content);
				hist->SetBinError  (ibinX,ibinY,error);

			}//<<>>for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++)
		}//<<>>for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++)

	}//<<>>void scaleHist(TH2F *& hist, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY)

//----------------------------------------------------------------------
// Divide a TH2D numerator histogram by a denominator histogram with a
// denominator-content threshold.
//
// For each X,Y bin, the numerator bin content is divided by the
// corresponding denominator bin content only if the denominator content
// is greater than the supplied threshold.  Bins failing the threshold are
// set to zero content and zero error.
//
// The result is written back into the numerator histogram in place.
//
// This function is intended for ratio-style maps where poorly populated
// denominator bins should be suppressed.  Underflow and overflow bins are
// not modified.
//
// The input histograms are assumed to have compatible X- and Y-axis
// binning.  If either histogram pointer is null, the function prints an
// error and returns without modifying anything.
//----------------------------------------------------------------------

	void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

        if( !numi || !denom ){ std::cerr << "ERROR: null histogram pointer passed to thresDivTH2D" << std::endl; return;}
		std::cout << "Threshold Division - " << " hist: " << numi->GetName() << std::endl;

		const auto nXbins = numi->GetNbinsX();
		const auto nYbins = numi->GetNbinsY();

		for (auto ibinX = 1; ibinX <= nXbins; ibinX++){
			for (auto ibinY = 1; ibinY <= nYbins; ibinY++){

				// get content/error
				auto ncontent = numi->GetBinContent(ibinX,ibinY);
				auto nerror   = numi->GetBinError  (ibinX,ibinY);
				auto dcontent = denom->GetBinContent(ibinX,ibinY);
				auto derror   = denom->GetBinError  (ibinX,ibinY);
				// set new contents
				double content = 0.0;
				double error   = 0.0;
				if( dcontent > thres ){
    				content = ncontent/dcontent;
    				error   = std::sqrt( sq2(nerror/dcontent) + sq2(ncontent*derror/sq2(dcontent)) );
				}//<<>>if( dcontent > thres )
				numi->SetBinContent(ibinX,ibinY,content);
				numi->SetBinError  (ibinX,ibinY,error);

			}//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
		}//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

	}//<<>>void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

//----------------------------------------------------------------------
// Compute summary statistics for a vector of floating-point values, with
// optional per-value weights.
//
// If no weight vector is provided, unit weights are used for all values.
// If a weight vector is provided, it must have the same size as the value
// vector.  The temporary ROOT histogram is filled with the supplied
// weights and with the requested histogram bounds.
//
// Returned values are stored in the following order:
//
//   [0]  histogram weighted mean
//   [1]  histogram mean error
//   [2]  histogram median from ROOT quantiles
//   [3]  histogram median error estimate, 1.2533*meanError
//   [4]  histogram RMS
//   [5]  histogram skewness
//   [6]  direct weighted mean
//   [7]  direct weighted standard deviation
//   [8]  direct weighted RMS
//   [9]  direct weighted mean error
//   [10] direct unweighted median
//   [11] direct median error estimate, 1.2533*meanError
//   [12] total event weight
//
// The direct weighted mean, standard deviation, and RMS are computed from
// the values and weights directly.  The effective number of entries used
// for the weighted mean error is:
//
//     Neff = (sum w)^2 / sum(w^2)
//
// The direct median is the ordinary unweighted median of the input values.
// The factor 1.2533 is the large-sample conversion between the standard
// error of the mean and the standard error of the median for an
// approximately normal distribution.
//
// If the input is empty, the weights are inconsistent, or the total weight
// is zero, sentinel values are returned and a warning is printed.
//----------------------------------------------------------------------


    std::vector<float> getDistStats( std::vector<float> values, std::vector<float> wts = {}, 
                                     const float hLow = -25.0, const float hHigh = 25.0 ){

        std::vector<float> results;

        const int size   = values.size();
        const int wtsize = wts.size();

        if( size == 0 ){

            std::cout << "value group is empty" << std::endl;

            results.push_back(-29.25);//0 hist mean
            results.push_back(99.9);  //1 hist mean error
            results.push_back(-29.25);//2 hist median
            results.push_back(99.9);  //3 hist median error
            results.push_back(99.9);  //4 hist RMS
            results.push_back(9.9);   //5 hist skewness
            results.push_back(-29.25);//6 direct weighted mean
            results.push_back(99.9);  //7 direct weighted stdev
            results.push_back(99.9);  //8 direct weighted RMS
            results.push_back(99.9);  //9 direct weighted mean error
            results.push_back(-29.25);//10 direct median
            results.push_back(99.9);  //11 direct median error
            results.push_back(-10.0); //12 total weight

            return results;

        }//<<>>if( size == 0 )

        const bool useUnitWeights = ( wtsize == 0 );

        if( !useUnitWeights && size != wtsize ){

            std::cout << "value & weight groups not same size" << std::endl;

            results.push_back(-29.25);//0 hist mean
            results.push_back(99.9);  //1 hist mean error
            results.push_back(-29.25);//2 hist median
            results.push_back(99.9);  //3 hist median error
            results.push_back(99.9);  //4 hist RMS
            results.push_back(9.9);   //5 hist skewness
            results.push_back(-29.25);//6 direct weighted mean
            results.push_back(99.9);  //7 direct weighted stdev
            results.push_back(99.9);  //8 direct weighted RMS
            results.push_back(99.9);  //9 direct weighted mean error
            results.push_back(-29.25);//10 direct median
            results.push_back(99.9);  //11 direct median error
            results.push_back(-10.0); //12 total weight

            return results;

        }//<<>>if( !useUnitWeights && size != wtsize )

        TH1D valueDist( "temp_valueDist", "temp_valueDist", 500, hLow, hHigh );
        valueDist.SetDirectory(nullptr);

        double wtot  = 0.0;
        double wtot2 = 0.0;
        double sumwx = 0.0;
        double sumwx2 = 0.0;

        for( int it = 0; it < size; it++ ){

            const double value = values[it];
            const double wt    = useUnitWeights ? 1.0 : wts[it];

            valueDist.Fill( value, wt );

            wtot   += wt;
            wtot2  += wt*wt;
            sumwx  += wt*value;
            sumwx2 += wt*value*value;

        }//<<>>for( int it = 0; it < size; it++ )

        if( wtot == 0.0 ){

            std::cout << "total weight is zero" << std::endl;

            results.push_back(-29.25);//0 hist mean
            results.push_back(99.9);  //1 hist mean error
            results.push_back(-29.25);//2 hist median
            results.push_back(99.9);  //3 hist median error
            results.push_back(99.9);  //4 hist RMS
            results.push_back(9.9);   //5 hist skewness
            results.push_back(-29.25);//6 direct weighted mean
            results.push_back(99.9);  //7 direct weighted stdev
            results.push_back(99.9);  //8 direct weighted RMS
            results.push_back(99.9);  //9 direct weighted mean error
            results.push_back(-29.25);//10 direct median
            results.push_back(99.9);  //11 direct median error
            results.push_back(-10.0); //12 total weight

            return results;

        }//<<>>if( wtot == 0.0 )

        double hmedvalue = -9.9;
        const double nEntries = valueDist.GetEntries();

        if( nEntries > 0.0 ){

            double quant = 0.5;
            valueDist.ComputeIntegral();
            valueDist.GetQuantiles( 1, &hmedvalue, &quant );

        }//<<>>if( nEntries > 0.0 )

        const double herror = valueDist.GetMeanError();

        results.push_back( valueDist.GetMean() );      //0 histogram weighted mean
        results.push_back( herror );                   //1 histogram mean error
        results.push_back( hmedvalue );                //2 histogram median
        results.push_back( 1.2533*herror );            //3 histogram median error, normal approx
        results.push_back( valueDist.GetRMS() );       //4 histogram RMS
        results.push_back( valueDist.GetSkewness() );  //5 histogram skewness

        const double rmu = sumwx/wtot;
        results.push_back( rmu );                      //6 direct weighted mean

        double var = (sumwx2/wtot) - rmu*rmu;
        if( var < 0.0 ) var = 0.0;

        const double rsd = std::sqrt(var);
        results.push_back( rsd );                      //7 direct weighted stdev

        const double wrms = std::sqrt(sumwx2/wtot);
        results.push_back( wrms );                     //8 direct weighted RMS

        double neff = 0.0;
        if( wtot2 > 0.0 ) neff = wtot*wtot/wtot2;

        double err = 99.9;
        if( neff > 0.0 ) err = rsd/std::sqrt(neff);

        results.push_back( err );                      //9 direct weighted mean error

        std::sort( values.begin(), values.end() );

        double median = 0.0;
        if( size%2 == 0 ){
            median = 0.5*( values[(size/2)-1] + values[size/2] );
        } else {
            median = values[size/2];
        }//<<>>if( size%2 == 0 )

        results.push_back( median );                   //10 direct unweighted median
        results.push_back( 1.2533*err );               //11 direct median error, normal approx
        results.push_back( wtot );                     //12 total weight

        return results;

    }//>>>> std::vector<float> getDistStats( std::vector<float> values, std::vector<float> wts = {}, const float hLow = -25.0, const float hHigh = 25.0 )


//----------------------------------------------------------------------
// Return the index of the largest eigenvalue in a 2D eigenvalue vector.
//----------------------------------------------------------------------
    int getMaxEigenIndex2D( const TVectorD& eigenVal ){

        int index = 1;

        if( eigenVal(0) >= eigenVal(1) ){
            index = 0;
        }//<<>>if( eigenVal(0) >= eigenVal(1) )

        return index;

    }//<<>>int getMaxEigenIndex2D( const TVectorD& eigenVal )


//----------------------------------------------------------------------
// Return the indices of the largest, middle, and smallest eigenvalues in
// a 3D eigenvalue vector.
//
// The returned vector contains:
//
//   [0] index of largest eigenvalue
//   [1] index of middle eigenvalue
//   [2] index of smallest eigenvalue
//----------------------------------------------------------------------
    std::vector<int> getEigenOrder3D( const TVectorD& eigenVal ){

        std::vector<int> order;
        order.push_back(0);
        order.push_back(1);
        order.push_back(2);

        std::sort( order.begin(), order.end(),
                   [&]( const int a, const int b ){ return eigenVal(a) > eigenVal(b); } );

        return order;

    }//<<>>std::vector<int> getEigenOrder3D( const TVectorD& eigenVal )

//----------------------------------------------------------------------
// Compute the dominant weighted angular eigen-axis for a rechit group.
//
// This overload builds a 2D angular shape matrix using weighted sin/cos
// moments of the input coordinate:
//
//   [ sin^2   sin*cos ]
//   [ sin*cos cos^2   ]
//
// The matrix is diagonalized and the eigenvector corresponding to the
// largest eigenvalue is returned.
//
// Returns:
//   [0] dominant eigenvector x component
//   [1] dominant eigenvector y component
//   [2] dominant eigenvalue
//
// The input coordinate and weight vectors are required to have the same
// nonzero size.
//----------------------------------------------------------------------
    std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> wts ){
    // spherical / angular

        std::vector<float> results;

        const int size   = xs.size();
        const int wtsize = wts.size();

        if( size == 0 || size != wtsize ){

            std::cerr << "ERROR: invalid input vectors passed to getRhGrpEigen(xs,wts)"
                      << " xs.size()=" << size
                      << " wts.size()=" << wtsize
                      << std::endl;

            results.push_back(0.0);//0 ex
            results.push_back(0.0);//1 ey
            results.push_back(0.0);//2 eigenvalue

            return results;

        }//<<>>if( size == 0 || size != wtsize )

        const double swts = vfsum(wts);

        if( swts == 0.0 ){

            std::cerr << "ERROR: zero total weight passed to getRhGrpEigen(xs,wts)" << std::endl;

            results.push_back(0.0);//0 ex
            results.push_back(0.0);//1 ey
            results.push_back(0.0);//2 eigenvalue

            return results;

        }//<<>>if( swts == 0.0 )

        const double ts2 = wsin2( xs, wts );
        const double tc2 = wcos2( xs, wts );
        const double tsc = wsincos( xs, wts );

        double array[] = { ts2, tsc,
                           tsc, tc2 };

        TMatrixDSym mat( 2, array );
        TMatrixDSymEigen eigen(mat);

        const TVectorD& eigenVal = eigen.GetEigenValues();
        const TMatrixD& eigenVec = eigen.GetEigenVectors();

        const int index = getMaxEigenIndex2D( eigenVal );

        const double ex = eigenVec(index,0);
        const double ey = eigenVec(index,1);
        const double ev = eigenVal(index);

        results.push_back(ex);//0
        results.push_back(ey);//1
        results.push_back(ev);//2

        return results;

    }//<<>>std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> wts )

//----------------------------------------------------------------------
// Compute the dominant weighted 2D covariance eigen-axis for a rechit
// group.
//
// This overload builds a weighted x-y covariance matrix:
//
//   [ var_x   cov_xy ]
//   [ cov_xy  var_y  ]
//
// The matrix is diagonalized and the eigenvector corresponding to the
// largest eigenvalue is returned.
//
// Returns:
//   [0] dominant eigenvector x component
//   [1] dominant eigenvector y component
//   [2] dominant eigenvalue
//
// The input coordinate and weight vectors are required to have the same
// nonzero size and nonzero total weight.
//----------------------------------------------------------------------
    std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> ys, std::vector<float> wts ){

        std::vector<float> results;

        const int size   = xs.size();
        const int ysize  = ys.size();
        const int wtsize = wts.size();

        if( size == 0 || size != ysize || size != wtsize ){

            std::cerr << "ERROR: invalid input vectors passed to getRhGrpEigen(xs,ys,wts)"
                      << " xs.size()=" << size
                      << " ys.size()=" << ysize
                      << " wts.size()=" << wtsize
                      << std::endl;

            results.push_back(0.0);//0 ex
            results.push_back(0.0);//1 ey
            results.push_back(0.0);//2 eigenvalue

            return results;

        }//<<>>if( size == 0 || size != ysize || size != wtsize )

        const double swts = vfsum(wts);

        if( swts == 0.0 ){

            std::cerr << "ERROR: zero total weight passed to getRhGrpEigen(xs,ys,wts)" << std::endl;

            results.push_back(0.0);//0 ex
            results.push_back(0.0);//1 ey
            results.push_back(0.0);//2 eigenvalue

            return results;

        }//<<>>if( swts == 0.0 )

        const double mean_x = mean( xs, wts );
        const double mean_y = mean( ys, wts );

        const double var_x  = var ( xs, mean_x, wts, swts );
        const double var_y  = var ( ys, mean_y, wts, swts );
        const double var_xy = cvar( xs, mean_x, ys, mean_y, wts, swts );

        double array[] = { var_x,  var_xy,
                           var_xy, var_y  };

        TMatrixDSym mat( 2, array );
        TMatrixDSymEigen eigen(mat);

        const TVectorD& eigenVal = eigen.GetEigenValues();
        const TMatrixD& eigenVec = eigen.GetEigenVectors();

        const int index = getMaxEigenIndex2D( eigenVal );

        const double ex = eigenVec(index,0);
        const double ey = eigenVec(index,1);
        const double ev = eigenVal(index);

        results.push_back(ex);//0
        results.push_back(ey);//1
        results.push_back(ev);//2

        return results;

    }//<<>>std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> ys, std::vector<float> wts )

//----------------------------------------------------------------------
// Compute the dominant weighted 3D covariance eigen-axis for a rechit
// group.
//
// This overload builds a weighted x-y-z covariance matrix:
//
//   [ var_x   cov_xy  cov_xz ]
//   [ cov_xy  var_y   cov_yz ]
//   [ cov_xz  cov_yz  var_z  ]
//
// The matrix is diagonalized and the eigenvector corresponding to the
// largest eigenvalue is returned.  The final returned value is a
// normalized dominant-eigenvalue shape measure, rather than the raw
// largest eigenvalue.
//
// Returns:
//   [0] dominant eigenvector x component
//   [1] dominant eigenvector y component
//   [2] dominant eigenvector z component
//   [3] normalized dominant-eigenvalue shape measure
//
// The input coordinate and weight vectors are required to have the same
// nonzero size and nonzero total weight.
//----------------------------------------------------------------------
    std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> ys, std::vector<float> zs, std::vector<float> wts ){
    // ieipt

        std::vector<float> results;

        const int size   = xs.size();
        const int ysize  = ys.size();
        const int zsize  = zs.size();
        const int wtsize = wts.size();

        if( size == 0 || size != ysize || size != zsize || size != wtsize ){

            std::cerr << "ERROR: invalid input vectors passed to getRhGrpEigen(xs,ys,zs,wts)"
                      << " xs.size()=" << size
                      << " ys.size()=" << ysize
                      << " zs.size()=" << zsize
                      << " wts.size()=" << wtsize
                      << std::endl;

            results.push_back(0.0);//0 ex
            results.push_back(0.0);//1 ey
            results.push_back(0.0);//2 ez
            results.push_back(0.0);//3 shape eigenvalue

            return results;

        }//<<>>if( size == 0 || size != ysize || size != zsize || size != wtsize )

        const double swts = vfsum(wts);

        if( swts == 0.0 ){

            std::cerr << "ERROR: zero total weight passed to getRhGrpEigen(xs,ys,zs,wts)" << std::endl;

            results.push_back(0.0);//0 ex
            results.push_back(0.0);//1 ey
            results.push_back(0.0);//2 ez
            results.push_back(0.0);//3 shape eigenvalue

            return results;

        }//<<>>if( swts == 0.0 )

        const double mean_x = mean( xs, wts );
        const double mean_y = mean( ys, wts );
        const double mean_z = mean( zs, wts );

        const double var_x  = var ( xs, mean_x, wts, swts );
        const double var_y  = var ( ys, mean_y, wts, swts );
        const double var_z  = var ( zs, mean_z, wts, swts );
        const double var_xy = cvar( xs, mean_x, ys, mean_y, wts, swts );
        const double var_xz = cvar( xs, mean_x, zs, mean_z, wts, swts );
        const double var_yz = cvar( ys, mean_y, zs, mean_z, wts, swts );

        double array[] = { var_x,  var_xy, var_xz,
                           var_xy, var_y,  var_yz,
                           var_xz, var_yz, var_z  };

        TMatrixDSym mat( 3, array );
        TMatrixDSymEigen eigen(mat);

        const TVectorD& eigenVal = eigen.GetEigenValues();
        const TMatrixD& eigenVec = eigen.GetEigenVectors();

        std::vector<int> order = getEigenOrder3D( eigenVal );

        const int index  = order[0];
        const int index2 = order[1];
        const int index3 = order[2];

        if( eigenVal(index) == eigenVal(index2) && eigenVal(index) == eigenVal(index3) ){
            std::cout << " -- rhGrp is Spherical" << std::endl;
        } else if( eigenVal(index2) == eigenVal(index3) ){
            std::cout << " -- rhGrp is a Flattened Sphere" << std::endl;
        }//<<>>if( eigenVal(index) == eigenVal(index2) && eigenVal(index) == eigenVal(index3) )

        const double ev23hypo = hypo( eigenVal(index2), eigenVal(index3) );
        const double ev1Shypo = hypo( eigenVal(index), ev23hypo );

        double speigenval = 0.0;
        if( ev1Shypo > 0.0 ){
            speigenval = sq2( eigenVal(index)/ev1Shypo );
        }//<<>>if( ev1Shypo > 0.0 )

        const double ex = eigenVec(index,0);
        const double ey = eigenVec(index,1);
        const double ez = eigenVec(index,2);
        const double ev = speigenval;

        results.push_back(ex);//0
        results.push_back(ey);//1
        results.push_back(ez);//2
        results.push_back(ev);//3

        return results;

    }//<<>>std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> ys, std::vector<float> zs, std::vector<float> wts )
    
//----------------------------------------------------------------------
// Check whether a string ends with a requested suffix.
//
// Returns true if the input string ends with the supplied suffix and
// false otherwise.  This is used for simple filename/path filtering, such
// as selecting files ending in ".root".
//
// If the suffix is longer than the input string, the function returns
// false.
//----------------------------------------------------------------------
    bool endsWith( const std::string& str, const std::string& suffix ){

        if( str.size() < suffix.size() ){
            return false;
        }//<<>>if( str.size() < suffix.size() )

        return std::equal( suffix.rbegin(), suffix.rend(), str.rbegin() );

    }//<<>>bool endsWith( const std::string& str, const std::string& suffix )

//----------------------------------------------------------------------
// Run a shell command and collect non-empty stdout lines.
//
// The command is executed with popen(), and each non-empty output line is
// returned as one string in the output vector.  Trailing newline and
// carriage-return characters are removed from each line.
//
// If the command cannot be started, an empty vector is returned and an
// error is printed.  If the command exits with a nonzero status, the
// collected output is still returned, but a warning is printed.
//
// Important: popen() runs the command through the shell.  Do not pass
// untrusted or unsanitized user input into this helper.
//----------------------------------------------------------------------
    std::vector<std::string> runCommandLines( const std::string& cmd ){

        std::vector<std::string> lines;

        FILE* pipe = popen( cmd.c_str(), "r" );

        if( !pipe ){
            std::cerr << "ERROR: could not run command: " << cmd << std::endl;
            return lines;
        }//<<>>if( !pipe )

        char buffer[4096];

        while( fgets( buffer, sizeof(buffer), pipe ) != nullptr ){

            std::string line(buffer);

            while( !line.empty() && ( line.back() == '\n' || line.back() == '\r' ) ){
                line.pop_back();
            }//<<>>while( !line.empty() && ( line.back() == '\n' || line.back() == '\r' ) )

            if( !line.empty() ){
                lines.push_back(line);
            }//<<>>if( !line.empty() )

        }//<<>>while( fgets( buffer, sizeof(buffer), pipe ) != nullptr )

        const int status = pclose(pipe);

        if( status != 0 ){
            std::cerr << "WARNING: command exited nonzero: " << cmd
                      << " status=" << status
                      << std::endl;
        }//<<>>if( status != 0 )

        return lines;

    }//<<>>std::vector<std::string> runCommandLines( const std::string& cmd )
    
    // ------------------------------------------------------------
    // Recursively scan EOS using xrdfs.
    // eosh example: "cmseos.fnal.gov"
    // eosDir example: "/store/user/lpcsusylep/jaking/KUCMSNtuple/someDir"
    // matchStr example: "somename.root"
    //   - if matchString is empty, all .root files are accepted
    //   - if matchString is non-empty, filename/path must contain it
    // ------------------------------------------------------------
    std::vector<std::string> findEOSRootFiles( const std::string& eosh, const std::string& eosDir, const std::string& matchStr = ""){
 
       std::vector<std::string> rootFiles;
    
        std::string cmd = "xrdfs " + eosh + " ls -R " + eosDir;
        std::cout << "Scanning EOS with command:\n  " << cmd << std::endl;
    
        auto lines = runCommandLines(cmd);
    
        for (const auto& path : lines) {
    
            // xrdfs returns full EOS paths like:
            // /store/user/...
            if (!endsWith(path, ".root")) continue;
    
            if (!matchStr.empty()) {
                if (path.find(matchStr) == std::string::npos) continue;
            }
    
            std::string fullXrdPath = "root://" + eosh + "/" + path;
    
            // This produces:
            // root://cmseos.fnal.gov//store/user/...
            rootFiles.push_back(fullXrdPath);
        }
    
        return rootFiles;
    }//<<>>std::vector<std::string> findEOSRootFiles( const std::string& eosh, const std::string& eosDir, const std::string& matchStr = "")

};//<<>> class KUCMSRootHelperBaseClass : KUCMSHelperBaseClass

#endif
//----------------------------------------------------------------------------------------------------------------------



