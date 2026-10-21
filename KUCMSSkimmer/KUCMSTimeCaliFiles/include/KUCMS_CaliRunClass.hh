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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#ifndef KUCMS_CaliRunClass
#define KUCMS_CaliRunClass

// Per-run calibration state is accumulated in these records and later converted
// into crystal or trigger-tower mean, error, occupancy, and resolution maps.

struct kucms_sumCnt {

  kucms_sumCnt(){}
  kucms_sumCnt( double s, double s2, double c ) : sum(s), sumsqr(s2), cnt(c) {}

  double sum;
  double sumsqr;
  double cnt;

};//<<>>struct sumCnt

// Represents one calibration run range and owns the ROOT objects used to build
// its maps.  The class also tracks whether the range is complete and whether
// additional gain or resolution information still needs to be produced.
class kucms_CaliRunClass : public KUCMS_RootHelperBaseClass {

  public:

  kucms_CaliRunClass(){}
  //kucms_CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi );
  //~CaliRunClass();

  std::string histMapName;
  int startRun;
  int endRun;
  int lastRun;
  float lumi;

  float noise;
  float stoch;
  float stant;

  bool isNew;
  bool updated;
  bool hasGID2;
  bool doGID2;
  bool has2DResMap;
  bool hasResParams;
  bool isExternal;

  std::map<uInt,kucms_sumCnt> sumCntMap;
  std::map<uInt,double> meanMap;
  std::map<uInt,double> errMap;
  std::map<uInt,TH1F*> detIdHists;
  //detIdHists[cmsswId] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",500,-5,5);

  //void fillSumCnt( uInt detID, float val, int cnt = 1 );
  //void makeMeanMap( bool filter = false );


  //inline kucms_CaliRunClass::kucms_CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi )
  kucms_CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi )
  : histMapName(tmpxtalmap), startRun(tstart), endRun(tend), lastRun(last), lumi(tlumi)
  { isNew = true; isExternal = false; updated = false; hasGID2 = false; doGID2 = false;
    has2DResMap = false; hasResParams = false; noise = 0; stoch = 0; stant = 0; }

  //inline void kucms_CaliRunClass::makeMeanMap( bool filter ){
  bool makeMeanMap( bool filter = false ){

    meanMap.clear();
    errMap.clear();
    bool valid = true;
    unsigned int invalidBins = 0;
    //std::cout << " make mean maps: " << std::endl;
    for( auto& entry : sumCntMap ){

      double cnt = double( entry.second.cnt );
      double mean = 0;
      double err = 9;
      bool binValid = cnt > 0 && std::isfinite(cnt) && std::isfinite(entry.second.sum) &&
      std::isfinite(entry.second.sumsqr);
      if( binValid ){
        mean = entry.second.sum / cnt;
        double secondMoment = entry.second.sumsqr/cnt;
        double meanSquared = mean*mean;
        double variance = secondMoment - meanSquared;
        // Sum and sum-of-squares are persisted independently in TH2F bins. Their
        // float rounding can make a truly zero variance slightly negative after reload.
        const double momentScale = std::max({1.0,std::abs(secondMoment),std::abs(meanSquared)});
        const double doubleTolerance = 1.e-12*momentScale;
        const double floatStorageTolerance =
        8.0*std::numeric_limits<float>::epsilon()*momentScale;
        const double varianceTolerance = std::max(doubleTolerance,floatStorageTolerance);
        if( variance < 0 && variance >= -varianceTolerance ) variance = 0;
        binValid = std::isfinite(mean) && mean > -29.0 && std::isfinite(variance) && variance >= 0;
        if( binValid ){
          err = std::sqrt(variance/cnt);
          binValid = std::isfinite(err);
        }//<<>>if( binValid )
      }//<<>>if( binValid )
      if( not binValid ){
        mean = 0;
        err = 9;
        valid = false;
        invalidBins++;
      }//<<>>if( not binValid )
      //std::cout << " - calc: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl;
      auto hist = detIdHists.find(entry.first);
      if( binValid && filter && ( cnt > 10 ) && hist != detIdHists.end() && hist->second ){
        double high = mean + 3*err;
        double low = mean - 3*err;
        auto fitFunc  = new TF1("gfit","gaus",low,high);
        auto thefit = hist->second->Fit("gfit","QNRL");
        double fmean = fitFunc->GetParameter(1);
        double ferr = fitFunc->GetParError(1);
        //std::cout << " --- fit: result " << thefit << " mean " << fmean << " +/- " << ferr << std::endl;
        if( thefit != 4 && ferr != 0 && ferr < err ){ mean = fmean; err = ferr; }
        delete fitFunc;
      }//<<>>if( filter )
      if( hist != detIdHists.end() ) delete hist->second;
      meanMap[entry.first] = mean;
      errMap[entry.first] = err;
      //std::cout << " - mmm: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl;
    }//<<>>for( auto& entry : sumCntMap )
    detIdHists.clear();
    if( invalidBins > 0 ){
      std::cerr << "ERROR: makeMeanMap found " << invalidBins << " invalid calibration bins." << std::endl;
    }//<<>>if( invalidBins > 0 )
    return valid;

  }//<<>>bool makeMeanMap( bool filter = false )

  //inline void kucms_CaliRunClass::fillSumCnt( uInt detid, float val, int cnt ){
  bool fillSumCnt( uInt detid, double val, double cnt = 1, bool fillHist = false ){

    //std::cout << "Cecking " << detid << " with " << endRun << " " << lastRun <<  " " << hasGID2 << " " << doGID2  << std::endl;
    if( endRun == lastRun ){
      if( hasGID2 == true ) return true;
      if( doGID2 == false ) return true;
    }//<<>>if( ( endRun == lastRun )
    if( not std::isfinite(val) || val <= -29.0 || not std::isfinite(cnt) || cnt <= 0 ){
      std::cerr << "ERROR: refusing invalid calibration value for detID " << detid
      << ": value=" << val << " weight=" << cnt << std::endl;
      return false;
    }//<<>>if( invalid value or weight )
    //std::cout << "Filling " << detid << " with " << val << " " << cnt << std::endl;
    updated = true;
    double wval = val*cnt;
    double wvalsqr = val*val*cnt;
    auto result = sumCntMap.try_emplace(detid, 0.0, 0.0, 0.0);
    auto& sums = result.first->second;
    sums.sum += wval;
    sums.sumsqr += wvalsqr;
    sums.cnt += cnt;

    //std::cout << " -- hist Fill " << detid << " with " << val << std::endl; }
    if( fillHist && detIdHists.find(detid) != detIdHists.end() ){ detIdHists[detid]->Fill(val,cnt); }
    else if( fillHist ){
      //std::cout << " -- hist Make " << detid << " with " << val << std::endl;
      std::string histname = histMapName + "_SumCntHist_" + std::to_string( startRun ) + "_" + std::to_string(detid) ;
      detIdHists[detid] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",200,-5,5);
      detIdHists[detid]->Sumw2();
      detIdHists[detid]->Fill(val,cnt);
    }//<<>>if( detIdHists.find(detid) != detIdHists.end() )
    return true;

  }//<<>>bool CaliRunClass::fillSumCnt( uInt detID, float sum, int cnt )

  bool mergeSumCnt( uInt detid, double sum, double sumsqr, double cnt ){
    if( not std::isfinite(sum) || not std::isfinite(sumsqr) || not std::isfinite(cnt) || cnt <= 0 ){
      std::cerr << "ERROR: refusing invalid aggregate for detID " << detid << std::endl;
      return false;
    }//<<>>if( invalid aggregate )
    updated = true;
    auto result = sumCntMap.try_emplace(detid, 0.0, 0.0, 0.0);
    auto& totals = result.first->second;
    totals.sum += sum;
    totals.sumsqr += sumsqr;
    totals.cnt += cnt;
    return true;
  }//<<>>bool mergeSumCnt( uInt detID, double sum, double sumsqr, double cnt )

};//<<>>TimeCaliTagStruct

#endif
