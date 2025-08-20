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
 
#ifndef KUCMSCaliRunClass 
#define KUCMSCaliRunClass 

struct sumCnt { 

    sumCnt(){} 
    sumCnt( float s, float s2, int c ) : sum(s), sumsqr(s2), cnt(c) {}  
 
    float sum; 
    float sumsqr; 
    int cnt; 
 
};//<<>>struct sumCnt 

class CaliRunClass : public KUCMSRootHelperBaseClass { 
     
    public:  
 
    CaliRunClass(){} 
    CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi ); 
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
    bool has2DResMap; 
    bool hasResParams; 
    bool isExternal; 
 
    std::map<uInt,sumCnt> sumCntMap; 
    std::map<uInt,float> meanMap; 
    std::map<uInt,float> errMap; 
    std::map<uInt,TH1F*> detIdHists; 
    //detIdHists[cmsswId] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",500,-5,5);     
 
    void fillSumCnt( uInt detID, float val, int cnt = 1 ); 
    void makeMeanMap( bool filter = false ); 

};//<<>>TimeCaliTagStruct 
 
inline CaliRunClass::CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi ) 
    : histMapName(tmpxtalmap), startRun(tstart), endRun(tend), lastRun(last), lumi(tlumi)  
    { isNew = true; isExternal = false; updated = false; has2DResMap = false; hasResParams = false; noise = 0; stoch = 0; stant = 0; } 

inline void CaliRunClass::makeMeanMap( bool filter ){ 
 
    meanMap.clear(); 
    errMap.clear(); 
    //std::cout << " make mean maps: " << std::endl; 
    for( auto& entry : sumCntMap ){ 
 
        float cnt = float( entry.second.cnt ); 
        if( cnt == 0 ) cnt = 1; 
        float mean = entry.second.sum / cnt;  
        float err = sqrt( (entry.second.sumsqr/cnt - mean*mean)/cnt ); 
        if( err == 0 ){ err = 9; mean = 0; } 
        //std::cout << " - calc: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl; 
        if( filter && cnt > 100 ){ 
            auto fitFunc  = new TF1("gfit","gaus",-3.0,3.0); 
            auto thefit = detIdHists[entry.first]->Fit("gfit","QNRL"); 
            float fmean = fitFunc->GetParameter(1); 
            float ferr = fitFunc->GetParError(1); 
            //std::cout << " --- fit: result " << thefit << " mean " << fmean << " +/- " << ferr << std::endl; 
            if( thefit != 4 && ferr != 0 && ferr < err ){ mean = fmean; err = ferr; } 
            delete fitFunc; 
        }//<<>>if( filter ) 
        delete detIdHists[entry.first]; 
        meanMap[entry.first] = mean; 
        errMap[entry.first] = err;  
        //std::cout << " - mmm: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl; 
    }//<<>>for( auto& entry : sumCntMap ) 
 
}//<<>>void CaliRunClass::makeMeanMap() 
 
inline void CaliRunClass::fillSumCnt( uInt detid, float val, int cnt ){ 
 
    if( endRun == lastRun ) return; 
    //std::cout << "Filling " << detid << " with " << val << " " << cnt << std::endl; 
    updated = true; 
    if( sumCntMap.find(detid) != sumCntMap.end() ){  
        sumCntMap[detid].sum += val; 
        sumCntMap[detid].sumsqr += val*val;  
        sumCntMap[detid].cnt += cnt;  
    } else { 
        sumCntMap[detid] = { val, val*val, cnt };  
    }//<<>>if( sumCntMap.find(detid) != sumCntMap.end() ) 
 
    if( detIdHists.find(detid) != detIdHists.end() ){ detIdHists[detid]->Fill(val); } //std::cout << " -- hist Fill " << detid << " with " << val << std::endl; } 
    else { 
        //std::cout << " -- hist Make " << detid << " with " << val << std::endl; 
        std::string histname = histMapName + "_SumCntHist_" + std::to_string( startRun ) + "_" + std::to_string(detid) ; 
        detIdHists[detid] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",200,-5,5); 
        detIdHists[detid]->Sumw2(); 
        detIdHists[detid]->Fill(val); 
    }//<<>>if( detIdHists.find(detid) != detIdHists.end() ) 

}//<<>>void CaliRunClass::fillSumCnt( uInt detID, float sum, int cnt ) 

#endif
