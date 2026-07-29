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
 
#ifndef KUCMS_CaliRunClass 
#define KUCMS_CaliRunClass 

struct kucms_sumCnt { 

    kucms_sumCnt(){} 
    kucms_sumCnt( double s, double s2, double c ) : sum(s), sumsqr(s2), cnt(c) {}  
 
    double sum; 
    double sumsqr; 
    double cnt; 
 
};//<<>>struct sumCnt 

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
	void makeMeanMap( bool filter = false ){
     
        meanMap.clear(); 
        errMap.clear(); 
        //std::cout << " make mean maps: " << std::endl; 
        for( auto& entry : sumCntMap ){ 
     
            double cnt = double( entry.second.cnt ); 
            if( cnt == 0 ) cnt = 1; 
            double mean = entry.second.sum / cnt;  
            double err = sqrt( (entry.second.sumsqr/cnt - mean*mean)/cnt ); 
            if( err == 0 ){ err = 9; mean = 0; } 
            //std::cout << " - calc: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl; 
            if( filter && ( cnt > 10 ) ){ 
				double high = mean + 3*err;
                double low = mean - 3*err;
                auto fitFunc  = new TF1("gfit","gaus",low,high); 
                auto thefit = detIdHists[entry.first]->Fit("gfit","QNRL"); 
                double fmean = fitFunc->GetParameter(1); 
                double ferr = fitFunc->GetParError(1); 
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
     
    //inline void kucms_CaliRunClass::fillSumCnt( uInt detid, float val, int cnt ){ 
	void fillSumCnt( uInt detid, double val, double cnt = 1 ){
    
		//std::cout << "Cecking " << detid << " with " << endRun << " " << lastRun <<  " " << hasGID2 << " " << doGID2  << std::endl; 
        if( endRun == lastRun ){ 
			if( hasGID2 == true ) return; 
			if( doGID2 == false ) return;
		}//<<>>if( ( endRun == lastRun )
        //std::cout << "Filling " << detid << " with " << val << " " << cnt << std::endl; 
        updated = true; 
        double wval = val*cnt;
        if( sumCntMap.find(detid) != sumCntMap.end() ){ 
            sumCntMap[detid].sum += wval; 
            sumCntMap[detid].sumsqr += wval*wval;  
            sumCntMap[detid].cnt += cnt;  
        } else { 
            sumCntMap[detid] = { wval, wval*wval, cnt };  
        }//<<>>if( sumCntMap.find(detid) != sumCntMap.end() ) 
    
    	//std::cout << " -- hist Fill " << detid << " with " << val << std::endl; } 
        if( detIdHists.find(detid) != detIdHists.end() ){ detIdHists[detid]->Fill(val,cnt); } 
        else { 
            //std::cout << " -- hist Make " << detid << " with " << val << std::endl; 
            std::string histname = histMapName + "_SumCntHist_" + std::to_string( startRun ) + "_" + std::to_string(detid) ; 
            detIdHists[detid] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",200,-5,5); 
            detIdHists[detid]->Sumw2(); 
            detIdHists[detid]->Fill(val,cnt); 
        }//<<>>if( detIdHists.find(detid) != detIdHists.end() ) 
    
    }//<<>>void CaliRunClass::fillSumCnt( uInt detID, float sum, int cnt ) 

};//<<>>TimeCaliTagStruct 

#endif
