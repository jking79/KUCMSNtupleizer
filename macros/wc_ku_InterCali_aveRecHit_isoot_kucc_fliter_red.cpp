// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TText.h"
#include "TChain.h"

// STL includes
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <sys/stat.h>

#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"

using namespace std;

//----------------------------------------------------------------------------------------------------

void wc_ku_InterCali_aveRecHit_mini( string indir, string infilelistname, string outfilename ){

    //const int  nAlgos = 1; // Mini, MfootCCStc
    //const double offset = 0.0;
    const int bin_offset = 86;
	const float minRhEnergy = 5.0;

	//const bool debug = true;
    const bool debug = false;
	const bool useEnergy = true;

    const string treename("tree/llpgtree");

    std::cout << "Opening Outfile : " << outfilename << std::endl;
    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    TH2F * IcMapEB;
    TH2F * IcMapErrEB;
    TH2F * IcMapOccEB;
    TH1F * IcDistEB;
    TH1F * IcDistErrEB;
    TH2F * IcMapEP;
    TH2F * IcMapEM;
    TH2F * IcDistMeanEBTTMap;
    TH2F * IcDistMeanEBTTMapOcc;
    TH2F * IcDistMeanEBMDMap;
    TH2F * IcDistMeanEBMDMapOcc;
    TH2F * IcDistMeanEBSMMap;
    TH2F * IcDistMeanEBSMMapOcc;

    TH1F * IcDistMeanEBEta;
    TH1F * IcDistMeanErrEBEta;
    TH1F * IcDistMeanEBPhi;
    TH1F * IcDistMeanErrEBPhi;
    TH1F * IcMapETEta;
    TH1F * IcMapETPhi;
    TH1F * IcMapETErrEta;
    TH1F * IcMapETErrPhi;
    TH1F * IcDistMeanEBTT;
    TH1F * IcDistMeanErrEBTT;
    TH1F * IcDistMeanEBMD;
    TH1F * IcDistMeanErrEBMD;
    TH1F * IcDistMeanEBSM;
    TH1F * IcDistMeanErrEBSM;

    //string algostring[1] = { "" };
    //for( auto i = 0; i < nAlgos; i++){
        string hnameEB( "AveXtalRecTimeEBMap");
        string htitleEB( "AveXtal RecTimeEBMap EB ");
        IcMapEB = new TH2F(hnameEB.c_str(),htitleEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
        IcMapEB->Sumw2();
        string hnameEP( "AveXtalRecTimeEPMap");
        string htitleEP( "AveXtal RecTimeEPMap EP ");
        IcMapEP = new TH2F(hnameEP.c_str(),htitleEP.c_str(),100,0.5,100.5,100,0.5,100.5);
        IcMapEP->Sumw2();
        string hnameEM( "AveXtalRecTimeEMMap");
        string htitleEM( "AveXtal RecTimeEBMap EM ");
        IcMapEM = new TH2F(hnameEM.c_str(),htitleEM.c_str(),100,0.5,100.5,100,0.5,100.5);
        IcMapEM->Sumw2();
        string hnameOccEB( "AveXtalOccEBMap");
        string htitleOccEB( "AveXtal OccEBMap EB ");
        IcMapOccEB = new TH2F(hnameOccEB.c_str(),htitleOccEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
        IcMapOccEB->Sumw2();
        string hnameErrEB( "AveXtalRecTimeErrEBMap");
        string htitleErrEB( "AveXtal RecTimeErrEBMap EB ");
        IcMapErrEB = new TH2F(hnameErrEB.c_str(),htitleErrEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
        IcMapErrEB->Sumw2();
        string hnameDistEB( "AveXtalEBdist");
        string htitleDistEB( "AveXtal EBDist EB ");
        IcDistEB = new TH1F(hnameDistEB.c_str(),htitleDistEB.c_str(),320,-4,4);
        IcDistEB->Sumw2();
        string hnameErrDistEB( "AveXtalEBErrdist");
        string htitleErrDistEB( "AveXtal EBErrDist EB ");
        IcDistErrEB = new TH1F(hnameErrDistEB.c_str(),htitleErrDistEB.c_str(),1000,0,1.0);
        IcDistErrEB->Sumw2();
    //}//<<>>for( auto i = 0; i < nAlgos; i++)

    IcDistMeanEBEta = new TH1F("AveEtaRecTimeEta","AveEtaRecTimeEta;iEta;MeanTime [ns]",171,-85.5,85.5);
    IcDistMeanErrEBEta = new TH1F("AvePhiRecTimeEtaDist","AvePhiRecTimeEtaDist;MeanTime [ns]",160,-1,1);
    IcDistMeanEBEta->Sumw2();
    IcDistMeanErrEBEta->Sumw2();

    IcDistMeanEBPhi = new TH1F("AveEtaRecTimePhi","AveEtaRecTimePhi;iPhi;MeanTime [ns]",360,0.5,360.5);
    IcDistMeanErrEBPhi = new TH1F("AveEtaRecTimePhiDist","AveEtaRecTimePhiDist;MeanTime [ns]",160,-1,1);
    IcDistMeanEBPhi->Sumw2();
    IcDistMeanErrEBPhi->Sumw2();

    IcMapETEta = new TH1F("AveEtaRecMMTimeEta","AveEtaRecMeanMeanTimeEta;iEta;Mean MeanTime [ns]",171,-85.5,85.5);
    IcMapETPhi = new TH1F("AvePhiRecMMTimePhi","AvePhiRecMeanMeanTimePhi;iPhi;Mean MeanTime [ns]",360,0.5,360.5);
    IcMapETPhi->Sumw2();
    IcMapETEta->Sumw2();

    IcMapETErrEta = new TH1F("AveEtaRecMMTimeEtaErr","AveEtaRecMMTimeEtaErrt;Mean MeanTime [ns]",160,-1,1);
    IcMapETErrPhi = new TH1F("AvePhiRecMMTimePhiErr","AvePhiRecMMTimePhiErr;Mean MeanTime [ns]",160,-1,1);
    IcMapETErrEta->Sumw2();
    IcMapETErrPhi->Sumw2();

    IcDistMeanEBTTMap = new TH2F("AveTTRecTimeMap","AveTTRecTimeMap",34,0,34,72,0,72);
    IcDistMeanEBTT = new TH1F("AveTTRecTimeDist","AveTTRecTimeDist;MeanTime [ns]",160,-4,4);
    IcDistMeanErrEBTT = new TH1F("AveTTRecTimeErr","AveTTRecTimeErr;MeanTime Error [ns]",200,0,0.2);
    IcDistMeanEBTTMapOcc = new TH2F("AveTTRecTimeOccMap","AveTTRecTimeOccMap",34,0,34,72,0,72);
    IcDistMeanEBTTMap->Sumw2();
    IcDistMeanEBTT->Sumw2();
    IcDistMeanErrEBTT->Sumw2();
    IcDistMeanEBTTMapOcc->Sumw2();

    IcDistMeanEBMDMap = new TH2F("AveMDRecTimeMap","AveMDRecTimeMap",8,0,8,18,0,18);
    IcDistMeanEBMD = new TH1F("AveMDRecTimeDist","AveMDRecTimeDist;MeanTime [ns]",120,-3,3);
    IcDistMeanErrEBMD = new TH1F("AveMDRecTimeErr","AveMDRecTimeErr;MeanTime Error [ns]",400,0,0.04);
    IcDistMeanEBMDMapOcc = new TH2F("AveMDRecTimeOccMap","AveMDRecTimeOccMap",8,0,8,18,0,18);
    IcDistMeanEBMDMap->Sumw2();
    IcDistMeanEBMD->Sumw2();
    IcDistMeanErrEBMD->Sumw2();
    IcDistMeanEBMDMapOcc->Sumw2();

    IcDistMeanEBSMMap = new TH2F("AveSMRecTimeMap","AveSMRecTimeMap",2,0,2,18,0,18);
    IcDistMeanEBSM = new TH1F("AveSMRecTimeDist","AveSMRecTimeDist;MeanTime [ns]",80,-2,2);
    IcDistMeanErrEBSM = new TH1F("AveSMRecTimeErr","AveSMRecTimeErr;MeanTime Error [ns]",200,0,0.02);
    IcDistMeanEBSMMapOcc = new TH2F("AveSMRecTimeOccMap","AveSMRecTimeOccMap",2,0,2,18,0,18);
    IcDistMeanEBSMMap->Sumw2();
    IcDistMeanEBSM->Sumw2();
    IcDistMeanErrEBSM->Sumw2();
    IcDistMeanEBSMMapOcc->Sumw2();

    int tdiv = 240;
    float tstart = -6.0;
    float tend = 6.0;
    string locsname("SRO_");
    string locdname("DRO_");
    string globname("ZEE_");
    string histnames("Data_SHist");
    auto lochists = locsname+histnames;
    auto glohists = globname+histnames;
    auto theSHistLS = new TH2F(lochists.c_str(),lochists.c_str(), 240, 0, 1200, tdiv, tstart, tend);
    auto theSHistGS = new TH2F(glohists.c_str(),glohists.c_str(), 240, 0, 1200, tdiv, tstart, tend);

	auto dist0TOF = new TH1F("dist0tof","cms 000 TOF",300,0,15);
    auto distPvTOF = new TH1F("distpvtof","pv TOF",300,0,15);
    auto distRhTime = new TH1F("distrhtime","rh Time",2000,-50,50);
    auto distRhTimeTof = new TH1F("distrhtimetof","rh |Time|<25 + TOF",1000,-25,25);

    std::cout << "Setting up DetIDs." << std::endl;
    std::map<UInt_t,DetIDStruct> DetIDMap;
    SetupDetIDsEB( DetIDMap );
    SetupDetIDsEE( DetIDMap );
	//for( auto iter : DetIDMap ){ std::cout << iter.first << " " << iter.second.i1 << " " << iter.second.i2 << std::endl; }

    std::map<UInt_t,Float_t> sumXtalMiniRecTime;
    std::map<int,Float_t> sumXtalEtaRecTime;
    std::map<int,Float_t> sumXtalPhiRecTime;
    std::map<UInt_t,Float_t> sumXtal2MiniRecTime;
    std::map<int,Float_t> sumXtal2EtaRecTime;
    std::map<int,Float_t> sumXtal2PhiRecTime;
    std::map<UInt_t,UInt_t> numXtalMiniRecTime;
	std::map<int,UInt_t> numXtalEtaRecTime;
    std::map<int,UInt_t> numXtalPhiRecTime;

	std::map<UInt_t,Float_t> xtalMeanRtRecTime;

    std::map<int,Float_t> sumXtalEtaRecMMTime;
    std::map<int,Float_t> sumXtalPhiRecMMTime;
    std::map<int,Float_t> sumXtal2EtaRecMMTime;
    std::map<int,Float_t> sumXtal2PhiRecMMTime;
    std::map<int,UInt_t> numXtalEtaRecMMTime;
    std::map<int,UInt_t> numXtalPhiRecMMTime;

    std::map<int,Float_t> sumXtalTTRecTime;
    std::map<int,Float_t> sumXtal2TTRecTime;
    std::map<int,UInt_t> numXtalTTRecTime;
    std::map<int,Float_t> sumXtalMDRecTime;
    std::map<int,Float_t> sumXtal2MDRecTime;
    std::map<int,UInt_t> numXtalMDRecTime;
    std::map<int,Float_t> sumXtalSMRecTime;
    std::map<int,Float_t> sumXtal2SMRecTime;
    std::map<int,UInt_t> numXtalSMRecTime;


	// Declaration of leaf types for kuntuple
    UInt_t          Evt_run;
    UInt_t          Evt_event;
    vector<vector<unsigned int> > *Photon_gloResRhId;
    vector<vector<unsigned int> > *Photon_locResRhId;
    vector<float>   *ECALRecHit_energy;
    vector<unsigned int> *ECALRecHit_ID;
    vector<float>   *ECALRecHit_0TOF;
    vector<float>   *ECALRecHit_pvTOF;
    vector<float>   *ECALRecHit_time;
    vector<float>   *ECALRecHit_amplitude;
    vector<float>   *ECALRecHit_ampres;

    // List of branches for kuntuple
    TBranch        *b_Evt_run;   //!
    TBranch        *b_Evt_event;   //!
    TBranch        *b_Photon_gloResRhId;   //!
    TBranch        *b_Photon_locResRhId;   //!
    TBranch        *b_ECALRecHit_energy;   //!
    TBranch        *b_ECALRecHit_ID;   //!
    TBranch        *b_ECALRecHit_0TOF;   //!
    TBranch        *b_ECALRecHit_pvTOF;   //!
    TBranch        *b_ECALRecHit_time;   //!
    TBranch        *b_ECALRecHit_amplitude;   //!
    TBranch        *b_ECALRecHit_ampres;   //!

    std::ifstream infilelist(infilelistname);
    std::string infilestr;
    while (std::getline(infilelist,infilestr)){

    	std::stringstream ss(infilestr);
        std::string infilename;
        std::string srunstr;
        std::string erunstr;
        ss >> infilename >> srunstr >> erunstr;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "For Run " << srunstr << " to Run " << erunstr << std::endl;
        auto srun = std::stoi(srunstr);
        auto erun = std::stoi(erunstr);

    	std::ifstream infile(infilename);
    	std::string instr;
        auto fInTree = new TChain(treename.c_str());
        std::cout << "Adding files to TChain." << std::endl;
        const std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");
		//const std::string eosdir("/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/test/");
		//const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
		//const std::string eosdir("root://cmseos.fnal.gov//store/user/");	
		//int fcnt = 0;
        while (std::getline(infile,instr)){
			//if( fcnt > 10 ) continue;
			auto tfilename = eosdir + indir + instr;
         	//auto tfilename = indir + "/" + str;
         	//std::cout << "--  adding file: " << tfilename << std::endl;
			std::cout << "-";
         	fInTree->Add(tfilename.c_str());
			//fcnt++;
        }//<<>>while (std::getline(infile,str))
		std::cout << std::endl;

   		Photon_gloResRhId = 0;
   		Photon_locResRhId = 0;
   		ECALRecHit_energy = 0;
   		ECALRecHit_ID = 0;
   		ECALRecHit_0TOF = 0;
   		ECALRecHit_pvTOF = 0;
   		ECALRecHit_time = 0;
   		ECALRecHit_amplitude = 0;
   		ECALRecHit_ampres = 0;
	
   		fInTree->SetBranchAddress("Photon_gloResRhId", &Photon_gloResRhId, &b_Photon_gloResRhId);
   		fInTree->SetBranchAddress("Photon_locResRhId", &Photon_locResRhId, &b_Photon_locResRhId);
   		fInTree->SetBranchAddress("ECALRecHit_energy", &ECALRecHit_energy, &b_ECALRecHit_energy);
   		fInTree->SetBranchAddress("ECALRecHit_ID", &ECALRecHit_ID, &b_ECALRecHit_ID);
   		fInTree->SetBranchAddress("ECALRecHit_0TOF", &ECALRecHit_0TOF, &b_ECALRecHit_0TOF);
   		fInTree->SetBranchAddress("ECALRecHit_pvTOF", &ECALRecHit_pvTOF, &b_ECALRecHit_pvTOF);
   		fInTree->SetBranchAddress("ECALRecHit_time", &ECALRecHit_time, &b_ECALRecHit_time);
   		fInTree->SetBranchAddress("ECALRecHit_amplitude", &ECALRecHit_amplitude, &b_ECALRecHit_amplitude);
   		fInTree->SetBranchAddress("ECALRecHit_ampres", &ECALRecHit_ampres, &b_ECALRecHit_ampres);
   		fInTree->SetBranchAddress("Evt_run", &Evt_run, &b_Evt_run);
   		fInTree->SetBranchAddress("Evt_event", &Evt_event, &b_Evt_event);

         // >> calcs  <<
     
        auto nEntries = fInTree->GetEntries();
        std::cout << "Starting entry loops for " << nEntries << " Entries" <<  std::endl;
		//nEntries = 1000;
        if( debug ) nEntries = 100;
        for (Long64_t centry = 0; centry < nEntries; centry++){
			  
     		if( centry%100000 == 0 or centry == 0){ 
				std::cout << "Proccessed " << centry << " of " << nEntries; 
                std::cout << " (" << static_cast<float>((10000*centry)/nEntries)/(100) << "%)" << std::endl;
			}//<<>>if( centry%1000000 == 0 or centry == 0)

			auto entry = fInTree->LoadTree(centry);

            b_Photon_gloResRhId->GetEntry(entry);   //!
            b_Photon_locResRhId->GetEntry(entry);   //!

            b_ECALRecHit_energy->GetEntry(entry);   //!
            b_ECALRecHit_ID->GetEntry(entry);   //!
            b_ECALRecHit_0TOF->GetEntry(entry);   //!
            b_ECALRecHit_pvTOF->GetEntry(entry);   //!
            b_ECALRecHit_time->GetEntry(entry);   //!
            b_ECALRecHit_amplitude->GetEntry(entry);   //!
            b_ECALRecHit_ampres->GetEntry(entry);   //!
            b_Evt_run->GetEntry(entry);   //!
            b_Evt_event->GetEntry(entry);   //!

			if( Evt_run < srun || Evt_run > erun ) continue;

			const auto locResRhIds = (*Photon_locResRhId)[0];
            const auto nLocResRhIds = locResRhIds.size();
            const auto nRecHits = ECALRecHit_ID->size(); //(cluster[ipho0])->size();
			//std::cout << " -- " << nLocResRhIds << " - " << nRecHits << std::endl;
			for( int idx = 0; idx < nLocResRhIds; idx += 2 ){

				int rhidx1 = -1;
                int rhidx2 = -1;
				uInt locRhId1 = locResRhIds[idx];
                uInt locRhId2 = locResRhIds[idx+1];
                const auto & idinfo1 = DetIDMap[locRhId1];
                const auto & idinfo2 = DetIDMap[locRhId2];
				for( int i = 0; i < nRecHits; i++ ){ 
					if( (*ECALRecHit_ID)[i] == locRhId1 ) rhidx1 = i;
                    if( (*ECALRecHit_ID)[i] == locRhId2 ) rhidx2 = i;
				}//<<>>for( int i = 0; i < nRecHits; i++ )

                if( rhidx1 < 0 || rhidx2 < 0 ) continue;
	            if( idinfo1.ecal != ECAL::EB || (*ECALRecHit_time)[rhidx1] == 0.0 ) continue;
                if( idinfo2.ecal != ECAL::EB || (*ECALRecHit_time)[rhidx2] == 0.0 ) continue;
				if( (*ECALRecHit_energy)[rhidx1] < 0.5 || (*ECALRecHit_energy)[rhidx2] < 0.5 ) continue;				

				float loct1 = (*ECALRecHit_time)[rhidx1] + (*ECALRecHit_0TOF)[rhidx1] - (*ECALRecHit_pvTOF)[rhidx1];
				float loct2 = (*ECALRecHit_time)[rhidx2] + (*ECALRecHit_0TOF)[rhidx2] - (*ECALRecHit_pvTOF)[rhidx2];
				float dt = loct1 - loct2;
				float effampnum = (*ECALRecHit_ampres)[rhidx1] * (*ECALRecHit_ampres)[rhidx2];
		        float effampden = std::sqrt(sq2((*ECALRecHit_ampres)[rhidx1])+sq2((*ECALRecHit_ampres)[rhidx2]));
				float effamp = effampnum/effampden;
				theSHistLS->Fill(effamp,dt);

			}//<<>>for( auto idx = 0; idx < nLocResRhIds; idx += 2 )

            const auto gloResRhIds = (*Photon_gloResRhId)[0];
            const auto nGloResRhIds = gloResRhIds.size();
            //const auto nRecHits = ECALRecHit_ID->size(); //(cluster[ipho0])->size();
            //std::cout << " -- " << nLocResRhIds << " - " << nRecHits << std::endl;
            for( int idx = 0; idx < nGloResRhIds; idx += 2 ){

                int rhidx1 = -1;
                int rhidx2 = -1;
                uInt gloRhId1 = gloResRhIds[idx];
                uInt gloRhId2 = gloResRhIds[idx+1];
				const auto & idinfo1 = DetIDMap[gloRhId1];
				const auto & idinfo2 = DetIDMap[gloRhId2];
                for( int i = 0; i < nRecHits; i++ ){
                    if( (*ECALRecHit_ID)[i] == gloRhId1 ) rhidx1 = i;
                    if( (*ECALRecHit_ID)[i] == gloRhId2 ) rhidx2 = i;
                }//<<>>for( int i = 0; i < nRecHits; i++ )

                if( rhidx1 < 0 || rhidx2 < 0 ) continue;
                if( idinfo1.ecal != ECAL::EB || (*ECALRecHit_time)[rhidx1] == 0.0 ) continue;
                if( idinfo2.ecal != ECAL::EB || (*ECALRecHit_time)[rhidx2] == 0.0 ) continue;
                if( (*ECALRecHit_energy)[rhidx1] < 0.5 || (*ECALRecHit_energy)[rhidx2] < 0.5 ) continue;

                float glot1 = (*ECALRecHit_time)[rhidx1] + (*ECALRecHit_0TOF)[rhidx1] - (*ECALRecHit_pvTOF)[rhidx1];
                float glot2 = (*ECALRecHit_time)[rhidx2] + (*ECALRecHit_0TOF)[rhidx2] - (*ECALRecHit_pvTOF)[rhidx2];
                float dt = glot1 - glot2;
                float effampnum = (*ECALRecHit_ampres)[rhidx1] * (*ECALRecHit_ampres)[rhidx2];
                float effampden = std::sqrt(sq2((*ECALRecHit_ampres)[rhidx1])+sq2((*ECALRecHit_ampres)[rhidx2]));
                float effamp = effampnum/effampden;
                theSHistGS->Fill(effamp,dt);

            }//<<>>for( auto idx = 0; idx < nLocResRhIds; idx += 2 )

            //const auto nRecHits1 = ECALRecHit_ID->size(); //(cluster[ipho0])->size();
            if( debug ) std::cout << "Looping over first recHits"  << std::endl;
            for ( auto rh_i = 0U; rh_i < nRecHits; rh_i++ ){

				//auto rhe = useEnergy ? (*ECALRecHit_energy)[rh_i] : 999;	
				//if( rhe < minRhEnergy ) continue;
                bool underMinEnergy = (*ECALRecHit_energy)[rh_i] < 5.0;
                bool rhTimeZero = (*ECALRecHit_time )[rh_i] == 0.0;
				bool timeOutOfRange = std::abs((*ECALRecHit_time )[rh_i]) > 25.0;
				if( underMinEnergy or rhTimeZero ) continue;
                auto id_i = (*ECALRecHit_ID)[rh_i];
                const auto & fill_idinfo = DetIDMap[id_i];
                auto iEta = fill_idinfo.i2;
                auto iPhi = fill_idinfo.i1;
				bool isEB = (fill_idinfo.ecal == ECAL::EB);
                auto Mini_t_i = (*ECALRecHit_time )[rh_i] + (*ECALRecHit_0TOF)[rh_i] - (*ECALRecHit_pvTOF)[rh_i];
                if( isEB ) dist0TOF->Fill((*ECALRecHit_0TOF)[rh_i]);
                if( isEB ) distPvTOF->Fill((*ECALRecHit_pvTOF)[rh_i]);
                if( isEB ) distRhTime->Fill((*ECALRecHit_time)[rh_i]);
                if( timeOutOfRange ) continue;
                if( isEB ) distRhTimeTof->Fill(Mini_t_i);
     	     	if( debug ) std::cout << "Getting maps " << std::endl;
                sumXtalMiniRecTime[id_i] += Mini_t_i; 
				numXtalMiniRecTime[id_i] += 1;
                sumXtal2MiniRecTime[id_i] += Mini_t_i*Mini_t_i;
				int ttphi = (iPhi-1)/5;
                int tteta = (std::abs(iEta)-1)/5;
				int ttidx = ( iEta < 0 ) ? (ttphi+(tteta+1)*100) * -1 : (ttphi+tteta*100);
                int mdphi = (iPhi-1)/20;
				int mdeta = ( std::abs(iEta) < 26 ) ? 0 : 1 + (std::abs(iEta)-26)/20;  
                int mdidx = ( iEta < 0 ) ? (mdphi+(mdeta+1)*100) * -1 : (mdphi+mdeta*100);
                int smphi = (iPhi-1)/20;
                int smidx = ( iEta < 0 ) ? (smphi+1) * -1 : smphi;
				//std::cout << " TT idx: " << ttidx << std::endl;
                if( isEB ){

                	sumXtalEtaRecTime[iEta] += Mini_t_i;
                  	numXtalEtaRecTime[iEta] += 1;
                  	sumXtal2EtaRecTime[iEta] += Mini_t_i*Mini_t_i;

                  	sumXtalPhiRecTime[iPhi] += Mini_t_i;
                  	numXtalPhiRecTime[iPhi] += 1;
                  	sumXtal2PhiRecTime[iPhi] += Mini_t_i*Mini_t_i;

					sumXtalTTRecTime[ttidx] += Mini_t_i;
					sumXtal2TTRecTime[ttidx] += Mini_t_i*Mini_t_i;
                    numXtalTTRecTime[ttidx] += 1;

                    sumXtalMDRecTime[mdidx] += Mini_t_i;
                    sumXtal2MDRecTime[mdidx] += Mini_t_i*Mini_t_i;
                    numXtalMDRecTime[mdidx] += 1;

                    sumXtalSMRecTime[smidx] += Mini_t_i;
                    sumXtal2SMRecTime[smidx] += Mini_t_i*Mini_t_i;
                    numXtalSMRecTime[smidx] += 1;

                }//<<>>if( fill_idinfo.ecal == ECAL::EB )

             }//<<>>for (auto i = 0U; i < nRecHits1; i++) // end loop over rechits
             if( debug ) std::cout << "RecHits Loop done "<< std::endl;
         }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)
		delete fInTree;

    } // while (std::getline(infilelist,infiles))

	std::cout << "Filling Calibration Maps" <<  std::endl;
    //std::map<UInt_t,Float_t> *  icmaps[nAlgos] = {&sumXtalMiniRecTime};
    //std::map<UInt_t,UInt_t> *  nicmaps[nAlgos] = {&numXtalMiniRecTime};
    //std::map<UInt_t,Float_t> *  ic2maps[nAlgos] = {&sumXtal2MiniRecTime};
	//std::map<UInt_t,Float_t> *  meanMaps[nAlgos] = {&xtalMeanRtRecTime};
    //for( auto ai = 0; ai < nAlgos; ai++ ){
         for( std::map<UInt_t,Float_t>::iterator it=sumXtalMiniRecTime.begin(); it!=sumXtalMiniRecTime.end(); ++it){
            const auto & fill_idinfo = DetIDMap[it->first];
            const auto & map_time = ((sumXtalMiniRecTime[it->first])/(numXtalMiniRecTime[it->first]));// + offset; 
				// - (drift/(icmaps[ai]->size()))) + offset;
            const auto & map_occ = numXtalMiniRecTime[it->first];
            const auto & map_err = sqrt(((sumXtal2MiniRecTime[it->first])/map_occ - map_time*map_time)/map_occ);
			if( debug ) std::cout << "Fill hist for Algo at " << fill_idinfo.i2; 
			if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << " for iter " << std::endl;
            if( fill_idinfo.ecal == ECAL::EB ){
		   		if( debug ) std::cout << "Fill EB hist for Algo at " << fill_idinfo.i2 << " "; 
                if( debug ) std::cout << fill_idinfo.i1 << " with " << map_time << std::endl;
            	IcMapEB->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
                IcMapOccEB->Fill( fill_idinfo.i2, fill_idinfo.i1, map_occ );
                IcMapErrEB->Fill( fill_idinfo.i2, fill_idinfo.i1, map_err );
                IcDistEB->Fill(map_time);
                IcDistErrEB->Fill(map_err);
				xtalMeanRtRecTime[it->first] = map_time;
                if( map_time != 0.0 ){
                    sumXtalEtaRecMMTime[fill_idinfo.i2] += map_time;
                    sumXtalPhiRecMMTime[fill_idinfo.i1] += map_time;
                    sumXtal2EtaRecMMTime[fill_idinfo.i2] += map_time*map_time;
                    sumXtal2PhiRecMMTime[fill_idinfo.i1] += map_time*map_time;
                    numXtalEtaRecMMTime[fill_idinfo.i2] += 1;
                    numXtalPhiRecMMTime[fill_idinfo.i1] += 1;
                }//<<>>if( map_time != 0.0 )
            } else if( fill_idinfo.ecal == ECAL::EP ){
                if( debug ) std::cout << "Fill EP hist for Algo  at " << fill_idinfo.i2; 
                if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                IcMapEP->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            } else if( fill_idinfo.ecal == ECAL::EM ){
                if( debug ) std::cout << "Fill EM hist for Algo  at " << fill_idinfo.i2;
                if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                IcMapEM->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            }//<<>>if( fill_idinfo.ecal == ECAL::EB )
         }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)
    //}//<<>>for( auto ai = 0; ai < nAlgos; ai++ )

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecTime.begin(); it!=sumXtalEtaRecTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecTime[iEta]/numXtalEtaRecTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcDistMeanEBEta->SetBinContent(iEta+86,map_time);
           IcDistMeanEBEta->SetBinError(iEta+86,map_err);
           IcDistMeanErrEBEta->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalPhiRecTime.begin(); it!=sumXtalPhiRecTime.end(); ++it){
       const auto iPhi = it->first;
       const auto & map_time = sumXtalPhiRecTime[iPhi]/numXtalPhiRecTime[iPhi]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalPhiRecTime[iPhi];
       const auto & map_err = sqrt((sumXtal2PhiRecTime[iPhi]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcDistMeanEBPhi->SetBinContent(iPhi,map_time);
           IcDistMeanEBPhi->SetBinError(iPhi,map_err);
           IcDistMeanErrEBPhi->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecMMTime.begin(); it!=sumXtalEtaRecMMTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecMMTime[iEta]/numXtalEtaRecMMTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecMMTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecMMTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           //IcMapETEta->Fill( iEta, map_time );
           IcMapETErrEta->Fill( map_time );
           IcMapETEta->SetBinContent(iEta+86,map_time);
           IcMapETEta->SetBinError(iEta+86,map_err);
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalPhiRecMMTime.begin(); it!=sumXtalPhiRecMMTime.end(); ++it){
       const auto iPhi = it->first;
       const auto & map_time = sumXtalPhiRecMMTime[iPhi]/numXtalPhiRecMMTime[iPhi]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalPhiRecMMTime[iPhi];
       const auto & map_err = sqrt((sumXtal2PhiRecMMTime[iPhi]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcMapETErrPhi->Fill( map_time );
           IcMapETPhi->SetBinContent(iPhi,map_time);
           IcMapETPhi->SetBinError(iPhi,map_err);
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    std::cout << "TT : map fills " << std::endl;
    for( std::map<int,Float_t>::iterator it=sumXtalTTRecTime.begin(); it!=sumXtalTTRecTime.end(); ++it){
       const int iTT = it->first;
	   const int eta = iTT/100;
	   const int iPhiTT = std::abs(iTT) - std::abs(eta)*100;
       const int iEtaTT = eta + 17;
	   //std::cout << "TT : iTT = " << iTT << " phi: " << iPhiTT << " eta: " << iEtaTT << std::endl;
       const float & map_time = sumXtalTTRecTime[iTT]/numXtalTTRecTime[iTT]; // - (drift/(icmaps[ai]->size()))) + offset;
       const int & map_occ = numXtalTTRecTime[iTT];
       const float & map_err = sqrt((sumXtal2TTRecTime[iTT]/map_occ - map_time*map_time)/map_occ);
	   //std::cout << "TT vals : time: " << map_time << " occ: " << map_occ << " err: " << map_err << std::endl;
       //if( fill_idinfo.ecal == ECAL::EB ){
           (IcDistMeanEBTTMap)->Fill(iEtaTT,iPhiTT,map_time);
           (IcDistMeanEBTTMapOcc)->Fill(iEtaTT,iPhiTT,map_occ);
           IcDistMeanEBTT->Fill( map_time );
           IcDistMeanErrEBTT->Fill( map_err );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    std::cout << "MD : map fills " << std::endl;
    for( std::map<int,Float_t>::iterator it=sumXtalMDRecTime.begin(); it!=sumXtalMDRecTime.end(); ++it){
       const int iMD = it->first;
       const int eta = iMD/100;
       const int iPhiMD = std::abs(iMD) - std::abs(eta)*100;
       const int iEtaMD = eta + 4;
       //std::cout << "MD : iMD = " << iMD << " phi: " << iPhiMD << " eta: " << iEtaMD << std::endl;
       const float & map_time = sumXtalMDRecTime[iMD]/numXtalMDRecTime[iMD]; // - (drift/(icmaps[ai]->size()))) + offset;
       const int & map_occ = numXtalMDRecTime[iMD];
       const float & map_err = sqrt((sumXtal2MDRecTime[iMD]/map_occ - map_time*map_time)/map_occ);
       //std::cout << "MD vals : time: " << map_time << " occ: " << map_occ << " err: " << map_err << std::endl;
       //if( fill_idinfo.ecal == ECAL::EB ){
           (IcDistMeanEBMDMap)->Fill(iEtaMD,iPhiMD,map_time);
           (IcDistMeanEBMDMapOcc)->Fill(iEtaMD,iPhiMD,map_occ);
           IcDistMeanEBMD->Fill( map_time );
           IcDistMeanErrEBMD->Fill( map_err );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    std::cout << "SM : map fills " << std::endl;
    for( std::map<int,Float_t>::iterator it=sumXtalSMRecTime.begin(); it!=sumXtalSMRecTime.end(); ++it){
       const int iSM = it->first;
       //const int eta = iSM/100;
       const int iPhiSM = ( iSM < 0 ) ? std::abs(iSM) - 1 : std::abs(iSM);
       const int iEtaSM = ( iSM < 0 ) ? 0 : 1;
       //std::cout << "SM : iSM = " << iSM << " phi: " << iPhiSM << " eta: " << iEtaSM << std::endl;
       const float & map_time = sumXtalSMRecTime[iSM]/numXtalSMRecTime[iSM]; // - (drift/(icmaps[ai]->size()))) + offset;
       const int & map_occ = numXtalSMRecTime[iSM];
       const float & map_err = sqrt((sumXtal2SMRecTime[iSM]/map_occ - map_time*map_time)/map_occ);
       //std::cout << "SM vals : time: " << map_time << " occ: " << map_occ << " err: " << map_err << std::endl;
       //if( fill_idinfo.ecal == ECAL::EB ){
           (IcDistMeanEBSMMap)->Fill(iEtaSM,iPhiSM,map_time);
           (IcDistMeanEBSMMapOcc)->Fill(iEtaSM,iPhiSM,map_occ);
           IcDistMeanEBSM->Fill( map_time );
           IcDistMeanErrEBSM->Fill( map_err );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    //for( auto ai = 0; ai < nAlgos; ai++ ){ (*icmaps[ai]).clear(); (*nicmaps[ai]).clear(); (*meanMaps[ai]).clear(); }

    fOutFile->cd();

    std::cout << "Write AveXtal Rechit Time Maps" << std::endl;

    theSHistLS->Write();
    theSHistGS->Write();
	delete theSHistLS;
    delete theSHistGS;

	dist0TOF->Write();
	distPvTOF->Write();
	distRhTime->Write();
    distRhTimeTof->Write();    
	delete dist0TOF;
    delete distPvTOF;
    delete distRhTime;
    delete distRhTimeTof;

    IcMapETErrPhi->Write();
    IcMapETPhi->Write();
    IcMapETErrEta->Write();
    IcMapETEta->Write();
    IcDistMeanErrEBPhi->Write();
    IcDistMeanEBPhi->Write();
    IcDistMeanErrEBEta->Write();
    IcDistMeanEBEta->Write();

	IcDistMeanEBTTMap->Write();
	IcDistMeanEBTTMapOcc->Write();
	IcDistMeanEBTT->Write();
	IcDistMeanErrEBTT->Write();

    IcDistMeanEBMDMap->Write();
    IcDistMeanEBMDMapOcc->Write();
    IcDistMeanEBMD->Write();
    IcDistMeanErrEBMD->Write();

    IcDistMeanEBSMMap->Write();
    IcDistMeanEBSMMapOcc->Write();
    IcDistMeanEBSM->Write();
    IcDistMeanErrEBSM->Write();

    delete IcMapETErrPhi;
    delete IcMapETPhi;
    delete IcMapETErrEta;
    delete IcMapETEta;
    delete IcDistMeanErrEBPhi;
    delete IcDistMeanEBPhi;
    delete IcDistMeanErrEBEta;
    delete IcDistMeanEBEta;

    delete IcDistMeanEBTTMap;
    delete IcDistMeanEBTTMapOcc;
    delete IcDistMeanEBTT;
    delete IcDistMeanErrEBTT;

    delete IcDistMeanEBMDMap;
    delete IcDistMeanEBMDMapOcc;
    delete IcDistMeanEBMD;
    delete IcDistMeanErrEBMD;

    delete IcDistMeanEBSMMap;
    delete IcDistMeanEBSMMapOcc;
    delete IcDistMeanEBSM;
    delete IcDistMeanErrEBSM;

    //for( auto i = 0; i < nAlgos; i++){

         IcMapEB->Write();
         IcMapErrEB->Write();
         IcMapOccEB->Write();
         IcDistEB->Write();
         IcDistErrEB->Write();
         IcMapEP->Write();
         IcMapEM->Write();

         delete IcMapEB;
         delete IcMapErrEB;
         delete IcMapOccEB;
         delete IcDistEB;
         delete IcDistErrEB;
         delete IcMapEP;
         delete IcMapEM;

    //}//<<>>for( auto i = 0; i < nAlgos; i++)

    delete fOutFile;

}//<<>>void wc_ku_InterCali_aveRecHit_mini( string indir, string infilelistname, string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {

        //auto indir = "";
        //auto indir = "jaking/ecalTiming/gammares_tt_kucc_126_v3/EGamma/";
		//auto indir = "ecalTiming/gammares_tt_kucc_126_v7_diag/EGamma/";
		//auto indir = "ecalTiming/gammares_ttcc_140_v11_diag/EGamma1/";
        //auto indir = "ecalTiming/gammares_ttcc_140_v11_diag_ebsf_ccgt/EGamma1/";
        //auto indir = "ecalTiming/gammares_mc/ZprimeToEE_M-6000_TuneCP5_13p6TeV_pythia8/";
        auto indir = "KUCMSNtuple/kucmsntuple_MET_R17E_MET75_v20/";
        //auto indir = "KUCMSNtuple/kucmsntuple_GMSB_R17_MET75_v20/";
        //auto indir = "KUCMSNtuple/kucmsntuple_GJETS_R17_MET75_v20/";
        //auto indir = "KUCMSNtuple/kucmsntuple_QCD_R17_MET75_v20/"; 

		//auto infilelist = "test_califilelist.txt";
        //auto infilelist = "tt_run3_2023B_Prompt_366323_367065_1301_gammares_v11_reso_califilelist.txt";
        //auto infilelist = "egamma1_run3_2024B_Prompt_378794_379190_1404_gammares_v11_reso_califilelist.txt";
        //auto infilelist = "list_files/egamma1_run3_2024C_Prompt_379415-380238_1404_gammares_v11_reso_califilelist.txt";
        //auto infilelist = "list_files/egamma1_run3_2024D_Prompt_380066-380947_1404_gammares_v11_reso_califilelist.txt";
		//auto infilelist = "list_files/egamma1_run3_2023D_Prompt_369844_369999_1404_gammares_v11_reso_califilelist.txt";
        //auto infilelist = "list_files/egamma1_run3_2023D_Prompt_370293-370580_1404_ccgt_gammares_v11_reso_califilelist.txt";
		//auto infilelist = "list_files/egamma_mc_ZprimeToEE_M-6000_gammares_v11_reso_califilelist.txt";
        auto infilelist = "cali_list_files/KUCMS_MET_Met75_R17_v20_califilelist.txt";
        //auto infilelist = "cali_list_files/KUCMS_GMSB_Met75_R17_v20_L300_califilelist.txt";
        //auto infilelist = "cali_list_files/KUCMS_GJets_Met75_R17_v20_califilelist.txt";
        //auto infilelist = "cali_list_files/KUCMS_QCD_Met75_R17_v20_califilelist.txt";

        auto outfilename = "kuntuple_MET_R17_MET75_v20_v8_Cali.root";
        //auto outfilename = "kuntuple_GJets_R17_MET75_v20_Cali.root";
        //auto outfilename = "kuntuple_QCD_R17_MET75_v20_v7_Cali.root";

        wc_ku_InterCali_aveRecHit_mini( indir, infilelist, outfilename );

    //}//<<>>if( argc != 4 )
    return 1;

}//<<>>int main ( int argc, char *argv[] )

