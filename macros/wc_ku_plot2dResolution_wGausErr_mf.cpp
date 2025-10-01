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
#include <TRandom.h>

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
#include "wc_ku_timefitter_wErr_func.cpp"

using namespace std;

void plot2dResolution( std::string indir, std::string infilelistname, std::string outfilename, std::string tvarname, std::string calimapname, std::string isd_type, bool useAmp, std::string xbinstr ){

	bool debug = false;
	//bool debug = true;

	TRandom* getRandom = new TRandom();
	getRandom->SetSeed(0);
	unsigned int nRand(100);

    std::cout << "opening output file" << std::endl;
    string histoutfilename(outfilename+".root");
    TFile* fOutFile = new TFile( histoutfilename.c_str(), "RECREATE" );
    std::cout << "fOutFile : " << fOutFile << std::endl;

    float phoseedtimeCaliIc_0(0.0);
    float phoseedtimeCaliIcErr_0(0.0);
    float phoseedtimeCaliIc_1(0.0);
    float phoseedtimeCaliIcErr_1(0.0);

    // Declaration of leaf types
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

    // List of branches
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

    // >> calcs  <<
    std::cout << "Setting up 2D plot" << std::endl;

	string locsname("SRO_");
    string locdname("DRO_");
    string globname("ZEE_");

    string histname("Data_Hist");
    string histnameg("Data_GHist");
    string histnames("Data_SHist");
    string histname0("Data_OthHist");
    string histname1("td_Hist");
    string histname2("tdc_Hist");
    string histname3("teffa_Hist");
    string histname4("tdfc_Hist");
    string histname5("td_eta_Hist");
    string histname6("td_phi_Hist");
    string fTitle("#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)");
    string fTitle0("#Delta(Photon Seed Time) [ns] vs. A_{0}/#sigma_{n} (EBEB)");
    string fTitle1("Photon Seed Time [ns]");
    string fTitle2("Photon Seed Time Calibrated [ns]");
    string fTitle3("Amplitude 0 vs Amplitude 1");
    string fTitle4("Photon Seed Time Filtered Calibrated [ns]");
    string fTitle5("#Delta(Photon Seed Time) [ns] vs iEta");
    string fTitle6("#Delta(Photon Seed Time) [ns] vs iPhi");
    string fXTitle("A_{eff}/#sigma_{n} (EBEB)");
    std::vector<Double_t> fXBins;
    Bool_t fXVarBins = false;//dummy not used
    string fYTitle("#Delta(Photon Seed Time) [ns] (EBEB)");
    std::vector<Double_t> fYBins;
    Bool_t fYVarBins = false;//dummy not used;
    //string ybinstr("CONSTANT 1920 -3 3");
    //string ybinstr("CONSTANT 600 -3 3");
    string ybinstr("CONSTANT 240 -3 3");
    string fZTitle("");

    setBins(xbinstr,fXBins,fXVarBins);
    setBins(ybinstr,fYBins,fYVarBins);

    const auto xbins = &fXBins[0];
    const auto ybins = &fYBins[0];
    int nMyBins = fXBins.size()-1;

    int tdiv = 240;
    float tstart = -3.0;
    float tend = 3.0;

	int ediv = 720;
	float eend = 3600;

    std::string tehistname = "timeErrHist";
	auto timeErrBinHist = new TH1F(tehistname.c_str(), tehistname.c_str(),1000,0,1);
    TH1F * caliErrBinHist[nMyBins+1];
    std::string eahistname = "effAHist";
    caliErrBinHist[0] = new TH1F(eahistname.c_str(),eahistname.c_str(),ediv,0,eend);
    caliErrBinHist[0]->Sumw2();
    for( auto ibin = 1; ibin < nMyBins+1; ibin++ ){
		auto binhistname = "bin" + to_string(ibin) + "ErrHist";
		caliErrBinHist[ibin] = new TH1F(binhistname.c_str(),binhistname.c_str(),500,0,0.5);
        caliErrBinHist[ibin]->Sumw2();
    }//<<>>for( auto ibin = 1; ibin < nMyBins+1; ibin++ )

	// local same
	auto lochist = locsname+histname;
    auto theHistLS = new TH2F(lochist.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto lochistg = locsname+histnameg;
    auto theGHistLS = new TH2F(lochistg.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto lochists = locsname+histnames;
    auto theSHistLS = new TH2F(lochists.c_str(),fTitle.c_str(), ediv, 0, eend, tdiv, tstart, tend);
    auto lochist0 = locsname+histname0;
    auto theOthHistLS = new TH2F(lochist0.c_str(),fTitle0.c_str(), ediv, 0, eend, tdiv, tstart, tend);
    auto lochist1 = locsname+histname1;
    auto thetdHistLS = new TH1F(lochist1.c_str(),fTitle1.c_str(),tdiv,tstart,tend);
    auto lochist2 = locsname+histname2;
    auto thetdcHistLS = new TH1F(lochist2.c_str(),fTitle2.c_str(),tdiv,tstart,tend);

    //auto thetdfcHist = new TH1F(histname4.c_str(),fTitle4.c_str(),tdiv,tstart,tend);
    auto lochist3 = locsname+histname3;
    auto theEffaHistLS = new TH2F(lochist3.c_str(),fTitle3.c_str(),ediv,0,eend,ediv,0,eend);
    auto lochist5 = locsname+histname5;
    auto theetaHistLS = new TH2F(lochist5.c_str(),fTitle5.c_str(),182,-90.5,90.5,tdiv,tstart,tend);
    auto lochist6 = locsname+histname6;
    auto thephiHistLS = new TH2F(lochist6.c_str(),fTitle6.c_str(),363,-1.5,361.5,tdiv,tstart,tend);

    theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
    theGHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    theGHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    theGHistLS->GetZaxis()->SetTitle(fZTitle.c_str());

	// local diffrnent
    auto locdhist = locdname+histname;
    auto theHistLD = new TH2F(locdhist.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());

	// global
    auto globhist = globname+histname;
    auto theHistGB = new TH2F(globhist.c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());

    std::cout << "Setting up DetIDs." << std::endl;
    std::map<UInt_t,DetIDStruct> DetIDMap;
    SetupDetIDsEB( DetIDMap );
    SetupDetIDsEE( DetIDMap );

    float goodlev(0);
    float goodlin(0);
	float goodgev(0);
    float goodgin(0);
    float gevents(0);

    std::cout << "open input files list : " << infilelistname << std::endl;
    string disphotreename("tree/llpgtree");

    std::ifstream infilelist(infilelistname);
    std::string infiles;
    while (std::getline(infilelist,infiles)){

        std::stringstream ss(infiles);
        std::string infilename;
        std::string califilename;
        std::string srunstr;
        std::string erunstr;
        std::string subdir;
        ss >> infilename >> califilename >> srunstr >> erunstr >> subdir;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "open input cali : " << califilename << std::endl;
		auto srun = std::stoi(srunstr);
        auto erun = std::stoi(erunstr);
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Subdir: " << subdir << std::endl;

        std::ifstream infile( infilename );
        std::string str;

        auto fInTree = new TChain(disphotreename.c_str());
        std::cout << "Adding files to TChain." << std::endl;
        while (std::getline(infile,str)){
        	const std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");
            //const std::string eosdir("root://cmseos.fnal.gov//store/user/");
            //const std::string eosdir("/uscms/home/jaking/nobackup/el8/ecaltiming/CMSSW_14_0_4/src/GammaResTool/GammaResTool/");
            auto tfilename = eosdir + indir + subdir + str;
            //auto tfilename = indir + "/" + str;
            //std::cout << "--  adding file: " << tfilename << std::endl;
			std::cout << "-";
        	fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
		std::cout << std::endl;

		//auto calidir = "/home/t3-ku/jaking/ecaltiming/skimmed_trees/local_chain/";
		string calidir = "cali_root_files/";
		TFile* fCaliFile(NULL);
		if( calimapname != "none" ){
			//string calidir = "";
			string califilein = calidir+califilename;
			std::cout << "Using caliFile : " << califilein << std::endl;
        	if( califilename != "none" ) fCaliFile = TFile::Open( califilein.c_str(), "read" );
		}//<<>>if( calimapname != "none" ){
        std::cout << "fInFile : " << infilename  << " fInTree : " << fInTree << " fCaliFile : " << fCaliFile << std::endl;

        std::cout << "set branches to get from fInFile : fInTree" << std::endl;

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

		vector<TH2F*> calimaps{0,0,0,0,0,0,0,0};
        std::cout << "get maps from fCaliFile" << std::endl;
        if( fCaliFile ){
			fCaliFile->cd();
        	string cmbsrt("AveXtalRecTime");
        	string itcnt("");
        	string ebmaprt(cmbsrt+"EBMap"+itcnt);
        	string eberrmaprt(cmbsrt+"ErrEBMap"+itcnt);
        	string epmaprt(cmbsrt+"EPMap"+itcnt);
        	string emmaprt(cmbsrt+"EMMap"+itcnt);
        	calimaps[0] = (TH2F*)fCaliFile->Get(ebmaprt.c_str());
        	calimaps[1] = (TH2F*)fCaliFile->Get(eberrmaprt.c_str());
        	calimaps[2] = (TH2F*)fCaliFile->Get(epmaprt.c_str());
        	calimaps[3] = (TH2F*)fCaliFile->Get(emmaprt.c_str());
		}//<<>>if( calimapname != "none" )

        std::cout << "Getting calibration values and plotting" << std::endl;

        auto nEntries = fInTree->GetEntries();
        //auto nEntries = 10000;
		int loopcnt = 250000;
        //int loopcnt = 1000000;
		if( debug ) nEntries = ( nEntries < 100000 ) ? nEntries : 100000;
		if( debug ) std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
        for (auto centry = 0U; centry < nEntries; centry++){

			if( centry%loopcnt == 0 or centry == 0){ 
				std::cout << "Proccessed " << centry << " of " << nEntries;
				std::cout << " " << (1000.0*static_cast<float>(centry)/static_cast<float>(nEntries))/10.0 << "%" << std::endl;
			}//<<>>if( centry%10000000 == 0 or centry == 0)

			auto entry = fInTree->LoadTree(centry);

			if(debug) std::cout << " - Start loop " << std::endl;

			gevents++;

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

			//int didx = 0;
			//for( int didx = debug?0:4; didx < 4; didx++ ){
			//	if(debug) std::cout << "Run " << run << " id " << (*resRhID)[didx] << " Amp " << (*resAmp)[didx] << " E " << (*resE)[didx];
            //	if(debug) std::cout << " Rt " << (*resRtTime)[didx]  << " TOF " << (*resTOF)[didx] << std::endl;
			//}//<<>>for( int didx = 0; didx < 4; didx++ ){
			if( srun != 0 && ( Evt_run < srun || Evt_run > erun ) ) continue;

			if(debug) std::cout << " - Finshed Get Entry " << std::endl;


            /////////////////////  loop over reshits
			const auto locResRhIds = (*Photon_locResRhId)[0];
            const auto nLocResRhIds = locResRhIds.size();
            const auto nRecHits = ECALRecHit_ID->size(); //(cluster[ipho0])->size();
            //std::cout << " -- " << nLocResRhIds << " - " << nRecHits << std::endl;
            const auto gloResRhIds = (*Photon_gloResRhId)[0];
            const auto nGloResRhIds = gloResRhIds.size();
			const auto nResRhIds = ( nLocResRhIds > nGloResRhIds ) ? nLocResRhIds : nGloResRhIds;
            for( int idx = 0; idx < nResRhIds; idx += 2 ){

				const bool useLoc = idx < nLocResRhIds;
	            const bool useGlo = idx < nGloResRhIds;

				Int_t i1L0 = -99;
                Int_t i1L1 = -99;
                Int_t i2L0 = -99;
                Int_t i2L1 = -99;
				Int_t L0EB = -99;
                Int_t L1EB = -99;
                Int_t i1G0 = -99;
                Int_t i1G1 = -99;
                Int_t i2G0 = -99;
                Int_t i2G1 = -99;
                Int_t G0EB = -99;
                Int_t G1EB = -99;

                int lrhidx1 = -1;
                int lrhidx2 = -1;
                int grhidx1 = -1;
                int grhidx2 = -1;

                uInt locRhId1 = ( useLoc ) ? locResRhIds[idx] : 0;
                uInt locRhId2 = ( useLoc ) ? locResRhIds[idx+1] : 0;
                uInt gloRhId1 = ( useGlo ) ? gloResRhIds[idx] : 0;
                uInt gloRhId2 = ( useGlo ) ? gloResRhIds[idx+1] : 0;

                for( int i = 0; i < nRecHits; i++ ){
                    if( (*ECALRecHit_ID)[i] == locRhId1 ) lrhidx1 = i;
                    if( (*ECALRecHit_ID)[i] == locRhId2 ) lrhidx2 = i;
                    if( (*ECALRecHit_ID)[i] == gloRhId1 ) grhidx1 = i;
                    if( (*ECALRecHit_ID)[i] == gloRhId2 ) grhidx2 = i;
                }//<<>>for( int i = 0; i < nRecHits; i++ )	

                bool isd_cut = true; //idinfoL0.TT == idinfoL1.TT; // true = same, false = different
				if( useLoc ){

					auto idinfoL0 = DetIDMap[(*ECALRecHit_ID)[lrhidx1]];
	            	auto idinfoL1 = DetIDMap[(*ECALRecHit_ID)[lrhidx2]];
	
	            	i1L0 = idinfoL0.i1;
	            	i1L1 = idinfoL1.i1;
	            	i2L0 = idinfoL0.i2;     
	            	i2L1 = idinfoL1.i2;
	
					L0EB = idinfoL0.ecal;
	            	L1EB = idinfoL1.ecal;

					isd_cut = (idinfoL0.TT == idinfoL1.TT); // true = same, false = different

				}//<<>>if( useLoc ) 	

                if( useGlo ){

	            	auto idinfoG0 = DetIDMap[(*ECALRecHit_ID)[grhidx1]];
	            	auto idinfoG1 = DetIDMap[(*ECALRecHit_ID)[grhidx2]];
	
	            	i1G0 = idinfoG0.i1;
	            	i1G1 = idinfoG1.i1;
	            	i2G0 = idinfoG0.i2;
	            	i2G1 = idinfoG1.i2;
	
	            	G0EB = idinfoG0.ecal;
	            	G1EB = idinfoG1.ecal;

				}//<<>>if( useGlo )	

	            if(debug) std::cout << " - Finshed Setting IDInfo " << std::endl;
	
		        int bin_offset = 86;
		        int adjust = 0.0;
	
				vector<float> seedTimeIC0{0,0};
	            vector<float> seedTimeICE0{0,0};
	            vector<float> seedTimeIC1{0,0};
	            vector<float> seedTimeICE1{0,0};
	
		        if( fCaliFile ){ //and califilename != "none" ){
	
					if ( useLoc ){
	            		if ( L0EB == ECAL::EB ){
	            			seedTimeIC0[0] = calimaps[0]->GetBinContent( i2L0 + bin_offset, i1L0 ) - adjust;
							seedTimeICE0[0] = calimaps[1]->GetBinContent( i2L0 + bin_offset, i1L0 );
	            		}else if ( L0EB == ECAL::EP ){
	            	    	seedTimeIC0[0] = calimaps[2]->GetBinContent( i2L0, i1L0 ) - adjust;
	            		}else if (L0EB == ECAL::EM ){
	            	    	seedTimeIC0[0] = calimaps[3]->GetBinContent( i2L0, i1L0 ) - adjust;
	            		}//<<>>if ( L0EB == ECAL::EB )
					}//<<>>if ( (*resRhID)[0] )
	            
	                if ( useLoc ){
	        			if ( L1EB == ECAL::EB ){
	            			seedTimeIC1[0] = calimaps[0]->GetBinContent( i2L1 + bin_offset, i1L1 ) - adjust;
	                		seedTimeICE1[0] = calimaps[1]->GetBinContent( i2L1 + bin_offset, i1L1 );
	            		}else if ( L1EB == ECAL::EP ){
	                		seedTimeIC1[0] = calimaps[2]->GetBinContent( i2L1, i1L1 ) - adjust;
	            		}else if (L1EB == ECAL::EM ){
	                		seedTimeIC1[0] = calimaps[3]->GetBinContent( i2L1, i1L1 ) - adjust;
						}//<<>>if ( L1EB == ECAL::EB )
	                }//<<>>if ( (*resRhID)[1] )
	
	                if ( useGlo ){
	                	if ( G0EB == ECAL::EB ){
	                    	seedTimeIC0[1] = calimaps[0]->GetBinContent( i2G0 + bin_offset, i1G0 ) - adjust;
	                    	seedTimeICE0[1] = calimaps[1]->GetBinContent( i2G0 + bin_offset, i1G0 );
	                	}else if ( G0EB == ECAL::EP ){
	                    	seedTimeIC0[1] = calimaps[2]->GetBinContent( i2G0, i1G0 ) - adjust;
	                	}else if (G0EB == ECAL::EM ){
	                    	seedTimeIC0[1] = calimaps[3]->GetBinContent( i2G0, i1G0 ) - adjust;
	                	}//<<>>if ( G0EB == ECAL::EB )
	                }//<<>>if ( (*resRhID)[2] )
	
	                if ( useGlo ){
	                	if ( G1EB == ECAL::EB ){
	                    	seedTimeIC1[1] = calimaps[0]->GetBinContent( i2G1 + bin_offset, i1G1 ) - adjust;
	                    	seedTimeICE1[1] = calimaps[1]->GetBinContent( i2G1 + bin_offset, i1G1 );
	                	}else if ( G1EB == ECAL::EP ){
	                    	seedTimeIC1[1] = calimaps[2]->GetBinContent( i2G1, i1G1 ) - adjust;
	                	}else if (G1EB == ECAL::EM ){
	                    	seedTimeIC1[1] = calimaps[3]->GetBinContent( i2G1, i1G1 ) - adjust;
	                	}//<<>>if ( G1EB == ECAL::EB )
	                }//<<>>if ( (*resRhID)[3] )
	
		        }//<<>>if( calimapname != "none" )
	
	
	//-------------------set for local, repo calcs for global --------------------------------
	
	            if(debug) std::cout << " - Calc 2D Hist" << std::endl;
	
				float lyfill = -99;
				float restime0 = -99;
                float restime1 = -99;
                float crestime0 = -99;
                float crestime1 = -99;
				if( useLoc ){
	            	//double lyfill = ((*resRtTime)[0]-seedTimeIC0[0])-((*resRtTime)[1]-seedTimeIC1[0])+(*resTOF)[0]-(*resTOF)[1];
                	float loct1 = (*ECALRecHit_time)[lrhidx1] + (*ECALRecHit_0TOF)[lrhidx1] - (*ECALRecHit_pvTOF)[lrhidx1] - seedTimeIC0[0];
                	float loct2 = (*ECALRecHit_time)[lrhidx2] + (*ECALRecHit_0TOF)[lrhidx2] - (*ECALRecHit_pvTOF)[lrhidx2] - seedTimeIC1[0];
                	lyfill = loct1 - loct2;
					restime0 = (*ECALRecHit_time)[lrhidx1];
                    restime1 = (*ECALRecHit_time)[lrhidx2];
                    crestime0 = (*ECALRecHit_time)[lrhidx1] - seedTimeIC0[0];
                    crestime1 = (*ECALRecHit_time)[lrhidx2] - seedTimeIC1[0];
				}//<<>>if( useLoc )

                float gyfill = -99;
                if( useGlo ){
	            	//double gyfill = ((*resRtTime)[2]-seedTimeIC0[1])-((*resRtTime)[3]-seedTimeIC1[1])+(*resTOF)[2]-(*resTOF)[3];
                	float glot1 = (*ECALRecHit_time)[grhidx1] + (*ECALRecHit_0TOF)[grhidx1] - (*ECALRecHit_pvTOF)[grhidx1] - seedTimeIC0[1];
                	float glot2 = (*ECALRecHit_time)[grhidx2] + (*ECALRecHit_0TOF)[grhidx2] - (*ECALRecHit_pvTOF)[grhidx2] - seedTimeIC1[1];
                	gyfill = glot1 - glot2;
                }//<<>>if( useGlo )

                float lxfill = -99;
                float effenergy = 0;
				float leffa0 = 0;
                float leffa1 = 0;
                if( useLoc ){
	            	float effa0 = useAmp ? (*ECALRecHit_ampres)[lrhidx1] : (*ECALRecHit_energy)[lrhidx1]; 
	            	float effa1 = useAmp ? (*ECALRecHit_ampres)[lrhidx2] : (*ECALRecHit_energy)[lrhidx2]; 
                	float effampnum = effa0 * effa1;
                	float effampden = std::sqrt(sq2(effa0)+sq2(effa1));
                	lxfill = effampnum/effampden;
					leffa0 = (*ECALRecHit_ampres)[lrhidx1];
					leffa1 = (*ECALRecHit_ampres)[lrhidx2];
					float effenergyn = (*ECALRecHit_energy)[lrhidx1] * (*ECALRecHit_energy)[lrhidx2];
					float effenergyd = std::sqrt(sq2((*ECALRecHit_energy)[lrhidx1])+sq2((*ECALRecHit_energy)[lrhidx2]));
					effenergy = effenergyn/effenergyd;
                }//<<>>if( useLoc )

                float gxfill = -99;
                if( useGlo ){
                    float effa0 = useAmp ? (*ECALRecHit_ampres)[grhidx1] : (*ECALRecHit_energy)[grhidx1];
                    float effa1 = useAmp ? (*ECALRecHit_ampres)[grhidx2] : (*ECALRecHit_energy)[grhidx2];
                    float effampnum = effa0 * effa1;
                    float effampden = std::sqrt(sq2(effa0)+sq2(effa1));
                    gxfill = effampnum/effampden;
                }//<<>>if( useGlo )

	            float timeErr0 = 0; //phoseedtimeErr_0/25.0;
	            float timeErr1 = 0; //phoseedtimeErr_1/25.0;
				float ldtserr0 = timeErr0*timeErr0+seedTimeICE0[0]*seedTimeICE0[0];
	            float ldtserr1 = timeErr1*timeErr1+seedTimeICE1[0]*seedTimeICE1[0];
				float ldterr = sqrt(ldtserr0+ldtserr1);
	            float gdtserr0 = timeErr0*timeErr0+seedTimeICE0[1]*seedTimeICE0[1];
	            float gdtserr1 = timeErr1*timeErr1+seedTimeICE1[1]*seedTimeICE1[1];
	            float gdterr = sqrt(gdtserr0+gdtserr1);

				float resE0 = useLoc ? (*ECALRecHit_energy)[lrhidx1] : -1;
                float resE1 = useLoc ? (*ECALRecHit_energy)[lrhidx2] : -1;
                float resE2 = useGlo ? (*ECALRecHit_energy)[grhidx1] : -1;
                float resE3 = useGlo ? (*ECALRecHit_energy)[grhidx2] : -1;	
			    bool goodLocEnergy = (resE0>=10)&&(resE0<=120)&&(resE1>=10)&&(resE1<=120);
	            bool goodGloEnergy = (resE2>=10)&&(resE2<=120)&&(resE3>=10)&&(resE3<=120);

		        bool isLocEB = (L0EB == ECAL::EB)&&(L1EB == ECAL::EB);
	            bool isGloEB = (G0EB == ECAL::EB)&&(G1EB == ECAL::EB);

                bool goodLocRHs = lrhidx1 >= 0 && lrhidx2 >= 0;
                bool goodGloRHs = grhidx1 >= 0 && grhidx2 >= 0;
	            bool goodLocTime = useLoc ? (*ECALRecHit_time)[lrhidx1] != 0 && (*ECALRecHit_time)[lrhidx2] != 0 : false; 
	            bool goodGloTime = useGlo ? (*ECALRecHit_time)[grhidx1] != 0 && (*ECALRecHit_time)[grhidx2] != 0 : false; 
	
				//auto isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, fasle = different
		        bool levent_good = goodLocEnergy && isLocEB && goodLocRHs && goodLocTime;
	            bool gevent_good = goodGloEnergy && isGloEB && goodGloRHs && goodGloTime;
	
	            if( levent_good ) goodlev++;
	            if( goodLocRHs ) goodlin++;
	            if( gevent_good ) goodgev++;
	            if( goodGloRHs ) goodgin++;
	
				if(debug) std::cout << " - lxfill : " << lxfill << " lyfill : " << lyfill;
				if(debug) std::cout  << " flag : " << goodLocEnergy << " " << isLocEB << " " << goodLocRHs << std::endl; 
	            if(debug) std::cout << " - gxfill : " << gxfill << " gyfill : " << gyfill;
	            if(debug) std::cout  << " flag : " << goodGloEnergy << " " << isGloEB << " " << goodGloRHs << std::endl;
	
	//-----------------   set up fills for local same diffrent and global -----------------------
	
				if(debug) std::cout << " - Fill 2D Hist" << std::endl;
		        if( levent_good && isd_cut ){ // isd_cut : true = same, fasle = different
					theHistLS->Fill(lxfill,lyfill);
					theSHistLS->Fill(lxfill,lyfill);
					theOthHistLS->Fill(effenergy,lyfill);
					//for( auto i = 0; i < nRand; i++ ){ 
					//	//virtual Double_t	Gaus(Double_t mean = 0, Double_t sigma = 1)
					//	auto gsfill = getRandom->Gaus(lyfill,ldterr);
	                //  	theGHistLS->Fill(lxfill,gsfill);
					//}//<<>>for( int i = 0; i < 100; i++ )	
					if(debug) std::cout << " - Fill effA dist Hist" << std::endl;
					//string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250");
                    if( false ){//turn following section off
					caliErrBinHist[0]->Fill(lxfill);
					timeErrBinHist->Fill(timeErr0); timeErrBinHist->Fill(timeErr1);
	               	if(debug) std::cout << " - Fill bin Hists" << std::endl;
					if( lxfill < 75 ) caliErrBinHist[1]->Fill(ldterr);
					else if( lxfill < 100 ) caliErrBinHist[2]->Fill(ldterr);
	               	else if( lxfill < 125 ) caliErrBinHist[3]->Fill(ldterr);
	               	else if( lxfill < 150 ) caliErrBinHist[4]->Fill(ldterr);
	               	else if( lxfill < 175 ) caliErrBinHist[5]->Fill(ldterr);
	               	else if( lxfill < 225 ) caliErrBinHist[6]->Fill(ldterr);
	               	else if( lxfill < 275 ) caliErrBinHist[7]->Fill(ldterr);
	               	else if( lxfill < 325 ) caliErrBinHist[8]->Fill(ldterr);
	               	else if( lxfill < 375 ) caliErrBinHist[9]->Fill(ldterr);
	               	else if( lxfill < 475 ) caliErrBinHist[10]->Fill(ldterr);
	               	else if( lxfill < 600 ) caliErrBinHist[11]->Fill(ldterr);
	               	else if( lxfill < 950 ) caliErrBinHist[12]->Fill(ldterr);
	               	else if( lxfill < 2250 ) caliErrBinHist[13]->Fill(ldterr);
					else std::cout << "Over 2250" << std::endl;
					}//turn above section off
	            }//<<>>if( levent_good )
	            if(debug) std::cout << " - Fill 1D Hist" << std::endl;
	            if( levent_good ){
					theEffaHistLS->Fill(leffa0,leffa1); 
					thetdHistLS->Fill(restime0); thetdHistLS->Fill(restime1);
	                thetdcHistLS->Fill(crestime0); thetdcHistLS->Fill(crestime1);
					theetaHistLS->Fill(i2L0,lyfill); theetaHistLS->Fill(i2L1,lyfill);
	                thephiHistLS->Fill(i1L0,lyfill); thephiHistLS->Fill(i1L1,lyfill);
				}//<<>>if( levent_good )
				if( levent_good && not isd_cut ){ theHistLD->Fill(lxfill,lyfill); }
	            if( gevent_good ){ theHistGB->Fill(gxfill,gyfill); }
	
				if(debug) std::cout << " - Fill hists done" << std::endl;
	
	        }//<<>>for( int idx = 0; idx < nResRhIds; idx += 2 )
	        //////////////  loop over res

        } // for (auto entry = 0U; entry < nEntries; entry++)
	 	//delete fInFile;

        //delete fCaliFile;  <<<<<<<<<<<<<<   delete califiles ????????????????
		if(debug) std::cout << " -------- Next Input file " << std::endl;

    } // while (std::getline(infilelist,infiles))

// --------------  process new histos ----------------------------------------------

	if(debug) std::cout << " - Scale 2D Hist" << std::endl;
    scaleHist(theHistLS,false,fXVarBins,fYVarBins);
    scaleHist(theGHistLS,false,fXVarBins,fYVarBins);
    scaleHist(theHistLD,false,fXVarBins,fYVarBins);
    scaleHist(theHistGB,false,fXVarBins,fYVarBins);
    normTH2F(theetaHistLS);
    normTH2F(thephiHistLS);

    fOutFile->cd();
    theHistLS->Write();
    theGHistLS->Write();
    theSHistLS->Write();
    theOthHistLS->Write();
    theEffaHistLS->Write();
    thetdHistLS->Write();
    thetdcHistLS->Write();
    theetaHistLS->Write();
    thephiHistLS->Write();
    //timeErrBinHist->Write();
    //for( auto ibin = 0; ibin < nMyBins; ibin++ ){ caliErrBinHist[ibin]->Write(); }
	theHistLD->Write();
    theHistGB->Write();

    delete theHistLS;
    delete theGHistLS;
    delete theSHistLS;
    delete theOthHistLS;
    delete theEffaHistLS;
    delete thetdHistLS;
    delete thetdcHistLS;
    delete theetaHistLS;
    delete thephiHistLS;
    delete timeErrBinHist;
    for( auto ibin = 0; ibin < nMyBins; ibin++ ){ delete caliErrBinHist[ibin]; }
    delete theHistLD;
    delete theHistGB;

    delete fOutFile;
	delete getRandom;

	auto passlev = 100.0*goodlev/gevents;
    auto passlin = 100.0*goodlin/gevents;
    auto passgev = 100.0*goodgev/gevents;
    auto passgin = 100.0*goodgin/gevents;

    if(gevents) std::cout << "Processed " << gevents << " with %" << passlev << " from %" << passlin << " of Local"; 
	if(gevents) std::cout << " and %" << passgev << " from %" << passgin << " of Global" << std::endl;
    std::cout << "Thats all Folks!" << std::endl;

}//<<>> void plot2dResolution( std::string indir, std::string infilelistname, 

int main ( int argc, char *argv[] ){

    	//std::string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250 9000"); 
   		//std::string xbinstr("VARIABLE 0 75 100 125 150 175 200 225 250 275 325 375 450 550 725 925 1325 1700 2250"); 
    	//std::string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 750 950 1275 1700 2250");
        //std::string xbinstr("VARIABLE 60 80 100 120 150 180 240 320 400 480 560 640 800 1600 3200");
        std::string xbinstr("VARIABLE 60 120 180 240 300 420 510 600 780 960 1260 1740 2220 3200");
		bool useAmp(true);
        //std::string xbinstr("VARIABLE 0 10 12 14 16 18 20 22 24 26 30 34 40 48 62 78 120");
        //bool useAmp(false);

        //if( argc != 7 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {

            //auto indir = "jaking/ecalTiming/";
            //auto indir = "jaking/ecalTiming/gammares_ttcc_140_v11_diag_ebsf_ccgt/";
            //auto indir = "jaking/ecalTiming/gammares_mc/";
        	//auto indir = "KUCMSNtuple/kucmsntuple_MET_R17E_MET75_v20/";
            auto indir = "KUCMSNtuple/kucmsntuple_GJETS_R17_MET75_v20/";
            //auto indir = "KUCMSNtuple/kucmsntuple_QCD_R17_MET75_v20/";

			//auto infilelistname = "list_files/ku_24D_eg01_cali_140_gammares_v12_plotfilelist.txt";
            //auto infilelistname = "list_files/ku_24CD_eg01_mixed_cali_140_gammares_v12_plotfilelist.txt";
            //auto infilelistname = "list_files/ku_24CE_test_140_gammares_v12_plotfilelist.txt";
            //auto infilelistname = "list_files/ku_23D_eg1_sfeb_ccgt_140_gammares_v12_plotfilelist.txt";
            //auto infilelistname = "list_files/ku_W24_dy_eb_st_14011_gammares_v12_plotfilelist.txt";
            //auto infilelistname = "list_files/KUCMS_MET_Met75_R17_v20_plotfile.txt";
            auto infilelistname = "list_files/KUCMS_GJets_Met75_R17_v20_plotfilelist.txt";
            //auto infilelistname = "list_files/KUCMS_QCD_Met75_R17_v20_plotfilelist.txt";

            //std::string outfilename = "egres_Run2024C_goldjson_eg01_379367_379543_140_v10_resplots";
            //std::string outfilename = "egres_Run2024C_goldjson_eg01_380066_380238_140_v10_resplots";
            //std::string outfilename = "egres_Run2024D_nojson_eg01_Full_cali_140_v10_resplots";
            //std::string outfilename = "egres_Run2024CD_mixed_eg01_Full_cali_140_v10_resplots";
            //std::string outfilename = "egres_Run2023D_SF_140_v10_resplots";
            //std::string outfilename = "egres_Winter24_DY_14011_v10_resplots";
            //std::string outfilename = "kuntuple_egres_MET_Met75_R17_v2_v20_resplots";
            //std::string outfilename = "kuntuple_egres_MET_Met75_R17_ecut2_v20_Cali_resplots";
            std::string outfilename = "kuntuple_egres_GJets_Met75_R17_v3_v20_resplots";
            //std::string outfilename = "kuntuple_egres_QCD_Met75_R17_ecut2_v20_resplots";
            //std::string outfilename = "kuntuple_egres_QCD_Met75_R17_ecut2_v20_Cali_resplots";

            auto tvarname = ""; //argv[4];
			auto calimapname = "infile"; //argv[5];
            auto isd_type = "yes"; //argv[6];
            //auto brun = std::stoi(argv[7]);
            //auto erun = std::stoi(argv[8]);
            //auto leta = std::stoi(argv[9]);
            //auto heta = std::stoi(argv[10]);
      		plot2dResolution( indir, infilelistname, outfilename, tvarname, calimapname, isd_type, useAmp, xbinstr );
			std::string fitInFile = outfilename + ".root";
			runTimeFitter( fitInFile, "", "", "", "", outfilename, "SRO_Data_Hist", xbinstr );
            runTimeFitter( fitInFile, "", "", "", "", outfilename, "DRO_Data_Hist", xbinstr );
            runTimeFitter( fitInFile, "", "", "", "", outfilename, "ZEE_Data_Hist", xbinstr );

        //}//<<>>//if( argc != 7 )
        return 1;

}//<<>>int main ( int argc, char *argv[] )
