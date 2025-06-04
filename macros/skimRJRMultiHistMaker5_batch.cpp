//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "skimHistMaker.hh"

//#define DEBUG true
#define DEBUG false
#define doEBEEmaps false

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// HistMaker class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle, 
		int nphos, int njets, std::vector<float> m_vec, std::vector<float> r_vec, std::vector<float> rv_vec ){

	//listdir, infilenameWJLNu, outfilenameswjlnu, htitleswjlnu, jrjtype, modtype, nphos, rjrcut, rjrcut2

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");
	const std::string ofnending = "_RjrSkim_v24_cfcatHists.root";

    cutselection = 1; // jrj type -> pho visible
    fl = true; // first event loop

	nRjrPhos = nphos;
	nRjrJets = njets;

	scale = 400;

    std::ifstream infile(infilelist);
    std::string str;
    while (std::getline(infile,str)){
        //std::cout << "--  adding from: " << str << std::endl;
        if( str[0] == '#' ) continue;
        if( str == " " ) continue;
        if( str == "" ) continue;
        auto instrs = splitString( str, " " );
        std::string inpath = instrs[0];
        std::string intitle = instrs[1];
		std::string batch = instrs[2];
        std::string tfilename = eosdir + indir + inpath;
        std::cout << "--  adding file: " << tfilename << std::endl;

    	std::cout << "Producing Histograms for : " << outfilename << std::endl;
        //std::ifstream infile(indir+infilelist);
        //std::ifstream infile(infilelist);
        auto fInTree = new TChain(disphotreename.c_str());
        auto fConfigTree = new TChain(configtreename.c_str());
        //std::cout << "Adding files to TChain." << std::endl;
        //std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
        //std::string str;
    	int cnt = 1;
		if( batch == "yes" ){
			std::ifstream inpathfile(inpath);
			std::string tempstr;
        	while( std::getline(inpathfile,tempstr) ){
    			//std::cout << "--  adding from: " << str << std::endl;
            	if( tempstr[0] == '#' ) continue;
            	if( tempstr == " " ) continue;
    			if( tempstr == "" ) continue;
            	auto batchtfilename = eosdir + indir + tempstr;
            	std::cout << "--  adding file: " << batchtfilename << std::endl;
            	fInTree->Add(batchtfilename.c_str());
            	fConfigTree->Add(batchtfilename.c_str());
            	cnt++;
			}//<<>>while( std::getline(inpathfile,tempstr) )
		} else { // <<>> if( batch == "yes" )
            fInTree->Add(tfilename.c_str());
            fConfigTree->Add(tfilename.c_str());
    		cnt++;
        }//<<>>if( batch == "yes" ) else
    
        std::cout << "Setting up For Main Loop." << std::endl;
    
    	Init(fInTree);
    	initHists(htitle+"_"+intitle);
    
        std::cout << "Filling Config Map." << std::endl;
    
        UInt_t          nEvents;
        UInt_t          nSelectedEvents;
        string          *sKey;
        Float_t         sCrossSection;
        Float_t         sGMSBGravMass;
        Float_t         sGMSBChi1Mass;
        Float_t         sMCWgt;
        Int_t           sMCType;
        Float_t         sumEvtGenWgt;
    
        Int_t           cf_m_gt2jets;
        Int_t           cf_met150;
        Int_t           cf_mj_gt1phos;
        Int_t           cf_mjp_leadPhoPt30;
        Int_t           cf_nFltrdEvts;
        Int_t           cf_nTotEvts;	
    	Int_t           cf_sel_j;
       	Int_t           cf_sel_m;
       	Int_t           cf_sel_p;
       	Int_t           cf_sel_ppt;
       	Int_t           nFltrdEvts;
       	Int_t           nTotEvts;
       	Float_t         sumEvtWgt;
       	Float_t         sumFltrdEvtWgt;
    
        TBranch        *b_nEvents;   //!
        TBranch        *b_nSelectedEvents;   //!
        TBranch        *b_sKey;   //!
        TBranch        *b_sCrossSection;   //!
        TBranch        *b_sGMSBGravMass;   //!
        TBranch        *b_sGMSBChi1Mass;   //!
        TBranch        *b_sMCWgt;   //!
        TBranch        *b_sMCType;   //!
        TBranch        *b_sumEvtGenWgt;
    
       	TBranch        *b_cf_m_gt2jets;   //!
       	TBranch        *b_cf_met150;   //!
       	TBranch        *b_cf_mj_gt1phos;   //!
       	TBranch        *b_cf_mjp_leadPhoPt30;   //!
       	TBranch        *b_cf_nFltrdEvts;   //!
       	TBranch        *b_cf_nTotEvts;   //!
       	TBranch        *b_cf_sel_j;   //!
       	TBranch        *b_cf_sel_m;   //!
       	TBranch        *b_cf_sel_p;   //!
       	TBranch        *b_cf_sel_ppt;   //!
       	TBranch        *b_nFltrdEvts;   //!
       	TBranch        *b_nTotEvts;   //!
       	TBranch        *b_sumEvtWgt;   //!
       	TBranch        *b_sumFltrdEvtWgt;   //!
    
    
        sKey = 0;
    
        fConfigTree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
        fConfigTree->SetBranchAddress("nSelectedEvents", &nSelectedEvents, &b_nSelectedEvents);
        fConfigTree->SetBranchAddress("sKey", &sKey, &b_sKey);
        fConfigTree->SetBranchAddress("sCrossSection", &sCrossSection, &b_sCrossSection);
        fConfigTree->SetBranchAddress("sGMSBGravMass", &sGMSBGravMass, &b_sGMSBGravMass);
        fConfigTree->SetBranchAddress("sGMSBChi1Mass", &sGMSBChi1Mass, &b_sGMSBChi1Mass);
        fConfigTree->SetBranchAddress("sMCWgt", &sMCWgt, &b_sMCWgt);
        fConfigTree->SetBranchAddress("sMCType", &sMCType, &b_sMCType);
        fConfigTree->SetBranchAddress("sumEvtGenWgt", &sumEvtGenWgt, &b_sumEvtGenWgt);
    
        fConfigTree->SetBranchAddress("cf_m_gt2jets", &cf_m_gt2jets, &b_cf_m_gt2jets);
        fConfigTree->SetBranchAddress("cf_met150", &cf_met150, &b_cf_met150);
        fConfigTree->SetBranchAddress("cf_mj_gt1phos", &cf_mj_gt1phos, &b_cf_mj_gt1phos);
        fConfigTree->SetBranchAddress("cf_mjp_leadPhoPt30", &cf_mjp_leadPhoPt30, &b_cf_mjp_leadPhoPt30);
        fConfigTree->SetBranchAddress("cf_nFltrdEvts", &cf_nFltrdEvts, &b_cf_nFltrdEvts);
        fConfigTree->SetBranchAddress("cf_nTotEvts", &cf_nTotEvts, &b_cf_nTotEvts);
        fConfigTree->SetBranchAddress("cf_sel_j", &cf_sel_j, &b_cf_sel_j);
        fConfigTree->SetBranchAddress("cf_sel_m", &cf_sel_m, &b_cf_sel_m);
        fConfigTree->SetBranchAddress("cf_sel_p", &cf_sel_p, &b_cf_sel_p);
        fConfigTree->SetBranchAddress("cf_sel_ppt", &cf_sel_ppt, &b_cf_sel_ppt);
        fConfigTree->SetBranchAddress("nFltrdEvts", &nFltrdEvts, &b_nFltrdEvts);
        fConfigTree->SetBranchAddress("nTotEvts", &nTotEvts, &b_nTotEvts);
        fConfigTree->SetBranchAddress("sumEvtWgt", &sumEvtWgt, &b_sumEvtWgt);
        fConfigTree->SetBranchAddress("sumFltrdEvtWgt", &sumFltrdEvtWgt, &b_sumFltrdEvtWgt);
    
        auto nConfigEntries = fConfigTree->GetEntries();
        std::cout << "Proccessing " << nConfigEntries << " config entries." << std::endl;
    	configInfo.clear();
    	cutflowInfo.clear();
        for (Long64_t centry = 0; centry < nConfigEntries; centry++){
    		
    		auto entry = fConfigTree->LoadTree(centry);
    
            if(debug) std::cout << " - Getting Branches. " << std::endl;
        	b_nEvents->GetEntry(entry);   //!
        	b_nSelectedEvents->GetEntry(entry);   //!
        	b_sKey->GetEntry(entry);   //!  
        	b_sCrossSection->GetEntry(entry);   //!
        	b_sGMSBGravMass->GetEntry(entry);   //!
        	b_sGMSBChi1Mass->GetEntry(entry);   //!
        	b_sMCWgt->GetEntry(entry);   //!
        	b_sMCType->GetEntry(entry);   //!
            std::string configKey(*sKey);
    
            b_sumEvtGenWgt->GetEntry(entry);
            b_cf_m_gt2jets->GetEntry(entry);   //!
            b_cf_met150->GetEntry(entry);   //!
            b_cf_mj_gt1phos->GetEntry(entry);   //!
            b_cf_mjp_leadPhoPt30->GetEntry(entry);   //!
            b_cf_nFltrdEvts->GetEntry(entry);   //!
            b_cf_nTotEvts->GetEntry(entry);   //!
            b_cf_sel_j->GetEntry(entry);   //!
            b_cf_sel_m->GetEntry(entry);   //!
            b_cf_sel_p->GetEntry(entry);   //!
            b_cf_sel_ppt->GetEntry(entry);   //!
            b_nFltrdEvts->GetEntry(entry);   //!
            b_nTotEvts->GetEntry(entry);   //!
            b_sumEvtWgt->GetEntry(entry);   //!
            b_sumFltrdEvtWgt->GetEntry(entry);   //!
    
    		if( not configInfo.count(configKey) ){
            	if(debug) std::cout << " - Filling Cutflow/Wieghts. " << std::endl;
            	//fConfigTree->Add(tfilename.c_str());			
    			std::map< std::string, float > configValues;
            	configValues["nEvents"] = nEvents;
            	configValues["nSelectedEvents"] = nSelectedEvents;
                configValues["sumEvtGenWgt"] = sumEvtGenWgt;
                configValues["sumEvtWgt"] = sumEvtWgt;
                if(debug) std::cout << " - Filling configValues. " << std::endl;
            	configValues["sCrossSection"] = sCrossSection;
            	configValues["sGMSBGravMass"] = sGMSBGravMass;
            	configValues["sGMSBChi1Mass"] = sGMSBChi1Mass;
            	configValues["sMCWgt"] = sMCWgt;
            	configValues["sMCType"] = sMCType;
                configValues["nTotEvts"] = nTotEvts;
            	if(debug) std::cout << " - Filling configInfo. " << std::endl;
            	configInfo[configKey] = configValues;
    		} else {
    			auto & configValues = configInfo[configKey];
    			configValues["nEvents"] += nEvents;
                configValues["nTotEvts"] += nTotEvts;
    			configValues["nSelectedEvents"] += nSelectedEvents;
                configValues["sumEvtGenWgt"] += sumEvtGenWgt;
                configValues["sumEvtWgt"] += sumEvtWgt;
    		}//<<>>if( not configInfo.count(configKey) )
    
    		//float fillwt = scale * ( sCrossSection * 1000 ) * ( 1 / sumEvtWgt );
            float fillwt = scale * ( sCrossSection * 1000 ) * ( 1.0 / float( nTotEvts) );
			std::cout << " fillwt : " << scale <<  " * ( " << sCrossSection << " * 1000 ) * ( " << " 1 / " << nTotEvts << " ) = " 
						<< fillwt << std::endl;
            if( not cutflowInfo.count("nTotEvts") ){
                if(debug) std::cout << " - Filling Cutflow/Wieghts. " << std::endl;
                cutflowInfo["cf_m_gt2jets"] = cf_m_gt2jets*fillwt;
                cutflowInfo["cf_met150"] = cf_met150*fillwt;
                cutflowInfo["cf_mj_gt1phos"] = cf_mj_gt1phos*fillwt;
                cutflowInfo["cf_mjp_leadPhoPt30"] = cf_mjp_leadPhoPt30*fillwt;
                cutflowInfo["cf_nFltrdEvts"] = cf_nFltrdEvts*fillwt;
                cutflowInfo["cf_nTotEvts"] = float(cf_nTotEvts)*fillwt;
                cutflowInfo["cf_sel_j"] = cf_sel_j*fillwt;
                cutflowInfo["cf_sel_m"] = cf_sel_m*fillwt;
                cutflowInfo["cf_sel_p"] = cf_sel_p*fillwt;
                cutflowInfo["cf_sel_ppt"] = cf_sel_ppt*fillwt;
                cutflowInfo["nFltrdEvts"] = nFltrdEvts*fillwt;
                cutflowInfo["nTotEvts"] = float(nTotEvts)*fillwt;
                //cutflowWgts["sumEvtWgt"] = sumEvtWgt;
                //cutflowWgts["sumFltrdEvtWgt"] = sumFltrdEvtWgt;
            } else {
                cutflowInfo["cf_m_gt2jets"] += cf_m_gt2jets*fillwt;
                cutflowInfo["cf_met150"] += cf_met150*fillwt;
                cutflowInfo["cf_mj_gt1phos"] += cf_mj_gt1phos*fillwt;
                cutflowInfo["cf_mjp_leadPhoPt30"] += cf_mjp_leadPhoPt30*fillwt;
                cutflowInfo["cf_nFltrdEvts"] += cf_nFltrdEvts*fillwt;
                cutflowInfo["cf_nTotEvts"] += float(cf_nTotEvts)*fillwt;
                cutflowInfo["cf_sel_j"] += cf_sel_j*fillwt;
                cutflowInfo["cf_sel_m"] += cf_sel_m*fillwt;
                cutflowInfo["cf_sel_p"] += cf_sel_p*fillwt;
                cutflowInfo["cf_sel_ppt"] += cf_sel_ppt*fillwt;
                cutflowInfo["nFltrdEvts"] += nFltrdEvts*fillwt;
                cutflowInfo["nTotEvts"] += float(nTotEvts)*fillwt;
                //cutflowWgts["sumEvtWgt"] += sumEvtWgt;
                //cutflowWgts["sumFltrdEvtWgt"] += sumFltrdEvtWgt;
            }//<<>>if( not configInfo.count(configKey) )
    
    	}//<<>>for (Long64_t centry = 0; centry < nConfigEntries; centry++)
    
        for( auto item : configInfo ){ 
    		std::cout << item.first << " ( ";  
    		for( auto line : item.second ){ 
    			std::cout << line.first <<  " " << line.second << " ";
    		}//<<>>for( auto line : item.second )
    		std::cout << ")" << std::endl;
    	}//<<>>for( auto item : configInfo )
    
        std::cout << " CutFlow Info ( ";
        for( auto line : cutflowInfo ){
            std::cout << line.first <<  " " << line.second << " ";
        }//<<>>for( auto line : item.second )
        std::cout << ")" << std::endl;   

 
        std::cout << "<<<<<<<< Processing Event Loop <<<<<<<<<<<<<< " << std::endl;
    	int loopCounter(100000);
        auto nEntries = fInTree->GetEntries();
        if(debug){ nEntries = 10; loopCounter = 1; }
        std::cout << "Proccessing " << nEntries << " entries." << std::endl;
        for (Long64_t centry = 0; centry < nEntries; centry++){
            if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
            if(debug) std::cout << "*****************************************************************************" << std::endl;
            auto entry = fInTree->LoadTree(centry);
    		if(debug) std::cout << " - getBranches " << std::endl;
    		getBranches(entry);
    		if(debug) std::cout << " - eventLoop " << std::endl;
    		eventLoop(entry, m_vec, r_vec, rv_vec );
        }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
       
        if(debug) std::cout << " - Creating output file " << std::endl;
    	std::string fulloutfilename =  outfilename + intitle + ofnending;
        TFile* fOutFile = new TFile( fulloutfilename.c_str(), "RECREATE" );
        fOutFile->cd();
    
        std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;
    
    	endJobs();
    	for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
        for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
        for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }
    
    	nMaps = 0;
    	if( doEBEEmaps ){ for( int it = 0; it < nEBEEMaps; it++ ){ 
    		ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
    		ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
    		ebeeMapR[it]->Write(); delete ebeeMapR[it];
    	}}//<<>>for( int it = 0; it < nEBEEMaps; it++ )
    
        fOutFile->Close();

	}//<<>>while (std::getline(infile,str))

    std::cout << "histMaker : Thats all Folks!!" << std::endl;
}//<<>>void kucmsSkimmer
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------- event loop -------------------------------------------------------------------------
void HistMaker::eventLoop( Long64_t entry, std::vector<float> m_vec, std::vector<float> r_vec, std::vector<float> rv_vec ){

	int cs = ( cutselection == 0 ) ? 0 : 1;

    //cutv = cut;  n rjr jets per side
    //cutva = va;  n rjr phos in event
    //cutvb = vb;  n sel phos in event
    //cutvc = vc;
    //cutvd = vd;

// ( minJetsPer, nphos, maxSelPhos, rjrcut2, metcut )
    //auto dskey  = *DataSetKey
    //float evtgwt = evtGenWgt;
    float evtgwt = 1;
	float cfbin = 1;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    float segwt = (configInfo[configKey])["nTotEvts"];
    //float segwt = (configInfo[configKey])["sumEvtWgt"];
    float fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
    if( DEBUG ) std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;

    //nRjrPhos = nphos;
    //nRjrJets = njets;
    //hist1d[2]->Fill(7,fillwt);// first eventloop cutflow -> should be same as n-1 cutflow

	for( int i = 0; i < 1; i++ ){ // continue loop - can use continue to skip
	/////////////////////////////////////////////////////////////////////////
	//int nRjrPhotons = rjrNRjrPhotons->at(cs);
	//if( nRjrPhotons == 0 ) continue;
    //hist1d[2]->Fill(7,fillwt);

	//!!!!! nRjrPhotons Cut
    //bool nRjrPhoCut( nRjrPhotons != nRjrPhos );
    //if( cutva == 0 ) nRjrPhoCut = false;
	//cfbin++;
    //if( nRjrPhos == 1 ) hist1d[2]->Fill(cfbin,fillwt);
    //if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"rph1");
    //cfbin++;
    //if( nRjrPhos == 2 ) hist1d[2]->Fill(cfbin,fillwt);
    //if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"rph2");
    if( rjrNRjrPhotons->at(cs) != nRjrPhos ) continue;
    //cfbin++;
    //hist1d[2]->Fill(cfbin,fillwt);
    //if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"rphc");

	//if( DEBUG ) std::cout << "RJR cut at " << cutvc << " with value " << rjrX2GMean->at(0) << std::endl;

	//!!!!! MET Cut
    auto metCPt = hypo(selCMetPx,selCMetPy);
    auto metPt = hypo(selMetPx,selMetPy);
	if( metCPt < 150 ) continue;
    //hist1d[2]->Fill(12,fillwt);

// ( minJetsPer, nphos, maxSelPhos, rjrcut2, metcut )
	//exploritory cut flows // nPhotons
	int numPhos = int( nSelPhotons );
	hist1d[100]->Fill(numPhos,fillwt);
	//if( numPhos > 0 ) hist1d[2]->Fill(10,fillwt);
    //hist1d[2]->GetXaxis()->SetBinLabel(10,"ph0+");
    //if( numPhos > 1 ) hist1d[2]->Fill(11,fillwt);
    //hist1d[2]->GetXaxis()->SetBinLabel(11,"ph1+");
    //bool XPhoCut = numPhos > cutva;
    //cfbin++;
    //if( numPhos > cutva ) hist1d[2]->Fill(cfbin,fillwt);
    //if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"xtph");
	//if(  ) continue;
	//hist1d[2]->Fill(13,fillwt);
    //hist1d[2]->GetXaxis()->SetBinLabel(13,"cut");

//------------------------------------------------------------------------------------

    //cfbin++;
    //if( (*rjrNJetsJa)[cs] > 0 && (*rjrNJetsJb)[cs] > 0 ) hist1d[2]->Fill(cfbin,fillwt);
    //if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"jab0");
    //cfbin++;
    //if( (*rjrNJetsJa)[cs] > 1 && (*rjrNJetsJb)[cs] > 1 ) hist1d[2]->Fill(cfbin,fillwt);
    //if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"jab1");
    //if( (*rjrNJetsJa)[cs] > 3 && (*rjrNJetsJb)[cs] > 3 ) hist1d[2]->Fill(cfbin++,fillwt);

	// Min # jets/side
	if( (*rjrNJetsJa)[cs] < nRjrJets || (*rjrNJetsJb)[cs] < nRjrJets ) continue;

    cfbin++;
    hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"PSel");


	// rjr var cuts
	//cutvc = vc; ax2nq
    //cutvd = vd; nvsum
	// ve	asmass
	// vf vdiff

    float nmass = (*rjrASMass)[cs];
    float rjrVDiff = (*rjrMVDiff)[cs]; //((*rjrMVa)[cs]-(*rjrMVb)[cs])/((*rjrMVa)[cs]+(*rjrMVb)[cs]); 
    float rjrVSum = (*rjrMVSum)[cs]; //std::sqrt((sq2((*rjrMVa)[cs])+sq2((*rjrMVb)[cs]))/2);
    float rjrNVSum = 2*rjrVSum/nmass;

    float NormMBetaEql = (*selPhoMBetaEql)[0]/nmass;
    float NormMBetaPmt = (*selPhoMBetaPmt)[0]/nmass;

	//(*rjrASMass)[cs] -- M
	//(*rjrAX2NQSum)[cs] -- R
	//rjrNVSum -- Rv

	bool MT = (*rjrASMass)[cs] > m_vec[2];
    bool MM = (*rjrASMass)[cs] > m_vec[1];
    bool ML = (*rjrASMass)[cs] > m_vec[0];

	bool RT = (*rjrAX2NQSum)[cs] > r_vec[2];
    bool RM = (*rjrAX2NQSum)[cs] > r_vec[1];
    bool RL = (*rjrAX2NQSum)[cs] > r_vec[0];

    bool RvT = rjrNVSum > rv_vec[2];
    bool RvM = rjrNVSum > rv_vec[1];
    bool RvL = rjrNVSum > rv_vec[0];

    cfbin++;
    cfbin++;
    if( ML && RL && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LLL");
    cfbin++;
    if( ML && RL && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LLM");
    cfbin++;
    if( ML && RL && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LLT");
    //cfbin++;
    cfbin++;
    if( ML && RM && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LML");
    cfbin++;
    if( ML && RM && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LMM");
    cfbin++;
    if( ML && RM && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LMT");
    //cfbin++;
    cfbin++;
    if( ML && RT && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LTL");
    cfbin++;
    if( ML && RT && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LTM");
    cfbin++;
    if( ML && RT && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"LTT");


    cfbin++;
    cfbin++;
    if( MM && RL && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MLL");
    cfbin++;
    if( MM && RL && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MLM");
    cfbin++;
    if( MM && RL && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MLT");
    //cfbin++;
    cfbin++;
    if( MM && RM && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MML");
    cfbin++;
    if( MM && RM && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MMM");
    cfbin++;
    if( MM && RM && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MMT");
    //cfbin++;
    cfbin++;
    if( MM && RT && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MTL");
    cfbin++;
    if( MM && RT && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MTM");
    cfbin++;
    if( MM && RT && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"MTT");

    cfbin++;
    cfbin++;
    if( MT && RL && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TLL");
    cfbin++;
    if( MT && RL && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TLM");
    cfbin++;
    if( MT && RL && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TLT");
    //cfbin++;
    cfbin++;
    if( MT && RM && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TML");
    cfbin++;
    if( MT && RM && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TMM");
    cfbin++;
    if( MT && RM && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TMT");
    //cfbin++;
    cfbin++;
    if( MT && RT && RvL ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TTL");
    cfbin++;
    if( MT && RT && RvM ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TTM");
    cfbin++;
    if( MT && RT && RvT ) hist1d[2]->Fill(cfbin,fillwt);
    if( fl ) hist1d[2]->GetXaxis()->SetBinLabel(cfbin,"TTT");

	//var hist fill

    hist1d[402]->Fill( rjrMET->at(cs), fillwt );
    hist1d[403]->Fill( metCPt, fillwt );

	if( DEBUG ) std::cout << " -- Filling Histograms set 1" << std::endl;
    if( true ){

	    hist1d[255]->Fill((*rjrNJetsJa)[cs], fillwt ); 
	    hist1d[256]->Fill((*rjrNJetsJb)[cs], fillwt ); 

    if( DEBUG ) std::cout << " -- Filling Histograms set 12" << std::endl;

        hist1d[268]->Fill((*selPhoMBetaEql)[0], fillwt );
        hist1d[269]->Fill((*selPhoMBetaPmt)[0], fillwt );

        hist1d[272]->Fill(rjrVDiff, fillwt );
        hist1d[273]->Fill(rjrNVSum, fillwt );

        hist1d[274]->Fill(NormMBetaEql, fillwt );
        hist1d[275]->Fill(NormMBetaPmt, fillwt );

    if( DEBUG ) std::cout << " -- Filling Histograms set 2" << std::endl;

	    hist1d[351]->Fill((*rjrASMass)[cs], fillwt ); //("SMass", mkht(ht,"SMass"), 500, 0, 5000);
	    hist1d[366]->Fill((*rjrAX2NQSum)[cs], fillwt );

    if( DEBUG ) std::cout << " -- Filling Histograms set 3" << std::endl;

		hist2d[0]->Fill( (*rjrASMass)[cs], (*rjrAX2NQSum)[cs], fillwt );
    	hist2d[2]->Fill( (*rjrASMass)[cs], rjrNVSum, fillwt ); 
    	hist2d[3]->Fill( (*rjrASMass)[cs], rjrVDiff, fillwt ); 
    	hist2d[4]->Fill( (*rjrAX2NQSum)[cs], rjrNVSum, fillwt ); 
    	hist2d[5]->Fill( (*rjrAX2NQSum)[cs], rjrVDiff, fillwt ); 
    	hist2d[9]->Fill( rjrNVSum, rjrVDiff, fillwt ); 
		hist2d[10]->Fill( (*rjrNJetsJa)[cs], (*rjrNJetsJb)[cs], fillwt );


	}//<<>>if( rjrAX2NGMean->at(cs) < cutvc && rjrAX2GMean->at(cs) < cutvd )
	}//<<>>for( int i = 0; i < 1; i++ ) -- continue loop

	//fl = false;
	
}//<<>>void HistMaker::eventLoop(` Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

	std::cout << " -- Filling cutflow bins 1-4 " << cutflowInfo["cf_nTotEvts"] << " : " << cutflowInfo["cf_m_gt2jets"] << std::endl;
    hist1d[2]->SetBinContent(1,cutflowInfo["cf_nTotEvts"]);
    //hist1d[2]->SetBinError(1,std::sqrt(cutflowInfo["cf_nTotEvts"]));
	hist1d[2]->GetXaxis()->SetBinLabel(1,"Evts");
    //hist1d[2]->SetBinContent(2,cutflowInfo["cf_nFltrdEvts"]);
    //hist1d[2]->SetBinError(2,std::sqrt(cutflowInfo["cf_nFltrdEvts"]));
    //hist1d[2]->GetXaxis()->SetBinLabel(2,"Fltr");
    //hist1d[2]->SetBinContent(3,cutflowInfo["cf_met150"]);
    //hist1d[2]->SetBinError(3,std::sqrt(cutflowInfo["cf_met150"]));
    //hist1d[2]->GetXaxis()->SetBinLabel(3,"M150");
    //hist1d[2]->SetBinContent(4,cutflowInfo["cf_m_gt2jets"]);
    //hist1d[2]->SetBinError(4,std::sqrt(cutflowInfo["cf_m_gt2jets"]));
    //hist1d[2]->GetXaxis()->SetBinLabel(4,"gt2J");
    //hist1d[2]->SetBinContent(5,cutflowInfo["cf_mj_gt1phos"]);
    //hist1d[2]->SetBinError(5,std::sqrt(cutflowInfo["cf_mj_gt1phos"]));
    //hist1d[2]->SetBinContent(6,cutflowInfo["cf_mjp_leadPhoPt30"]);
    //hist1d[2]->SetBinError(6,std::sqrt(cutflowInfo["cf_mjp_leadPhoPt30"]));
    //hist1d[2]->SetBinContent(2,cutflowInfo["cf_mjp_leadPhoPt30"]);
    //hist1d[2]->GetXaxis()->SetBinLabel(2,"PSel");

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

    std::cout << " title test : " << addstr(ht,"_varTitle") << std::endl;

    hist1d[2] = new TH1D("cutflow", addstr(ht,"_cutFlow").c_str(), 32, 0.5, 32.5);
	hist1d[2]->GetXaxis()->SetLabelSize(.02);

	hist1d[100] = new TH1D("nSelPho", addstr(ht,"nSelPho").c_str(), 20, 0, 20);

    hist1d[255] = new TH1D("rjrNJetsJa", addstr(ht,"rjrNJetsJa").c_str(), 20, 0, 20);
    hist1d[256] = new TH1D("rjrNJetsJb", addstr(ht,"rjrNJetsJb").c_str(), 20, 0, 20);

    hist1d[268] = new TH1D("selPhoMBetaEql", addstr(ht,"selPhoBetaMass").c_str(), 60, 0, 1200);//----------- BetaMass
    hist1d[269] = new TH1D("selPhoMBetaPmt", addstr(ht,"selPhoMBetaPmt").c_str(), 60, 0, 1200);

    hist1d[272] = new TH1D("rjrVDiff", addstr(ht,"NVisMassDiff").c_str(), 25, -1, 1);//---------------VisMassDiff !!!!!!!!!!
    hist1d[273] = new TH1D("rjrVSum", addstr(ht,"NVisMass").c_str(), 25, 0, 1.0);//--------------VisMassSum !!!!!!!!!!

    hist1d[274] = new TH1D("rjrPhoNMBetaEql", addstr(ht,"rjrPhoNBetaEql").c_str(), 100, 0, 2);//----------- BetaMass
    hist1d[275] = new TH1D("rjrPhoNMBetaPmt", addstr(ht,"rjrPhoNMBetaPmt").c_str(), 100, 0, 2);

    hist1d[351] = new TH1D("ASMass", addstr(ht,"ASMass").c_str(), 48, 0, 12000); //----------------------ASMass !!!!!!!!!!!!!!!!!!1
    hist1d[366] = new TH1D("rjrAX2NQSum", addstr(ht,"X2MassRatio").c_str(), 25, 0, 1.0);//----------X2MassRatio !!!!!!!!!!!!!!!!!


	//-------- event vars 400 - 450
    hist1d[402] = new TH1D("evtRjrMetPt", addstr(ht,"evtRjrMetPt ").c_str(), 100, 0, 1500 );
    hist1d[403] = new TH1D("evtCorrMetPt", addstr(ht,"evtCorrMetPt ").c_str(), 100, 0, 1500 );

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    hist2d[0] = new TH2D("ASMvAX2NQSum", addstr(ht,"ASMvAX2NQSum;ASMass [GeV];X2ANQSum").c_str(), 48, 0, 12000, 25, 0, 1.0 );
    hist2d[2] = new TH2D("ASMvNVMass", addstr(ht,"ASMvNVSum;ASMass [GeV];NVSum").c_str(), 48, 0, 12000, 25, 0, 1.0 );
    hist2d[3] = new TH2D("ASMvVDiff", addstr(ht,"ASMvVDiff;ASMass [GeV];VDiff").c_str(), 48, 0, 12000, 25, -1.0, 1.0 );
    hist2d[4] = new TH2D("AX2NQSumvVMass", addstr(ht,"AX2NQSumvNVSum;AX2NQSum;NVSum").c_str(), 25, 0, 1.0, 25, 0, 1.0 );
    hist2d[5] = new TH2D("AX2NQSumvVDiff", addstr(ht,"AX2NQSumvVDiff;AX2NQSum;VDiff").c_str(), 25, 0, 1.0, 25, -1.0, 1.0 );
    hist2d[9] = new TH2D("VMassvVDiff", addstr(ht,"NVSumvVDiff;NVSum;VDiff").c_str(), 25, 0, 1.0, 25, -1.0, 1.0 );
    hist2d[10] = new TH2D("NJetsJavNJetsJb", addstr(ht,"NJetsJavNJetsJb;NJetsJa;NJetsJb").c_str(), 20, 0, 20,20, 0, 20 );

	//------- jets ( time ) 0-49 ------------------------------

	//---jet id stuff 50 - 99 ---------------------------------------------------

	//--- Photons 200 - 349 -------------------------------------------
	

	//60 - 63

	//--- rechit collections 350 - 399 -------------------------------------------------

    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------
	nMaps = 0;
    std::string label("baseHists");
	if( doEBEEmaps ){ for(int it=0; it<nEBEEMaps; it++){
		fMap[it] = false;
		std::string label(";iEta;iPhi");
        std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
        ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
        std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
        ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
		std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
        ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
	}}//<<>>for(int it=0; it<nEBEEMaps; it++)

}//<<>>void HistMaker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                std::string listdir = "/uscms/home/jaking/nobackup/llpana_skims/";
            
				std::string infilenameJ = "rjr_skim_files/KUCMS_RJR_GIGI_ootmet_Skim_List.txt";
				std::string infilenameBG = "rjr_skim_files/KUCMS_RJR_BG_ootmet_Skim_List.txt";

				std::string version = "_v24_";
				std::string sigtype = "llpana";
				std::string ofnstart = "KUCMS_";

                std::string htitleBG = "BG_"+sigtype+version;
				std::string htitleJ = "Sgg10_"+sigtype+version;

                HistMaker base;

                //float metcut = 150;  // cmet > 150
    			// rjr var cuts
    			//cutvc = vc; ax2nq
    			//cutvd = vd; nvsum
    			// ve   asmass
    			// vf vdiff
    
				//int nj = 1;
				//int np = 1;
                for( int np = 1; np < 3; np++ ){
                for( int nj = 1; nj < 3; nj++ ){

				std::string subdir = "cf_" + std::to_string(np) + "pho_" + std::to_string(nj) + "jet/";

				// 1-1
				//std::vector<float> m_vec{1000,2000,3000}; // L-M-T
				std::vector<float> m_vec{2000,2750,3500}; // L-M-T
                std::vector<float> r_vec{0.2,0.275,0.35}; // L-M-T
                std::vector<float> rv_vec{0.1,0.25,0.4}; // L-M-T
				//std::string outdir = "cf_1-2-3_2-275-35_1-25-40/" + subdir;
                std::string outdir = "cf_2-275-35_2-275-35_1-25-40/" + subdir;
				// 2-2   cf_1pho_4jet
                //std::vector<float> m_vec{1000,1500,2000}; // L-M-T
                //std::vector<float> r_vec{0.15,0.225,0.3}; // L-M-T
                //std::vector<float> rv_vec{0.1,0.175,0.25}; // L-M-T

				std::string isoline = "mj" + std::to_string( nj ) + "_";
				isoline += "rp" + std::to_string( np ) + "_" ;
                std::string outfilenameJ = outdir + ofnstart + htitleJ + isoline;
				std::string htitlefullJ =  htitleJ + isoline;
				//base.histMaker( listdir, infilenameJ, outfilenameJ, htitlefullJ, np, nj, m_vec, r_vec, rv_vec );
                std::string outfilenameBG = outdir + ofnstart + htitleBG + isoline;
                std::string htitlefullBG =  htitleBG + isoline;
                base.histMaker( listdir, infilenameBG, outfilenameBG, htitlefullBG, np, nj, m_vec, r_vec, rv_vec );

				}}


    return 1;


}//<<>>int main ( int argc, char *argv[] )

