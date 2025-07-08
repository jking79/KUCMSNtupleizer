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
	const std::string ofnending = "_RjrSkim_v24_batchHists.root";

    cutselection = 1; // jrj type -> pho visible
    fl = true; // first event loop

	nRjrPhos = nphos;
	nRjrJets = njets;

	scale = 400;
	float glumass = 0;
	datalist.clear();

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
		glumass = ( instrs.size() > 2 ) ? std::stof( instrs[3] ) : 0;
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
   
        std::cout << "Intilizing BatchVars" << std::endl;
 
		skimlabel = intitle;
		datalist.push_back( skimlabel );
        std::string bglumass = "glueMass_" + skimlabel;
        batchVars[bglumass] = glumass;

    	std::string bfillwgt = "sumFillWgt_" + skimlabel;
    	batchVars[bfillwgt] = 0;
        std::string psfillwgt = "sumPSFillWgt_" + skimlabel;
        batchVars[psfillwgt] = 0;

        std::string bmetfilter = "metFilter" + skimlabel;
        batchVars[bmetfilter] = 0;
        std::string bhaloFilter = "haloFilter" + skimlabel;
        batchVars[bhaloFilter] = 0;
        std::string bnSVHad = "nSVHad" + skimlabel;
        batchVars[bnSVHad] = 0;
        std::string bnSVLep = "nSVLep" + skimlabel;
        batchVars[bnSVLep] = 0;
        std::string bnEleSVM = "nEleSVM" + skimlabel;
        batchVars[bnEleSVM] = 0;
        std::string bnSelIsoEle = "nSelIsoEle" + skimlabel;
        batchVars[bnSelIsoEle] = 0;
        std::string bnLooseMuons = "nLooseMuons" + skimlabel;
        batchVars[bnLooseMuons] = 0;
        std::string bnMuonSVM = "nMuonSVM" + skimlabel;
        batchVars[bnMuonSVM] = 0;
    	std::string blepSVM = "SVLepVeto" + skimlabel;
        batchVars[blepSVM] = 0;
    	std::string bisolept = "LeptonVeto" + skimlabel;
        batchVars[bisolept] = 0;

        std::string bCutA = "CutA_" + skimlabel;        
        batchVars[bCutA] = 0;

        std::string bfinal = "final_" + skimlabel;        
        batchVars[bfinal] = 0;

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
        //Float_t         sumEvtGenWgt;
    
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
        //TBranch        *b_sumEvtGenWgt;
    
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
        //fConfigTree->SetBranchAddress("sumEvtGenWgt", &sumEvtGenWgt, &b_sumEvtGenWgt);
    
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
    
            //b_sumEvtGenWgt->GetEntry(entry);
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
                //configValues["sumEvtGenWgt"] = sumEvtGenWgt;
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
                //configValues["sumEvtGenWgt"] += sumEvtGenWgt;
                configValues["sumEvtWgt"] += sumEvtWgt;
    		}//<<>>if( not configInfo.count(configKey) )
    
    		float fillwt = scale * ( sCrossSection * 1000 ) * ( 1 / sumEvtWgt );
            //float fillwt = scale * ( sCrossSection * 1000 ) * ( 1.0 / float( nTotEvts) );
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
      
		endJobs();
 
        if(debug) std::cout << " - Creating output file " << std::endl;
    	//std::string fulloutfilename =  outfilename + intitle + ofnending;
        //TFile* fOutFile = new TFile( fulloutfilename.c_str(), "RECREATE" );
        //fOutFile->cd();
    
        std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;
    
    	endJobs();
    	//for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
        //for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
        //for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }
        //for( int it = 0; it < n1dHists; it++ ){ delete hist1d[it]; } }
        //for( int it = 0; it < n2dHists; it++ ){ delete hist2d[it]; } }
        //for( int it = 0; it < n3dHists; it++ ){ delete hist3d[it]; } }    


    	nMaps = 0;
    	if( doEBEEmaps ){ for( int it = 0; it < nEBEEMaps; it++ ){ 
    		//ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
    		//ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
    		//ebeeMapR[it]->Write(); delete ebeeMapR[it];
			delete ebeeMapP[it];
			delete ebeeMapT[it];
			delete ebeeMapR[it];
    	}}//<<>>for( int it = 0; it < nEBEEMaps; it++ )
    
        //fOutFile->Close();

	}//<<>>while (std::getline(infile,str))

	endBatchJobs();

	std::string boutfilename = "Batch_" + outfilename;
    std::string bintitle = skimlabel;
    std::string bofnending = ofnending;
    if(debug) std::cout << " - Creating batch output file " << std::endl;
    std::string bfulloutfilename =  boutfilename + bintitle + bofnending;
    TFile* fOutFile = new TFile( bfulloutfilename.c_str(), "RECREATE" );
    fOutFile->cd();

	//  save stuff here ......
    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }



    fOutFile->Close();

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

    std::string bfillwgt = "sumFillWgt_" + skimlabel;
	std::string psfillwgt = "sumPSFillWgt_" + skimlabel;

    std::string bmetfilter = "metFilter_" + skimlabel;
    std::string bhaloFilter = "haloFilter_" + skimlabel;
    std::string bnSVHad = "nSVHad_" + skimlabel;
    std::string bnSVLep = "nSVLep_" + skimlabel;
    std::string bnEleSVM = "nEleSVM_" + skimlabel;
    std::string bnMuonSVM = "nMuonSVM_" + skimlabel;
    std::string blepSVM = "SVLepVeto_" + skimlabel;
    std::string bnSelIsoEle = "nSelIsoEle_" + skimlabel;
    std::string bnLooseMuons = "nLooseMuons_" + skimlabel;
	std::string bisolept = "LeptonVeto_" + skimlabel;

    std::string bglumass = "glueMass_" + skimlabel;

    std::string bCutA = "CutA_" + skimlabel; 
    std::string bfinal = "final_" + skimlabel;


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
	batchVars[bfillwgt] += fillwt;

    //nRjrPhos = nphos;
    //nRjrJets = njets;
    //hist1d[2]->Fill(7,fillwt);// first eventloop cutflow -> should be same as n-1 cutflow

	for( int i = 0; i < 1; i++ ){ // continue loop - can use continue to skip
	/////////////////////////////////////////////////////////////////////////

	//!!!!! nRjrPhotons Cut
    if( ( nRjrPhos != 0 ) &&  rjrNRjrPhotons->at(cs) != nRjrPhos ) continue;
	if( ( nRjrPhos == 0 ) &&  ( ( SV_nHadronic < 1 ) && ( SV_nLeptonic < 1 ) ) ) continue;

	//if( DEBUG ) std::cout << "RJR cut at " << cutvc << " with value " << rjrX2GMean->at(0) << std::endl;

    if( DEBUG ) std::cout << "Met cut " << std::endl;
	//!!!!! MET Cut
    auto metCPt = hypo(selCMetPx,selCMetPy);
    auto metPt = hypo(selMetPx,selMetPy);
	if( metCPt < 150 ) continue;
    //hist1d[2]->Fill(12,fillwt);

	if( DEBUG ) std::cout << "RJR jet cut " << std::endl;
	// Min # jets/side
	if( (*rjrNJetsJa)[cs] < nRjrJets || (*rjrNJetsJb)[cs] < nRjrJets ) continue;

    batchVars[psfillwgt] += fillwt;

	bool lElectronF = false;
	for( auto ele : *eleIsLoose	){ if( ele ){ lElectronF = true; break; } }
    bool lMuonF = false;
    for( auto muon : *muonIsLoose ){ if( muon ){ lMuonF = true; break; } }
    //for( auto muon : *muonIsMedium ){ if( muon ){ lMuonF = true; break; } }
	bool hasLLepF = lElectronF || lMuonF;

    if( DEBUG ) std::cout << "Batch bools " << std::endl;
	//bool passmet = not MetFilter;
	//bool passhalo = not HaloFilter;
	bool svhardv = SV_nHadronic > 0;
	bool svlepv = SV_nLeptonic > 0;	

	bool esvmatv = nEleSVMatched < 1;
	bool msvmatv = nMuonSVMatched < 1;
	bool svlepmatv = esvmatv && msvmatv;

	bool isoev = nSelIsoElectrons < 1;
	bool isomv = nLooseMuons < 1;
	bool isolepv = isoev && isomv;

    if( DEBUG ) std::cout << "Batch fill " << std::endl;
	if( !MetFilter ) batchVars[bmetfilter] += fillwt;
    if( !HaloFilter ) batchVars[bhaloFilter] += fillwt;
    if( svhardv ) batchVars[bnSVHad] += fillwt;
    if( svlepv ) batchVars[bnSVLep] += fillwt;
	
    if( esvmatv ) batchVars[bnEleSVM] += fillwt;
    if( msvmatv ) batchVars[bnMuonSVM] += fillwt;
    if( svlepmatv ) batchVars[blepSVM] += fillwt;

    if( isoev ) batchVars[bnSelIsoEle] += fillwt;
    if( isomv ) batchVars[bnLooseMuons] += fillwt;
    if( isolepv ) batchVars[bisolept] += fillwt;

    //if( !lElectronF ) batchVars[bnSelIsoEle] += fillwt;
    //if( !lMuonF ) batchVars[bnLooseMuons] += fillwt;
    //if( !hasLLepF ) batchVars[bisolept] += fillwt;

	if( !MetFilter && svlepmatv ) batchVars[bCutA] += fillwt;
	if( !MetFilter && !hasLLepF ) batchVars[bfinal] += fillwt;

	// rjr var cuts
	//cutvc = vc; ax2nq
    //cutvd = vd; nvsum
	// ve	asmass
	// vf vdiff

    if( DEBUG ) std::cout << "RJR info " << std::endl;

    float mr = (*rjr_Mr)[cs];
    float rjrNVSum = (*rjr_Rv)[cs];
	float rs = (*rjr_R)[cs];

    //float NormMBetaEql = (*selPhoMBetaEql)[0]/mr;
    //float NormMBetaPmt = (*selPhoMBetaPmt)[0]/mr;

	//(*rjrASMass)[cs] -- M
	//(*rjrAX2NQSum)[cs] -- R
	//rjrNVSum -- Rv

	bool MT =  mr > m_vec[2];
    bool MM =  mr > m_vec[1];
    bool ML =  mr > m_vec[0];

	bool RT = rs > r_vec[2];
    bool RM = rs > r_vec[1];
    bool RL = rs > r_vec[0];

    bool RvT = rjrNVSum > rv_vec[2];
    bool RvM = rjrNVSum > rv_vec[1];
    bool RvL = rjrNVSum > rv_vec[0];

	//var hist fill

	//if( not( ML && RL && RvM ) ) continue;


	}//<<>>for( int i = 0; i < 1; i++ ) -- continue loop

	
}//<<>>void HistMaker::eventLoop(` Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){


}//<<>>void HistMaker::endJobs()

void HistMaker::endBatchJobs(){

  	//TH1D *sigHist[n1dHists];
    //TH1D *bkgHist[n1dHists];
    //TH1D *dataHist[n1dHists];
    //std::map< std::string, float > batchVars;
    //std::map< std::string, std::vector< float > > batchVecVars;
	//batchVars["glumass"]

    for( auto skimlabel : datalist ){

    	std::string bfillwgt = "sumFillWgt_" + skimlabel;
    	std::string psfillwgt = "sumPSFillWgt_" + skimlabel;

    	std::string bmetfilter = "metFilter_" + skimlabel;
    	std::string bhaloFilter = "haloFilter_" + skimlabel;
    	std::string bnSVHad = "nSVHad_" + skimlabel;
    	std::string bnSVLep = "nSVLep_" + skimlabel;
    	std::string bnEleSVM = "nEleSVM_" + skimlabel;
    	std::string bnMuonSVM = "nMuonSVM_" + skimlabel;
    	std::string blepSVM = "SVLepVeto_" + skimlabel;
    	std::string bnSelIsoEle = "nSelIsoEle_" + skimlabel;
    	std::string bnLooseMuons = "nLooseMuons_" + skimlabel;
    	std::string bisolept = "LeptonVeto_" + skimlabel;
	
    	std::string bglumass = "glueMass_" + skimlabel;
	
    	std::string bCutA = "CutA_" + skimlabel;
    	std::string bfinal = "final_" + skimlabel;

		std::cout << skimlabel << " :" << std::endl;
		std::cout << " PS " << round( batchVars[psfillwgt]/batchVars[bfillwgt] );
        std::cout << " | hSV " << round( batchVars[bnSVHad]/batchVars[psfillwgt] );
        std::cout << " : lSV " << round( batchVars[bnSVLep]/batchVars[psfillwgt] );
        std::cout << " | eSVm " << round( batchVars[bnEleSVM]/batchVars[psfillwgt] );
        std::cout << " : mSVm " << round( batchVars[bnMuonSVM]/batchVars[psfillwgt] );
        std::cout << " = lSVm " << round( batchVars[blepSVM]/batchVars[psfillwgt] );
        std::cout << " | mF " << round( batchVars[bmetfilter]/batchVars[psfillwgt] );
        std::cout << " : hF " << round( batchVars[bhaloFilter]/batchVars[psfillwgt] );
      	std::cout << " | eV " << round( batchVars[bnSelIsoEle]/batchVars[psfillwgt] );
        std::cout << " - mV " << round( batchVars[bnLooseMuons]/batchVars[psfillwgt] );
        std::cout << " = lV " << round( batchVars[bisolept]/batchVars[psfillwgt] );
        //std::cout << " | met+lSVm " << round( batchVars[bCutA]/batchVars[bfillwgt] );
        std::cout << " | met+lV " << round( batchVars[bfinal]/batchVars[psfillwgt] );
		
		std::cout << std::endl;		

    }//<<>>for( auto line : item.second )


}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

    std::cout << " title test : " << addstr(ht,"varTitle") << std::endl;


    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------


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
            
				//std::string infilenameJ = "rjr_skim_files/KUCMS_RJR_GIGI_ootmet_Skim_List.txt";
                std::string infilenameJ = "rjr_skim_files/KUCMS_RJR_SMS_ootmet_Skim_List.txt";
				//std::string infilenameBG = "rjr_skim_files/KUCMS_RJR_BG_ootmet_Skim_List.txt";
				std::string infilenameBG = "rjr_skim_files/KUCMS_RJR_BG_v28_ootmet_Skim_List.txt";

				std::string version = "_v24_";
				std::string sigtype = "llpana_v28nt";
				std::string ofnstart = "KUCMS_";

                std::string htitleBG = "BG_"+sigtype+version;
                //std::string htitleBG = "BG_sv_"+sigtype+version;
				std::string htitleJ = sigtype+version;

                HistMaker base;

                //float metcut = 150;  // cmet > 150
    			// rjr var cuts
    			//cutvc = vc; ax2nq
    			//cutvd = vd; nvsum
    			// ve   asmass
    			// vf vdiff
    
				//int nj = 1;
				//int np = 1;
                for( int np = 1; np < 2; np++ ){
                for( int nj = 1; nj < 2; nj++ ){

				//std::string subdir = "cf_" + std::to_string(np) + "pho_" + std::to_string(nj) + "jet/";
				std::string subdir = "";

				// 1-1
				//std::vector<float> m_vec{1000,2000,3000}; // L-M-T
				std::vector<float> m_vec{2000,2750,3500}; // L-M-T
				//std::vector<float> m_vec{0,0,0}; // L-M-T
                std::vector<float> r_vec{0.2,0.275,0.35}; // L-M-T
                ////std::vector<float> rv_vec{0.1,0.25,0.4}; // L-M-T
				//std::vector<float> r_vec{0,0,0};
				std::vector<float> rv_vec{0.0,0.15,0.3}; // L-M-T
				//std::vector<float> rv_vec{0,0,0};
				//std::string outdir = "cf_1-2-3_2-275-35_1-25-40/" + subdir;
                //std::string outdir = "cf_2-275-35_2-275-35_1-25-40/" + subdir;
                //std::string outdir = "cf_2-275-35_2-275-35_0-15-3/" + subdir;
                //std::string outdir = "cf_0_0_0/" + subdir;
                std::string outdir = "";
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

