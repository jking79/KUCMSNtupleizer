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

	bool hasGenInfo = true;

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");
	const std::string ofnending = "_RjrSkim_Hists.root";

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
        if( debug ) std::cout << "--  adding from: " << str << std::endl;
        if( str[0] == '#' ) continue;
        if( str == " " ) continue;
        if( str == "" ) continue;
        auto instrs = splitString( str, " " );
        std::string inpath = instrs[0];
        std::string intitle = instrs[1];
		std::string batch = instrs[2];
		glumass = ( instrs.size() > 3 ) ? std::stof( instrs[3] ) : 0;
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
  
		if( glumass < 0 ){ glumass = 0; scale = 0; }
 
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
    
    	Init( fInTree, hasGenInfo );
    	initHists(htitle+"_"+intitle);
    
        std::cout << "Filling Config Map." << std::endl;
    
        UInt_t          nEvents;
        UInt_t          nSelectedEvents;
        std::string     *sKey;
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
			if( scale == 0 ) fillwt = 1;
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
        //int loopCounter(1000);
        auto nEntries = fInTree->GetEntries();
        if(debug){ nEntries = 10; loopCounter = 1; }
        std::cout << "Proccessing " << nEntries << " entries." << std::endl;
        for (Long64_t centry = 0; centry < nEntries; centry++){
            if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
            if(debug) std::cout << "*****************************************************************************" << std::endl;
            long long entry = fInTree->LoadTree(centry);
    		if(debug) std::cout << " - getBranches " << std::endl;
    		getBranches( entry, hasGenInfo );
    		if(debug) std::cout << " - eventLoop " << std::endl;
    		eventLoop(entry, m_vec, r_vec, rv_vec );
        }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
      
		//endJobs();
 
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
    		//ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
    		//ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
    		//ebeeMapR[it]->Write(); delete ebeeMapR[it];
			delete ebeeMapP[it];
			delete ebeeMapT[it];
			delete ebeeMapR[it];
    	}}//<<>>for( int it = 0; it < nEBEEMaps; it++ )
    
        fOutFile->Close();

	}//<<>>while (std::getline(infile,str))

	//endBatchJobs();

	//std::string boutfilename = "Batch_" + outfilename;
    //std::string bintitle = skimlabel;
    //std::string bofnending = ofnending;
    //if(debug) std::cout << " - Creating batch output file " << std::endl;
    //std::string bfulloutfilename =  boutfilename + bintitle + bofnending;
    //TFile* fOutFile = new TFile( bfulloutfilename.c_str(), "RECREATE" );
    //fOutFile->cd();

	//  save stuff here ......
    //for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    //for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    //for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }



    //fOutFile->Close();

    std::cout << "histMaker : Thats all Folks!!" << std::endl;
}//<<>>void kucmsSkimmer
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------- event loop -------------------------------------------------------------------------
void HistMaker::eventLoop( Long64_t entry, std::vector<float> m_vec, std::vector<float> r_vec, std::vector<float> rv_vec ){

	//int cs = ( cutselection == 0 ) ? 0 : 1;
	int cs = 0;

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
    float evtgwt = evtGenWgt;
    //float evtgwt = 1;
	float cfbin = 1;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    //float segwt = (configInfo[configKey])["nTotEvts"];
    float segwt = (configInfo[configKey])["sumEvtWgt"];
    float fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
	if( scale == 0 ) fillwt = 1;
    if( DEBUG ) std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
	batchVars[bfillwgt] += fillwt;

    //nRjrPhos = nphos;
    //nRjrJets = njets;
    //hist1d[2]->Fill(7,fillwt);// first eventloop cutflow -> should be same as n-1 cutflow

	for( int i = 0; i < 1; i++ ){ // continue loop - can use continue to skip
	/////////////////////////////////////////////////////////////////////////

	//!!!!! nRjrPhotons Cut
    //if( ( nRjrPhos != 0 ) &&  ( rjrNPhotons->at(cs) != nRjrPhos ) ) continue;
	//if( ( nRjrPhos == 0 ) &&  ( ( SV_nHadronic < 1 ) && ( SV_nLeptonic < 1 ) ) ) continue;

	//if( DEBUG ) std::cout << "RJR cut at " << cutvc << " with value " << rjrX2GMean->at(0) << std::endl;

    if( DEBUG ) std::cout << "Met cut " << std::endl;
	//!!!!! MET Cut
    //auto metCPt = hypo(selCMetPx,selCMetPy);
    //auto metPt = hypo(selMetPx,selMetPy);
	//if( metPt < 150 ) continue;
    //hist1d[2]->Fill(12,fillwt);

	//if( DEBUG ) std::cout << "RJR jet cut " << std::endl;
	// Min # jets/side
	//if( (*rjrNJetsJa)[cs] < nRjrJets || (*rjrNJetsJb)[cs] < nRjrJets ) continue;

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

    bool metFlag = Flag_BadChargedCandidateFilter && Flag_BadPFMuonDzFilter && Flag_BadPFMuonFilter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter;
    metFlag = metFlag && Flag_ecalBadCalibFilter && Flag_eeBadScFilter && Flag_goodVertices && Flag_hfNoisyHitsFilter;
	//metFlag = metFlag && Flag_EcalDeadCellTriggerPrimitiveFilter; // not working 

    bool haloFilter = Flag_globalSuperTightHalo2016Filter;

	bool trigger = Trigger_PFMET120_PFMHT120_IDTight or Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight; 
	trigger = trigger or Trigger_PFMET120_PFMHT120_IDTight_PFHT60 or Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;

	bool hemVeto = Flag_hemVeto;


    if( DEBUG ) std::cout << "Batch fill " << std::endl;
	//if( !MetFilter ) batchVars[bmetfilter] += fillwt;
    //if( !HaloFilter ) batchVars[bhaloFilter] += fillwt;
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

	//if( !MetFilter && svlepmatv ) batchVars[bCutA] += fillwt;
	//if( !MetFilter && !hasLLepF ) batchVars[bfinal] += fillwt;

	// rjr var cuts
	//cutvc = vc; ax2nq
    //cutvd = vd; nvsum
	// ve	asmass
	// vf vdiff

    if( DEBUG ) std::cout << "RJR info " << std::endl;

    int nJets = nSelJets;
    float jht = 0;
    TLorentzVector sJet;
    for( int jidx = 0; jidx < nJets; jidx++ ){
		if( (*selJetNeEmEF)[jidx] > 0.9 ) continue;
		if( (*selJetChEmEF)[jidx] > 0.8 ) continue;
		if( std::abs((*selJetEta)[jidx]) > 1.4 ) continue;
        sJet.SetPtEtaPhiM( (*selJetPt)[jidx], (*selJetEta)[jidx], (*selJetPhi)[jidx], (*selJetMass)[jidx] );
        jht += (*selJetPt)[jidx];
    }//<<>>for( int jidx = 0; jidx < nJets; jidx++ )
    float bJetPhi = sJet.Phi();
	float bJetPt = sJet.Pt();

	int nPho = (*rjrNPhotons)[cs];
    int nAJets = (*rjrNJetsJa)[cs];
    int nBJets = (*rjrNJetsJb)[cs];
	int nVis = (*rjrNVisObjects)[cs];

	if( ( nPho + nAJets + nBJets ) != nVis ) std::cout << "Vis Object counting Wrong !!!!!!!!!!" << std::endl;

    float mr = (*rjr_Mr)[cs];
    float rv = (*rjr_Rv)[cs];
	float r = (*rjr_R)[cs];

	float pts = (*rjrPTS)[cs];
	float dphisi = (*rjrSdphiBV)[cs];

    float ms = (*rjr_Ms)[cs]; // pHs41; NEW !!!!!!!
    float msq = (*rjr_pHs21)[cs]; // Ms for sqsq?
    float mva = (*rjr_Mva)[cs]; // float m_MVa = Ja->GetMass();
    float mvb = (*rjr_Mvb)[cs]; // float m_MVb = Jb->GetMass();

    float rm = (*rjr_Rm)[cs]; // float m_MVab = std::sqrt(sq2(m_MVa)+sq2(m_MVb))/phs40;  ? good in less compressed 
    float rs = (*rjr_Rs)[cs]; // float Rs = (phxa11+phxb11)/pHs22;   NEW !!!!!!!!
    float rx = (*rjr_Rx)[cs]/std::sqrt(2); // float Rx = hypo(phxa11/phxa21,phxb11/phxb21); NEW !!!!!!!!
    float rx0a = (*rjr_pHxa10)[cs]/(*rjr_pHxa20)[cs]; //(*rjr_Rxa)[cs]; // float Rxa = phxa10/phxa20;
    float rx0b = (*rjr_pHxb10)[cs]/(*rjr_pHxb20)[cs]; //(*rjr_Rxb)[cs]; // float Rxb = phxb10/phxb20;
    float rxa = (*rjr_pHxa11)[cs]/(*rjr_pHxa21)[cs];
    float rxb = (*rjr_pHxb11)[cs]/(*rjr_pHxb21)[cs];


	float rx0 = std::sqrt( rx0a*rx0a + rx0b*rx0b)/std::sqrt(2); 
	float rxmin = std::min( rxb, rxa );

	float phs2040 = (*rjr_pHs20)[cs]/(*rjr_pHs40)[cs];
    float phs2141 = (*rjr_pHs21)[cs]/(*rjr_pHs41)[cs];
    float phts2040 = (*rjr_pHts20)[cs]/(*rjr_pHts40)[cs];
    float phts2141 = (*rjr_pHts21)[cs]/(*rjr_pHts41)[cs];

	float pho1pts = (*rjr_p1Pts)[cs];
	float pho1pxa = (*rjr_p1Ptxa)[cs];
    float pho1pxb = (*rjr_p1Ptxb)[cs];
    float pho2pts = (*rjr_p2Pts)[cs];
    float pho2pxa = (*rjr_p2Ptxa)[cs];
    float pho2pxb = (*rjr_p2Ptxb)[cs];

	float pho1PtsH41 = pho1pts/(*rjr_pHs41)[cs];
    float pho2PtsH41 = pho2pts/(*rjr_pHs41)[cs];
    float pho1PtsH21 = pho1pts/(*rjr_pHs21)[cs];
    float pho2PtsH21 = pho2pts/(*rjr_pHs21)[cs];
	float minPhoPtsH41 = ( pho2PtsH41 > 0 ) ? std::min( pho1PtsH41, pho2PtsH41 ) : pho1PtsH41;
    float minPhoPtsH21 = ( pho2PtsH21 > 0 ) ? std::min( pho1PtsH21, pho2PtsH21 ) : pho1PtsH21;
	float sumPhoPtsH41 = ( pho1PtsH41 + pho2PtsH41 );
    float sumPhoPtsH21 = ( pho1PtsH21 + pho2PtsH21 );

	float pho1PxaH21 = pho1pxa/(*rjr_pHtxa21)[cs];
    float pho2PxaH21 = pho2pxa/(*rjr_pHtxa21)[cs];
    float pho1PxbH21 = pho1pxb/(*rjr_pHtxb21)[cs];
    float pho2PxbH21 = pho2pxb/(*rjr_pHtxb21)[cs];
	float pho1PxH21 = pho1PxaH21 + pho1PxbH21;
    float pho2PxH21 = pho2PxaH21 + pho2PxbH21;
	float minPhoPxH21 = ( pho2PxH21 > 0 ) ? std::min( pho1PxH21, pho2PxH21 ) : pho1PxH21;

    //float NormMBetaEql = (*selPhoMBetaEql)[0]/mr;
    //float NormMBetaPmt = (*selPhoMBetaPmt)[0]/mr;

	//(*rjrASMass)[cs] -- M
	//(*rjrAX2NQSum)[cs] -- R
	//rjrNVSum -- Rv


    fillwt = fillwt * m_vec[0];

	//var hist fill
	//bool timeCtrReg( ( gtime < -0.5 ) || ( sgtime < -0.5 ) );
    //bool p1tv( gtime < -25.0 || gtime > 25.0 );
    //bool p2tv( sgtime < -25.0 || sgtime > 25.0 );

	//bool rxabcust( rxa == 1 || rxb == 1 );
    //bool rxab7cust( rxa == 1 || ( rxb == 1 && rxa < 0.7 ) );
	bool jab1cust( nAJets < 1 && nBJets < 1 );
    bool jab2cust( nAJets < 2 && nBJets < 2 );

	bool nphocust( rjrNPhotons->at(cs) < 1 );
    bool n1phocust( rjrNPhotons->at(cs) != 1 );
    bool n2phocust( rjrNPhotons->at(cs) < 2 );
	bool timecust = false;
	//if( not nphocust && p1tv ) timecust = true;
    //if( not n2phocust && p2tv ) timecust = true;

	if( nphocust ) continue;
	if( timecust )  continue;

	if( not metFlag ) continue;
    if( not haloFilter ) continue;
	//if( not trigger ) continue;
    //if( hemVeto ) continue;

	// GGG cut sets
	//if( nRjrPhos < 20 && nphocust ) continue;  
    if( nRjrPhos == 0 ){
        //if( ms < 1000 ) continue;
        //if( rs < 0.25 ) continue;
	}//<<>>if( nRjrPhos == 0 )
	if( nRjrPhos == 1 ){
		//if( pho1PtsH21 < 0.1  ) continue;
		if( ms > 1000 ) continue;
        //if( rs > 0.25 ) continue;
	}//<<>>if( nRjrPhos == 1 )
    if( nRjrPhos == 2 ){
        if( ms < 1200 ) continue;
    }//<<>>if( nRjrPhos == 2 )
	if( nRjrPhos == 3 ){
		if( ms < 2000 ) continue;
	}//<<>>if( nRjrPhos == 3 )
    if( nRjrPhos == 4 ){
        if( ms < 3000 ) continue;
    }//<<>>if( nRjrPhos == 4 )
    if( nRjrPhos == 5 ){
        if( rs < 0.3 ) continue;
    }//<<>>if( nRjrPhos == 5 )
    if( nRjrPhos == 6 ){
        if( rs < 0.35 ) continue;
    }//<<>>if( nRjrPhos == 6 )
    if( nRjrPhos == 7 ){
        if( rs < 0.45 ) continue;
    }//<<>>if( nRjrPhos == 7 )
    if( nRjrPhos == 8 ){
        if( ms < 3000 ) continue;
        if( rs < 0.30 ) continue;
    }//<<>>if( nRjrPhos == 8 )
    if( nRjrPhos == 9 ){
        if( ms < 2000 ) continue;
        if( rs < 0.35 ) continue;
    }//<<>>if( nRjrPhos == 9 )
    if( nRjrPhos == 10 ){
        if( ms < 1200 ) continue;
        if( rs < 0.45 ) continue;
    }//<<>>if( nRjrPhos == 10 )


    if( ( nRjrPhos > 19 && nRjrPhos < 30 ) && ( nphocust ) ) continue;
    if( nRjrPhos == 25 ){
        if( ms < 1000 ) continue;
    }//<<>>if( nRjrPhos == 4 )
    if( nRjrPhos == 26 ){
        if( ms < 2200 ) continue;
    }//<<>>if( nRjrPhos == 4 )
    if( nRjrPhos == 27 ){
        if( ms < 3400 ) continue;
    }//<<>>if( nRjrPhos == 4 )


	//if( nphocust ) continue;
    if( nRjrPhos == 30 ){
        if( ms < 1000 ) continue;
        if( rs < 0.25 ) continue;
    }//<<>>if( nRjrPhos == 10 )
    if( nRjrPhos == 31 ){
        if( ms < 3000 ) continue;
        if( rs < 0.30 ) continue;
		//if( jabcust ) continue;
    }//<<>>if( nRjrPhos == 11 )
    if( nRjrPhos == 32 ){
    }//<<>>if( nRjrPhos == 12 )
    if( nRjrPhos == 33 ){
    }//<<>>if( nRjrPhos == 13 )


	hist1d[0]->Fill(mr,fillwt);
    hist1d[1]->Fill(ms,fillwt);
    hist1d[2]->Fill(mva,fillwt);
    hist1d[3]->Fill(mvb,fillwt);
    hist1d[8]->Fill(msq,fillwt);

    hist1d[10]->Fill(r,fillwt);
    hist1d[11]->Fill(rm,fillwt);
    hist1d[12]->Fill(rs,fillwt);
    hist1d[13]->Fill(rv,fillwt);
    hist1d[14]->Fill(rx,fillwt);
    hist1d[15]->Fill(rxa,fillwt);
    hist1d[16]->Fill(rxb,fillwt);
	if( rxb == 1 ) hist1d[17]->Fill(rxa,fillwt);

    hist2d[0]->Fill(ms,rs,fillwt);
    hist2d[1]->Fill(ms,rx,fillwt);
    hist2d[2]->Fill(ms,rm,fillwt);
    hist2d[3]->Fill(mva,mvb,fillwt);
    hist2d[4]->Fill(dphisi,pts,fillwt); //= new TH2D("dPhiSI_v_PtS"

    hist2d[10]->Fill(rs,rx,fillwt);
    hist2d[11]->Fill(rs,rm,fillwt);
    hist2d[12]->Fill(rx,rm,fillwt);
	hist2d[13]->Fill(rxa,rxb,fillwt);

    hist1d[18]->Fill(rx0,fillwt); // = new TH1D("Rx0", addstr(ht,"Rx0").c_str(), 120, 0, 1.2);
    hist1d[19]->Fill(rxmin,fillwt); // = new TH1D("Rxmin", addstr(ht,"Rxmin").c_str(), 120, 0, 1.2);

    hist1d[20]->Fill(rx0a,fillwt); //  = new TH1D("Rx0a", addstr(ht,"Rx0a").c_str(), 120, 0, 1.2);
    hist1d[21]->Fill(rx0b,fillwt); //  = new TH1D("Rx0b", addstr(ht,"Rx0b").c_str(), 120, 0, 1.2);

    hist1d[22]->Fill(phs2040,fillwt); //  = new TH1D("phs2040", addstr(ht,"phs2040").c_str(), 120, 0, 1.2);
    hist1d[23]->Fill(phs2141,fillwt); //  = new TH1D("phts2141", addstr(ht,"phts2141").c_str(), 120, 0, 1.2);
    hist1d[24]->Fill(phts2040,fillwt); //  = new TH1D("phs2040", addstr(ht,"phs2040").c_str(), 120, 0, 1.2);
    hist1d[25]->Fill(phts2141,fillwt); //  = new TH1D("phts2141", addstr(ht,"phts2141").c_str(), 120, 0, 1.2);

    hist2d[14]->Fill(rx0a,rx0b,fillwt);
    hist2d[15]->Fill(nAJets,nBJets,fillwt);
    hist2d[16]->Fill(nAJets,rxa,fillwt);
    hist2d[17]->Fill(nBJets,rxb,fillwt);

	//if( not( ML && RL && RvM ) ) continue;    
	
	hist1d[31]->Fill(pho1PtsH41,fillwt);
    hist1d[32]->Fill(pho2PtsH41,fillwt);
    hist1d[33]->Fill(pho1PtsH21,fillwt);
    hist1d[34]->Fill(pho2PtsH21,fillwt);
    hist1d[35]->Fill(minPhoPtsH41,fillwt);
    hist1d[36]->Fill(minPhoPtsH21,fillwt);
    hist1d[40]->Fill(sumPhoPtsH41,fillwt);
    hist1d[41]->Fill(sumPhoPtsH21,fillwt);

    hist1d[37]->Fill(pho1PxH21,fillwt);
    hist1d[38]->Fill(pho2PxH21,fillwt);
    hist1d[39]->Fill(minPhoPxH21,fillwt);

    hist2d[20]->Fill(minPhoPtsH41,ms,fillwt);
    hist2d[21]->Fill(minPhoPtsH41,rs,fillwt);
    hist2d[22]->Fill(sumPhoPtsH41,ms,fillwt);
    hist2d[23]->Fill(sumPhoPtsH41,rs,fillwt);

    hist2d[24]->Fill(pho1PxH21,pho2PxH21,fillwt);
    hist2d[25]->Fill(pho1pxa+pho1pxb,pho2pxa+pho2pxb,fillwt);


	// Time Sig

    float gtime = -99;
    float gpt = -99;
    float gphi = -99;
    float geta = -99;
    float gstsig = -99;
    float gltsig = -99;
    float gwtsig = -99;
    float gstime = -99;
    float gltime = -99;
    float gwtime = -99;
    float gstres = -99;
    float gltres = -99;
    float gwtres = -99;
    float genergy = -99;
    float gcenergy = -99;
    float glsx = 99;
    float gssx = 99;
    float gcrn = 99;
    float gr9 = 99;
    float gs4 = 99;
    float gsmj = 99;
    float gsmn = 99;
    float gsieie = 99;
    bool ggs = false;

    float sr2 = std::sqrt(1);
    float sgtime = 99;
    float sgpt = -99;
    if( selPhoTime->size() > 0 ){

        gtime = (*selPhoTime)[0];
        gpt = (*selPhoPt)[0];
        gphi = (*selPhoPhi)[0];
        geta = (*selPhoEta)[0];
        gstsig = (*selPhoSTimeSig)[0]/sr2;
        gwtsig = (*selPhoWTimeSig)[0]/sr2;
        gltsig = (*selPhoLTimeSig)[0]/sr2;
        gstime = (*selPhoSTime)[0];
        gwtime = (*selPhoWTime)[0];
        gltime = (*selPhoLTime)[0];
        gstres = (*selPhoSTimeRes)[0]/sr2;
        gltres = (*selPhoWTimeRes)[0]/sr2;
        gwtres = (*selPhoLTimeRes)[0]/sr2;
        genergy = (*selPhoEnergy)[0];
        gcenergy = (*selPhoCorEnergy)[0];

        glsx = (*selPhoLSCross)[0];
        gssx = (*selPhoSSCross)[0];
        ggs = (*selPhoShasGS)[0];

        gcrn = (*selPhoClstrRn)[0];
        gr9 = (*selPhoR9)[0];
        gs4 = (*selPhoS4)[0];
        gsmj = (*selPhoSMaj)[0];
        gsmn = (*selPhoSMin)[0];
        gsieie = (*selPhoSigmaIEtaIEta)[0];


    }//<<>>if( selPhoTime->size() > 0 )

    bool isPJphi = dPhi( bJetPhi, gphi ) >= ( PI - 0.3 );
    bool pjeb = ( std::min( bJetPt, gpt ) / std::max(  bJetPt, gpt ) ) > 0.6;

    if( std::abs(geta) > 1.4 ) continue;
    if( not isPJphi || not pjeb ) continue;

	//if( not ggs ) continue;

    hist1d[4]->Fill(gtime,fillwt);
    hist1d[5]->Fill(sgtime,fillwt);
    hist1d[6]->Fill(gpt,fillwt);
    hist1d[7]->Fill(sgpt,fillwt);

    hist1d[51]->Fill(gstsig,fillwt);
    hist1d[52]->Fill(gwtsig,fillwt);
    hist1d[53]->Fill(gltsig,fillwt);
    hist1d[54]->Fill(gwtime,fillwt);
    hist1d[55]->Fill(genergy,fillwt);
    hist1d[56]->Fill(gcenergy,fillwt);

    hist2d[40]->Fill(genergy,gcenergy,fillwt);
    hist2d[41]->Fill(gwtsig,gcenergy,fillwt);

    hist2d[51]->Fill(gwtsig,genergy,fillwt);
    hist2d[52]->Fill(gwtsig,gstime,fillwt);
    hist2d[53]->Fill(gwtsig,gltime,fillwt);
    hist2d[54]->Fill(gwtsig,gwtime,fillwt);

    hist2d[56]->Fill(gwtsig,glsx,fillwt);
    hist2d[57]->Fill(gwtsig,gcrn,fillwt);
    hist2d[58]->Fill(gwtsig,gr9,fillwt);
    hist2d[59]->Fill(gwtsig,gs4,fillwt);
    hist2d[60]->Fill(gwtsig,gsmj,fillwt);
    hist2d[61]->Fill(gwtsig,gsmn,fillwt);
    hist2d[62]->Fill(gwtsig,gsieie,fillwt);
    hist2d[63]->Fill(gwtsig,gssx,fillwt);

    hist2d[80]->Fill(gwtime,genergy,fillwt);
    hist2d[81]->Fill(gwtime,glsx,fillwt);
    hist2d[82]->Fill(gwtime,gstime,fillwt);
    hist2d[83]->Fill(gwtime,gltime,fillwt);
    hist2d[84]->Fill(gwtime,gcrn,fillwt);
    hist2d[85]->Fill(gwtime,gr9,fillwt);
    hist2d[86]->Fill(gwtime,gs4,fillwt);
    hist2d[87]->Fill(gwtime,gsmj,fillwt);
    hist2d[88]->Fill(gwtime,gsmn,fillwt);
    hist2d[89]->Fill(gwtime,gsieie,fillwt);
    hist2d[90]->Fill(gwtime,gssx,fillwt);

    hist2d[100]->Fill(gltime,genergy,fillwt);
    hist2d[101]->Fill(gltime,glsx,fillwt);
    hist2d[102]->Fill(gltime,gstime,fillwt);
    hist2d[103]->Fill(gltime,gltime,fillwt);
    hist2d[104]->Fill(gltime,gcrn,fillwt);
    hist2d[105]->Fill(gltime,gr9,fillwt);
    hist2d[106]->Fill(gltime,gs4,fillwt);
    hist2d[107]->Fill(gltime,gsmj,fillwt);
    hist2d[108]->Fill(gltime,gsmn,fillwt);
    hist2d[109]->Fill(gltime,gsieie,fillwt);
    hist2d[110]->Fill(gltime,gssx,fillwt);

    hist2d[200]->Fill(gstime,genergy,fillwt);
    hist2d[201]->Fill(gstime,glsx,fillwt);
    hist2d[202]->Fill(gstime,gstime,fillwt);
    hist2d[203]->Fill(gstime,gltime,fillwt);
    hist2d[204]->Fill(gstime,gcrn,fillwt);
    hist2d[205]->Fill(gstime,gr9,fillwt);
    hist2d[206]->Fill(gstime,gs4,fillwt);
    hist2d[207]->Fill(gstime,gsmj,fillwt);
    hist2d[208]->Fill(gstime,gsmn,fillwt);
    hist2d[209]->Fill(gstime,gsieie,fillwt);
    hist2d[210]->Fill(gstime,gssx,fillwt);


	}//<<>>for( int i = 0; i < 1; i++ ) -- continue loop

}//<<>>void HistMaker::eventLoop(` Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

    for( int it = 10; it < 30; it++ ){ if(hist1d[it]) smoothTH1D_v2(hist1d[it]);}

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

    hist1d[0] = new TH1D("Mr", addstr(ht," Mr;Mr [GeV]").c_str(), 120, 0, 12000);
    hist1d[1] = new TH1D("Ms", addstr(ht," Ms;Ms [GeV]").c_str(), 120, 0, 12000);
    hist1d[2] = new TH1D("Mva", addstr(ht," Mva;Mva [GeV]").c_str(), 80, 0, 4000);
    hist1d[3] = new TH1D("Mvb", addstr(ht," Mvb;Mvb [GeV]").c_str(), 80, 0, 4000);
    hist1d[4] = new TH1D("photime", addstr(ht," Lead Photon Time;Lead Photon Time [ns]").c_str(), 800, -2.0, 8.0);
    hist1d[5] = new TH1D("sphotime", addstr(ht," Second Photon Time;Second Photon Time [ns]").c_str(), 250, -25.0, 25.0);
    hist1d[6] = new TH1D("phopt", addstr(ht," Lead Photon Pt;Lead Photon Pt [GeV]").c_str(), 800, 0, 4000);
    hist1d[7] = new TH1D("sphopt", addstr(ht," Second Photon Pt;Second Photon Pt [GeV]").c_str(), 800, 0, 4000);
    hist1d[8] = new TH1D("Msq", addstr(ht," Msq;Msq [GeV]").c_str(), 120, 0, 12000);

    hist1d[10] = new TH1D("R", addstr(ht," R;R").c_str(), 120, 0, 1.2);
    hist1d[13] = new TH1D("Rv", addstr(ht," Rv'Rv").c_str(), 120, 0, 1.2);

    hist1d[11] = new TH1D("Rm", addstr(ht," Rm;Rm").c_str(), 120, 0, 1.2);
    hist1d[12] = new TH1D("Rs", addstr(ht," Rs;Rs").c_str(), 120, 0, 1.2);
    hist1d[14] = new TH1D("Rx", addstr(ht," Rx;Rx").c_str(), 120, 0, 1.2);

    hist1d[15] = new TH1D("Rxa", addstr(ht," Rxa;Rxa").c_str(), 120, 0, 1.2);
    hist1d[16] = new TH1D("Rxb", addstr(ht," Rxb;Rxb").c_str(), 120, 0, 1.2);
    hist1d[17] = new TH1D("Rxa_Rxb1", addstr(ht," Rxa : Rxb = 1;Rxa").c_str(), 120, 0, 1.2);
	
    hist1d[18] = new TH1D("Rx0", addstr(ht," Rx0;Rx0").c_str(), 120, 0, 1.2);
    hist1d[19] = new TH1D("Rxmin", addstr(ht," Rxmin;Rxmin").c_str(), 120, 0, 1.2);

    hist1d[20] = new TH1D("Rx0a", addstr(ht," Rx0a;Rx0a").c_str(), 120, 0, 1.2);
    hist1d[21] = new TH1D("Rx0b", addstr(ht," Rx0b;Rx0b").c_str(), 120, 0, 1.2);

    hist1d[22] = new TH1D("phs2040", addstr(ht," phs2040;phs2040").c_str(), 120, 0, 1.2);
    hist1d[23] = new TH1D("phs2141", addstr(ht," phts2141;phts2141").c_str(), 120, 0, 1.2);
    hist1d[24] = new TH1D("phts2040", addstr(ht," phs2040;phs2040").c_str(), 120, 0, 1.2);
    hist1d[25] = new TH1D("phts2141", addstr(ht," phts2141;phts2141").c_str(), 120, 0, 1.2);



    hist1d[31] = new TH1D("pho1PtsH41", addstr(ht," pho1PtsH41;pho1PtsH41").c_str(), 200, 0, 2.0);
    hist1d[32] = new TH1D("pho2PtsH41", addstr(ht," pho2PtsH41;pho2PtsH41").c_str(), 200, 0, 2.0);
    hist1d[33] = new TH1D("pho1PtsH21", addstr(ht," pho1PtsH21;pho1PtsH21").c_str(), 200, 0, 2.0);
    hist1d[34] = new TH1D("pho2PtsH21", addstr(ht," pho2PtsH21;pho2PtsH21").c_str(), 200, 0, 2.0);
    hist1d[35] = new TH1D("minPhoPtsH41", addstr(ht," minPhoPtsH41;minPhoPtsH41").c_str(), 200, 0, 2.0);
    hist1d[36] = new TH1D("minPhoPtsH21", addstr(ht," minPhoPtsH21;minPhoPtsH21").c_str(), 200, 0, 2.0);

    hist1d[37] = new TH1D("pho1PtxH21", addstr(ht," pho1PtxH21;pho1PtxH21").c_str(), 200, 0, 2.0);
    hist1d[38] = new TH1D("pho2PtxH21", addstr(ht," pho2PtxH21;pho2PtxH21").c_str(), 200, 0, 2.0);
    hist1d[39] = new TH1D("minPhoPtxH21", addstr(ht," minPhoPtxH21;minPhoPtxH21").c_str(), 200, 0, 2.0);

    hist1d[40] = new TH1D("sumPhoPtsH41", addstr(ht," sumPhoPtsH41;sumPhoPtsH41").c_str(), 200, 0, 2.0);
    hist1d[41] = new TH1D("sumPhoPtsH21", addstr(ht," sumPhoPtsH21;sumPhoPtsH21").c_str(), 200, 0, 2.0);


	hist1d[51] = new TH1D("selPhoSTimeSig", addstr(ht," selPhoSTimeSig;selPhoSTimeSig").c_str(), 1200, -60, 60);
    hist1d[52] = new TH1D("selPhoWTimeSig", addstr(ht," selPhoWTimeSig;selPhoWTimeSig").c_str(), 1200, -60, 60);
    hist1d[53] = new TH1D("selPhoLTimeSig", addstr(ht," selPhoLTimeSig;selPhoLTimeSig").c_str(), 1200, -60, 60);
    hist1d[54] = new TH1D("selPhoTOFTime", addstr(ht," selPhoTOFTime;selPhoTOFTime").c_str(), 500, -25, 25);
    hist1d[55] = new TH1D("selPhoEnergy", addstr(ht," selPhoEnergy;selPhoEnergy").c_str(), 300, 0, 3000);
    hist1d[56] = new TH1D("selPhoCEnergy", addstr(ht," selPhoCEnergy;selPhoCEnergy").c_str(), 300, 0, 3000);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    hist2d[0] = new TH2D("Ms_v_Rs", addstr(ht," Ms_v_Rs;M_{s} [GeV];R_{s}").c_str(), 120, 0, 12000, 120, 0, 1.2 );
    hist2d[1] = new TH2D("Ms_v_Rx", addstr(ht," Ms_v_Rx;M_{s} [GeV];R_{x}").c_str(), 120, 0, 12000, 120, 0, 1.2 );
    hist2d[2] = new TH2D("Ms_v_Rm", addstr(ht," Ms_v_Rm;M_{s} [GeV];R_{m}").c_str(), 120, 0, 12000, 120, 0, 1.2 );
    hist2d[3] = new TH2D("Mva_v_Mvb", addstr(ht," Mva_v_Mvb;M_{va} [GeV];M_{vb} [GeV]").c_str(), 80, 0, 4000, 80, 0, 4000 );
    hist2d[4] = new TH2D("dPhiSI_v_PtS", addstr(ht," dPhiSI_v_PtS;sPhi_{SI};P_{t}^{S}").c_str(), 320, 0, 3.2, 120, 0, 1200 );

    hist2d[10] = new TH2D("Rs_v_Rx", addstr(ht," Rs_v_Rx;R_{s};R_{x}").c_str(), 120, 0, 1.2, 120, 0, 1.2 );
    hist2d[11] = new TH2D("Rs_v_Rm", addstr(ht," Rs_v_Rm;R_{s};R_{m}").c_str(), 120, 0, 1.2, 120, 0, 1.2 );
    hist2d[12] = new TH2D("Rx_v_Rm", addstr(ht," Rx_v_Rm;R_{x};R_{m}").c_str(), 120, 0, 1.2, 120, 0, 1.2 );
    hist2d[13] = new TH2D("Rxa_v_Rxb", addstr(ht," Rxa_v_Rxb;R_{xa};R_{xb}").c_str(), 120, 0, 1.2, 120, 0, 1.2 );
    hist2d[14] = new TH2D("Rx0a_v_Rx0b", addstr(ht," Rx0a_v_Rx0b;R_{x0a};R_{x0b}").c_str(), 120, 0, 1.2, 120, 0, 1.2 );
    hist2d[15] = new TH2D("nAJets_v_nBJets", addstr(ht," nAJets_v_nBJets;nAJets;nBJets").c_str(), 20, 0, 20, 20, 0, 20 );
    hist2d[16] = new TH2D("nAJets_v_Rxa", addstr(ht," nAJets_v_Rxa;nAJets;Rxa").c_str(), 20, 0, 20, 120, 0, 1.2 );
    hist2d[17] = new TH2D("nBJets_v_Rxb", addstr(ht," nAJets_v_Rxb;nBJets;Rxb").c_str(), 20, 0, 20, 120, 0, 1.2 );

    hist2d[20] = new TH2D("gpts_v_ms", addstr(ht," gpts_v_ms;minPhoPtsH41;Ms [GeV]").c_str(), 200, 0, 2.0, 120, 0, 12000 );
    hist2d[21] = new TH2D("gpts_v_rs", addstr(ht," gpts_v_rs;minPhoPtsH41;Rs").c_str(), 200, 0, 2.0, 120, 0, 1.2 );
    hist2d[22] = new TH2D("gspts_v_ms", addstr(ht," gspts_v_ms;sumPhoPtsH41;Ms [GeV]").c_str(), 200, 0, 2.0, 120, 0, 12000 );
    hist2d[23] = new TH2D("gspts_v_rs", addstr(ht," gspts_v_rs;sumPhoPtsH41;Rs").c_str(), 200, 0, 2.0, 120, 0, 1.2 );

    hist2d[24] = new TH2D("pho1PxH21_v_pho2PxH21", addstr(ht," pho1PxH21_v_pho2PxH21;pho1PxH21;pho2PxH21").c_str(), 500, 0, 5.0, 500, 0, 5.0 );
    hist2d[25] = new TH2D("pho1px_v_pho2px", addstr(ht," pho1px_v_pho2px;pho1px;pho2px").c_str(), 200, 0, 2000, 200, 0, 2000 );


    hist2d[40] = new TH2D("Energy_v_CEnergy", addstr(ht," Energy_v_CEnergy;selPhoEnergy;selPhoCEnergy").c_str(), 200, 0, 2000, 200, 0, 2000);
    hist2d[41] = new TH2D("WSig_v_CEnergy", addstr(ht," WSig_v_CEnergy;selPhoSTimeSig;selPhoCEnergy").c_str(), 400, -80, 80, 200, 0, 2000);


    hist2d[51] = new TH2D("WSig_v_PhoEnergy", addstr(ht," WSig_v_Energy;selPhoWTimeSig;selPhoEnergy").c_str(), 400, -80, 80, 200, 0, 2000);
    hist2d[52] = new TH2D("WSig_v_STime", addstr(ht," WSig_v_STime;selPhoWTimeSig;selPhoSTime").c_str(), 400, -80, 80, 400, -20, 20);
    hist2d[53] = new TH2D("WSig_v_LTime", addstr(ht," WSig_v_LTime;selPhoWTimeSig;selPhoLTime").c_str(), 400, -80, 80, 400, -20, 20);
    hist2d[54] = new TH2D("WSig_v_WTime", addstr(ht," WSig_v_WTime;selPhoWTimeSig;selPhoWTime").c_str(), 400, -80, 80, 400, -20, 20);
    hist2d[56] = new TH2D("WSig_v_LSCross", addstr(ht," WSig_v_LSCross;selPhoWTimeSig;selPhoLSCross").c_str(), 400, -80, 80, 200, -1, 1);
    hist2d[57] = new TH2D("WSig_v_CRn", addstr(ht," WSig_v_VRn;selPhoWTimeSig;selPhoClstrRn").c_str(), 400, -80, 80, 100, 0, 1);
    hist2d[58] = new TH2D("WSig_v_R9", addstr(ht," WSig_v_R9;selPhoWTimeSig;selPhoR9").c_str(), 400, -80, 80, 100, 0, 1);
    hist2d[59] = new TH2D("WSig_v_S4", addstr(ht," WSig_v_S4;selPhoWTimeSig;selPhoS4").c_str(), 400, -80, 80, 100, 0, 1);
    hist2d[60] = new TH2D("WSig_v_Smj", addstr(ht," WSig_v_Smj;selPhoWTimeSig;selPhoSMaj").c_str(), 400, -80, 80, 500, 0, 5);
    hist2d[61] = new TH2D("WSig_v_Smn", addstr(ht," WSig_v_Smn;selPhoWTimeSig;selPhoSMin").c_str(), 400, -80, 80, 200, 0, 2);
    hist2d[62] = new TH2D("WSig_v_Sieie", addstr(ht," WSig_v_Sieie;selPhoWTimeSig;selPhoSigmaIEtaIEta").c_str(), 400, -80, 80, 100, 0, 0.025);
    hist2d[63] = new TH2D("WSig_v_SSCross", addstr(ht," WSig_v_SSCross;selPhoWTimeSig;selPhoSSCross").c_str(), 400, -80, 80, 200, -1, 1);

    hist2d[80] = new TH2D("WTime_v_PhoEnergy", addstr(ht," WTime_v_PhoEnergy;selPhoWTime;selPhoEnergy").c_str(), 400, -20, 20, 200, 0, 2000);
    hist2d[81] = new TH2D("WTime_v_LSCross", addstr(ht," WTime_v_LSCross;selPhoWTime;selPhoLSCross").c_str(), 400, -20, 20, 200, -1, 1);
    hist2d[82] = new TH2D("WTime_v_STime", addstr(ht," WTime_v_STime;selPhoWTime;selPhoSTime").c_str(), 400, -20, 20, 400, -20, 20);
    hist2d[83] = new TH2D("WTime_v_LTime", addstr(ht," WTime_v_LTime;selPhoWTime;selPhoLTime").c_str(), 400, -20, 20, 400, -20, 20);
    hist2d[84] = new TH2D("WTime_v_CRn", addstr(ht," WTime_v_VRn;selPhoWTime;selPhoClstrRn").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[85] = new TH2D("WTime_v_R9", addstr(ht," WTime_v_R9;selPhoWTime;selPhoR9").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[86] = new TH2D("WTime_v_S4", addstr(ht," WTime_v_S4;selPhoWTime;selPhoS4").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[87] = new TH2D("WTime_v_Smj", addstr(ht," WTime_v_Smj;selPhoWTime;selPhoSMaj").c_str(), 400, -20, 20, 500, 0, 5);
    hist2d[88] = new TH2D("WTime_v_Smn", addstr(ht," WTime_v_Smn;selPhoWTime;selPhoSMin").c_str(), 400, -20, 20, 200, 0, 2);
    hist2d[89] = new TH2D("WTime_v_Sieie", addstr(ht," WTime_v_Sieie;selPhoWTime;selPhoSigmaIEtaIEta").c_str(), 400, -20, 20, 100, 0, 0.025);
    hist2d[90] = new TH2D("WTime_v_SSCross", addstr(ht," WTime_v_SSCross;selPhoWTime;selPhoSSCross").c_str(), 400, -20, 20, 200, -1, 1);

    hist2d[100] = new TH2D("LTime_v_PhoEnergy", addstr(ht," LTime_v_PhoEnergy;selPhoLTime;selPhoEnergy").c_str(), 400, -20, 20, 200, 0, 2000);
    hist2d[101] = new TH2D("LTime_v_LSCross", addstr(ht," LTime_v_LSCross;selPhoLTime;selPhoLSCross").c_str(), 400, -20, 20, 200, -1, 1);
    hist2d[102] = new TH2D("LTime_v_STime", addstr(ht," LTime_v_STime;selPhoLTime;selPhoSTime").c_str(), 400, -20, 20, 400, -20, 20);
    hist2d[103] = new TH2D("LTime_v_LTime", addstr(ht," LTime_v_LTime;selPhoLTime;selPhoLTime").c_str(), 400, -20, 20, 400, -20, 20);
    hist2d[104] = new TH2D("LTime_v_CRn", addstr(ht," LTime_v_VRn;selPhoLTime;selPhoClstrRn").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[105] = new TH2D("LTime_v_R9", addstr(ht," LTime_v_R9;selPhoLTime;selPhoR9").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[106] = new TH2D("LTime_v_S4", addstr(ht," LTime_v_S4;selPhoLTime;selPhoS4").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[107] = new TH2D("LTime_v_Smj", addstr(ht," LTime_v_Smj;selPhoLTime;selPhoSMaj").c_str(), 400, -20, 20, 500, 0, 5);
    hist2d[108] = new TH2D("LTime_v_Smn", addstr(ht," LTime_v_Smn;selPhoLTime;selPhoSMin").c_str(), 400, -20, 20, 200, 0, 2);
    hist2d[109] = new TH2D("LTime_v_Sieie", addstr(ht," LTime_v_Sieie;selPhoLTime;selPhoSigmaIEtaIEta").c_str(), 400, -20, 20, 100, 0, 0.025);
    hist2d[110] = new TH2D("LTime_v_SSCross", addstr(ht," LTime_v_SSCross;selPhoLTime;selPhoSSCross").c_str(), 400, -20, 20, 200, -1, 1);

    hist2d[200] = new TH2D("STime_v_PhoEnergy", addstr(ht," STime_v_PhoEnergy;selPhoSTime;selPhoEnergy").c_str(), 400, -20, 20, 200, 0, 2000);
    hist2d[201] = new TH2D("STime_v_LSCross", addstr(ht," STime_v_LSCross;selPhoSTime;selPhoLSCross").c_str(), 400, -20, 20, 200, -1, 1);
    hist2d[202] = new TH2D("STime_v_STime", addstr(ht," STime_v_STime;selPhoSTime;selPhoSTime").c_str(), 400, -20, 20, 400, -20, 20);
    hist2d[203] = new TH2D("STime_v_LTime", addstr(ht," STime_v_LTime;selPhoSTime;selPhoLTime").c_str(), 400, -20, 20, 400, -20, 20);
    hist2d[204] = new TH2D("STime_v_CRn", addstr(ht," STime_v_VRn;selPhoSTime;selPhoClstrRn").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[205] = new TH2D("STime_v_R9", addstr(ht," STime_v_R9;selPhoSTime;selPhoR9").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[206] = new TH2D("STime_v_S4", addstr(ht," STime_v_S4;selPhoSTime;selPhoS4").c_str(), 400, -20, 20, 100, 0, 1);
    hist2d[207] = new TH2D("STime_v_Smj", addstr(ht," STime_v_Smj;selPhoSTime;selPhoSMaj").c_str(), 400, -20, 20, 500, 0, 5);
    hist2d[208] = new TH2D("STime_v_Smn", addstr(ht," STime_v_Smn;selPhoSTime;selPhoSMin").c_str(), 400, -20, 20, 200, 0, 2);
    hist2d[209] = new TH2D("STime_v_Sieie", addstr(ht," STime_v_Sieie;selPhoSTime;selPhoSigmaIEtaIEta").c_str(), 400, -20, 20, 100, 0, 0.025);
    hist2d[210] = new TH2D("STime_v_SSCross", addstr(ht," STime_v_SSCross;selPhoSTime;selPhoSSCross").c_str(), 400, -20, 20, 200, -1, 1);

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

}//<<>>void Hist`Maker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                //std::string listdir = "/uscms/home/jaking/nobackup/llpana_skims/";
				//std::string listdir = "/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer/tsig_skims/"; 
                std::string listdir = "/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer/";

				//std::string infilenameJ = "rjr_skim_files/KUCMS_RJR_GIGI_ootmet_Skim_List.txt";
                //std::string infilenameJ = "rjr_skim_files/KUCMS_RJR_SMS_ootmet_Skim_List.txt";
				//std::string infilenameBG = "rjr_skim_files/KUCMS_RJR_BG_ootmet_Skim_List.txt";
				//std::string infilenameBG = "rjr_skim_files/KUCMS_RJR_BG_v28_ootmet_Skim_List.txt";

                std::string infilenameJ = "rjr_skim_files/KUCMS_RJR_SMS_v39_Skim_List.txt";
                std::string infilenameBG = "rjr_skim_files/KUCMS_RJR_BG_v39_Skim_List.txt";
                std::string infilenameD = "rjr_skim_files/KUCMS_RJR_DATA_v40_Skim_List.txt";

				std::string version = "v40_";
				std::string sigtype = "llpana_";
				std::string ofnstart = "KUCMS_";

                std::string htitleBG = "BG_"+sigtype+version;
                //std::string htitleBG = "BG_sv_"+sigtype+version;
				std::string htitleJ = sigtype+version;
                std::string htitleD = "Data_"+sigtype+version;

                HistMaker base;

				// new rjr ken studies w/ H vats  -----------------------------------------------------------------------------

				//int nj = 1;
				//int np = 1; : 2,7,10
                for( int np = 0; np < 1; np++ ){
                for( int nj = 0; nj < 1; nj++ ){

				//std::string subdir = "cf_" + std::to_string(np) + "pho_" + std::to_string(nj) + "jet/";
				std::string subdir = "";

                std::vector<float> r_vec{1.0}; // BG scaling
				std::vector<float> rv_vec{1.0}; // Sig scaling
                std::string outdir = "";

                std::string isoline = "P1MfHa_RJR0_TSig_EGJ_"; // EGJ = EGamma + Jets selecton
                //std::string isoline = "P1TrMfHeHa_RJR0_"; // Tv = time valid, Mf = metfilters, Ha = halofilters, Tr = trigger, He = hemfilter;
				//std::string isoline = "P1TrMfHeHa_TSig_GS1_"; // Tv = time valid, Mf = metfilters, Ha = halofilters, Tr = trigger, He = hemfilter;
				isoline += "cv" + std::to_string( np ) + "_";
                std::string outfilenameJ = outdir + ofnstart + htitleJ + isoline;
				std::string htitlefullJ =  htitleJ + isoline;
				//base.histMaker( listdir, infilenameJ, outfilenameJ, htitlefullJ, np, nj, rv_vec, r_vec, rv_vec );
                std::string outfilenameBG = outdir + ofnstart + htitleBG + isoline;
                std::string htitlefullBG =  htitleBG + isoline;
                //base.histMaker( listdir, infilenameBG, outfilenameBG, htitlefullBG, np, nj, r_vec, r_vec, rv_vec );
                std::string outfilenameD = outdir + ofnstart + htitleD + isoline;
                std::string htitlefullD =  htitleD + isoline;
                base.histMaker( listdir, infilenameD, outfilenameD, htitlefullD, np, nj, r_vec, r_vec, rv_vec );

				}}


    return 1;

}//<<>>int main ( int argc, char *argv[] )

