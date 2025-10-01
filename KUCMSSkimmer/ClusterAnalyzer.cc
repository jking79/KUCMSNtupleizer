#include "ClusterAnalyzer.hh"

ClusterAnalyzer::ClusterAnalyzer(){
	_algo = nullptr;
	_gev = 1;
	_radius = 1.29;
}

//sets transfer factor for energy weighting in clustering
void ClusterAnalyzer::SetTransferFactor(double g){ _gev = g; }

//add rechit to list of rechits to be clustered
//can be run in a for-loop, looping over all rechits
void ClusterAnalyzer::AddRecHit(double rhx, double rhy, double rhz, double rhE, double rht, int rhId){
	//TODO - maybe set the (0,0,0) point with a separate function s.t. it can be set from beamspot
	//take center point to be (0,0,0) for now to account for TOF to detector
	double dx = rhx
	double dy = rhy
	double dz = rhz
	double d_rh = sqrt(dx*dx + dy*dy + dz*dz)/_SOL;
	rht += d_rh;	

	JetPoint rh(rhx, rhy, rhz, rht);
	rh.SetEnergy(rhE);
	rh.SetWeight(rhE*_gev);
	
	Jet jet(rh);
	_rhs.push_back(jet);
}

void ClusterAnalyzer::ClearRecHitList(){
	_rhs.clear();
}


//should be run after all rechits for clustering have been added
ClusterObj ClusterAnalyzer::RunClustering(){
	_algo = new BayesCluster(_rhs);	
	
	//hard coding parameters that won't change
	double cell = acos(-1)/180;
	//time resolution parameters - set by detector time resolution
	double tresCte = 0.1727;//times given in ns//0.133913 * 1e-9;
	double tresStoch = 0.5109;//1.60666 * 1e-9; 
	double tresNoise = 2.106;//0.00691415 * 1e-9;
	_algo->SetMeasErrParams(cell, tresCte, tresStoch*_gev, tresNoise*_gev); 
	_//set time resolution smearing
	_//if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n);
	double thresh = 0.1;
	_algo->SetThresh(thresh);
	double alpha = 1e-300;
	_algo->SetAlpha(alpha);
	double emAlpha = 1e-5;
	_algo->SetSubclusterAlpha(_emAlpha);
	_algo->SetVerbosity(0);
	map<string, Matrix> prior_params;
	//beta
	prior_params["scale"] = Matrix(1e-3);
	//nu
	prior_params["dof"] = Matrix(3);
	//W
	Matrix W(3,3);
	W.SetEntry(0.013333,0,0);
	W.SetEntry(0.013333,1,1);
	W.SetEntry(33.33333,2,2);
	prior_params["scalemat"] = W;
	//m
	prior_params["mean"] = Matrix(3,1);
	_algo->SetPriorParameters(_prior_params);


	//do hierarchical clustering for subcluster constraints
	vector<node*> trees = _algo->NlnNCluster();
	vector<Jet> objs;
	_treesToObjs(trees, objs);	
	return objs[0];

}


void ClusterAnalyzer::_treesToObjs(vector<node*>& trees, vector<ClusterObj>& objs){
	objs.clear();
	double x, y, z, eta, phi, t, theta, px, py, pz;
	int njets_tot = 0;
	for(int i = 0; i < trees.size(); i++){
		//get points from tree
		PointCollection* pc = trees[i]->points;
		//at least 2 points (rhs)
		if(pc->GetNPoints() < 2) continue;
		rhs.clear();
		Jet predJet(trees[i]->model, _PV, _gev, _radius);
		//add Jet to jets	
		objs.push_back(ClusterObj(predJet));	
	}
}
