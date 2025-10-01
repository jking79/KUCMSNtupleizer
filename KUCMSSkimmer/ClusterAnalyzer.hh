#ifndef CLUSTERANALYZER_HH
#define CLUSTERANALYZER_HH

#include "BayesianClustering/BayesianClustering.hh"
#include "fdeep/fdeep.hpp"


//wrapper class for clustering
class ClusterAnalyzer{
	public:
		ClusterAnalyzer();
		virtual ~ClusterAnalyzer(){
			delete _algo;
		};

		//setting options for clustering algorithm
		void SetTransferFactor(double g);
		//needs PV info for geometric corrections and correct momentum calculations of clustered elements
		void SetPV(double pvx, double pvy, double pvz){ _PV = BayesPoint({pvx, pvy, pvz});}

		//adding data to clustering algorithm
		void AddRecHit(double rhx, double rhy, double rhz, double rhE, double rht, double rhId);
		void ClearRecHitList();

		//TODO - may not want to declare detidmap here because it already exists in the skimmer
		std::map<UInt_t,DetIDStruct> _detIDmap;
		void SetDetIdMap(std::map<UInt_t,DetIDStruct> detidmap){ _detIDmap = detidmap; }
		//functions for calculating PV time, cluster observables, etc
		struct ClusterObj{
			Jet _jet;
			ClusterObj(Jet jet){ _jet = Jet; }

			vector<bool> _PUscores;
			//will fully remove PU clusters - ie set all w_nk = 0 for PU subcluster k
			void CleanOutPU(){ _jet.CleanOutPU(true, _PUscores);}
			
			//per-object information
			double _pvTime = -999; //time at PV
			double _detTime = -999;
			double GetObjTime_PV(){ return _pvTime; }
			double GetObjTime_Det(){ return _detTime; }
			void CalculateObjTimes(){ 
				//obj time = sum_i w_i * t_i / w_i
				double t = 0;
				double norm = 0;
				PointCollection locs;
				for(int k = 0; k < _jet.GetNConstituents(); k++){
					Jet subcl = _jet.GetConstituent(k);
					vector<JetPoint> rhs = subcl.GetJetPoints();	
					for(int r = 0; r < rhs.size(); r++){
						double w = rhs[r].GetWeight(); //weight = r_nk*isPU_k
						BayesPoint loc({rhs[r].x(), rhs[r].y(), rhs[r].z(), rhs[r].t()});
						loc.SetWeight(w);
						locs += loc;	
					}
				}
				_detTime = locs.Centroid(3);
				//include geo factor to PV for PV time
				BayesPoint center;
				center.SetValue(locs.Centroid(0),0);	
				center.SetValue(locs.Centroid(1),1);	
				center.SetValue(locs.Centroid(2),2);	
				
				double dx = rhs[r].x() - _PV.at(0);
				double dy = rhs[r].y() - _PV.at(1);
				double dz = rhs[r].z() - _PV.at(2);
				double d_rh_pv = sqrt(dx*dx + dy*dy + dz*dz)/_SOL;
				_pvTime = _detTime - d_rh_pv;
			}

			double _timeRes = -999;
			double GetObjTimeSig(){
				if(_timeRes == -999){
					cout << "Time resolution not calculated. Please run CalculateObjTimeSig(res), where res is a map of rh ids to their time resolution." << endl;
					return -999;
				} 
				if(_detTime == -999){
					cout << "Object time at detector not calculated. Please run CalculateObjTimes()." << endl;
					return -999;
				}
				return _detTime / _timeRes;
			} 
			void CalculateObjTimeSig(map<int, double> res){ //calculate time significance - needs to be called after CalculateObjTime
				//match resolutions via rh ids
				double restot = 0;
				double norm = 0;
				for(int k = 0; k < _jet.GetNConstituents(); k++){
					Jet subcl = _jet.GetConstituent(k);
					vector<JetPoint> rhs = subcl.GetJetPoints();	
					for(int r = 0; r < rhs.size(); r++){
						double w = rhs[r].GetWeight(); //weight = r_nk*isPU_k
						restot += w*res[rhs[r].GetId()];
						norm += w;
					}
				}
				restot /= norm;
				_timeRes = restot; 
			}


		
			//per-subcluster information - vector of length k for k subclusters
			void GetPUScores(vector<double>& scores){
				if(_PUscores.size() == 0){ cout << "PU scores not calculated. Please run CalculatePUSCores()" << endl; return; }
				scores.clear(); scores = _PUscores;
			}
			void CalculatePUScores(){
				if(_PUscores.size() == 0) _jet.CleanOutPU(true, _PUscores);
			} 
			
			vector<pair<int, double>> _detBkgScore; //is a vector of pairs s.t. _detBkgScore[k] = pair(max_class, max_score) 
			GetDetBkgScore(vector<pair<int,double>>& detBkgScores){ detBkgScores.clear(); detBkgScores = _detBkgScores; }
			fdeep::model _nnmodel = fdeep::load_model("json/small3CNN_EMultr.json");
			vector<string> _nnfeatures = {"Er"};
			int _ngrid = 7;
			void SetCNNModel(string model){ _nnmodel = fdeep::load_model(model); }
			void GetCenterXtal(JetPoint& center){
				//get center of pts in ieta, iphi -> max E point
				vector<JetPoint> rhs = _jet.GetJetPoints();
				double maxE = 0;
				for(int r = 0; r < rhs.size(); r++){
					if(rhs[r].E() > maxE){
						maxE = rhs[r].E();
						center = rhs[r];
					}
				}
			}
			//void MakeCNNInputGrid(BasePDFMixture* model, int k, vector<JetPoint>& rhs, JetPoint& center, map<string,double>& mapobs){
			void MakeCNNInputGrid(int k, JetPoint& center, map<string,double>& mapobs){
				mapobs.clear();
				map<pair<int,int>, vector<double>> grid;
				//make sure ngrid is odd to include center crystal
				if(_ngrid % 2 == 0)
					_ngrid++;

				int ngrid_boundary = (_ngrid-1)/2;
				//set default channel values to 0
				//{E, t, r, Er}
				for(int i = -ngrid_boundary; i < ngrid_boundary+1; i++)
					for(int j = -ngrid_boundary; j < ngrid_boundary+1; j++)
						grid[make_pair(i,j)] = {0., 0., 0.};


				//get ngrid x ngrid around center point 
				int ieta, iphi;
				int rh_ieta = _detIDmap[center.rhId()].i2;
				int rh_iphi = _detIDmap[center.rhId()].i1;
				int deta, dphi;
				Jet subcl = _jet.GetConstituent(k);
				vector<JetPoint> rhs = subcl.GetJetPoints();
				for(int j = 0; j < rhs.size(); j++){
					ieta = _detIDmap[rhs[j].rhId()].i2;
					iphi = _detIDmap[rhs[j].rhId()].i1;
					//do eta flip
					if(rh_ieta < 0)
						deta = -(ieta - rh_ieta);
					else
						deta = ieta - rh_ieta;
					dphi = iphi - rh_iphi; 
					//needs wraparound
					if(dphi > 180)
						dphi = 360 - dphi;
					else if(dphi < -180)
						dphi = -(360 + dphi);
					if(fabs(deta) <= ngrid_boundary && fabs(dphi) <= ngrid_boundary){
						//posterior is weighted s.t. sum_k post_nk = w_n = E*_gev, need to just have unweighted probs since E is already here
						//r_nk = post_nk/w_n s.t. sum_k (post_nk/w_n) = w_n/w_n = 1
						grid[make_pair(deta, dphi)] = {rhs[j].E()}; //subcluster jets already have rh energies "projected in" to their subclusters

					}
				}
				pair<int, int> icoords_grid;
				for(int i = -(_ngrid-1)/2.; i < (_ngrid-1)/2+1; i++){
					for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
						icoords_grid = make_pair(i,j);
						mapobs["CNNgrid_Er_cell"+to_string(i)+"_"+to_string(j)] = grid[icoords_grid][0];
					}
				}
			}
			void CalculateDetBkgScores(){
				int nchan = _nnfeatures.size();
				fdeep::tensor_shape tensor_shape(_ngrid, _ngrid, nchan);
				fdeep::tensor input_tensor(tensor_shape, 0.0f);
				//make grid for each subclusters in jet
				JetPoint center;
				GetCenterXtal(center);
				for(int k = 0; k < _jet.GetNConstituents(); k++){
					map<string, double> inputs;
					MakeCNNGrid(k, center, inputs);
					//transform grid to input_sample
					for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++){
						for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
							for(int c = 0; c < _nnfeatures.size(); c++){
								double val = obs["CNNgrid_"+_nnfeatures[c]+"_cell"+to_string(i)+"_"+to_string(j)];
								input_tensor.set(fdeep::tensor_pos(i+(_ngrid-1)/2, j+(_ngrid-1)/2, c), val);
							}
						}							
					}
					//predict_class returns predicted class number and value of max output neuron
					pair<size_t, double> result = _nnmodel.predict_class_with_confidence({input_tensor});
					_detBkgScores.push_back(make_pair((int)result.first, result.second))
				}
			}







				

		}



	private:
		BayesCluster* _algo;
		vector<Jet> _rhs;
		double _gev; 
		BayesPoint _PV;
		double _radius; //detector radius in m
		double _SOL = 29.9792458; //in cm/ns

		void _treesToObjs(vector<node*>& trees, vector<ClusterObj>& objs);




};
#endif
