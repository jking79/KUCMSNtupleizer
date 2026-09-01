#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

using std::string;
enum Syst{
	dn = -1,
	up = 1,
	nom = 0
};

struct EffHTBin {
    double htLow;
    double htHigh;
    double mu;
    double sigma;
};

using HTBins   = std::vector<EffHTBin>;
using SystMap = std::unordered_map<int, HTBins>;
using YearMap = std::unordered_map<int, SystMap>;
using EffMap  = std::unordered_map<std::string, YearMap>;

class METTriggerTool{
	public:
		METTriggerTool(string csvfile){
			_parseCSV(csvfile);
		}
		
		virtual ~METTriggerTool(){ };

		void SetTrigger(string trigger){
			try {
    			    // Throws std::out_of_range if trigger is not found
    			    YearMap value = _effparams.at(trigger); 
    			} catch (const std::out_of_range& e) {
    			    std::cout << "Trigger not found in efficiency map: " << e.what() << '\n';
			}
			_trigger = trigger;
		}
		double GetEffWeight(double met, bool data, int year, double ht, Syst sys = nom);

	private:
		void _parseCSV(const std::string& filename); 
		string _trigger;
		EffMap _effparams;
		Syst parseSyst(const std::string& syst)
		{
		    if (syst == "up")
		        return up;
		    else if (syst == "dn")
		        return dn;
		    else if (syst == "nom")
		        return nom;
		
		    throw std::runtime_error("Unknown systematic: " + syst);
		}

};
