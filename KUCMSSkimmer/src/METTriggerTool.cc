#include "METTriggerTool.hh"
#include <cmath>

double METTriggerTool::GetEffWeight(double met, bool data, int year, double ht, Syst sys){
	if(_trigger == ""){
        	std::cerr << "Error: trigger not defined. Use SetTrigger(trigname) before getting weight." << std::endl;
        	return -1;
	}
	if(data)
		return 1.;
	const auto& bins = _effparams
	    .at(_trigger)
	    .at(year)
	    .at(sys);
	for (const auto& bin : bins) {
	    if (ht >= bin.htLow && ht < bin.htHigh) {
		//std::cout << "mu " << bin.mu << " sigma " << bin.sigma << " x " << met << std::endl;
		return 0.5 * (1.0 + std::erf((met - bin.mu) / (bin.sigma * std::sqrt(2.0))));
	    }
	}
	std::cerr << "HT " << ht << " not defined in bins for trigger " << _trigger << " year " << year << " sys " << sys << std::endl;
	return -1;
}

// Function to parse the CSV file
void METTriggerTool::_parseCSV(const std::string& filename) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filename);

    // Verify if the file opened successfully
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return;
    }

    std::string line;
    // Read the file line by line
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;

        // Split the line using the comma ',' delimiter
        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        // Add the row to our 2D data structure
        data.push_back(row);
    }

    file.close();

    //build map
    // Skip row 0 if it contains the CSV header
    for (size_t i = 1; i < data.size(); ++i) {

        const auto& row = data[i];

        const std::string& trigger = row[0];
        int year                = std::stoi(row[1]);
        double mu                = std::stod(row[2]);
        double sigma             = std::stod(row[3]);
        Syst syst  		 = parseSyst(row[4]);
        double htLow             = std::stod(row[5]);
        double htHigh            = std::stod(row[6]);

        EffHTBin bin{
            htLow,
            htHigh,
            mu,
            sigma
        };

        _effparams[trigger][year][syst].push_back(bin);
    }

}


