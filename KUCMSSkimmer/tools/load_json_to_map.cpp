#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <nlohmann_json.hpp>

using json = nlohmann::json;

int main() {
    // Define your target data structure
    std::map<int, std::vector<std::map<int, int>>> lumiData;

    // Open JSON file
    std::ifstream file("Cert_Example.json");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open JSON file.\n";
        return 1;
    }

    // Parse JSON
    json j;
    file >> j;

    // Convert JSON structure into C++ map
    for (auto& [runStr, lumiRanges] : j.items()) {
        int run = std::stoi(runStr);
        std::vector<std::map<int, int>> rangeList;

        for (auto& pair : lumiRanges) {
            if (pair.is_array() && pair.size() == 2) {
                std::map<int, int> range;
                range[pair[0].get<int>()] = pair[1].get<int>();
                rangeList.push_back(range);
            }
        }

        lumiData[run] = rangeList;
    }

    // ✅ Test output
    for (const auto& [run, ranges] : lumiData) {
        std::cout << "Run " << run << ":\n";
        for (const auto& range : ranges) {
            for (const auto& [start, end] : range) {
                std::cout << "  Lumi " << start << " → " << end << "\n";
            }
        }
    }

    return 0;
}

