#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main() {
    // -------------------------------------------------------------------------
    // 1️⃣ Example data structure (you can fill this however you like)
    // -------------------------------------------------------------------------
    std::map<int, std::vector<std::map<int, int>>> lumiData;

    // Example: Run 315257 → two lumi ranges {1→45, 50→90}
    lumiData[315257] = { {{1, 45}}, {{50, 90}} };
    // Example: Run 315258 → one lumi range {1→120}
    lumiData[315258] = { {{1, 120}} };

    // -------------------------------------------------------------------------
    // 2️⃣ Convert to JSON
    // -------------------------------------------------------------------------
    json j;

    for (const auto& [run, ranges] : lumiData) {
        json runRanges = json::array();

        for (const auto& range : ranges) {
            for (const auto& [start, end] : range) {
                runRanges.push_back({start, end});  // store as [start, end]
            }
        }

        j[std::to_string(run)] = runRanges;
    }

    // -------------------------------------------------------------------------
    // 3️⃣ Write JSON to file (pretty printed)
    // -------------------------------------------------------------------------
    std::ofstream out("output_lumi.json");
    if (!out.is_open()) {
        std::cerr << "❌ Error: Could not open output_lumi.json for writing\n";
        return 1;
    }

    out << std::setw(2) << j << std::endl;
    out.close();

    std::cout << "✅ JSON file written: output_lumi.json\n";

    // -------------------------------------------------------------------------
    // 4️⃣ (Optional) Print to console
    // -------------------------------------------------------------------------
    std::cout << j.dump(2) << "\n";

    return 0;
}

