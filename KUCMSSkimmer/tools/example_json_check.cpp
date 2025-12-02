#include <iostream>
#include <map>
#include <vector>

// Function to check if a given run and lumi are valid
bool isValidLumisection(
    int run,
    int lumi,
    const std::map<int, std::vector<std::map<int, int>>>& lumiData)
{
    // Look for the run
    auto runIt = lumiData.find(run);
    if (runIt == lumiData.end())
        return false;  // Run not found

    // Iterate over all lumi ranges for this run
    for (const auto& range : runIt->second) {
        for (const auto& [start, end] : range) {
            if (lumi >= start && lumi <= end)
                return true;  // Lumi within a valid range
        }
    }

    return false;  // Not within any range
}
