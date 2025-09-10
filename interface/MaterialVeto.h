#pragma once

#include <memory>
#include <string>
#include <iostream>
#include "TH2D.h"
#include "TFile.h"

class MaterialVeto {
public:
    MaterialVeto() = default;

    explicit MaterialVeto(const std::string& path) {
        Load(path);
    }

    // Load a veto map from file
    bool Load(const std::string& path) {
        TFile* file = TFile::Open(path.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open ROOT file: " << path << std::endl;
            return false;
        }

        TH2D* material_mask_raw = nullptr;
        file->GetObject("material_mask", material_mask_raw);
        if (!material_mask_raw) {
            std::cerr << "Error: Could not find histogram 'material_mask' in file!" << std::endl;
            file->Close();
            delete file;
            return false;
        }

        material_mask_raw->SetDirectory(0); // detach from file
        vetoMap_.reset(material_mask_raw);

        file->Close();
        delete file;

        return true;
    }

    // Check if a point passes the veto
    bool PassVeto(double x, double y) const {
        if (!vetoMap_) {
            std::cerr << "Error: Material mask is not loaded!" << std::endl;
            return false;
        }

        int bin_x = vetoMap_->GetXaxis()->FindBin(x);
        int bin_y = vetoMap_->GetYaxis()->FindBin(y);

        if (bin_x < 1 || bin_x > vetoMap_->GetNbinsX() ||
            bin_y < 1 || bin_y > vetoMap_->GetNbinsY()) {
            return true; // out of range = pass
        }

        double bin_content = vetoMap_->GetBinContent(bin_x, bin_y);
        return bin_content == 0;
    }

    // Access the underlying histogram if needed
    TH2D* GetMap() const { return vetoMap_.get(); }

private:
    std::unique_ptr<TH2D> vetoMap_{nullptr};
};
