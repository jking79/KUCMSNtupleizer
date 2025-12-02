// -*- C++ -*-
//
// Original Author:  Jack W King III
// Created: Wed, 27 Jan 2021 19:19:35 GMT
//
// Header-only version: fully inline, safe for templates
//

#ifndef KUCMSBRANCHMANAGER2_HH
#define KUCMSBRANCHMANAGER2_HH

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <iostream>

// KUCMS includes
#include "KUCMSItemManager.hh"
#include "KUCMSBranch2.hh"

// ROOT includes
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////////////////////

class KUCMSBranchManager {

public:

    KUCMSBranchManager() {}
    
    // destructor must be inline in header for header-only
    ~KUCMSBranchManager() {
        for (auto &b : theBranches) {
            delete b.second;
        }
    }

    // ------------ Branch creation ------------

    inline void makeBranch(std::string key, std::string name, BType type, std::string doc = "") {
        
        switch (type) {

            case VVFLOAT: theBranches[key] = new KUCMSBranch<std::vector<float>>(name, type, doc); break;
            case VVUINT:  theBranches[key] = new KUCMSBranch<std::vector<uInt>>(name, type, doc); break;

            case VUINT:   theBranches[key] = new KUCMSBranch<uInt>(name, type, doc); break;
            case VINT:    theBranches[key] = new KUCMSBranch<int>(name, type, doc); break;
            case VFLOAT:  theBranches[key] = new KUCMSBranch<float>(name, type, doc); break;
            case VSTR:    theBranches[key] = new KUCMSBranch<std::string>(name, type, doc); break;
            case VBOOL:   theBranches[key] = new KUCMSBranch<bool>(name, type, doc); break;

            case UINT:    theBranches[key] = new KUCMSBranch<uInt>(name, type, doc); break;
            case INT:     theBranches[key] = new KUCMSBranch<int>(name, type, doc); break;
            case FLOAT:   theBranches[key] = new KUCMSBranch<float>(name, type, doc); break;
            case STR:     theBranches[key] = new KUCMSBranch<std::string>(name, type, doc); break;
            case BOOL:    theBranches[key] = new KUCMSBranch<bool>(name, type, doc); break;

            default:
                std::cout << " -- KUCMSBranch " << name
                          << " Error : BranchType error in makeBranch!!!! "
                          << std::endl;
        }
    }

    inline void makeBranch(std::string name, BType type, std::string doc = "") {
        makeBranch(name, name, type, doc);
    }

    // ------------ ROOT tree interface ------------

    inline void attachBranches(TTree* fOutTree) {
        for (auto &b : theBranches) {
            b.second->attachBranch(fOutTree);
        }
    }

    inline void clearBranches() {
        for (auto &b : theBranches) {
            b.second->clear();
        }
    }

    // ------------ Fillers ------------

    inline void fillBranch(std::string key, std::vector<float>   v) { if (valid(key)) theBranches[key]->fill(v); }
    inline void fillBranch(std::string key, std::vector<uInt>    v) { if (valid(key)) theBranches[key]->fill(v); }
    
    inline void fillBranch(std::string key, uInt        v) { if (valid(key)) theBranches[key]->fill(v); }
    inline void fillBranch(std::string key, int         v) { if (valid(key)) theBranches[key]->fill(v); }
    inline void fillBranch(std::string key, float       v) { if (valid(key)) theBranches[key]->fill(v); }
    inline void fillBranch(std::string key, std::string v) { if (valid(key)) theBranches[key]->fill(v); }
    inline void fillBranch(std::string key, bool        v) { if (valid(key)) theBranches[key]->fill(v); }

private:

    std::map<std::string, KUCMSBranchBase*> theBranches;

    inline bool valid(std::string key) {
        if (theBranches.find(key) == theBranches.end()) {
            std::cout << " -- BM Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!"
                      << std::endl;
            return false;
        }
        return true;
    }

};

#endif

