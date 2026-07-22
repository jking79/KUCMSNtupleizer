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
#include <set>
#include <string>
#include <cmath>
#include <iostream>
#include <fnmatch.h>

// KUCMS includes
#include "KUCMSItemManager.hh"
#include "KUCMSBranch2.hh"

// ROOT includes
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////////////////////

class KUCMSBranchManager {

public:

  KUCMSBranchManager(){}
    
  // destructor must be inline in header for header-only
  ~KUCMSBranchManager(){ for (auto &b : theBranches){ delete b.second; b.second = nullptr; } theBranches.clear(); }

  KUCMSBranchManager(const KUCMSBranchManager&) = delete;
  KUCMSBranchManager& operator=(const KUCMSBranchManager&) = delete;
  KUCMSBranchManager(KUCMSBranchManager&&) = delete;
  KUCMSBranchManager& operator=(KUCMSBranchManager&&) = delete;

  // ------------ Branch creation ------------

  inline void setBranchMask( const std::vector<std::string>& patterns ){
    branchMaskPatterns = patterns;
    branchMaskActive = !branchMaskPatterns.empty();
    maskedBranches.clear();
    maskedOutputBranches.clear();
    bookedBranches.clear();
    matchedPatterns.clear();
  }

  inline void clearBranchMask(){
    branchMaskPatterns.clear();
    maskedBranches.clear();
    maskedOutputBranches.clear();
    bookedBranches.clear();
    matchedPatterns.clear();
    branchMaskActive = false;
  }

  inline bool hasBranchMask() const { return branchMaskActive; }

  inline void makeBranch( std::string key, std::string name, BType type, std::string doc = "" ){

    auto it = theBranches.find( key );
    if( it != theBranches.end() ){ delete it->second; it->second = nullptr; theBranches.erase(it); }
    if( branchMaskActive && !branchEnabled( name ) && !branchEnabled( key ) ){
      maskedBranches.insert( key );
      maskedBranches.insert( name );
      maskedOutputBranches.insert( name );
      return;
    }
    bookedBranches.insert( name );

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
      std::cout << " -- KUCMSBranch " << name << " Error : BranchType error in makeBranch!!!! " << std::endl;
    }//<<>>switch (type) 
  }//<<>>inline void makeBranch(std::string key, std::string name, BType type, std::string doc = "")

  inline void makeBranch(std::string name, BType type, std::string doc = "") { makeBranch(name, name, type, doc); }

  // ------------ ROOT tree interface ------------

  inline void attachBranches(TTree* fOutTree) { for (auto &b : theBranches) { b.second->attachBranch(fOutTree); } }
  inline void clearBranches() { for (auto &b : theBranches) { b.second->clear(); } }

  inline void printBranchMaskSummary( const std::string& label ) const {
    if( !branchMaskActive ) return;
    std::cout << " -- Branch mask summary [" << label << "] : booked "
	      << bookedBranches.size() << " masked " << maskedOutputBranches.size() << std::endl;
  }

  inline void collectMatchedPatterns( std::set<std::string>& matched ) const {
    if( !branchMaskActive ) return;
    for( const auto& pattern : matchedPatterns ) matched.insert(pattern);
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
  std::vector<std::string> branchMaskPatterns;
  std::set<std::string> maskedBranches;
  std::set<std::string> maskedOutputBranches;
  std::set<std::string> bookedBranches;
  mutable std::set<std::string> matchedPatterns;
  bool branchMaskActive = false;

  inline bool matchPattern( const std::string& text, const std::string& pattern ) const {
    return fnmatch(pattern.c_str(), text.c_str(), 0) == 0;
  }

  inline bool branchEnabled( const std::string& name ) const {
    if( !branchMaskActive ) return true;
    for( const auto& pattern : branchMaskPatterns ){
      if( matchPattern( name, pattern ) ){
	matchedPatterns.insert(pattern);
	return true;
      }
    }
    return false;
  }

  inline bool valid( std::string key ){
    if( theBranches.find(key) == theBranches.end() ){
      if( maskedBranches.find(key) != maskedBranches.end() ) return false;
      std::cout << " -- BM Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      return false;
    }//<<>>if (theBranches.find(key) == theBranches.end())
    return true;
  }//<<>>inline bool valid(std::string key)

};

#endif
