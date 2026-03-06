// -*- C++ -*-
//
// Original Author:  Jack W King III
// Created: Wed, 27 Jan 2021 19:19:35 GMT
// Header-only version: fully inline
//

#ifndef KUCMSBRANCHMANAGER3_HH
#define KUCMSBRANCHMANAGER3_HH

#include <map>
#include <memory>
#include <string>
#include <vector>

// KUCMS includes
#include "KUCMSBranch3.hh"   // brings in KUCMSBranchBase, KUCMSBranch<T>, BType, uInt

// Forward declare to avoid ROOT headers here (KUCMSBranch2.hh may include TTree.h anyway)
class TTree;

// Optional logging hook (disabled by default for CMSSW cleanliness)
// Define KUCMS_BM_LOG(msg) before including this header if you want messages.
// Example (standalone): #define KUCMS_BM_LOG(x) std::cerr << x << std::endl
#ifndef KUCMS_BM_LOG
#define KUCMS_BM_LOG(x) do {} while (0)
#endif

////////////////////////////////////////////////////////////////////////////////////////////////

class KUCMSBranchManager {
public:
  KUCMSBranchManager() = default;
  ~KUCMSBranchManager() = default;

  KUCMSBranchManager(const KUCMSBranchManager&) = delete;
  KUCMSBranchManager& operator=(const KUCMSBranchManager&) = delete;
  KUCMSBranchManager(KUCMSBranchManager&&) = delete;
  KUCMSBranchManager& operator=(KUCMSBranchManager&&) = delete;

  // ------------ Branch creation ------------

  inline void makeBranch(const std::string& key,
                         const std::string& name,
                         BType type,
                         const std::string& doc = "") {
    // Replace if exists
    theBranches.erase(key);

    switch (type) {
      case VVFLOAT: theBranches[key] = std::make_unique<KUCMSBranch<std::vector<float>>>(name, type, doc); break;
      case VVUINT:  theBranches[key] = std::make_unique<KUCMSBranch<std::vector<uInt>>>(name, type, doc); break;

      case VUINT:   theBranches[key] = std::make_unique<KUCMSBranch<uInt>>(name, type, doc); break;
      case VINT:    theBranches[key] = std::make_unique<KUCMSBranch<int>>(name, type, doc); break;
      case VFLOAT:  theBranches[key] = std::make_unique<KUCMSBranch<float>>(name, type, doc); break;
      case VSTR:    theBranches[key] = std::make_unique<KUCMSBranch<std::string>>(name, type, doc); break;
      case VBOOL:   theBranches[key] = std::make_unique<KUCMSBranch<bool>>(name, type, doc); break;

      case UINT:    theBranches[key] = std::make_unique<KUCMSBranch<uInt>>(name, type, doc); break;
      case INT:     theBranches[key] = std::make_unique<KUCMSBranch<int>>(name, type, doc); break;
      case FLOAT:   theBranches[key] = std::make_unique<KUCMSBranch<float>>(name, type, doc); break;
      case STR:     theBranches[key] = std::make_unique<KUCMSBranch<std::string>>(name, type, doc); break;
      case BOOL:    theBranches[key] = std::make_unique<KUCMSBranch<bool>>(name, type, doc); break;

      default:
        KUCMS_BM_LOG("KUCMSBranchManager::makeBranch error: unknown BType for branch '" << name << "'");
        break;
    }
  }

  inline void makeBranch(const std::string& name, BType type, const std::string& doc = "") {
    makeBranch(name, name, type, doc);
  }

  // ------------ ROOT tree interface ------------

  inline void attachBranches(TTree* outTree) {
    for (auto& kv : theBranches) {
      kv.second->attachBranch(outTree);
    }
  }

  inline void clearBranches() {
    for (auto& kv : theBranches) {
      kv.second->clear();
    }
  }

  // ------------ Fillers (no extra lookups / no accidental insertion) ------------

  inline void fillBranch(const std::string& key, const std::vector<float>& v) {
    if (auto* b = find_(key)) b->fill(v);
  }

  inline void fillBranch(const std::string& key, const std::vector<uInt>& v) {
    if (auto* b = find_(key)) b->fill(v);
  }

  inline void fillBranch(const std::string& key, uInt v) {
    if (auto* b = find_(key)) b->fill(v);
  }

  inline void fillBranch(const std::string& key, int v) {
    if (auto* b = find_(key)) b->fill(v);
  }

  inline void fillBranch(const std::string& key, float v) {
    if (auto* b = find_(key)) b->fill(v);
  }

  inline void fillBranch(const std::string& key, const std::string& v) {
    if (auto* b = find_(key)) b->fill(v);
  }

  inline void fillBranch(const std::string& key, bool v) {
    if (auto* b = find_(key)) b->fill(v);
  }

private:
  using BranchPtr = std::unique_ptr<KUCMSBranchBase>;
  std::map<std::string, BranchPtr> theBranches;

  inline KUCMSBranchBase* find_(const std::string& key) {
    auto it = theBranches.find(key);
    if (it == theBranches.end()) {
      // silent by default (CMSSW-friendly)
      // KUCMS_BM_LOG("KUCMSBranchManager: missing key '" << key << "'");
      return nullptr;
    }
    return it->second.get();
  }
};

#endif
