#pragma once
// -*- C++ -*-
// KUCMSBranchManager2.hh (header-only)

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "KUCMSBranch2.hh"

namespace kucms {

class KUCMSBranchManager {
public:
  KUCMSBranchManager() = default;
  ~KUCMSBranchManager() = default;

  KUCMSBranchManager(const KUCMSBranchManager&) = delete;
  KUCMSBranchManager& operator=(const KUCMSBranchManager&) = delete;

  // -------- Preferred: type-safe booking --------
  template <typename T>
  void makeBranch(const std::string& key, const std::string& name, const std::string& doc = "") {
    branches_[key] = make_typed_branch<T>(name, doc);
  }

  template <typename T>
  void makeBranch(const std::string& name, const std::string& doc = "") {
    makeBranch<T>(name, name, doc);
  }

  // -------- Convenience: runtime booking via BType --------
  void makeBranch(const std::string& key, const std::string& name, BType type, const std::string& doc = "") {
    using BT = BType;
    switch (type) {
      case BT::BOOL:    branches_[key] = make_typed_branch<bool>(name, doc); break;
      case BT::I32:     branches_[key] = make_typed_branch<int32_t>(name, doc); break;
      case BT::U32:     branches_[key] = make_typed_branch<uint32_t>(name, doc); break;
      case BT::I64:     branches_[key] = make_typed_branch<int64_t>(name, doc); break;
      case BT::U64:     branches_[key] = make_typed_branch<uint64_t>(name, doc); break;
      case BT::F32:     branches_[key] = make_typed_branch<float>(name, doc); break;
      case BT::F64:     branches_[key] = make_typed_branch<double>(name, doc); break;
      case BT::STR:     branches_[key] = make_typed_branch<std::string>(name, doc); break;

      case BT::V_BOOL8: branches_[key] = make_vector_branch<uint8_t>(name, doc); break;
      case BT::V_I32:   branches_[key] = make_vector_branch<int32_t>(name, doc); break;
      case BT::V_U32:   branches_[key] = make_vector_branch<uint32_t>(name, doc); break;
      case BT::V_I64:   branches_[key] = make_vector_branch<int64_t>(name, doc); break;
      case BT::V_U64:   branches_[key] = make_vector_branch<uint64_t>(name, doc); break;
      case BT::V_F32:   branches_[key] = make_vector_branch<float>(name, doc); break;
      case BT::V_F64:   branches_[key] = make_vector_branch<double>(name, doc); break;
      case BT::V_STR:   branches_[key] = make_vector_branch<std::string>(name, doc); break;

      case BT::VV_I32:  branches_[key] = make_vv_branch<int32_t>(name, doc); break;
      case BT::VV_U32:  branches_[key] = make_vv_branch<uint32_t>(name, doc); break;
      case BT::VV_F32:  branches_[key] = make_vv_branch<float>(name, doc); break;
      case BT::VV_F64:  branches_[key] = make_vv_branch<double>(name, doc); break;
    }
  }

  void makeBranch(const std::string& name, BType type, const std::string& doc = "") {
    makeBranch(name, name, type, doc);
  }

  // -------- Tree ops --------
  void attachBranches(TTree* outTree) {
    for (auto& kv : branches_) kv.second->attachBranch(outTree);
  }

  void clearBranches() {
    for (auto& kv : branches_) kv.second->clear();
  }

  // -------- Fillers (typed overloads) --------
  void fillBranch(const std::string& key, bool v)                { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, int32_t v)             { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, uint32_t v)            { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, int64_t v)             { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, uint64_t v)            { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, float v)               { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, double v)              { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::string& v)  { if (auto* b = find_(key)) b->fill(v); }

  void fillBranch(const std::string& key, const std::vector<uint8_t>& v)        { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<int32_t>& v)        { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<uint32_t>& v)       { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<int64_t>& v)        { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<uint64_t>& v)       { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<float>& v)          { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<double>& v)         { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<std::string>& v)    { if (auto* b = find_(key)) b->fill(v); }

  void fillBranch(const std::string& key, const std::vector<std::vector<int32_t>>& v)  { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<std::vector<uint32_t>>& v) { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<std::vector<float>>& v)    { if (auto* b = find_(key)) b->fill(v); }
  void fillBranch(const std::string& key, const std::vector<std::vector<double>>& v)   { if (auto* b = find_(key)) b->fill(v); }

private:
  KUCMSBranchBase* find_(const std::string& key) {
    auto it = branches_.find(key);
    return (it == branches_.end()) ? nullptr : it->second.get();
  }

  std::map<std::string, std::unique_ptr<KUCMSBranchBase>> branches_;
};

}  // namespace kucms
