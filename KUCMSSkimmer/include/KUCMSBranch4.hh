#pragma once
// -*- C++ -*-
// KUCMSBranch2.hh (header-only; requires ROOT)

#include <cstdint>
#include <string>
#include <vector>

#include "KUCMSItemManager.hh"
#include "TTree.h"

namespace kucms {

template <class T>
class ScalarBranch : public Item<T> {
public:
  void attachBranch(TTree* outTree) {
    outTree->Branch(this->iName.c_str(), &this->iValue)->SetTitle(this->iDoc.c_str());
  }
};

template <class T>
class VectorBranch : public VectorItem<T> {
public:
  void attachBranch(TTree* outTree) {
    outTree->Branch(this->iName.c_str(), &this->iVector)->SetTitle(this->iDoc.c_str());
  }
};

// More explicit than relying on enum ordering
enum class BType {
  // scalars
  BOOL,
  I32,
  U32,
  I64,
  U64,
  F32,
  F64,
  STR,

  // vectors
  V_BOOL8,   // std::vector<uint8_t> recommended for flags
  V_I32,
  V_U32,
  V_I64,
  V_U64,
  V_F32,
  V_F64,
  V_STR,

  // nested vectors
  VV_I32,
  VV_U32,
  VV_F32,
  VV_F64
};

class KUCMSBranchBase {
public:
  virtual ~KUCMSBranchBase() = default;
  virtual void attachBranch(TTree* outTree) = 0;
  virtual void clear() = 0;

  // Fill overloads (base no-ops, derived typed branch will provide its own fill(T))
  virtual void fill(bool) {}
  virtual void fill(int32_t) {}
  virtual void fill(uint32_t) {}
  virtual void fill(int64_t) {}
  virtual void fill(uint64_t) {}
  virtual void fill(float) {}
  virtual void fill(double) {}
  virtual void fill(std::string) {}

  virtual void fill(std::vector<uint8_t>) {}
  virtual void fill(std::vector<int32_t>) {}
  virtual void fill(std::vector<uint32_t>) {}
  virtual void fill(std::vector<int64_t>) {}
  virtual void fill(std::vector<uint64_t>) {}
  virtual void fill(std::vector<float>) {}
  virtual void fill(std::vector<double>) {}
  virtual void fill(std::vector<std::string>) {}

  virtual void fill(std::vector<std::vector<int32_t>>) {}
  virtual void fill(std::vector<std::vector<uint32_t>>) {}
  virtual void fill(std::vector<std::vector<float>>) {}
  virtual void fill(std::vector<std::vector<double>>) {}
};

template <class T>
class KUCMSBranch : public KUCMSBranchBase {
public:
  KUCMSBranch(const std::string& name, const std::string& doc)
      : is_scalar_(true) {
    bscalar_.make(name, doc);
  }

  // vector branch
  template <class U = T, typename std::enable_if_t<std::is_same_v<U, std::vector<typename U::value_type>>, int> = 0>
  KUCMSBranch(const std::string& name, const std::string& doc, int /*tag*/)
      : is_scalar_(false) {
    bvector_.make(name, doc);
  }

  void attachBranch(TTree* outTree) override {
    if (is_scalar_) bscalar_.attachBranch(outTree);
    else bvector_.attachBranch(outTree);
  }

  void clear() override {
    if (is_scalar_) bscalar_.clear();
    else bvector_.clear();
  }

  void fill(T val) { // matches your existing pattern
    if (is_scalar_) bscalar_.fill(val);
    else bvector_.fill(val);
  }

private:
  bool is_scalar_{true};
  ScalarBranch<T> bscalar_{};
  VectorBranch<typename T::value_type> bvector_{}; // only used if T is std::vector<...>
};

// Helper factories for scalar vs vector vs nested-vector
template <typename T>
inline std::unique_ptr<KUCMSBranchBase> make_typed_branch(const std::string& name, const std::string& doc) {
  // scalar
  return std::make_unique<KUCMSBranch<T>>(name, doc);
}

template <typename T>
inline std::unique_ptr<KUCMSBranchBase> make_vector_branch(const std::string& name, const std::string& doc) {
  // vector uses the special ctor tag to select vector storage
  return std::unique_ptr<KUCMSBranchBase>(new KUCMSBranch<std::vector<T>>(name, doc, 0));
}

template <typename T>
inline std::unique_ptr<KUCMSBranchBase> make_vv_branch(const std::string& name, const std::string& doc) {
  return std::unique_ptr<KUCMSBranchBase>(new KUCMSBranch<std::vector<std::vector<T>>>(name, doc, 0));
}

}  // namespace kucms
