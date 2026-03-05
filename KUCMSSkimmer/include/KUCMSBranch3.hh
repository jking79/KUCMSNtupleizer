#pragma once
// -*- C++ -*-
//
// Original Author: Jack W King III
// Created: Wed, 27 Jan 2021 19:19:35 GMT
//

#include <string>
#include <vector>

// KUCMS includes
#include "KUCMSItemManager3.hh"

// ROOT includes
#include "TTree.h"

using uInt = unsigned int;

//------------------------------------------------------------------------------
// ScalarBranch / VectorBranch
//------------------------------------------------------------------------------

template <class T>
class ScalarBranch : public Item<T> {
public:
  void attachBranch(TTree* outTree) {
    // Item<T>::iName/iDoc/iValue are expected to be accessible (public or protected)
    outTree->Branch(this->iName.c_str(), &this->iValue)->SetTitle(this->iDoc.c_str());
  }
};

template <class T>
class VectorBranch : public VectorItem<T> {
public:
  void attachBranch(TTree* outTree) {
    // VectorItem<T>::iName/iDoc/iVector are expected to be accessible to derived
    outTree->Branch(this->iName.c_str(), &this->iVector)->SetTitle(this->iDoc.c_str());
  }
};

// NOTE: enum ordering is relied upon (>= UINT means scalar)
enum BType { VVFLOAT, VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };

//------------------------------------------------------------------------------
// Base interface
//------------------------------------------------------------------------------

class KUCMSBranchBase {
public:
  virtual ~KUCMSBranchBase() = default;

  virtual void attachBranch(TTree* outTree) = 0;
  virtual void clear() = 0;

  // Optional no-op overloads (kept for backward compatibility with existing call sites)
  virtual void fill(std::vector<float>) {}
  virtual void fill(std::vector<uInt>) {}
  virtual void fill(uInt) {}
  virtual void fill(int) {}
  virtual void fill(float) {}
  virtual void fill(std::string) {}
  virtual void fill(bool) {}
};

//------------------------------------------------------------------------------
// KUCMSBranch<T>
//------------------------------------------------------------------------------

template <class T>
class KUCMSBranch : public KUCMSBranchBase {
public:
  KUCMSBranch(const std::string& name, BType type, const std::string& doc)
      : vartype(type) {
    if (isScalar()) bscalar.make(name, doc);
    else           bvector.make(name, doc);
  }

  void attachBranch(TTree* outTree) override {
    if (isScalar()) bscalar.attachBranch(outTree);
    else           bvector.attachBranch(outTree);
  }

  void clear() override {
    if (isScalar()) bscalar.clear();
    else           bvector.clear();
  }

  // type-safe fill for this T
  void fill(T val) { // (kept by-value to match your original signature)
    if (isScalar()) bscalar.fill(val);
    else           bvector.fill(val);
  }

private:
  bool isScalar() const { return vartype >= UINT; }

  BType vartype;
  ScalarBranch<T> bscalar;
  VectorBranch<T> bvector;
};
