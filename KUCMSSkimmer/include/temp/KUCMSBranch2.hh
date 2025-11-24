// -*- C++ -*-
//
// KUCMSBranch2 – Combined Header-Only Template Implementation
// Original Author: Jack W King III
//

#ifndef KUCMSBranchHeader_hh
#define KUCMSBranchHeader_hh

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>

// KUCMS includes
#include "KUCMSItemManager.hh"

// Root includes
#include "TTree.h"

typedef unsigned int uInt;

// ============================================================================
// ScalarBranch<T>
// ============================================================================

template <class T>
class ScalarBranch : public Item<T> {

public:

    ScalarBranch() : Item<T>("") {}
    ScalarBranch(std::string name, std::string doc = "") : Item<T>(name, doc) {}

    // Create the TTree branch for scalar data
    void attachBranch(TTree* fOutTree) {
        fOutTree->Branch(Item<T>::iName.c_str(), &Item<T>::iValue)
                 ->SetTitle(Item<T>::iDoc.c_str());
    }
};

// ============================================================================
// VectorBranch<T>
// ============================================================================

template <class T>
class VectorBranch : public VectorItem<T> {

public:

    VectorBranch() : VectorItem<T>("") {}
    VectorBranch(std::string name, std::string doc = "") : VectorItem<T>(name, doc) {}

    // Create the TTree branch for vector data
    void attachBranch(TTree* fOutTree) {
        fOutTree->Branch(VectorItem<T>::iName.c_str(), &VectorItem<T>::iVector)
                 ->SetTitle(VectorItem<T>::iDoc.c_str());
    }
};

// ============================================================================
// BType enum
// ============================================================================

//enum BType{ VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };
enum BType{ VVFLOAT, VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };

// ============================================================================
// Base class for branch wrappers
// ============================================================================

class KUCMSBranchBase {

public:

    KUCMSBranchBase() {}
    virtual ~KUCMSBranchBase() {}

    virtual void attachBranch(TTree* fOutTree) {}
    virtual void clear() {}

    virtual void fill(std::vector<float> val) {}
    virtual void fill(std::vector<uInt> val) {}
    virtual void fill(uInt val) {}
    virtual void fill(int val) {}
    virtual void fill(float val) {}
    virtual void fill(std::string val) {}
    virtual void fill(bool val) {}
};


// ============================================================================
// KUCMSBranch<T>
// One branch instance (scalar or vector) depending on BType
// ============================================================================

template <class T>
class KUCMSBranch : public KUCMSBranchBase {

public:

    KUCMSBranch(std::string name, BType type, std::string doc);
    ~KUCMSBranch() {}

    void attachBranch(TTree* fOutTree);
    void clear();
    void fill(T val);

private:

    BType vartype;
    ScalarBranch<T> bscalar;
    VectorBranch<T> bvector;
};


// ============================================================================
// TEMPLATE IMPLEMENTATIONS
// ============================================================================

template <class T>
KUCMSBranch<T>::KUCMSBranch(std::string name, BType type, std::string doc)
{
    vartype = type;

    if (vartype >= UINT)
        bscalar = ScalarBranch<T>(name, doc);
    else
        bvector = VectorBranch<T>(name, doc);
}


template <class T>
void KUCMSBranch<T>::attachBranch(TTree* fOutTree)
{
    (vartype >= UINT) ?
        bscalar.attachBranch(fOutTree) :
        bvector.attachBranch(fOutTree);
}


template <class T>
void KUCMSBranch<T>::clear()
{
    (vartype >= UINT) ?
        bscalar.clear() :
        bvector.clear();
}


template <class T>
void KUCMSBranch<T>::fill(T val)
{
    (vartype >= UINT) ?
        bscalar.fill(val) :
        bvector.fill(val);
}

#endif // KUCMSBranchHeader_hh

