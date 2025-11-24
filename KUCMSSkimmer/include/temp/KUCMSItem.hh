// -*- C++ -*-
//
// Original Author:  Jack W King III
// Created: Wed, 27 Jan 2021 19:19:35 GMT
//

#ifndef KUCMSItemHeader_hh
#define KUCMSItemHeader_hh

// Standard includes
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

// ROOT
#include "TTree.h"

// --------------------------------------------------------------------------------------------
// Base class for all items
// --------------------------------------------------------------------------------------------

class ItemBase {
public:
    ItemBase() {}
    virtual ~ItemBase() {}

    virtual void clear() = 0;
    virtual void attachBranch(TTree* fOutTree) = 0;
};


// --------------------------------------------------------------------------------------------
// Template class Item<T>
// Holds a single value of type T
// --------------------------------------------------------------------------------------------

template <class T>
class Item : public ItemBase {
public:
    Item(std::string name) : theName(name), theVal(T()) {}
    virtual ~Item() {}

    inline void fill(T val) { theVal = val; }
    inline T value() const { return theVal; }

    void attachBranch(TTree* fOutTree) override {
        fOutTree->Branch(theName.c_str(), &theVal);
    }

    void clear() override { theVal = T(); }

private:
    std::string theName;
    T theVal;
};


// --------------------------------------------------------------------------------------------
// Template class VectorItem<T>
// Holds a vector<T>
// --------------------------------------------------------------------------------------------

template <class T>
class VectorItem : public ItemBase {
public:
    VectorItem(std::string name) : theName(name), theVals() {}
    virtual ~VectorItem() {}

    inline void fill(T v) { theVals.push_back(v); }
    inline std::vector<T>& value() { return theVals; }

    void attachBranch(TTree* fOutTree) override {
        fOutTree->Branch(theName.c_str(), &theVals);
    }

    void clear() override { theVals.clear(); }

private:
    std::string theName;
    std::vector<T> theVals;
};


// --------------------------------------------------------------------------------------------
// Explicit template specializations (if any were in your .cc)
// None were explicitly declared, but your code uses:
//   Item<bool>, Item<int>, Item<float>, Item<uInt>, Item<std::string>
//   VectorItem<float>, VectorItem<uInt>, etc.
// These are automatically handled because all implementations are here.
// --------------------------------------------------------------------------------------------


#endif // KUCMSItemHeader_hh

