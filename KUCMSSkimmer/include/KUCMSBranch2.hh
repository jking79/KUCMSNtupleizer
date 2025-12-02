// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>

//KUCMS includes
#include "KUCMSItemManager.hh"

//Root includes
#include "TTree.h"

#ifndef KUCMSBranchHeader
#define KUCMSBranchHeader

typedef unsigned int uInt;

template < class T >
class ScalarBranch : public Item<T> {

    public:
    
    void attachBranch( TTree* fOutTree ){ 
           fOutTree->Branch( Item<T>::iName.c_str(), &this->Item<T>::iValue )->SetTitle( Item<T>::iDoc.c_str() ); };

};//<<>> class Branch : Item<T> 


template < class T >
class VectorBranch : public VectorItem<T> {

    public:

    void attachBranch( TTree* fOutTree ){
           fOutTree->Branch( VectorItem<T>::iName.c_str(), &this->VectorItem<T>::iVector )->SetTitle( VectorItem<T>::iDoc.c_str() ); };

};//<<>> class Branch : Item<T> 

//enum BType{ VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };
enum BType{ VVFLOAT, VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };

class KUCMSBranchBase {

    public:

    KUCMSBranchBase(){};
    virtual ~KUCMSBranchBase(){};

    virtual void attachBranch( TTree* fOutTree ){};
    virtual void clear(){};

    virtual void fill( std::vector<float> val ){};
    virtual void fill( std::vector<uInt> val ){};
    virtual void fill( uInt val ){};
    virtual void fill( int val ){};
    virtual void fill( float val ){};
    virtual void fill( std::string val ){};
    virtual void fill( bool val ){};

};//<<>> class KUCMSBranchBase

template <class T>
class KUCMSBranch : public KUCMSBranchBase {

    public:

    KUCMSBranch( std::string name, BType type, std::string doc  );
    ~KUCMSBranch(){};

    void attachBranch( TTree* fOutTree );
    void clear();
    void fill( T val );
    
    private:

    BType vartype;
    ScalarBranch<T> bscalar;
    VectorBranch<T> bvector;
 
};//class Branch 

////////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
KUCMSBranch<T>::KUCMSBranch( std::string name, BType type, std::string doc  ){

    vartype = type;
    ( vartype >= UINT ) ? bscalar.make( name, doc ) : bvector.make( name, doc );

}//<<>>KUCMSBranch<T>::KUCMSBranch( std::string name, BType type, std::string doc  )

template <class T>
void KUCMSBranch<T>::attachBranch( TTree* fOutTree ){ 

    ( vartype >= UINT ) ? bscalar.attachBranch( fOutTree ) : bvector.attachBranch( fOutTree );

}//<<>>void KUCMSBranch<T>::attachBranch()

template <class T>
void KUCMSBranch<T>::clear(){

    ( vartype >= UINT ) ? bscalar.clear() : bvector.clear();

}//<<>>void KUCMSBranch<T>::clear()

template <class T>
void KUCMSBranch<T>::fill( T val ){

    ( vartype >= UINT ) ? bscalar.fill( val ) : bvector.fill( val );

}//<<>>void KUCMSBranch<T>::fill( T val )

////////////////////////////////////////////////////////////////////////////

#endif
