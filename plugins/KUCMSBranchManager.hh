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

//Root includes
#include "TTree.h"

#ifndef KUCMSBranchHeader
#define KUCMSBranchHeader

typedef unsigned int uInt;

    // BranchType List
    // --------------------------
    // vector<vector<uInt>  VVUINT
    // vector<uInt>         VUINT 
    // vector<int>          VINT
    // vector<float>        VFLOAT
    // vector<string>       VSTR
    // vector<bool>         VBOOL
    // uInt                 UINT
    // int                  INT
    // float                FLOAT
    // string               STR
    // bool                 BOOL

enum BType{ VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };


class KUCMSBranch {

	public:

	KUCMSBranch();
	KUCMSBranch( BType type, std::string name, std::string doc );
	//~KUCMSBranch();
	
	void initBranch( TTree* fOutTree );
	void clearBranch();

	void fillBranch( std::vector<uInt> val );
    void fillBranch( uInt val );
    void fillBranch( int val );
    void fillBranch( float val );
    void fillBranch( std::string val );
    void fillBranch( bool val );

	void getValue( std::vector<uInt>& val, uInt index );
    void getValue( uInt& val, uInt index );
    void getValue( int& val, uInt index );
    void getValue( float& val, uInt index );
    void getValue( std::string& val, uInt index );
    void getValue( bool& val, uInt index );

    std::vector<uInt> getVUINTValue( uInt index );
    uInt getUINTValue( uInt index );
    int getINTValue( uInt index );
    float getFLOATValue( uInt index );
    std::string getSTRValue( uInt index );
    bool getBOOLValue( uInt index );

	float getMaxValue();
	uInt getLeadIndex();
	uInt getSubLeadIndex();

    float max( std::vector<float> x );
    uInt leadIdx( std::vector<float> x );
    uInt subldIdx( std::vector<float> x, uInt ldx );

	private:

	BType BranchType;
	std::string BranchName;
	std::string BranchDoc;
	std::vector<std::vector<uInt>> VVUINTBranch;
	std::vector<uInt> VUINTBranch;
	std::vector<int> VINTBranch;
    std::vector<float> VFLOATBranch;
    std::vector<std::string> VSTRBranch;
    std::vector<bool> VBOOLBranch;
	uInt UINTBranch;
	int INTBranch;
	float FLOATBranch;
	std::string STRBranch;
	bool BOOLBranch;

};//<<>>class KUCMSBranch

KUCMSBranch::KUCMSBranch(){};

KUCMSBranch::KUCMSBranch( BType type, std::string name, std::string doc = "" ):

	BranchType(type),
	BranchName(name),
	BranchDoc(doc)

{}//<<>>KUCMSBranch::KUCMSBranch( BType type, std::string name, std::string doc )


void KUCMSBranch::initBranch( TTree* fOutTree ){

	switch( BranchType ){

		case VVUINT	: fOutTree->Branch( BranchName.c_str(), &this->VVUINTBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VUINT	: fOutTree->Branch( BranchName.c_str(), &this->VUINTBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VINT	: fOutTree->Branch( BranchName.c_str(), &this->VINTBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VFLOAT	: fOutTree->Branch( BranchName.c_str(), &this->VFLOATBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VSTR	: fOutTree->Branch( BranchName.c_str(), &this->VSTRBranch )->SetTitle( BranchDoc.c_str() ); break;
        case VBOOL  : fOutTree->Branch( BranchName.c_str(), &this->VBOOLBranch )->SetTitle( BranchDoc.c_str() ); break;
		case UINT	: fOutTree->Branch( BranchName.c_str(), &this->UINTBranch )->SetTitle( BranchDoc.c_str() ); break;
		case INT	: fOutTree->Branch( BranchName.c_str(), &this->INTBranch )->SetTitle( BranchDoc.c_str() ); break;
		case FLOAT	: fOutTree->Branch( BranchName.c_str(), &this->FLOATBranch )->SetTitle( BranchDoc.c_str() ); break;
		case STR	: fOutTree->Branch( BranchName.c_str(), &this->STRBranch )->SetTitle( BranchDoc.c_str() ); break;
        case BOOL	: fOutTree->Branch( BranchName.c_str(), &this->BOOLBranch )->SetTitle( BranchDoc.c_str() ); break;
		default	: std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in setBranch!!!! " << std::endl;

	}//<<>>switch( BranchType )

}//<<>>void setBranch( TTree* fOutTree )

void KUCMSBranch::clearBranch(){

    switch( BranchType ){

        case VVUINT	: VVUINTBranch.clear(); break;
        case VUINT	: VUINTBranch.clear(); break;
        case VINT	: VINTBranch.clear(); break;
        case VFLOAT	: VFLOATBranch.clear(); break;
        case VSTR	: VSTRBranch.clear(); break;
        case VBOOL  : VBOOLBranch.clear(); break;
        case UINT   : UINTBranch = std::numeric_limits<unsigned int>::max(); break;
        case INT    : INTBranch = std::numeric_limits<int>::max(); break;
        case FLOAT  : FLOATBranch = std::numeric_limits<float>::max(); break;
        case STR    : STRBranch = ""; break;
        case BOOL   : BOOLBranch = false; break;
		default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in clearBranch!!!! " << std::endl;    

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::clearBranch()

void KUCMSBranch::fillBranch( std::vector<uInt> val ){

    switch( BranchType ){

        case VVUINT  : VVUINTBranch.push_back(val); break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( std::vector<uInt> val )

void KUCMSBranch::fillBranch( uInt val ){

    switch( BranchType ){

        case VUINT  : VUINTBranch.push_back(val); break;
        case UINT   : UINTBranch = val; break;
   		default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;
	 
    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( uInt val )

void KUCMSBranch::fillBranch( int val ){

    switch( BranchType ){

        case VINT	: VINTBranch.push_back(val); break;
        case INT    : INTBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( int val )

void KUCMSBranch::fillBranch( float val ){

    switch( BranchType ){

        case VFLOAT : VFLOATBranch.push_back(val); break;
        case FLOAT  : FLOATBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( float val )

void KUCMSBranch::fillBranch( std::string val ){

    switch( BranchType ){

        case VSTR   : VSTRBranch.push_back(val); break;
        case STR    : STRBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( std::string val )

void KUCMSBranch::fillBranch( bool val ){

    switch( BranchType ){

        case VBOOL  : VBOOLBranch.push_back(val); break;
        case BOOL   : BOOLBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( bool val )

void KUCMSBranch::getValue( std::vector<uInt>& val, uInt index = 0 ){

    switch( BranchType ){

        case VVUINT  : val = VVUINTBranch[index]; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )	

}//<<>>void KUCMSBranch::getBranch( std::vector<uInt>& val, uInt index = 0 )

void KUCMSBranch::getValue( uInt& val, uInt index = 0 ){

    switch( BranchType ){

        case VUINT   : val = VUINTBranch[index]; break;
        case UINT    : val = UINTBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( uInt& val, uInt index = 0 )

void KUCMSBranch::getValue( int& val, uInt index = 0 ){

    switch( BranchType ){

        case VINT   : val = VINTBranch[index]; break;
        case INT    : val = INTBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( int& val, uInt index = 0 )

void KUCMSBranch::getValue( float& val, uInt index = 0 ){

    switch( BranchType ){

        case VFLOAT	: val = VFLOATBranch[index]; break;
        case FLOAT	: val = FLOATBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( float& val, uInt index = 0 )

void KUCMSBranch::getValue( std::string& val, uInt index = 0 ){

    switch( BranchType ){

        case VSTR	: val = VSTRBranch[index]; break;
        case STR    : val = STRBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( std::string& val, uInt index = 0 )

void KUCMSBranch::getValue( bool& val, uInt index = 0 ){

    switch( BranchType ){

        case VBOOL   : val = VBOOLBranch[index]; break;
        case BOOL   : val = BOOLBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( bool& val, uInt index = 0 )

std::vector<uInt> KUCMSBranch::getVUINTValue( uInt index = 0 ){ std::vector<uInt> val; getValue( val, index ); return val; }
uInt KUCMSBranch::getUINTValue( uInt index = 0 ){ uInt val; getValue( val, index ); return val; }
int KUCMSBranch::getINTValue( uInt index = 0 ){ int val; getValue( val, index ); return val; }
float KUCMSBranch::getFLOATValue( uInt index = 0 ){ float val; getValue( val, index ); return val; }
std::string KUCMSBranch::getSTRValue( uInt index = 0 ){ std::string val; getValue( val, index ); return val; }
bool KUCMSBranch::getBOOLValue( uInt index = 0 ){ bool val; getValue( val, index ); return val; }

float KUCMSBranch::max( std::vector<float> x ){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

uInt KUCMSBranch::leadIdx( std::vector<float> x ){

	float m(x[0]); uInt idx(0), it(0); 
	if( x.size() == 0 ) return -9; 
	if( x.size() == 1 ) return 0; 
	for(auto ix : x ){ if( ix > m ){ m = ix; idx = it; } it++; } 
	return idx;

}//<<>>uInt KUCMSBranch::leadIdx( std::vector<float> x )

uInt KUCMSBranch::subldIdx( std::vector<float> x, uInt ldx ){

	float m(x[0]); uInt idx(0), it(0); 
	if( x.size() <= 1 ) return -9; 
	if( ldx == 0 ){ m = x[1]; idx = 1;} 
	for(auto ix : x ){ if( ix > m && ix < x[ldx] ){ m = ix; idx = it; } it++; } 
	return idx;

}//<<>>uInt KUCMSBranch::subldIdx( std::vector<float> x, uInt ldx )


float KUCMSBranch::getMaxValue(){

	switch( BranchType ){

        case VFLOAT     : return max( VFLOATBranch );
        case FLOAT     : return FLOATBranch;
        default : 
			std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in getMaxValue!!!! " << std::endl;
			return -9999.f;
		
    }//<<>>switch( BranchType )

}//<<>>float KUCMSBranch::getMaxValue()

uInt KUCMSBranch::getLeadIndex(){

    switch( BranchType ){

        case VFLOAT     : return leadIdx( VFLOATBranch );
        default : 
            std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in getLeadIndex!!!! " << std::endl;
            return 9999;

    }//<<>>switch( BranchType )

}//<<>>uInt KUCMSBranch::getLeadIndex()

uInt KUCMSBranch::getSubLeadIndex(){

    switch( BranchType ){

        case VFLOAT     : return subldIdx( VFLOATBranch, leadIdx( VFLOATBranch ) );
        default :
            std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in getSubLeadIndex!!!! " << std::endl;
            return 9999;

    }//<<>>switch( BranchType )

}//<<>>uInt KUCMSBranch::getSubLeadIndex()

// -------------------------------------------------------------------------------------------------
// -----------  Manager for KUCMS Branches ---------------------------------------------------------
// -------------------------------------------------------------------------------------------------

class KUCMSBranchManager {

	public:

	void makeBranch( std::string key, std::string name, BType type, std::string doc );
    void makeBranch( std::string name, BType type, std::string doc );
    void attachBranches( TTree* fOutTree );
	void clearBranches();

    void fillBranch( std::string key, std::vector<uInt> val );
	void fillBranch( std::string key, uInt val );
    void fillBranch( std::string key, int val );
    void fillBranch( std::string key, float val );
    void fillBranch( std::string key, std::string val );
    void fillBranch( std::string key, bool val );

	void getBranchValue( std::string key, std::vector<uInt>& val, uInt index );
    void getBranchValue( std::string key, uInt& val, uInt index );
    void getBranchValue( std::string key, int& val, uInt index );
    void getBranchValue( std::string key, float& val, uInt index );
    void getBranchValue( std::string key, std::string& val, uInt index );
    void getBranchValue( std::string key, bool& val, uInt index );

    std::vector<uInt> getVUINTBranchValue( std::string key, uInt index );
    uInt getUINTBranchValue( std::string key, uInt index );
	int getINTBranchValue( std::string key, uInt index );
    float getFLOATBranchValue( std::string key, uInt index );
    std::string getSTRBranchValue( std::string key, uInt index );
    bool getBOOLBranchValue( std::string key, uInt index );

	float getBranchMaxValue( std::string key );
	uInt getBranchLeadIndex( std::string key );
	uInt getBranchSubLeadIndex( std::string key );

	private:

    std::map< std::string, KUCMSBranch > theBranches;
	bool valid( std::string key );

};//<<>>class KUCMSBranchManager 

void KUCMSBranchManager::makeBranch( std::string key, std::string name, BType type, std::string doc = "" ){

	KUCMSBranch newBranch( type, name, doc );
	theBranches[key] = newBranch;
	
}//<<>>void KUCMSBranchManager::makeBranch( std::string key, std::string name, KUCMSBranch::BT type, std::string doc )

void KUCMSBranchManager::makeBranch( std::string name, BType type, std::string doc = "" ){

    KUCMSBranch newBranch( type, name, doc );
    theBranches[name] = newBranch;

}//<<>>void KUCMSBranchManager::makeBranch( std::string name, KUCMSBranch::BT type, std::string doc )

bool KUCMSBranchManager::valid( std::string key ){ if( theBranches.find(key) == theBranches.end() ){ 
		std::cout << " -- Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; return false; } else return true;}

void KUCMSBranchManager::clearBranches(){ for( auto & branch : theBranches ){ (branch.second).clearBranch();}}
void KUCMSBranchManager::attachBranches( TTree* fOutTree ){ for( auto & branch : theBranches ){ branch.second.initBranch( fOutTree );}}

void KUCMSBranchManager::fillBranch( std::string key, std::vector<uInt> val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, uInt val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, int val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, float val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, std::string val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, bool val ){ if(valid(key)) theBranches[key].fillBranch( val );}

void KUCMSBranchManager::getBranchValue( std::string key, std::vector<uInt>& val, uInt index = 0 ){
							if(valid(key)) theBranches[key].getValue(val,index); else { std::vector<uInt> c{0}; val = c;}}
void KUCMSBranchManager::getBranchValue( std::string key, uInt& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); else val = 0;}
void KUCMSBranchManager::getBranchValue( std::string key, int& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); else val = 0;}
void KUCMSBranchManager::getBranchValue( std::string key, float& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); else val = 0;}
void KUCMSBranchManager::getBranchValue( std::string key, std::string& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); }
void KUCMSBranchManager::getBranchValue( std::string key, bool& val, uInt index = 0 ){
                            if(valid(key)) theBranches[key].getValue( val, index ); }
std::vector<uInt> KUCMSBranchManager::getVUINTBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getVUINTValue(index); else { std::vector<uInt> c{0}; return c; }}
uInt KUCMSBranchManager::getUINTBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getUINTValue(index); else return 0;}
int KUCMSBranchManager::getINTBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getINTValue(index); else return 0;}
float KUCMSBranchManager::getFLOATBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getFLOATValue(index); else return 0;}
std::string KUCMSBranchManager::getSTRBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getSTRValue(index); else return "";}
bool KUCMSBranchManager::getBOOLBranchValue( std::string key, uInt index = 0 ){
                            if(valid(key)) return theBranches[key].getBOOLValue(index); else return 0;}

float KUCMSBranchManager::getBranchMaxValue( std::string key ){ if(valid(key)) return theBranches[key].getMaxValue(); else return 0;}
uInt KUCMSBranchManager::getBranchLeadIndex( std::string key ){ if(valid(key)) return theBranches[key].getLeadIndex();else return 0;}
uInt KUCMSBranchManager::getBranchSubLeadIndex( std::string key ){ if(valid(key)) return theBranches[key].getSubLeadIndex(); else return 0;}

#endif
