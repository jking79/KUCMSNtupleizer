// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// system include files
#include <memory>

// basic C++ types
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <tuple>
#include <random>
#include <sys/stat.h>
#include <cmath>

#ifndef KUCMSHelperHeader
#define KUCMSHelperHeader

typedef unsigned int uInt;
typedef unsigned long int ulInt;
typedef const float CFlt;
typedef const double CDbl;
typedef const std::vector<float> CVFlt;
typedef const std::vector<double> CVDbl;

#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.14159265358979323846 // pie ...  

//
// Helper classes
//


template < class T* >
class hConstPtrWrap {

	public:  

    hConstPtrWrap() : use( NULL ) {};
	hConstPtrWrap( const T* it ) : use(it) {};

    const T* use;

};//<<>>class jwrap 

class v3fPoint {

	public:

	v3fPoint( float x, float y, float z ){ vx = x; vy = y; vz = z; };

	float x() const { return vx; };
    float y() const { return vy; };
    float z() const { return vz; };

	private: 

	float vx;
	float vy;
	float vz;

};//<<>>class v3fPoint 

//
// Helper functions ( single line function defs, mostly )
//

inline const auto sortByPt = [](auto & obj1, auto & obj2) {return obj1.pt() > obj2.pt();};

//
// The "crystalball" function for ROOT 5.x (mimics ROOT 6.x).
//
// Create the "crystalball" TF1 somewhere in your source code using:
// double xmin = 3., xmax = 8.; // whatever you need
// TF1 *crystalball = new TF1("crystalball", crystalball_function, xmin, xmax, 5);
// crystalball->SetParNames("Constant", "Mean", "Sigma", "Alpha", "N");
// crystalball->SetTitle("crystalball"); // not strictly necessary
//

// see math/mathcore/src/PdfFuncMathCore.cxx in ROOT 6.x
inline double crystalball_function(double x, double alpha, double n, double sigma, double mean) {

    // evaluate the crystal ball function
    if (sigma < 0.)     return 0.;
    double z = (x - mean)/sigma;
    if (alpha < 0) z = -z;
    double abs_alpha = std::abs(alpha);
    // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
    // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
    // double N = 1./(sigma*(C+D));
    if (z  > - abs_alpha) return std::exp(- 0.5 * z * z);
    else {
        //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
        double nDivAlpha = n/abs_alpha;
        double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
        double B = nDivAlpha -abs_alpha;
        double arg = nDivAlpha/(B-z);
        return AA * std::pow(arg,n);
    }//<<>> if (z  > - abs_alpha)

}//<<>>double crystalball_function(double x, double alpha, double n, double sigma, double mean)

inline double crystalball_function(const double *x, const double *p) {

    // if ((!x) || (!p)) return 0.; // just a precaution
    // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
    return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));

}//<<>>double crystalball_function(const double *x, const double *p)

inline double twosided_crystalball_function(double x, double alphalow, double alphahigh, double nlow, double nhigh, double sigma, double mean) {

    // evaluate the crystal ball function
    if(sigma < 0.) return 0.;
    auto z = (x - mean)/sigma;
    auto abs_alphalow = std::abs(alphalow);
    auto abs_alphahigh = std::abs(alphahigh);
    if( z  > abs_alphalow and z  < abs_alphahigh ) return std::exp(- 0.5 * z * z);
    else if( z  > abs_alphahigh ){
        double nDivAlpha = nhigh/abs_alphahigh;
        double AA =  std::exp(-0.5*abs_alphahigh*abs_alphahigh);
        double B = nDivAlpha -abs_alphahigh;
        double arg = nDivAlpha/(B-z);
        return AA * std::pow(arg,nhigh);
    } else {
        double nDivAlpha = nlow/abs_alphalow;
        double AA =  std::exp(-0.5*abs_alphalow*abs_alphalow);
        double B = nDivAlpha -abs_alphalow;
        double arg = nDivAlpha/(B-z);
        return AA * std::pow(arg,nlow);
    }//<<>>if( z  > abs_alphalow and z  < abs_alphahigh )

}//<<>>double double_crystalball_function(double x, ... double nlow, double nhigh, double sigma, double mean)

inline double twosided_crystalball_function(const double *x, const double *p) {

    return (p[0] * twosided_crystalball_function(x[0], p[3], p[4], p[5], p[6], p[2], p[1]));

}//<<>>double double_crystalball_function(const double *x, const double *p)

inline const float getMyAngle ( const float x, const float y){

    if( x == 0 && y == 0 ) return 6.39; //-999.0;
    if( x == 0 ){ if( y > 0 ) return PI/2; else return 3*PI/2; }
    if( y == 0 ){ if( x > 0 ) return 0; else return PI; }
    float m = std::sqrt(x*x+y*y);
    float a = std::asin(abs(y/m));
    if( x < 0 && y < 0 ) return PI+a;
    if( x < 0 ) return PI-a;
    if( y < 0 ) return 2*PI-a;
    return a;

}//<<>> const float getAngle (CFlt x, CFlt y)

inline const float getAngle ( const float x, const float y){

    if( x == 0 && y == 0) return 6.39;
    auto a = std::atan2(y,x);
    if( a < 0 ) a = 2*PI+a;
    return a;

}//<<>> const float getAngle (CFlt x, CFlt y) with atan2

inline const float getATan2 ( const float x, const float y){

    if( x == 0 && y == 0) return 6.39;
    else return std::atan2(y,x);

}//<<>> const float getAngle (CFlt x, CFlt y) with atan2

// math functions
inline const auto sq2      (CFlt x){return x*x;}
inline const auto sq2      (CDbl x){return x*x;}
inline const auto rad2     (CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
inline const auto hypo     (CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
inline const auto phi      (CFlt x, CFlt y){return std::atan2(y,x);}
inline const auto theta    (CFlt r, CFlt z){return std::atan2(r,z);}
inline const auto eta      (CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
inline const auto effMean  (CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}
inline const auto dIPhi    (CFlt x, CFlt y){auto dp(x-y); if( dp > 180 ){dp-=360.0;} else if( dp < -180 ){ dp+=360.0;} return dp;}
inline const auto vfsum    (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f);}
inline const auto max      (CVFlt x){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

// stats functions
inline const auto mean     (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f)/x.size();}
inline const auto mean     (CVFlt x, CFlt w){return std::accumulate(x.begin(),x.end(),0.0f)/w;}
inline const auto mean     (CVFlt x, CVFlt wv){
                        float sum(0.0), wt(0.0); int it(0); for( auto ix : x ){ sum+=ix*wv[it]; wt+=wv[it]; it++; } return sum/wt;}
inline const auto wnum     (CFlt it, CFlt w){return (((it-1)*w)/it);}
inline const auto stdev    (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return std::sqrt(sum/(x.size()-1));}
inline const auto var      (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return sum/(x.size()-1);}

inline const auto stdev    (CVFlt x, CFlt m, CVFlt wv, CFlt w){
                        float sum(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; }
                        return std::sqrt(sum/wnum(it,w));
                    }//inline const auto stdev

inline const auto var(CVFlt x, CFlt m, CVFlt wv, CFlt w){
    float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,w);}
inline const auto var(CVFlt x, CFlt m, CVFlt wv){
    float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,vfsum(wv));}

inline const auto cvar     (CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w){
                        float sum(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*(ix-mx)*(y[it]-my); it++; }
                        return sum/wnum(it,w);
                    }//<<>> inline const auto cvar(CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w)

inline const auto cvar     (CVFlt x, CFlt mx, CVFlt y, CFlt my){
                        float sum(0.0); int it(0);
                        for( auto ix : x ){ sum += (ix-mx)*(y[it]-my); it++; }
                        return sum/(x.size()-1);
                    }//inline const auto cvar

inline const auto rms      (CVFlt x){float sum(0.0); for(auto ix : x ){ sum += sq2(ix); } return std::sqrt(sum/x.size());}
inline const auto chisq    (CVFlt x, CFlt m, CFlt v){ float chi(0); for(auto ix : x ){ chi += sq2(ix-m)/v; } return chi; }

inline const auto meanIPhi (CVFlt x){
                        float sum(0.0);
                        auto maxphi = max(x);
                        for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360); else sum+=ix; }
                        auto rslt = sum/x.size();
                        if( rslt > 360 ) rslt-=360;
                        return rslt;
                    }//<<>> inline const auto meanPhi(CVFlt x)

inline const auto meanIPhi (CVFlt x, CVFlt wv){
                        float wt(0.0), sum(0.0); int it(0);
                        auto maxphi = max(x);
                        for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360)*wv[it]; else sum+=ix*wv[it]; wt+=wv[it]; it++; }
                        auto rslt = sum/wt;
                        if( rslt > 360 ) rslt-=360;
                        return rslt;
                     }//<<>> inline const auto meanPhi(CVFlt x, CVFlt wv)

inline const auto wsin2    (CVFlt x, CVFlt wv){
                        double sum(0.0), wt(0.0); int it(0);
                        for(auto ix : x ){
                            sum += wv[it]*sq2(sin(ix));
                            wt += wv[it];
                            it++;
                        }//for(auto ix : x )
                        return sum/wt;
                    }//inline const auto wsin2(CVFlt x, CVFlt wv)

inline const auto wcos2    (CVFlt x, CVFlt wv){
                        double sum(0.0), wt(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*sq2(cos(ix)); wt += wv[it]; it++;}
                        return sum/wt;
                    }//inline const auto wcos2

inline const auto wsincos  (CVFlt x, CVFlt wv){
                        double sum(0.0), wt(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*sin(ix)*cos(ix); wt += wv[it]; it++;}
                        return sum/wt;
                    }//const auto wsincos

#endif
//-------------------------------------------------------------------------------------------------------------------

