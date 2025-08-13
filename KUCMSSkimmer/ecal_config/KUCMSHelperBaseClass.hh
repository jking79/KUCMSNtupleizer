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
#include <cstdlib>
#include <utility>

// basic C++ types
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <cmath>
#include <limits>
#include <algorithm>
#include <sys/stat.h>
#include <cmath>

#include <algorithm>
#include <random>

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
// Helper functions ( single line function defs, mostly )
//

class KUCMSHelperBaseClass {

	public :

	//KUCMSHelperBaseClass(){};
	//~KUCMSHelperBaseClass(){};

	//auto sortByPt = [](auto & obj1, auto & obj2) {return obj1.pt() > obj2.pt();};

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
	double crystalball_function(double x, double alpha, double n, double sigma, double mean) {

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

	double crystalball_function(const double *x, const double *p) {

		// if ((!x) || (!p)) return 0.; // just a precaution
		// [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
		return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));

	}//<<>>double crystalball_function(const double *x, const double *p)

	double twosided_crystalball_function(double x, double alphalow, double alphahigh, double nlow, double nhigh, double sigma, double mean) {

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

	double twosided_crystalball_function(const double *x, const double *p) {

		return (p[0] * twosided_crystalball_function(x[0], p[3], p[4], p[5], p[6], p[2], p[1]));

	}//<<>>double double_crystalball_function(const double *x, const double *p)

	const float getMyAngle ( const float x, const float y){

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

	const float getAngle ( const float x, const float y){

		if( x == 0 && y == 0) return 6.39;
		auto a = std::atan2(y,x);
		if( a < 0 ) a = 2*PI+a;
		return a;

	}//<<>> const float getAngle (CFlt x, CFlt y) with atan2

	const float getATan2 ( const float x, const float y){

		if( x == 0 && y == 0) return 6.39;
		else return std::atan2(y,x);

	}//<<>> const float getAngle (CFlt x, CFlt y) with atan2

	// math functions
	const auto pow		(CFlt x, int p){ float r = x; for( int i = 1; i < p; i++ ){ r *= x; } return (p>0)?r:0; }
	const auto pow      (CDbl x, int p){ double r = x; for( int i = 1; i < p; i++ ){ r *= x; } return (p>0)?r:0; }
	const auto accum    (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f);}
	const auto sq2      (CFlt x){return x*x;}
	const auto sq2      (CDbl x){return x*x;}
	const auto rad2     (CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
	const auto hypo     (CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
	const auto phi      (CFlt x, CFlt y){return std::atan2(y,x);}
	const auto theta    (CFlt r, CFlt z){return std::atan2(r,z);}
	const auto eta      (CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
	const auto effMean  (CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}
	const auto dIPhi    (CFlt x, CFlt y){auto dp(x-y); if( dp > 180 ){dp-=360.0;} else if( dp < -180 ){ dp+=360.0;} return dp;}
	const auto dPhi     (CFlt x, CFlt y){auto dp(x-y); if( dp > PI ){dp-=2*PI;} else if( dp < -PI ){ dp+=2*PI;} return dp;}
	const auto vfsum    (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f);}
	const auto max      (CVFlt x){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}
	const auto dR1		(CFlt e1, CFlt p1, CFlt e2, CFlt p2){auto de = e1-e2; auto dp = dPhi(p1,p2); return hypo(de,dp);}

	// stats functions
	const auto mean     (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f)/x.size();}
	const auto mean     (CVFlt x, CFlt w){return std::accumulate(x.begin(),x.end(),0.0f)/w;}
	const auto mean     (CVFlt x, CVFlt wv){
							float sum(0.0), wt(0.0); int it(0); for( auto ix : x ){ sum+=ix*wv[it]; wt+=wv[it]; it++; } return sum/wt;}
	const auto wnum     (CFlt it, CFlt w){return (((it-1)*w)/it);}
	const auto stdev    (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return std::sqrt(sum/(x.size()-1));}
	const auto var      (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return sum/(x.size()-1);}

	const auto stdev    (CVFlt x, CFlt m, CVFlt wv, CFlt w){
							float sum(0.0); int it(0);
							for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; }
							return std::sqrt(sum/wnum(it,w));
						}//const auto stdev

	const auto var(CVFlt x, CFlt m, CVFlt wv, CFlt w){
		float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,w);}
	const auto var(CVFlt x, CFlt m, CVFlt wv){
		float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,vfsum(wv));}

	const auto cvar     (CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w){
							float sum(0.0); int it(0);
							for(auto ix : x ){ sum += wv[it]*(ix-mx)*(y[it]-my); it++; }
							return sum/wnum(it,w);
						}//<<>> const auto cvar(CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w)

	const auto cvar     (CVFlt x, CFlt mx, CVFlt y, CFlt my){
							float sum(0.0); int it(0);
							for( auto ix : x ){ sum += (ix-mx)*(y[it]-my); it++; }
							return sum/(x.size()-1);
						}//const auto cvar

	const auto rms      (CVFlt x){float sum(0.0); for(auto ix : x ){ sum += sq2(ix); } return std::sqrt(sum/x.size());}
	const auto chisq    (CVFlt x, CFlt m, CFlt v){ float chi(0); for(auto ix : x ){ chi += sq2(ix-m)/v; } return chi; }

	const auto meanIPhi (CVFlt x){
							float sum(0.0);
							auto maxphi = max(x);
							for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360); else sum+=ix; }
							auto rslt = sum/x.size();
							if( rslt > 360 ) rslt-=360;
							return rslt;
						}//<<>> const auto meanPhi(CVFlt x)

	const auto meanIPhi (CVFlt x, CVFlt wv){
							float wt(0.0), sum(0.0); int it(0);
							auto maxphi = max(x);
							for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360)*wv[it]; else sum+=ix*wv[it]; wt+=wv[it]; it++; }
							auto rslt = sum/wt;
							if( rslt > 360 ) rslt-=360;
							return rslt;
						 }//<<>> const auto meanPhi(CVFlt x, CVFlt wv)

	const auto wsin2    (CVFlt x, CVFlt wv){
							double sum(0.0), wt(0.0); int it(0);
							for(auto ix : x ){
								sum += wv[it]*sq2(sin(ix));
								wt += wv[it];
								it++;
							}//for(auto ix : x )
							return sum/wt;
						}//const auto wsin2(CVFlt x, CVFlt wv)

	const auto wcos2    (CVFlt x, CVFlt wv){
							double sum(0.0), wt(0.0); int it(0);
							for(auto ix : x ){ sum += wv[it]*sq2(cos(ix)); wt += wv[it]; it++;}
							return sum/wt;
						}//const auto wcos2

	const auto wsincos  (CVFlt x, CVFlt wv){
							double sum(0.0), wt(0.0); int it(0);
							for(auto ix : x ){ sum += wv[it]*sin(ix)*cos(ix); wt += wv[it]; it++;}
							return sum/wt;
						}//const auto wsincos

	//-----------  misc helper functions -------------------------------------------

	std::string addstr( std::string current, std::string input ){ return (current+input); }

	std::string RemoveDelim(std::string tmp, const std::string & delim){return tmp.erase(tmp.find(delim),delim.length());}

	const auto splitString( std::string str, const char* separator ) {

		std::vector < std::string > strings;
		int startIndex(0), endIndex(0);
		for ( int i = 0; i <= str.size(); i++ ){
			if ( str[i] == *separator || i == str.size() ){

				endIndex = i;
				std::string temp;
				temp.append(str, startIndex, endIndex - startIndex);
				strings.push_back(temp);
				startIndex = endIndex + 1;

			}//<<>>if (str[i] == separator || i == str.size())
		}//<<>>for (int i = 0; i <= str.size(); i++)
		return strings;

	}//<<>>const auto  splitString(string str, char separator)

	template <typename T> std::string toString(T value){

	   std::ostringstream os ;
	   os << value ;
	   return os.str() ;

	}//<<>>template <typename T> std::string to_string(T value)

	void setBins(std::string str, std::vector<double> & bins ){

		if(str.find("CONSTANT") != std::string::npos){

			bool var_bins = false;
			str = RemoveDelim(str,"CONSTANT");
			int nbins = 0; double low = 0.f, high = 0.f;
			std::stringstream ss(str);
			ss >> nbins >> low >> high;
			double bin_width = (high-low)/nbins;
			for (int ibin = 0; ibin <= nbins; ibin++){ bins.push_back(low+ibin*bin_width); }

		} else if(str.find("VARIABLE") != std::string::npos) {

			bool var_bins = true;
			str = RemoveDelim(str,"VARIABLE");
			double bin_edge;
			std::stringstream ss(str);
			std::cout << "Setting Var bins : ";
			while(ss >> bin_edge){
				std::cout << bin_edge << " ";
				bins.push_back(bin_edge);
			}//<<>>while (ss >> bin_edge)
			std::cout << std::endl;

		} else {

			std::cerr << "Aye... bins are either VARIABLE or CONSTANT! Exiting..." << std::endl;
			exit(1);

		}//<<>>if      (str.find("CONSTANT") != std::string::npos)

	}//<<>>void setBins(std::string & str, std::vector<Double_t> & bins, Bool_t & var_bins)

    void setBins(std::string str, std::vector<float> & bins ){

        if(str.find("CONSTANT") != std::string::npos){

            bool var_bins = false;
            str = RemoveDelim(str,"CONSTANT");
            int nbins = 0; float low = 0.f, high = 0.f;
            std::stringstream ss(str);
            ss >> nbins >> low >> high;
            float bin_width = (high-low)/nbins;
            for (int ibin = 0; ibin <= nbins; ibin++){ bins.push_back(low+ibin*bin_width); }

        } else if(str.find("VARIABLE") != std::string::npos) {

            bool var_bins = true;
            str = RemoveDelim(str,"VARIABLE");
            float bin_edge;
            std::stringstream ss(str);
            std::cout << "Setting Var bins : ";
            while(ss >> bin_edge){
                std::cout << bin_edge << " ";
                bins.push_back(bin_edge);
            }//<<>>while (ss >> bin_edge)
            std::cout << std::endl;

        } else {

            std::cerr << "Aye... bins are either VARIABLE or CONSTANT! Exiting..." << std::endl;
            exit(1);

        }//<<>>if      (str.find("CONSTANT") != std::string::npos)

    }//<<>>void setBins(std::string & str, std::vector<Double_t> & bins, Bool_t & var_bins)

};//<<>> class KUCMSHelperBaseClass 

#endif
//-------------------------------------------------------------------------------------------------------------------

