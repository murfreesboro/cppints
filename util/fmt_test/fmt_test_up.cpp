//
// This program is used to study how to calculate (SS|SS)^{m}
// for the recursive relation. Basically, it focus on the 
// fmt function:
// f_{m}(t) = \int^{1}_{0} u^{2m} e^{-tu^{2}} du 
// it needs the boost library to calcualte fmt
// fenglai liu
// Oct. 2013
//
// for reference, see Harris paper and the reference cited inside:
// Harris, Frank E
// Evaluation of GTO molecular integrals
// International Journal of Quantum Chemistry, 1983, Vol. 23, Page 1469--1478
//
// note:
// we found that the approximation of (2m-1)!!/R expression is 
// not numerical stable. Combined with down recursive relation,
// it introduces greater error compared with other methods. 
// Therefore we finally abadon it.
//
//
//

// common C head files used in the program
#include<cstdlib>
#include<cstdio>

// external math functions 
#include<cmath>

// common C++ head files may be used in the program
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
#include<list>
#include<set>
#include<map>
#include<iterator>
#include<algorithm>

#include <boost/math/special_functions/binomial.hpp>  // special functions in math
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/timer.hpp>
using namespace boost::math;
using namespace boost;
using namespace std;
#define PI              3.1415926535897932384626E0  

double fm(const double& a, const int& m); 
double doubleFac(const int& m, const int& n); 

double fm(const double& a, const int& m) {
	double x = 1.0E0;
	if (fabs(a)<1.0E-14) return (pow(x,2*m+1)/(2*m+1));
	double f12 = 1.0E0/2.0E0;
	double p   = 1.0E0/(2.0E0*pow(a,m+f12));
	double upperLimit = a*pow(x,2.0E0);
	return tgamma_lower(m+f12,upperLimit)*p;
}

double doubleFac(const int& m, const int& n) {
	// calculate (2m-1)!!/(2m+2n+1)!!
	double x = 1.0E0;
	for(int i=0; i<=n; i++) {
		x *= 1.0E0/(2*m+2*i+1); 
	}
	return x;
} 

int main() 
{
	///////////////////////////////////////////////////
	// global settings
	// step length and number of claculation for 
	// time performance
	// global error is the error range
	// if the difference is less than the error, we 
	// think the error could be omitted
	///////////////////////////////////////////////////
	double steplength  = 0.000001;
	double globalError = 1.0E-14;
	int M_limmit       = 8;
	double T_max_limit = 40.0E0;
	double T_min_limit = 1.8E0;
	int numSteps       = (T_max_limit-T_min_limit)/steplength;
	double largest     = 0.0E0;

	/////////////////////////////////////////////////////////////////////
	// this section testing the up recursice relation for M_limit = 8  //
	/////////////////////////////////////////////////////////////////////
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the up recursive relation for m = 1 to " << M_limmit << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from " << T_min_limit << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	for (int j = 0; j<numSteps; j++) {
		double T = T_min_limit+steplength*j;
		double sqT = sqrt(T);
		double fmt= erf(sqT)*sqrt(PI)/2.0E0;
		double et = exp(-T);
		double e0 = fmt/sqT;
		double oneO2T = 0.5E0/T;
		for(int i=1; i<=M_limmit; i++) {
			double s = fm(T,i);
			double e = oneO2T*((2*i-1)*e0-et);
			double d  = fabs(s-e);
			double d0 = s-e;
			if (d>globalError) {
				printf("m is  %d, T   is %-16.14f, difference is %-16.14f\n", i, T, d0);
				if (d>largest) largest = d;
			}
			e0 = e;
		}
	}
	printf("largest difference is %-16.14f\n", largest);

	/////////////////////////////////////////////////////////////////////
	// this section testing the up recursice relation for M_limit = 9  //
	/////////////////////////////////////////////////////////////////////
	cout << endl << endl;
	globalError = 1.0E-13;
	M_limmit    = 9;
	cout << "===========================================================" << endl;
	cout << "testing the up recursive relation for m = 1 to " << M_limmit << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from " << T_min_limit << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	largest = 0.0E0;
	for (int j = 0; j<numSteps; j++) {
		double T = T_min_limit+steplength*j;
		double sqT = sqrt(T);
		double fmt= erf(sqT)*sqrt(PI)/2.0E0;
		double et = exp(-T);
		double e0 = fmt/sqT;
		double oneO2T = 0.5E0/T;
		for(int i=1; i<=M_limmit; i++) {
			double s = fm(T,i);
			double e = oneO2T*((2*i-1)*e0-et);
			double d  = fabs(s-e);
			double d0 = s-e;
			if (d>globalError) {
				printf("m is  %d, T   is %-16.14f, difference is %-16.14f\n", i, T, d0);
				if (d>largest) largest = d;
			}
			e0 = e;
		}
	}
	printf("largest difference is %-16.14f\n", largest);

	/////////////////////////////////////////////////////////////////////
	// this section testing the up recursice relation for M_limit = 10 //
	/////////////////////////////////////////////////////////////////////
	cout << endl << endl;
	globalError = 1.0E-12;
	M_limmit    = 10;
	cout << "===========================================================" << endl;
	cout << "testing the up recursive relation for m = 1 to " << M_limmit << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from " << T_min_limit << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	largest = 0.0E0;
	for (int j = 0; j<numSteps; j++) {
		double T = T_min_limit+steplength*j;
		double sqT = sqrt(T);
		double fmt= erf(sqT)*sqrt(PI)/2.0E0;
		double et = exp(-T);
		double e0 = fmt/sqT;
		double oneO2T = 0.5E0/T;
		for(int i=1; i<=M_limmit; i++) {
			double s = fm(T,i);
			double e = oneO2T*((2*i-1)*e0-et);
			double d  = fabs(s-e);
			double d0 = s-e;
			if (d>globalError) {
				printf("m is  %d, T   is %-16.14f, difference is %-16.14f\n", i, T, d0);
				if (d>largest) largest = d;
			}
			e0 = e;
		}
	}
	printf("largest difference is %-16.14f\n", largest);

	return 0;
}
