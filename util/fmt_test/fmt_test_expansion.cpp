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
	int N = 10000000;
	double steplength  = 0.000001;
	double globalError = 1.0E-14;
	int M_limmit       = 40;
	double T_min_limit = 1.8E0;
   int nTerms = 17;
	float steplength_f  = 0.00001;
	float globalError_f = 1.0E-6;

	///////////////////////////////////////////////////
	// testing the power series
	// equation of 9 in Harris paper
	// We try to expand the power series 
	// here we test from m=1 to m=10
	// for T from 0 to 1.4
	///////////////////////////////////////////////////
	int numSteps = T_min_limit/steplength;
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the power series expansion with " << nTerms+1 << " terms" << endl;
	cout << "this is for double type variable" << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from 0 to " << T_min_limit << endl;
	cout << "===========================================================" << endl;
	for (int m = 1; m<=M_limmit; m++) {
		for (int j = 0; j<=numSteps; j++) {
			double T = steplength*j;
			double et = exp(-T);
			double x  = 1.0;
			double t2 = 2.0*T;
			for(int i=nTerms; i>=1; i--) {
				x = 1.0+t2/(2.0*(m+i)+1.0)*x;
			}
			x  = (1.0/(2*m+1))*x;
			x *= et;
			double result = fm(T,m);
			double d = fabs(result-x);
			if (d>globalError) {
				printf("fmt is %-16.14f\n", result);
				printf("for m=%d T=%f difference is %-16.14f\n", m, T, d);
			}
		}
	}

	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the power series expansion with " <<nTerms+1 << " terms" << endl;
	cout << "this is for float type variable" << endl;
	cout << "error range is " << globalError_f << endl;
	cout << "step length for T is " << steplength_f << endl;
	cout << "T is ranging from 0 to " << T_min_limit << endl;
	cout << "===========================================================" << endl;
	numSteps = T_min_limit/steplength_f;
	for (int m = 1; m<=M_limmit; m++) {
		for (int j = 0; j<=numSteps; j++) {
			float T  = steplength_f*j;
			float et = exp(-T);
			float x  = 1.0;
			float t2 = 2.0*T;
			for(int i=nTerms; i>=1; i--) {
				x = 1.0+t2/(2.0*(m+i)+1.0)*x;
			}
			x  = (1.0/(2*m+1))*x;
			x *= et;
			double result = fm(T,m);
			float d = fabs(result-x);
			if (d>globalError_f) {
				printf("fmt is %-16.14f\n", result);
				printf("for m=%d T=%f difference is %-16.14f\n", m, T, d);
			}
		}
	}

	return 0;
}
