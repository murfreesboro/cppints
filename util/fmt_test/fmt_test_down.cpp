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
	int top_M_limmit   = 40;
	double T_max_limit = 40.0E0;
	float steplength_f  = 0.00001;
	float globalError_f = 1.0E-6;

	///////////////////////////////////////////////////
	// this section testing the timing performance
	// between the different functions
	// so far, erf is 160 faster than the fmt function
	// in boost library
	///////////////////////////////////////////////////
	cout << "===========================================================" << endl;
	cout << "Time performance result: " << endl;
	cout << "Total sample number is " << N << endl;
	cout << "===========================================================" << endl;
	// testing the erf
	timer t1;
	for(int i=0; i<N; i++) {
		double a = i;
		double x = a/1000;
		double e = erf(x);
	}
	printf("erf time consuming: %-14.7f\n", t1.elapsed());

	// testing the sqrt
	timer t2;
	for(int i=0; i<N; i++) {
		double e = sqrt(i);
	}
	printf("sqrt time consuming: %-14.7f\n", t2.elapsed());
 
	// testing the exp
	timer t3;
	for(int i=0; i<N; i++) {
		double T = -(i+1)/N;
		double e = exp(T);
	}
	printf("exp(-T) time consuming: %-14.7f\n", t3.elapsed());

	// testing the fmt
	timer t4;
	for(int m=0; m<25; m++) {
		int n = N/25;
		for(int i=0; i<n; i++) {
			double T = 0.0001*(i+1);
			double e = fm(T,m);
		}
	}
	printf("fm(T,m) time consuming: %-14.7f\n", t4.elapsed());

	///////////////////////////////////////////////////
	// testing the fmt from down recursive relation
	///////////////////////////////////////////////////
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the down recursive relation for m = 1 to " << top_M_limmit << endl;
	cout << "this is for double type variable" << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from 0 " << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	int nSteps = T_max_limit/steplength;
	for (int n=1; n<=top_M_limmit; n++) {
		for (int j = 0; j<nSteps; j++) {
			double T = steplength*j;
			double fmt= fm(T,n);
			double et = exp(-T);
			double e0 = fmt;
			double t2 = 2.0*T;
			for(int i=n-1; i>=0; i--) {
				double s = fm(T,i);
				double e = (1.0/(2*i+1.0))*(t2*e0+et);
				double d = fabs(s-e);
				if (d>globalError) {
					printf("T is %-16.14f\n", T);
					printf("for m=%d difference is %-16.14f\n", n, d);
				}
				e0 = e;
			}
		}
	}

	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the down recursive relation for m = 1 to " << top_M_limmit << endl;
	cout << "this is for float type variable" << endl;
	cout << "error range is " << globalError_f << endl;
	cout << "step length for T is " << steplength_f << endl;
	cout << "T is ranging from 0 " << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	nSteps = T_max_limit/steplength_f;
	for (int n=1; n<=top_M_limmit; n++) {
		for (int j = 0; j<nSteps; j++) {
			float T = steplength_f*j;
			float fmt= fm(T,n);
			float et = exp(-T);
			float e0 = fmt;
			float t2 = 2.0*T;
			for(int i=n-1; i>=0; i--) {
				float s = fm(T,i);
				float e = (1.0/(2*i+1.0))*(t2*e0+et);
				float d = fabs(s-e);
				if (d>globalError_f) {
					printf("T is %-10.7f\n", T);
					printf("for m=%d difference is %-10.7f\n", n, d);
				}
				e0 = e;
			}
		}
	}

	return 0;
}
