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
	double globalError = 1.0E-10;
	int M_limmit       = 10;
	int top_M_limmit   = 30;
	double T_max_limit = 36.0E0;
	double T_min_limit = 1.1E0;

	float steplength_f  = 0.00001;
	float globalError_f = 1.0E-6;

	/*
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
	*/

	/*
	///////////////////////////////////////////////////
	// this section testing the up recursice relation
	// can we do it for m = 1 to 12?
	// 10 is perfect fine, however; 12 will have 
	// error around 1.0E-10
	// we note, that up recursive relation can not
	// go with float
	///////////////////////////////////////////////////
	int numSteps = T_max_limit/steplength;
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the up recursive relation for m = 1 to " << M_limmit << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from " << T_min_limit << " to " <<T_max_limit << endl;
	cout << "===========================================================" << endl;
	for (int j = 0; j<numSteps; j++) {
		double T = T_min_limit+steplength*(j+1);
		double sqT = sqrt(T);
		double fmt= erf(sqT)*sqrt(PI)/2.0E0;
		double et = exp(-T);
		double e0 = fmt/sqT;
		double oneO2T = 0.5E0/T;
		for(int i=1; i<=M_limmit; i++) {
			double s = fm(T,i);
			double e = oneO2T*((2*i-1)*e0-et);
			double d = fabs(s-e);
			if (d>globalError) {
				printf("fmt is %-14.12f\n", s);
				printf("T   is %-14.12f\n", T);
				printf("error/fmt is %-14.12f\n", d/fmt);
				printf("for m=%d difference is %-14.12f\n", i, d);
			}
			e0 = e;
		}
	}
	*/

	/*
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
			double T = steplength*(j+1);
			double fmt= fm(T,n);
			double et = exp(-T);
			double e0 = fmt;
			double t2 = 2.0*T;
			for(int i=n-1; i>=0; i--) {
				double s = fm(T,i);
				double e = (1.0/(2*i+1.0))*(t2*e0+et);
				double d = fabs(s-e);
				if (d>globalError) {
					printf("T is %-14.12f\n", T);
					printf("for m=%d difference is %-14.12f\n", n, d);
				}
				e0 = e;
			}
		}
	}
	*/

	/*
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
			float T = steplength_f*(j+1);
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
	*/

	/*
	///////////////////////////////////////////////////
	// testing the fmt from down recursive relation
	// this is not so accurate
	// abandoned, too
	///////////////////////////////////////////////////
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the down recursive relation for m = 1 to 10" << endl;
	cout << "the fm(t) is approximated by (2m-1)!!/(2t)^{m}*t^{1/2}" << endl;
	cout << "error range is 1.0E-10" << endl;
	cout << "step length for T is " << steplength << " from 30.0 to 50.0" << endl;
	cout << "===========================================================" << endl;
	int numSteps = (50-30)/steplength;
	double globalError = 1.0E-9;
	for (int n=1; n<=10; n++) {
		for (int j = 0; j<numSteps; j++) {

			// calculate fm(t)
			double T = 30.0E0+steplength*(j+1);
			double fmt= sqrt(PI)/(2*sqrt(T));
			for(int k=n; k>=1; k--) {
				fmt *= (2*k-1)/(2.0*T);
			}
			double s0 = fm(T,n);
			double d0 = fabs(s0-fmt);
			if (d0>globalError) {
				printf("T is %-14.12f\n", T);
				printf("error/fmt is %-14.12f\n", d0/fmt);
				printf("for m=%d difference is %-14.12f\n", n, d0);
			}

			// recursive relation
			double et = exp(-T);
			double e0 = fmt;
			double t2 = 2.0*T;
			for(int i=n-1; i>=0; i--) {
				double s = fm(T,i);
				double e = (1.0/(2*i+1.0))*(t2*e0+et);
				double d = fabs(s-e);
				if (d>globalError) {
					printf("T is %-14.12f\n", T);
					printf("error/fmt is %-14.12f\n", d/s);
					printf("for m=%d difference is %-14.12f\n", n, d);
				}
				e0 = e;
			}
		}
	}
	*/

	/*
	///////////////////////////////////////////////////
	// testing the power series
	// equation of 9 in Harris paper
	// We try to expand the power series up to 10 terms
	// here we test from m=1 to m=10
	// for T from 0 to 1.1
	///////////////////////////////////////////////////
	numSteps = T_min_limit/steplength;
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the power series expansion with 13 terms" << endl;
	cout << "this is for double type variable" << endl;
	cout << "error range is " << globalError << endl;
	cout << "step length for T is " << steplength << endl;
	cout << "T is ranging from 0 to " << T_min_limit << endl;
	cout << "===========================================================" << endl;
	for (int m = 1; m<=M_limmit; m++) {
		for (int j = 0; j<numSteps; j++) {
			double T = steplength*(j+1);
			double et = exp(-T);
			double x  = 1.0;
			double t2 = 2.0*T;
			for(int i=12; i>=1; i--) {
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
	*/

	/*
	cout << endl << endl;
	cout << "===========================================================" << endl;
	cout << "testing the power series expansion with 13 terms" << endl;
	cout << "this is for float type variable" << endl;
	cout << "error range is " << globalError_f << endl;
	cout << "step length for T is " << steplength_f << endl;
	cout << "T is ranging from 0 to " << T_min_limit << endl;
	cout << "===========================================================" << endl;
	numSteps = T_min_limit/steplength_f;
	for (int m = 1; m<=M_limmit; m++) {
		for (int j = 0; j<numSteps; j++) {
			float T  = steplength_f*(j+1);
			float et = exp(-T);
			float x  = 1.0;
			float t2 = 2.0*T;
			for(int i=12; i>=1; i--) {
				x = 1.0+t2/(2.0*(m+i)+1.0)*x;
			}
			x  = (1.0/(2*m+1))*x;
			x *= et;
			double result = fm(T,m);
			double d = fabs(result-x);
			if (d>globalError_f) {
				printf("fmt is %-16.14f\n", result);
				printf("for m=%d T=%f difference is %-16.14f\n", m, T, d);
			}
		}
	}
	*/

	////////////////////////////////////////////////////////
	// this is the continuous fraction implementation
	// of Harris equation, not correct yet, already abandoned
	// equation 24
	// testing the series
	// for m = 0, and T is from 0.01 to 25
	// test the convergence for the series
	////////////////////////////////////////////////////////
	/*
	int n = 20;
	for (int m=0; m<=30; m++) {
		for (int j = 0; j<2500; j++) {
			double T = 0.01*(j+1);
			double et = exp(-T);
			double minusn = 1.0E0;
			if ((n+1)%2 == 1) minusn = -1.0E0;
			double an = 2.0*(minusn*(m+n)-(m+1))/(4.0*(m+n)*(m+n)-1);
			double x  = 1+an*T;
			for(int i=n-1; i>1; i--) {
				double minusi = 1.0E0;
				if ((i+1)%2 == 1) minusi = -1.0E0;
				double ai = 2.0*(minusi*(m+i)-(m+1))/(4.0*(m+i)*(m+i)-1);
				x = (1+ai*T)/x;
			}
			double a1 = 2.0/(4.0*(m+1)*(m+1)-1);
			x = a1*T/x;
			x += 1/(2*m+1);
			x *= et;
			double result = fm(T,m);
			double d = fabs(result-x);
			if (d>0.00000000000001) {
				printf("for %d %f difference is %-16.14f\n", m, T, d);
			}
		}
	}
	*/

	return 0;
}
