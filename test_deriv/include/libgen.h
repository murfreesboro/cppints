#ifndef LIBGEN_H
#define LIBGEN_H
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
#include<iterator>
#include<algorithm>
using namespace std;

/**
 * re-define the types here
 * we do not use unsigned type for int, since the common int
 * should be good enough
 */
typedef int           Int;
typedef long long    LInt;
typedef size_t          UInt; 
#ifdef WITH_SINGLE_PRECISION
typedef float         Double;
#else
typedef double        Double;
#endif

// general constant used in the program
#define  MINUS_HALF   -0.5E0
#define  MINUS_ONE    -1.0E0
#define  MINUS_TWO    -2.0E0
#define  MINUS_THREE  -3.0E0
#define  MINUS_FOUR   -4.0E0
#define  MINUS_FIVE   -5.0E0
#define  MINUS_SIX    -6.0E0
#define  MINUS_SEVEN  -7.0E0
#define  MINUS_EIGHT  -8.0E0
#define  MINUS_NINE   -9.0E0
#define  MINUS_TEN    -10.0E0
#define  ZERO          0.0E0
#define  ONE           1.0E0
#define  TWO           2.0E0
#define  THREE         3.0E0
#define  FOUR          4.0E0
#define  FIVE          5.0E0
#define  SIX           6.0E0
#define  SEVEN         7.0E0
#define  EIGHT         8.0E0
#define  NINE          9.0E0
#define  TEN           10.0E0
#define  HALF          0.5E0
#define  PI            3.141592653589793239

#ifdef WITH_SINGLE_PRECISION
#define  THRESHOLD_MATH  1.0E-6
#define  THRESH          1.0E-5
#define  INTEGRAL_THRESH 1.0E-6
#else
#define  THRESHOLD_MATH  1.0E-14
#define  THRESH          1.0E-10
#define  INTEGRAL_THRESH 1.0E-10
#endif

// define the positions
#define  BRA1           1
#define  BRA2           2
#define  KET1           3
#define  KET2           4
#define  DERIV_X        5
#define  DERIV_Y        6
#define  DERIV_Z        7
#define  NO_DERIV       -1 

/**
 * \brief crash function is used to deal with the situation that something should not
 *        happen but happened, then we stop the whole program and give the information
 */
inline void crash(bool isCrash, string message) {
	if (isCrash) {
		cout<< "Fatal error occurs:" << endl;
		cout<< message << endl;
		exit(1);
	}
};

#endif

