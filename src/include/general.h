/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2012-2015 Fenglai Liu
 * This softare uses the MIT license as below:
 *
 *	Permission is hereby granted, free of charge, to any person obtaining 
 *	a copy of this software and associated documentation files (the "Software"), 
 *	to deal in the Software without restriction, including without limitation 
 *	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 *	and/or sell copies of the Software, and to permit persons to whom the Software 
 *	is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *						    
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 *	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
 *	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
 *	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 *	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * \file   general.h
 * \brief  primary head files for this program
 * \author Fenglai Liu 
 */
#ifndef GENERAL_H
#define GENERAL_H

// here to avoid the time const in the string class,
// we use an integer to express the "EUQAL", "LESS THAN",
// and "GREATER THAN" etc.
#define FAIL_COMPARE     -1
#define EQ                0
#define LT                1
#define GT                2
#define LE                3
#define GE                4

// define the bra1, bra2 etc.
// these values are defined different with the angular momentum
// changes 
#define BRA1              1000
#define BRA2              2000
#define KET1              3000
#define KET2              4000

//
// let's add in derivatives information
//
#define NO_DERIV          5000
#define DERIV_X           5001
#define DERIV_Y           5002
#define DERIV_Z           5003
#define DERIV_XX          5011
#define DERIV_XY          5012
#define DERIV_XZ          5013
#define DERIV_YY          5014
#define DERIV_YZ          5015
#define DERIV_ZZ          5016

// define the NULL position for shell and basis sets
// it's also used as null value
#define NULL_POS         -1

// define the operator
// MOM is the moment integrals
// esp is the electrostatic potential
// this group is within 100-900
#define TWOBODYOVERLAP    100
#define THREEBODYOVERLAP  200
#define KINETIC           300
#define NAI               400
#define ERI               500
#define TWOBODYERI        510
#define THREEBODYERI      520
#define MOM               600
#define THREEBODYKI       700
#define ESP               800

// define the VRR and HRR work type
// this group is within 10-90
#define OS                10
#define HRR               20

// define the BRA and KET
#define BRA               10000
#define KET               20000

// define the RR property 
// for given operator type
#define SQ_ON_L           100000
#define SQ_ON_L_M         200000
#define SQ_ON_L_OPER      300000

//
// the following status is used to characterize 
// the shell quartet in the RR process in terms
// of printing purpose.
//
// All of shell quartet could be divided into three 
// levels:
//
// 1  tmp result. The tmp result is just temporary 
// result generated. It will be used internally
// inside the VRR/HRR module to generate the 
// module result.
//
// 2  module result. Module result is the VRR/HRR
// output shell quartet and integrals. For example,
// the module results for VRR is input for HRR.
// we note, that here module result is particularly
// assigned to these shell quartets that they are 
// not the final results (so that we can distinguish
// with final results).
//
// 3  final result. final results is a special type 
// of module results (however, it is not included 
// into module results), that they are the final
// output for the whole integral file (they are those
// whose name is abcd).
//
// for how to determine the status of result sq, please
// refer to the sqStatusCheck function in rr.cpp.
//
//
#define TMP_RESULT        1000000
#define MODULE_RESULT     2000000
#define FINAL_RESULT      3000000

// common C head files used in the program
#include<cstdlib>
#include<cstdio>
#include<cmath>

// common C++ head files used in the program
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<list>
#include<set>
#include<algorithm>
using namespace std;

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

