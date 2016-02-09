/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2015 The State University of New York at Buffalo
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
// and "GREATER THAN" etc. (value starts from 1 to 9)
#define FAIL_COMPARE     -1
#define EQ                0
#define LT                1
#define GT                2
#define LE                3
#define GE                4

// define the VRR and HRR work type (value starts from 10 20)
#define OS                10
#define HRR               20

// define the possible code section names (from 30-99)
// DERIV                  :  derivatives code section
// DERIV_FUNC_STATEMENT   :  the functions declare for DERIV part
// NON_RR                 :  non-RR code section, like three body kinetic integral
// NON_RR_FUNC_STATEMENT  :  the functions declare for NON_RR part
// HRR2                   :  second part of HRR code section
// HRR1                   :  first  part of HRR code section
// HRR1_FUNC_STATEMENT    :  the HRR functions declare for HRR1 part
// HRR2_FUNC_STATEMENT    :  the HRR functions declare for HRR2 part
// VRR                    :  VRR code section
// VRR_FUNC_STATEMENT     :  the functions declare for VRR part
// VRR_HEAD               :  VRR head code section
// VRR_CONT               :  VRR contraction code section
// VRR_CONT_STATEMENT     :  VRR contraction code function statment section
#define DERIV                  30
#define DERIV_FUNC_STATEMENT   31
#define NON_RR                 40
#define NON_RR_FUNC_STATEMENT  41
#define HRR2                   50
#define HRR2_FUNC_STATEMENT    51
#define HRR1                   60
#define HRR1_FUNC_STATEMENT    61
#define VRR                    70
#define VRR_HEAD               71
#define VRR_CONT               72
#define VRR_FUNC_STATEMENT     73
#define VRR_CONT_STATEMENT     74

// define the operator (100-999)
// MOM is the moment integrals
// esp is the electrostatic potential
// this group is within 100-900
// here we note for expr12 operator:
// it's exp(-omega*r12^2)
#define TWOBODYOVERLAP    100
#define THREEBODYOVERLAP  200
#define KINETIC           300
#define NAI               400
#define ERI               500
#define TWOBODYERI        510
#define THREEBODYERI      520
#define EXPR12            530
#define MOM               600
#define THREEBODYKI       700
#define ESP               800

// define the BRA and KET

// define the bra1, bra2 etc.
// they are the symbol for the position (value between 1000-2000)
#define BRA1              1000
#define BRA2              1001
#define KET1              1002
#define KET2              1003
#define BRA               1010
#define KET               1100

//
// let's add in derivatives information
// basically NO_DERIV is same with NULL_POS (in range of 2000-3000)
//
#define NO_DERIV          -1
#define DERIV_X           2001
#define DERIV_Y           2002
#define DERIV_Z           2003

// define the NULL position for shell and basis sets
// it's also used as null value
#define NULL_POS          -1

// define the RR property 
// for given operator type (in range of 3000-4000)
#define SQ_ON_L           3001
#define SQ_ON_L_M         3002
#define SQ_ON_L_OPER      3003

//
// the following macros are used to define whether an arbitrary
// shell quartet is expressed as variable form or array form;
// or they are just the output of the whole cpp file.
//
// VARIABLE_SQ       :  shell quartets expressed as variable form
// ARRAY_SQ          :  shell quartets expressed as array form
// GLOBAL_RESULT_SQ  :  the output of final results, in name of "abcd"
//
// here it ranges from 4000 to 5000
//
#define VARIABLE_SQ            4001
#define ARRAY_SQ               4002
#define GLOBAL_RESULT_SQ       4003

//
// this is an constant used to specify the 
// maximum number of exponetial factors
// adding to the integrals/shell quartet 
//
#define MAX_EXP_FAC_LIST  4

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

/**
 * here we add some inline function to determine the VRR job
 */
inline bool isValidVRRJob(int vrrJob) {
	if (vrrJob == OS || vrrJob == VRR) return true;
	return false;
};

/**
 * here we add some inline function to determine the HRR job
 */
inline bool isValidHRRJob(int hrrJob) {
	if (hrrJob == HRR1 || hrrJob == HRR2 || hrrJob == HRR) return true;
	return false;
};

/**
 * whether this is global result
 */
inline bool isGlobalResult(int status) {
	if (status == GLOBAL_RESULT_SQ) return true;
	return false;
};

/**
 * whether this is shell quartet in array form
 */
inline bool inArrayStatus(int status) {
	if (status == ARRAY_SQ) return true;
	return false;
};

#endif

