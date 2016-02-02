//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2015 The State University of New York at Buffalo
// This softare uses the MIT license as below:
//
//	Permission is hereby granted, free of charge, to any person obtaining 
//	a copy of this software and associated documentation files (the "Software"), 
//	to deal in the Software without restriction, including without limitation 
//	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//	and/or sell copies of the Software, and to permit persons to whom the Software 
//	is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
//						    
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
//	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
//	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
//	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
#include "basis.h"
#include "shell.h"
#include "integral.h"
#include "shellquartet.h"
#include "rrbuild.h"
#include "inttype.h"
#include "derivinfor.h"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>
using namespace basis;
using namespace shell;
using namespace inttype;
using namespace derivinfor;
using namespace integral;
using namespace shellquartet;
using namespace rrbuild;
using namespace boost;

//  here below is a list of variables used in the program
//
//	 alpha is the bra1's exponent
//	 beta  is the bra2's exponent
//	 gamma is the ket1's exponent
//	 dleta is the ket2's exponent
//	 A is the nuclear center for bra1
//	 B is the nuclear center for bra2
//	 C is the nuclear center for ket1
//	 D is the nuclear center for ket2
//	 P is the new center after bra1 combined with bra2
//	 Q is the new center after ket1 combined with ket2
//	 W is the new center after P combined with Q
//
//	 variables:
//	 zeta      = alpha + beta
//	 eta       = gamma + delta
//  oned2z    = 1/(2*zeta)
//  oned2e    = 1/(2*eta)
//  onedz     = 1/zeta
//  onede     = 1/eta
//  kappa     = zeta + eta
//  onedk     = 1/kappa
//  oned2zeta = 1/(2*(alpha+beta+gamma))
//  xi        = alpha*beta*onedz
//  twoxi     = 2*alpha*beta*onedz
//  rho       = zeta*eta*onedk
//  rhod2zsq  = rho/(2*zeta*zeta)
//  rhod2esq  = rho/(2*eta*eta)
//  adz       = alpha*onedz
//  bdz       = beta*onedz
//  odorho    = omega/(rho+omega)
//  rhodorho  = rho/(rho+omega)
//  orhod2z2  = (rho/(2*zeta*zeta))*(omega/(rho+omega))
//  orhod2e2  = rho/(2*eta*eta)*(omega/(rho+omega))
//  od2k      = (1/(2*kappa))*(omega/(rho+omega))
//

/////////////////////////////////////////
//  ####     building RR               //
/////////////////////////////////////////
RRBuild::RRBuild(const int& rrType0, const int& oper0, 
		const int& pos0):rrType(rrType0),oper(oper0),position(pos0),
	firstDerivPos(NULL_POS),secondDerivPos(NULL_POS),
	firstDerivDir(NO_DERIV),secondDerivDir(NO_DERIV),
	derivOrder(0)
{
	// check the input information
	if(position != BRA1 && position != BRA2 && position != KET1 && position != KET2) {
		crash(true,"position information is invalid in RRBuild contructor for general RR");
	}
	if(! isValidVRRJob(rrType) && rrType != HRR) {
		crash(true,"RR type information is invalid in RRBuild contructor for general RR");
	}

	// get the result data length
	int nitems = 0;
	if (rrType == HRR) {
		nitems  = getHRRLenItems();
	}else{
		nitems  = getVRRLenItems();
	}

	// reserve space 
	coes.reserve(nitems);
	lenVars.reserve(nitems+1);
	lenVars.push_back(0); // for lenVars, the first value is always 0

	// does the integral type changes?
	if (intOperChanged(oper)){
		intTypes.reserve(nitems);
	}
	
	// does M value are defined?
	// only for VRR process
	if (isVRRWork() && hasMValueDefined(oper)) {
		mVals.reserve(nitems);
	}

	// for angs and vars, it actually depends on the variable number
	// the maximum variable number is just 2*nitems
	angs.reserve(2*nitems);
	vars.reserve(2*nitems);

	// fill in the data
	if (rrType == HRR) {
		buildGeneralHRR();
	}else{
		buildGeneralVRR();
	}
}

int RRBuild::getVRRLenItems() const
{
	int num = -1;
	if(rrType == OS) {
		if (oper == TWOBODYOVERLAP){
			num = 3;
		}else if (oper == THREEBODYOVERLAP){
			num = 4;
		}else if (oper == MOM){
			num = 4;
		}else if( oper == KINETIC){
			num = 5;
		}else if( oper == NAI){
			num = 6;
		}else if( oper == ESP){
			num = 6;
		}else if( oper == ERI){
			num = 8;
		}else if( oper == EXPR12) {
			num = 8;
		}else {
			cout << "operator is " << oper << endl;
			crash(true,"In getVRRLenItems operator is not supported");
		}
	}else{
		cout << "RR Type is given as " << rrType << endl;
		crash(true,"In getVRRLenItems rrType is not supported");
	}
	return num;
}

void RRBuild::buildGeneralVRR() 
{
	if(rrType == OS) {

		if(oper == TWOBODYOVERLAP){

			if(position == BRA1){

				// form the coefficients
				string k1 = "PA" + uni;
				string k2 = unangA + "*oned2z";  
				string k3 = unangB + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);

				// variable position
				lenVars.push_back(1);  // bra1_-1
				lenVars.push_back(2);  // bra1_-2
				lenVars.push_back(4);  // bra1_-1_bra2_-1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "PB" + uni;
				string k2 = unangA + "*oned2z";
				string k3 = unangB + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);

				// variable position
				lenVars.push_back(1);  // bra2_-1
				lenVars.push_back(3);  // bra1_-1_bra2_-1 
				lenVars.push_back(4);  // bra2_-2

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);

			}else {
				crash(true,"In two body overlap position is wrong");
			}

		}else if(oper == MOM) {

			if(position == BRA1){

				// form the coefficients
				string k1 = "PA" + uni;
				string k2 = unangA + "*oned2z";  
				string k3 = unangB + "*oned2z";
				string k4 = unangC + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);

				// variable position
				lenVars.push_back(1);  // bra1_-1
				lenVars.push_back(2);  // bra1_-2
				lenVars.push_back(4);  // bra1_-1_bra2_-1
				lenVars.push_back(6);  // bra1_-1_ket1_-1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "PB" + uni;
				string k2 = unangA + "*oned2z";
				string k3 = unangB + "*oned2z";
				string k4 = unangC + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);

				// variable position
				lenVars.push_back(1);  // bra2_-1
				lenVars.push_back(3);  // bra1_-1_bra2_-1 
				lenVars.push_back(4);  // bra2_-2
				lenVars.push_back(6);  // bra2_-1_ket1_-1

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);

			}else {
				crash(true,"In mometum integral position is wrong");
			}

		}else if (oper == THREEBODYOVERLAP){

			if (position == BRA1){

				// form the keys
				string k1 = "GA" + uni;
				string k2 = unangA + "*oned2zeta";
				string k3 = unangB + "*oned2zeta";
				string k4 = unangC + "*oned2zeta";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);

				// variable position
				lenVars.push_back(1); // bra1_-1
				lenVars.push_back(2); // bra1_-2
				lenVars.push_back(4); // bra1_-1_bra2_-1
				lenVars.push_back(6); // bra1_-1_ket1_-1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "GB" + uni;
				string k2 = unangA + "*oned2zeta";
				string k3 = unangB + "*oned2zeta";
				string k4 = unangC + "*oned2zeta";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);

				// variable position
				lenVars.push_back(1); // bra2_-1
				lenVars.push_back(3); // bra1_-1_bra2_-1
				lenVars.push_back(4); // bra2_-2
				lenVars.push_back(6); // bra2_-1_ket1_-1

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == KET1){

				// form the keys
				string k1 = "GC" + uni;
				string k2 = unangA + "*oned2zeta";
				string k3 = unangB + "*oned2zeta";
				string k4 = unangC + "*oned2zeta";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);

				// variable position
				lenVars.push_back(1); // ket1_-1
				lenVars.push_back(3); // bra1_-1_ket1_-1
				lenVars.push_back(5); // bra2_-1_ket1_-1
				lenVars.push_back(6); // ket1_-2

				// form the variables
				vars.push_back(KET1);
				vars.push_back(BRA1);
				vars.push_back(KET1);
				vars.push_back(BRA2);
				vars.push_back(KET1);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);

			}else {
				crash(true,"In three body overlap position is wrong");
			}

		}else if( oper == KINETIC){

			if (position == BRA1){

				// form the keys
				string k1 = "PA" + uni;
				string k2 = unangA + "*oned2z";
				string k3 = unangB + "*oned2z";
				string k4 = "twoxi";
				string k5 = "-" + unangA + "*bdz";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);

				// variable position
				lenVars.push_back(1); // bra1_-1
				lenVars.push_back(2); // bra1_-2
				lenVars.push_back(4); // bra1_-1_bra2_-1
				lenVars.push_back(5); // bra1_0_(TWOBODYOVERLAP)
				lenVars.push_back(6); // bra1_-2_(TWOBODYOVERLAP)

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(0);
				angs.push_back(-2);

				// change of the integral types
				intTypes.push_back(oper);
				intTypes.push_back(oper);
				intTypes.push_back(oper);
				intTypes.push_back(TWOBODYOVERLAP);
				intTypes.push_back(TWOBODYOVERLAP);

			}else if( position == BRA2) {

				// form the keys
				string k1 = "PB" + uni;
				string k2 = unangA + "*oned2z";
				string k3 = unangB + "*oned2z";
				string k4 = "twoxi";
				string k5 = "-" + unangB + "*adz";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);

				// variable position
				lenVars.push_back(1); // bra2_-1
				lenVars.push_back(3); // bra1_-1_bra2_-1
				lenVars.push_back(4); // bra2_-2
				lenVars.push_back(5); // bra2_0_(TWOBODYOVERLAP)
				lenVars.push_back(6); // bra2_-2_(TWOBODYOVERLAP)

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(0);
				angs.push_back(-2);

				// change of the integral types
				intTypes.push_back(oper);
				intTypes.push_back(oper);
				intTypes.push_back(oper);
				intTypes.push_back(TWOBODYOVERLAP);
				intTypes.push_back(TWOBODYOVERLAP);

			}else {
				crash(true,"In KINETIC position is wrong");
			}

		}else if( oper == NAI){

			if (position == BRA1){

				// form the keys
				string k1 = "PA" + uni;
				string k2 = "-PN" + uni;
				string k3 = unangA + "*oned2z";
				string k4 = "-" + unangA + "*oned2z";
				string k5 = unangB + "*oned2z";
				string k6 = "-" + unangB + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);

				// variable position
				lenVars.push_back(1); // bra1_-1_m+0
				lenVars.push_back(2); // bra1_-1_m+1
				lenVars.push_back(3); // bra1_-2_m+0
				lenVars.push_back(4); // bra1_-2_m+1
				lenVars.push_back(6); // bra1_-1_bra2_-1_m+0
				lenVars.push_back(8); // bra1_-1_bra2_-1_m+1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "PB" + uni;
				string k2 = "-PN" + uni;
				string k3 = unangB + "*oned2z";
				string k4 = "-" + unangB + "*oned2z";
				string k5 = unangA + "*oned2z";
				string k6 = "-" + unangA + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);

				// variable position
				lenVars.push_back(1); // bra2_-1_m+0
				lenVars.push_back(2); // bra2_-1_m+1
				lenVars.push_back(3); // bra2_-2_m+0
				lenVars.push_back(4); // bra2_-2_m+1
				lenVars.push_back(6); // bra1_-1_bra2_-1_m+0
				lenVars.push_back(8); // bra1_-1_bra2_-1_m+1

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);

			}else {
				crash(true,"In NAI position is wrong");
			}

		}else if( oper == ESP){

			if (position == BRA1){

				// form the keys
				string k1 = "PA" + uni;
				string k2 = "-PR" + uni;
				string k3 = unangA + "*oned2z";
				string k4 = "-" + unangA + "*oned2z";
				string k5 = unangB + "*oned2z";
				string k6 = "-" + unangB + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);

				// variable position
				lenVars.push_back(1); // bra1_-1_m+0
				lenVars.push_back(2); // bra1_-1_m+1
				lenVars.push_back(3); // bra1_-2_m+0
				lenVars.push_back(4); // bra1_-2_m+1
				lenVars.push_back(6); // bra1_-1_bra2_-1_m+0
				lenVars.push_back(8); // bra1_-1_bra2_-1_m+1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "PB" + uni;
				string k2 = "-PR" + uni;
				string k3 = unangB + "*oned2z";
				string k4 = "-" + unangB + "*oned2z";
				string k5 = unangA + "*oned2z";
				string k6 = "-" + unangA + "*oned2z";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);

				// variable position
				lenVars.push_back(1); // bra2_-1_m+0
				lenVars.push_back(2); // bra2_-1_m+1
				lenVars.push_back(3); // bra2_-2_m+0
				lenVars.push_back(4); // bra2_-2_m+1
				lenVars.push_back(6); // bra1_-1_bra2_-1_m+0
				lenVars.push_back(8); // bra1_-1_bra2_-1_m+1

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);

			}else {
				crash(true,"In ESP position is wrong");
			}

		}else if(oper == ERI){

			if (position == BRA1){

				// form the keys
				string k1 = "PA" + uni;
				string k2 = "WP" + uni;
				string k3 = unangA + "*oned2z";
				string k4 = "-" + unangA + "*rhod2zsq";
				string k5 = unangB + "*oned2z";
				string k6 = "-" + unangB + "*rhod2zsq";
				string k7 = unangC + "*oned2k";
				string k8 = unangD + "*oned2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // bra1_-1_m+0
				lenVars.push_back(2);  // bra1_-1_m+1
				lenVars.push_back(3);  // bra1_-2_m+0
				lenVars.push_back(4);  // bra1_-2_m+1
				lenVars.push_back(6);  // bra1_-1_bra2_-1_m+0
				lenVars.push_back(8);  // bra1_-1_bra2_-1_m+1
				lenVars.push_back(10); // bra1_-1_ket1_-1_m+1
				lenVars.push_back(12); // bra1_-1_ket2_-1_m+1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(KET1);
				vars.push_back(BRA1);
				vars.push_back(KET2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(1);
				mVals.push_back(1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "PB" + uni;
				string k2 = "WP" + uni;
				string k3 = unangA + "*oned2z";
				string k4 = "-" + unangA + "*rhod2zsq";
				string k5 = unangB + "*oned2z";
				string k6 = "-" + unangB + "*rhod2zsq";
				string k7 = unangC + "*oned2k";
				string k8 = unangD + "*oned2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // bra2_-1_m+0
				lenVars.push_back(2);  // bra2_-1_m+1
				lenVars.push_back(4);  // bra1_-1_bra2_-1_m+0
				lenVars.push_back(6);  // bra1_-1_bra2_-1_m+1
				lenVars.push_back(7);  // bra2_-2_m+0
				lenVars.push_back(8);  // bra2_-2_m+1
				lenVars.push_back(10); // bra2_-1_ket1_-1_m+1
				lenVars.push_back(12); // bra2_-1_ket2_-1_m+1

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(KET1);
				vars.push_back(BRA2);
				vars.push_back(KET2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(1);
				mVals.push_back(1);

			}else if( position == KET1){

				// form the keys
				string k1 = "QC" + uni;
				string k2 = "WQ" + uni;
				string k3 = unangC + "*oned2e";
				string k4 = "-" + unangC + "*rhod2esq";
				string k5 = unangD + "*oned2e";
				string k6 = "-" + unangD + "*rhod2esq";
				string k7 = unangA + "*oned2k";
				string k8 = unangB + "*oned2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // ket1_-1_m+0         
				lenVars.push_back(2);  // ket1_-1_m+1        
				lenVars.push_back(3);  // ket1_-2_m+0         
				lenVars.push_back(4);  // ket1_-2_m+1         
				lenVars.push_back(6);  // ket1_-1_ket2_-1_m+0 
				lenVars.push_back(8);  // ket1_-1_ket2_-1_m+1 
				lenVars.push_back(10); // ket1_-1_bra1_-1_m+1 
				lenVars.push_back(12); // ket1_-1_bra2_-1_m+1 

				// form the variables
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(BRA1);
				vars.push_back(KET1);
				vars.push_back(BRA2);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(1);
				mVals.push_back(1);

			}else if( position == KET2){

				// form the keys
				string k1 = "QD" + uni;
				string k2 = "WQ" + uni;
				string k3 = unangC + "*oned2e";
				string k4 = "-" + unangC + "*rhod2esq";
				string k5 = unangD + "*oned2e";
				string k6 = "-" + unangD + "*rhod2esq";
				string k7 = unangA + "*oned2k";
				string k8 = unangB + "*oned2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // ket2_-1_m+0        
				lenVars.push_back(2);  // ket2_-1_m+1        
				lenVars.push_back(4);  // ket1_-1_ket2_-1_m+0
				lenVars.push_back(6);  // ket1_-1_ket2_-1_m+1
				lenVars.push_back(7);  // ket2_-2_m+0        
				lenVars.push_back(8);  // ket2_-2_m+1        
				lenVars.push_back(10); // ket2_-1_bra1_-1_m+1
				lenVars.push_back(12); // ket2_-1_bra2_-1_m+1

				// form the variables
				vars.push_back(KET2);
				vars.push_back(KET2);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(KET2);
				vars.push_back(KET2);
				vars.push_back(BRA1);
				vars.push_back(KET2);
				vars.push_back(BRA2);
				vars.push_back(KET2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

				// m values change
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(0);
				mVals.push_back(1);
				mVals.push_back(1);
				mVals.push_back(1);

			}else {
				crash(true,"In ERI position is wrong");
			}

		}else if(oper == EXPR12) {

			if (position == BRA1){

				// form the keys
				string k1 = "PA" + uni;
				string k2 = "WP" + uni;
				string k3 = unangA + "*oned2z";
				string k4 = "-" + unangA + "*orhod2z2";
				string k5 = unangB + "*oned2z";
				string k6 = "-" + unangB + "*orhod2z2";
				string k7 = unangC + "*od2k";
				string k8 = unangD + "*od2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // bra1_-1
				lenVars.push_back(2);  // bra1_-1
				lenVars.push_back(3);  // bra1_-2
				lenVars.push_back(4);  // bra1_-2
				lenVars.push_back(6);  // bra1_-1_bra2_-1
				lenVars.push_back(8);  // bra1_-1_bra2_-1
				lenVars.push_back(10); // bra1_-1_ket1_-1
				lenVars.push_back(12); // bra1_-1_ket2_-1

				// form the variables
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(KET1);
				vars.push_back(BRA1);
				vars.push_back(KET2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == BRA2){

				// form the keys
				string k1 = "PB" + uni;
				string k2 = "WP" + uni;
				string k3 = unangA + "*oned2z";
				string k4 = "-" + unangA + "*orhod2z2";
				string k5 = unangB + "*oned2z";
				string k6 = "-" + unangB + "*orhod2z2";
				string k7 = unangC + "*od2k";
				string k8 = unangD + "*od2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // bra2_-1
				lenVars.push_back(2);  // bra2_-1
				lenVars.push_back(4);  // bra1_-1_bra2_-1
				lenVars.push_back(6);  // bra1_-1_bra2_-1
				lenVars.push_back(7);  // bra2_-2
				lenVars.push_back(8);  // bra2_-2
				lenVars.push_back(10); // bra2_-1_ket1_-1
				lenVars.push_back(12); // bra2_-1_ket2_-1

				// form the variables
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA1);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(BRA2);
				vars.push_back(KET1);
				vars.push_back(BRA2);
				vars.push_back(KET2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == KET1){

				// form the keys
				string k1 = "QC" + uni;
				string k2 = "WQ" + uni;
				string k3 = unangC + "*oned2e";
				string k4 = "-" + unangC + "*orhod2e2";
				string k5 = unangD + "*oned2e";
				string k6 = "-" + unangD + "*orhod2e2";
				string k7 = unangA + "*od2k";
				string k8 = unangB + "*od2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // ket1_-1
				lenVars.push_back(2);  // ket1_-1
				lenVars.push_back(3);  // ket1_-2
				lenVars.push_back(4);  // ket1_-2
				lenVars.push_back(6);  // ket1_-1_ket2_-1
				lenVars.push_back(8);  // ket1_-1_ket2_-1
				lenVars.push_back(10); // ket1_-1_bra1_-1
				lenVars.push_back(12); // ket1_-1_bra2_-1

				// form the variables
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(BRA1);
				vars.push_back(KET1);
				vars.push_back(BRA2);
				vars.push_back(KET1);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

			}else if( position == KET2){

				// form the keys
				string k1 = "QD" + uni;
				string k2 = "WQ" + uni;
				string k3 = unangC + "*oned2e";
				string k4 = "-" + unangC + "*orhod2e2";
				string k5 = unangD + "*oned2e";
				string k6 = "-" + unangD + "*orhod2e2";
				string k7 = unangA + "*od2k";
				string k8 = unangB + "*od2k";
				coes.push_back(k1);
				coes.push_back(k2);
				coes.push_back(k3);
				coes.push_back(k4);
				coes.push_back(k5);
				coes.push_back(k6);
				coes.push_back(k7);
				coes.push_back(k8);

				// variable position
				lenVars.push_back(1);  // ket2_-1
				lenVars.push_back(2);  // ket2_-1
				lenVars.push_back(4);  // ket1_-1_ket2_-1
				lenVars.push_back(6);  // ket1_-1_ket2_-1
				lenVars.push_back(7);  // ket2_-2
				lenVars.push_back(8);  // ket2_-2
				lenVars.push_back(10); // ket2_-1_bra1_-1
				lenVars.push_back(12); // ket2_-1_bra2_-1

				// form the variables
				vars.push_back(KET2);
				vars.push_back(KET2);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(KET1);
				vars.push_back(KET2);
				vars.push_back(KET2);
				vars.push_back(KET2);
				vars.push_back(BRA1);
				vars.push_back(KET2);
				vars.push_back(BRA2);
				vars.push_back(KET2);

				// angular momentum changes
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-2);
				angs.push_back(-2);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);
				angs.push_back(-1);

			}else {
				crash(true,"In EXPR12 position is wrong");
			}

		}else {
			crash(true,"In buildVRR operator is not supported");
		}

	}else{
		crash(true,"In buildVRR only OS now is supported");
	}
}

/*
 * we do not use this function anymore
 * the bottom integral is built in the sqints.cpp part
void RRBuild::buildSIntegral() 
{
	string expression;
	if (oper == TWOBODYOVERLAP) {
		expression = "ic*pow(pidz,1.5E0)";
	}else if (oper == THREEBODYOVERLAP) {
		expression = "ISSS";
	}else if (oper == KINETIC) {
		expression = "ic*pow(pidz,1.5E0)*xi*(3.0E0-twoxi*AB2)";
	}else if (oper == NAI || oper == ERI){
		expression = "prefactor";
		string fm = "fm(1.0E0,u," + unM + ")";
		expression = expression + "*" + fm;
	}else{
		crash(true,"Un-supported integral type in the S type of integral calculation");
	}
	zeroExpression = expression;
}
*/

void RRBuild::buildGeneralHRR() 
{

	// determine whether we do it to bra or ket
	int job = -1;
	if (position == BRA1 || position == BRA2) {
		job = 0;
	}else if (position == KET1 || position == KET2){
		job = 1;
	}else{
		cout << "Position is: " << position << endl;
		crash(true,"Un-supported position in the buildGeneralHRR");
	}

	// expansion is over the bra side
	if (job == 0) {
		if (position == BRA1){

			// form the coefficients
			string k1 = "1"; 
			string k2 = "-AB" + uni;
			coes.push_back(k1);
			coes.push_back(k2);

			// variable position
			lenVars.push_back(2);  // bra1_-1_bra2_+1
			lenVars.push_back(3);  // bra1_-1

			// form the variables
			vars.push_back(BRA1);
			vars.push_back(BRA2);
			vars.push_back(BRA1);

			// angular momentum changes
			angs.push_back(-1);
			angs.push_back(+1);
			angs.push_back(-1);

		}else{

			// form the coefficients
			string k1 = "1";
			string k2 = "AB" + uni; 
			coes.push_back(k1);
			coes.push_back(k2);

			// variable position
			lenVars.push_back(2);  // bra1_+1_bra2_-1
			lenVars.push_back(3);  // bra2_-1

			// form the variables
			vars.push_back(BRA1);
			vars.push_back(BRA2);
			vars.push_back(BRA2);

			// angular momentum changes
			angs.push_back(+1);
			angs.push_back(-1);
			angs.push_back(-1);

		}

	}else{

		if (position == KET1){

			// form the keys
			string k1 = "1"; 
			string k2 = "-CD" + uni;
			coes.push_back(k1);
			coes.push_back(k2);

			// variable position
			lenVars.push_back(2);  // ket1_-1_ket2_+1
			lenVars.push_back(3);  // ket1_-1

			// form the variables
			vars.push_back(KET1);
			vars.push_back(KET2);
			vars.push_back(KET1);

			// angular momentum changes
			angs.push_back(-1);
			angs.push_back(+1);
			angs.push_back(-1);

		}else{

			// form the keys
			string k1 = "1"; 
			string k2 = "CD" + uni; 
			coes.push_back(k1);
			coes.push_back(k2);

			// variable position
			lenVars.push_back(2);  // ket1_+1_ket2_-1
			lenVars.push_back(3);  // ket2_-1

			// form the variables
			vars.push_back(KET1);
			vars.push_back(KET2);
			vars.push_back(KET2);

			// angular momentum changes
			angs.push_back(+1);
			angs.push_back(-1);
			angs.push_back(-1);
		}
	}
}

/////////////////////////////////////////
//  ####     building derivatives      //
/////////////////////////////////////////
RRBuild::RRBuild(const int& oper0, const int& pos1, const int& pos2,
					const int& dir1, const int& dir2):rrType(NULL_POS),
	oper(oper0),position(NULL_POS),firstDerivPos(pos1),
	secondDerivPos(pos2),firstDerivDir(dir1),secondDerivDir(dir2),
	derivOrder(0)
{
	// check the input information
	// direction information could be NULL
	if(firstDerivPos != BRA1 && firstDerivPos != BRA2 && 
			firstDerivPos != KET1 && firstDerivPos != KET2) {
		crash(true,"first deriv position is invalid in RRBuild contructor");
	}
	if(secondDerivPos != BRA1 && secondDerivPos != BRA2 && 
			secondDerivPos != KET1 && secondDerivPos != KET2 && secondDerivPos != NULL_POS) {
		crash(true,"second deriv position is invalid in RRBuild contructor");
	}
	if(firstDerivDir != DERIV_X && firstDerivDir != DERIV_Y 
			&& firstDerivDir != DERIV_Z && firstDerivDir != NO_DERIV) {
		crash(true,"first deriv is invalid in RRBuild contructor");
	}
	if(secondDerivDir != DERIV_X && secondDerivDir != DERIV_Y 
			&& secondDerivDir != DERIV_Z && secondDerivDir != NO_DERIV) {
		crash(true,"second deriv is invalid in RRBuild contructor");
	}

	// now set the deriv order
	if (firstDerivPos>0) derivOrder++; 
	if (secondDerivPos>0) derivOrder++;; 

	// figure out the length of derivatives
	int nitems = getDerivLenItems();

	// let's initilize the arrays
	coes.reserve(nitems);
	lenVars.reserve(nitems+1);
	lenVars.push_back(0); // for lenVars, the first value is always 0

	// initialize angs and vars change for integral
	if (derivOrder == 1) {
		angs.reserve(nitems);
		vars.reserve(nitems);
	}else if (derivOrder == 2) {
		if (derivOnSamePosAndDir()) {
			angs.reserve(nitems);
			vars.reserve(nitems);
		}else{
			angs.reserve(2*nitems);
			vars.reserve(2*nitems);
		}
	}

	// initilize the exps
	exps.reserve(derivOrder*nitems);

	// finally, let's build the expression
	buildDerivExpression(); 
}

int RRBuild::getDerivLenItems() const
{
	int num = -1;
	if (derivOrder == 1) {
		num = 2;
	}else if (derivOrder == 2) {
		num = 4;
	}else{
		crash(true,"currently deriv order are only 1 or 2");
	}
	return num;
}

void RRBuild::buildDerivExpression() 
{
	if (derivOrder == 1) {

		// form the keys
		// assume that deriv pos is on bra1
		// this is what the comment means
		string k1 = "2.0E0";          // 2*alpha*I(l+1)
		string k2 = "-"+unang1;       // -l*I(l-1) 
		coes.push_back(k1);
		coes.push_back(k2);

		// variable position
		lenVars.push_back(1);         // I(l+1)
		lenVars.push_back(2);         // I(l-1)

		// form the variables
		vars.push_back(firstDerivPos);// I(l+1)
		vars.push_back(firstDerivPos);// I(l-1)

		// angular momentum changes
		angs.push_back(+1);
		angs.push_back(-1);
		
		// exponetial factors multiplied with integrals
		exps.push_back(firstDerivPos);
		exps.push_back(NULL_POS);    

	}else if (derivOrder == 2) {

		//
		// firstly, let's see for what kind of situations
		// we need to do in saperate way
		// else will go to general form
		//
		// note: you can not change the order of terms of 
		// k1, k2, k3 and k4. Actually because the second 
		// term may not appear (see the buildRRSQ function,
		// the exceptional case), therefore we will double
		// check the second term "-2*alpha*l*I(l+0)" in 
		// terms of the given shell quartet. Do not change 
		// the order of term below!!!!
		//
		if (derivOnSamePosAndDir()) {

			// form the keys
			// assume that deriv pos is on bra1 and bra1
			// this is what the comment means
			// l  indicates direction that first/second 
			// deriv pos working on, assume to be x
			string k1 = "4.0E0";                // 4*alpha*alpha*I(l+2)
			string k2 = "-2.0E0*"+unang1;       // -2*alpha*l*I(l+0), this term not exist for S shell 
			string k3 = "-2.0E0*"+unangp1;      // -2*alpha*(l+1)*I(l+0) 
			string k4 = unang1+"*"+unangm1;     // l*(l-1)*I(l-2) 
			coes.push_back(k1);
			coes.push_back(k2);
			coes.push_back(k3);
			coes.push_back(k4);

			// variable position
			lenVars.push_back(1);         // I(l+2)
			lenVars.push_back(2);         // I(l+0)
			lenVars.push_back(3);         // I(l+0)
			lenVars.push_back(4);         // I(l-2)

			// form the variables
			vars.push_back(firstDerivPos);   // I(l+2)
			vars.push_back(firstDerivPos);   // I(l+0)
			vars.push_back(firstDerivPos);   // I(l+0)
			vars.push_back(firstDerivPos);   // I(l-2)

			// angular momentum changes
			angs.push_back(+2);              // I(l+2)
			angs.push_back(0);               // I(l+0)
			angs.push_back(0);               // I(l+0)
			angs.push_back(-2);              // I(l-2)

			// exponetial factors multiplied with integrals
			exps.push_back(firstDerivPos);   // first in  I(l+2)
			exps.push_back(firstDerivPos);   // second in I(l+2)
			exps.push_back(firstDerivPos);   // first in  I(l+0)
			exps.push_back(NULL_POS);        // second in I(l+0)
			exps.push_back(firstDerivPos);   // first in  I(l+0)
			exps.push_back(NULL_POS);        // second in I(l+0)
			exps.push_back(NULL_POS);        // first in  I(l-2)
			exps.push_back(NULL_POS);        // second in I(l-2)

			//
			// note:
			//
			// for each term, the first 
			// changing variable is always associated with 
			// firstDerivPos/firstDerivDir, and the second variable
			// is always binding to secondDerivPos/secondDerivDir.
			// so please keep the first->second order in the expression
			// below
			//
			// We are using this feature in buildRRInt function. 
			// so be careful about it!!
			//
		}else{

			// form the keys
			// assume that deriv pos is on bra1 and bra2
			// this is what the comment means
			// l  indicates direction that first deriv pos working on, assume to be x
			// l' indicates direction that second deriv pos working on, assume to be x
			// this also applies to the situation that two deriv pos are same
			string k1 = "4.0E0";                // 4*alpha*beta*I(l+1,l'+1)
			string k2 = "-2.0E0*"+unang2;       // -2*alpha*l'*I(l+1,l'-1) 
			string k3 = "-2.0E0*"+unang1;       // -2*beta*l*I(l-1,l'+1) 
			string k4 = unang1+"*"+unang2;      // l*l'*I(l-1,l'-1) 
			coes.push_back(k1);
			coes.push_back(k2);
			coes.push_back(k3);
			coes.push_back(k4);

			// variable position
			lenVars.push_back(2);         // I(l+1,l'+1)
			lenVars.push_back(4);         // I(l+1,l'-1)
			lenVars.push_back(6);         // I(l-1,l'+1)
			lenVars.push_back(8);         // I(l-1,l'-1)

			// form the variables
			vars.push_back(firstDerivPos);   // I(l+1,l'+1)
			vars.push_back(secondDerivPos);  // I(l+1,l'+1) 
			vars.push_back(firstDerivPos);   // I(l+1,l'-1)
			vars.push_back(secondDerivPos);  // I(l+1,l'-1) 
			vars.push_back(firstDerivPos);   // I(l-1,l'+1)
			vars.push_back(secondDerivPos);  // I(l-1,l'+1) 
			vars.push_back(firstDerivPos);   // I(l-1,l'-1)
			vars.push_back(secondDerivPos);  // I(l-1,l'-1) 

			// angular momentum changes
			angs.push_back(+1);              // first in  I(l+1,l'+1)
			angs.push_back(+1);              // second in I(l+1,l'+1)
			angs.push_back(+1);              // first in  I(l+1,l'-1)
			angs.push_back(-1);              // second in I(l+1,l'-1)
			angs.push_back(-1);              // first in  I(l-1,l'+1)
			angs.push_back(+1);              // second in I(l-1,l'+1)
			angs.push_back(-1);              // first in  I(l-1,l'-1)
			angs.push_back(-1);              // second in I(l-1,l'-1)

			// exponetial factors multiplied with integrals
			exps.push_back(firstDerivPos);   // first in  I(l+1,l'+1)
			exps.push_back(secondDerivPos);  // second in I(l+1,l'+1)
			exps.push_back(firstDerivPos);   // first in  I(l+1,l'-1)
			exps.push_back(NULL_POS);        // second in I(l+1,l'-1)
			exps.push_back(NULL_POS);        // first in  I(l-1,l'+1)
			exps.push_back(secondDerivPos);  // second in I(l-1,l'+1)
			exps.push_back(NULL_POS);        // first in  I(l-1,l'-1)
			exps.push_back(NULL_POS);        // second in I(l-1,l'-1)
		}

	}else {
		crash(true,"deriv order now is not supported in buildDerivExpression");
	}
}

/////////////////////////////////////////////////////////
//  ####     building non-RR non-Deriv expression      //
/////////////////////////////////////////////////////////
RRBuild::RRBuild(const int& oper0):rrType(NULL_POS),
	oper(oper0),position(NULL_POS),firstDerivPos(NULL_POS),
	secondDerivPos(NULL_POS),firstDerivDir(NO_DERIV),
	secondDerivDir(NO_DERIV),derivOrder(0)
{
	// check the input information
	if(oper != THREEBODYKI) {
		crash(true,"oper is invalid in RRBuild contructor for non-RR non-Deriv work");
	}

	// finally, let's build the expression
	if(oper == THREEBODYKI) {
		buildThreeBodyKIExpression(); 
	}
}

void RRBuild::buildThreeBodyKIExpression() 
{
	// figure out the length terms
	int nitems = getNItems();

	// let's initilize the arrays
	coes.reserve(nitems);
	lenVars.reserve(nitems+1);
	lenVars.push_back(0); // for lenVars, the first value is always 0

	// initialize angs and vars etc.
	angs.reserve(2*nitems);
	vars.reserve(2*nitems);
	exps.reserve(2*nitems);
	intTypes.reserve(nitems);

	//
	// the expression below is for a single direction, for 
	// example; like x. The three body kinetic integrals
	// are sum over x, y and z. therefore, be careful
	// in building the integral exression you need to pass
	// all of x, y and z to get the full expression
	//
	// we note, that all of coefficients has been scaled with
	// 1/2, which appears in the kinetic operator
	//
	// NiA, NiB etc. are just the angular momentum value for
	// the LHS integrals, this is different from RR expression
	//

	// form the keys
	string k1 = "0.5E0*" + unangA + "*" + unangC;  // 0.5*NiA*NiC*I(bra1_-1_ket1_-1)
	string k2 = "-" + unangC;                      // -NiC*alpha*I(bra1_+1_ket1_-1)
	string k3 = "-" + unangA;                      // -NiA*gamma*I(bra1_-1_ket1_+1)
	string k4 = "2.0E0";                           // 2*alpha*gamma*I(bra1_+1_ket1_+1)
	string k5 = "0.5E0*" + unangB + "*" + unangC;  // 0.5*NiB*NiC*I(bra2_-1_ket1_-1)
	string k6 = "-" + unangC;                      // -NiC*beta*I(bra2_+1_ket1_-1)
	string k7 = "-" + unangB;                      // -NiB*gamma*I(bra2_-1_ket1_+1)
	string k8 = "2.0E0";                           // 2*beta*gamma*I(bra2_+1_ket1_+1)
	coes.push_back(k1);
	coes.push_back(k2);
	coes.push_back(k3);
	coes.push_back(k4);
	coes.push_back(k5);
	coes.push_back(k6);
	coes.push_back(k7);
	coes.push_back(k8);

	// variable position
	lenVars.push_back(2);  // bra1_-1_ket1_-1
	lenVars.push_back(4);  // bra1_+1_ket1_-1
	lenVars.push_back(6);  // bra1_-1_ket1_+1
	lenVars.push_back(8);  // bra1_+1_ket1_+1
	lenVars.push_back(10); // bra2_-1_ket1_-1
	lenVars.push_back(12); // bra2_+1_ket1_-1
	lenVars.push_back(14); // bra2_-1_ket1_+1
	lenVars.push_back(16); // bra2_+1_ket1_+1

	// form the variables
	// bra1 and ket1
	vars.push_back(BRA1);
	vars.push_back(KET1);
	vars.push_back(BRA1);
	vars.push_back(KET1);
	vars.push_back(BRA1);
	vars.push_back(KET1);
	vars.push_back(BRA1);
	vars.push_back(KET1);

	// bra2 and ket1
	vars.push_back(BRA2);
	vars.push_back(KET1);
	vars.push_back(BRA2);
	vars.push_back(KET1);
	vars.push_back(BRA2);
	vars.push_back(KET1);
	vars.push_back(BRA2);
	vars.push_back(KET1);

	// angular momentum changes
	// bra1 and ket1
	angs.push_back(-1);
	angs.push_back(-1);
	angs.push_back(+1);
	angs.push_back(-1);
	angs.push_back(-1);
	angs.push_back(+1);
	angs.push_back(+1);
	angs.push_back(+1);

	// bra2 and ket1
	angs.push_back(-1);
	angs.push_back(-1);
	angs.push_back(+1);
	angs.push_back(-1);
	angs.push_back(-1);
	angs.push_back(+1);
	angs.push_back(+1);
	angs.push_back(+1);

	// adding of the exponential factors
	exps.push_back(NULL_POS);
	exps.push_back(NULL_POS);
	exps.push_back(BRA1);
	exps.push_back(NULL_POS);
	exps.push_back(KET1);
	exps.push_back(NULL_POS);
	exps.push_back(BRA1);
	exps.push_back(KET1);
	exps.push_back(NULL_POS);
	exps.push_back(NULL_POS);
	exps.push_back(BRA2);
	exps.push_back(NULL_POS);
	exps.push_back(KET1);
	exps.push_back(NULL_POS);
	exps.push_back(BRA2);
	exps.push_back(KET1);

	// change of operator
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
	intTypes.push_back(THREEBODYOVERLAP);
}

/////////////////////////////////////////
//  #### general RRBuild functions     //
/////////////////////////////////////////
int RRBuild::getNItems() const 
{
	// is it RR work?
	int n = -1;
	if (isRRWork()) {
		if (rrType == HRR) {
			n = getHRRLenItems();
		}else{
			n = getVRRLenItems();
		}
	}else if (derivOrder>0) {
		n = getDerivLenItems();
	}else{
		if (oper == THREEBODYKI) {
			n = 8;
		}else{
			crash(true, "invalid situation in getNItems() in RRBuild");
		}
	}
	return n;
}

string RRBuild::determineDirection(const Basis& bas) const
{
	// get the original l,m,n
	int l,m,n;
	bas.getlmn(l,m,n);

	// we should not meet S type of basis set here!
	int L = l + m + n;
	if (L == 0) {
		crash(true,"No S type basis set in determineDirection!!!");
	}

	// find the smallest one
	int minL = max(l,m);
	minL     = max(minL,n);
	if (l > 0 && l < minL){
		minL = l;
	}
	if (m > 0 && m < minL){
		minL = m;
	}
	if (n > 0 && n < minL){
		minL = n;
	}

	// if we encounter the situation that minL = l and minL = n etc.,
	// we return it according to Z,Y,X order, that means, we try to
	// do it first in the Z direction, then Y, then X (I should tell you that
	// such decision is arbitrary....)
	if (minL == n){
		return "Z";
	}else if (minL == m){
		return "Y";
	}else{
		return "X";
	}
}

void RRBuild::determineNi(const Integral& I, const string& i,
		string& NiA, string& NiB, string& NiC, string& NiD) const 
{
	// firstly, if it's for derivatives calculation etc.; then we just
	// need natrual Ni value for the basis set inside the I
	// therefore let's check it first
	// here different from the RR process, we use the operator 
	// order to get the number of body
	if (isDerivWork() || isNonRRNonDerivWork()) {
		int nBody = getOperOrder(oper);
		for(int iBody=0; iBody<nBody; iBody++) {
			if (iBody == 0) {
				const Basis& b = I.getBasis(BRA1);
				NiA = Ni(b,i);
			} else if (iBody == 1) {
				const Basis& b = I.getBasis(BRA2);
				NiB = Ni(b,i);
			} else if (iBody == 2) {
				const Basis& b = I.getBasis(KET1);
				NiC = Ni(b,i);
			} else if (iBody == 3) {
				const Basis& b = I.getBasis(KET2);
				NiD = Ni(b,i);
			}
		}
		return;
	}

	// now below is the derivation for RR part
	// get information
	//
	// here we need to use the getRROrder
	// because the last dimension in MOM like 
	// integrals also involve into the RR process
	//
	// also pay attention to the VRR expression(right now
	// only OS scheme is involved), the RHS term use
	// the angular momentum for basMinus1, where I
	// is the LHS integral
	int order         = getRROrder(oper);
	const Basis& bas  = I.getBasis(position);
	Basis basMinus1   = getChangeBasisSet(bas,i,-1);

	// test that bas - l_{i} should not be None
	if (basMinus1.isnull()){
		crash(true,"Wrong basMinus1 in determineNi(), should not be None");
	}

	// get the basis set
	const Basis& bra1 = I.getBasis(BRA1);
	const Basis& bra2 = I.getBasis(BRA2);
	const Basis& ket1 = I.getBasis(KET1);
	const Basis& ket2 = I.getBasis(KET2);

	// for the RR work, the Ni value related to 
	// bas - l_{i}
	NiA = "0";
	NiB = "0";
	NiC = "0";
	NiD = "0";
	if (order == 2) {
		if (position == BRA1) {
			NiA = Ni(basMinus1,i);
			NiB = Ni(bra2,i);
		}else if (position == BRA2){
			NiB = Ni(basMinus1,i);
			NiA = Ni(bra1,i);
		}else{
			crash(true, "Un-supported position in determineNi()");
		}
	}else if (order == 3) {
		if (position == BRA1){
			NiA = Ni(basMinus1,i);
			NiB = Ni(bra2,i);
			NiC = Ni(ket1,i);
		}else if (position == BRA2) {
			NiA = Ni(bra1,i);
			NiB = Ni(basMinus1,i);
			NiC = Ni(ket1,i);
		}else if (position == KET1){
			NiA = Ni(bra1,i);
			NiB = Ni(bra2,i);
			NiC = Ni(basMinus1,i);
		}else{
			crash(true,"Un-supported position in determineNi()");
		}
	}else if (order == 4){
		if (position == BRA1) {
			NiA = Ni(basMinus1,i);
			NiB = Ni(bra2,i);
			NiC = Ni(ket1,i);
			NiD = Ni(ket2,i);
		}else if (position == BRA2){
			NiA = Ni(bra1,i);
			NiB = Ni(basMinus1,i);
			NiC = Ni(ket1,i);
			NiD = Ni(ket2,i);
		}else if (position == KET1){
			NiA = Ni(bra1,i);
			NiB = Ni(bra2,i);
			NiC = Ni(basMinus1,i);
			NiD = Ni(ket2,i);
		}else if (position == KET2){
			NiA = Ni(bra1,i);
			NiB = Ni(bra2,i);
			NiC = Ni(ket1,i);
			NiD = Ni(basMinus1,i);
		}else{
			crash(true,"Un-supported position in determineNi()");
		}
	}else{
		crash(true,"Un-supported integral order in determineNi()");
	}
}

Basis RRBuild::getChangeBasisSet(const Basis& bas, const string& i, const int& inc) const
{
	// test the basis set
	if (bas.isnull()) {
		crash(true,"It's wrong that in getChangeBasisSet the input bas is None");
	}

	int l,m,n; 
	bas.getlmn(l,m,n);
	if (i == "X"){
		if (l+inc >= 0){
			Basis b(l+inc,m,n);
			return b;
		}else{
			Basis b(NULL_POS,NULL_POS,NULL_POS);
			return b;
		}
	}else if (i == "Y"){
		if (m+inc >= 0){
			Basis b(l,m+inc,n);
			return b;
		}else{
			Basis b(NULL_POS,NULL_POS,NULL_POS);
			return b;
		}
	}else if (i == "Z"){
		if (n+inc >= 0){
			Basis b(l,m,n+inc);
			return b;
		}else{
			Basis b(NULL_POS,NULL_POS,NULL_POS);
			return b;
		}
	}else{
		crash(true,"Wrong i in the getChangeBasisSet");
		Basis b(NULL_POS,NULL_POS,NULL_POS);
		return b;
	}
}

string RRBuild::Ni(const Basis& bas,const string& state) const
{

	// the basis set could be none value, however; we should solved it in the
	// determineNi function, not here!!!
	if (bas.isnull()){
		crash(true, "Basis set should not be None in Ni()");
		return "0";
	} else{
		int l,m,n; 
		bas.getlmn(l,m,n);
		string ni = "0";
		if (state == "X") {
			ni = lexical_cast<string>(l);
		}else if (state == "Y"){
			ni = lexical_cast<string>(m);
		}else if (state == "Z"){
			ni = lexical_cast<string>(n);
		}else{
			crash(true,"Wrong state value in the Ni function");
		}
		return ni;
	}
}

void RRBuild::buildRRInt(const Integral& I, vector<string>& coeArray,
		vector<int>& indexArray) const
{
	//
	// originally buildRRInt only applies to the RR situation
	// then this function is extended to the derivative situation
	// 
	// we found that we can not directly use this function to
	// derive the non-RR non-deriv situation like three body
	// kinetic integrals, because the general expression only
	// involves one direction, but the result expression needs
	// all of x, y and z
	//
	// therefore, we have an independent function called
	// build3BodyKIInt for three body integral case
	//
	if (isNonRRNonDerivWork()) {
		crash(true,"the buildRRInt function does not apply to the non-RR non-deriv situation");
	}

	// firstly, let's take care of RR work(only VRR use it)
	// we determine the uncertainNiA etc.
	string NiA = "0";
	string NiB = "0";
	string NiC = "0";
	string NiD = "0";
	string xyz = "NONE";
	if (isRRWork()) {

		// according to the position, now we have to know
		// on each direction we do expanding, x, y or z?
		const Basis& b0 = I.getBasis(position);
		xyz = determineDirection(b0);

		// this is relying on the position and integral type
		// for RR, used only for VRR
		// HRR does not need this
		if (isVRRWork()) {
			determineNi(I,xyz,NiA,NiB,NiC,NiD);
		}
	}

	// now it's derivative expression
	// right now we have maximum order = 2
	// derivatives, this is for deriv order = 1
	// which is always needed for higher orders
	string NiA1 = "0";
	string NiB1 = "0";
	string NiC1 = "0";
	string NiD1 = "0";
	string xyz1 = "NONE";
	if (isDerivWork()) {

		// get the direction
		if (firstDerivDir == DERIV_X) {
			xyz1 = "X";
		}else if (firstDerivDir == DERIV_Y) {
			xyz1 = "Y";
		}else if (firstDerivDir == DERIV_Z) {
			xyz1 = "Z";
		}else{
			crash(true, "illegal firstDerivDir value in buildRRInt function");
		}

		// now determine the value
		determineNi(I,xyz1,NiA1,NiB1,NiC1,NiD1);
	}

	// here we also need L minus one value
	// for the derivatives on same position and same
	// direction, for order = 2
	string NiA1m = "0";
	string NiB1m = "0";
	string NiC1m = "0";
	string NiD1m = "0";
	string NiA1p = "0";
	string NiB1p = "0";
	string NiC1p = "0";
	string NiD1p = "0";
	if (derivOrder >= 2 && derivOnSamePosAndDir()) {
		
		// NiA-1
		int val = lexical_cast<int>(NiA1);
		NiA1p = lexical_cast<string>(val+1);
		if (val-1>0) {
			NiA1m = lexical_cast<string>(val-1);
		}

		// NiB-1
		val = lexical_cast<int>(NiB1);
		NiB1p = lexical_cast<string>(val+1);
		if (val-1>0) {
			NiB1m = lexical_cast<string>(val-1);
		}

		// NiC-1
		val = lexical_cast<int>(NiC1);
		NiC1p = lexical_cast<string>(val+1);
		if (val-1>0) {
			NiC1m = lexical_cast<string>(val-1);
		}

		// NiD-1
		val = lexical_cast<int>(NiD1);
		NiD1p = lexical_cast<string>(val+1);
		if (val-1>0) {
			NiD1m = lexical_cast<string>(val-1);
		}
	}

	// now do the second derivatives
	string NiA2 = "0";
	string NiB2 = "0";
	string NiC2 = "0";
	string NiD2 = "0";
	string xyz2 = "NONE";
	if (isDerivWork() && derivOrder >= 2) {

		// get the direction
		if (secondDerivDir == DERIV_X) {
			xyz2 = "X";
		}else if (secondDerivDir == DERIV_Y) {
			xyz2 = "Y";
		}else if (secondDerivDir == DERIV_Z) {
			xyz2 = "Z";
		}else{
			crash(true, "illegal secondDerivDir value in buildRRInt function");
		}

		// now determine the value
		determineNi(I,xyz2,NiA2,NiB2,NiC2,NiD2);
	}

	// loop over the coefficients and expressions in the RHS
	int nItems = getNItems();
	for(int item=0; item<nItems; item++) {

		// first anaylze the coefficient of this item
		string k = coes[item];
		string::size_type pos;
		if (isRRWork()) {
			if ((pos = k.find(uni)) != string::npos) {
				k.replace(pos,uni.size(),xyz);
			} 
			if ((pos = k.find(unangA)) != string::npos) {
				k.replace(pos,unangA.size(),NiA);
			} 
			if ((pos = k.find(unangB)) != string::npos) {
				k.replace(pos,unangB.size(),NiB);
			} 
			if ((pos = k.find(unangC)) != string::npos) {
				k.replace(pos,unangC.size(),NiC);
			} 
			if ((pos = k.find(unangD)) != string::npos) {
				k.replace(pos,unangD.size(),NiD);
			}
		} else if (isDerivWork()) {

			// consider 1st derivative
			if (derivOrder>=1) {
				if ((pos = k.find(unang1)) != string::npos) {
					if (firstDerivPos == BRA1) {
						k.replace(pos,unang1.size(),NiA1);
					}else if (firstDerivPos == BRA2) {
						k.replace(pos,unang1.size(),NiB1);
					}else if (firstDerivPos == KET1) {
						k.replace(pos,unang1.size(),NiC1);
					}else if (firstDerivPos == KET2) {
						k.replace(pos,unang1.size(),NiD1);
					}
				} 
			}

			// consider 2ed derivative
			if (derivOrder >=2) {
				if(derivOnSamePosAndDir()) {
					if ((pos = k.find(unangm1)) != string::npos) {
						if (secondDerivPos == BRA1) {
							k.replace(pos,unangm1.size(),NiA1m);
						}else if (secondDerivPos == BRA2) {
							k.replace(pos,unangm1.size(),NiB1m);
						}else if (secondDerivPos == KET1) {
							k.replace(pos,unangm1.size(),NiC1m);
						}else if (secondDerivPos == KET2) {
							k.replace(pos,unangm1.size(),NiD1m);
						}
					} 
					if ((pos = k.find(unangp1)) != string::npos) {
						if (secondDerivPos == BRA1) {
							k.replace(pos,unangp1.size(),NiA1p);
						}else if (secondDerivPos == BRA2) {
							k.replace(pos,unangp1.size(),NiB1p);
						}else if (secondDerivPos == KET1) {
							k.replace(pos,unangp1.size(),NiC1p);
						}else if (secondDerivPos == KET2) {
							k.replace(pos,unangp1.size(),NiD1p);
						}
					} 
				}else{
					if ((pos = k.find(unang2)) != string::npos) {
						if (secondDerivPos == BRA1) {
							k.replace(pos,unang2.size(),NiA2);
						}else if (secondDerivPos == BRA2) {
							k.replace(pos,unang2.size(),NiB2);
						}else if (secondDerivPos == KET1) {
							k.replace(pos,unang2.size(),NiC2);
						}else if (secondDerivPos == KET2) {
							k.replace(pos,unang2.size(),NiD2);
						}
					} 
				}
			}
		}

		//
		// check wether this item has zero coefficient
		// since the coefficients is merely a series of 
		// multiplication, therefore the 0 could appear
		// in the head, in the middle or in the last position.
		// based on this point, we have rules following:
		// 1  for the coefficient, the first element is 0
		//    and second element is *
		// 2  if it's not 1, then we check the following 
		//    situations:
		//    2.1   *0* 
		//    2.2   -0*
		//    2.3   +0*
		//    2.4   *0, 0 is in the last position
		//
		if (k[0] == '0' && k[1] == '*') {
			int index = NULL_POS;
			string coeff = " ";
			coeArray.push_back(coeff);
			indexArray.push_back(index);
			continue;
		}else {

			// consider the other cases
			bool bypass = false;

			// whether 0 is in the last position?
			int len = k.length();
			if (len >= 2 && k[len-1] == '0' && k[len-2] == '*') bypass = true;

			// now consider that 0 is in the middle
			string::size_type pos1 = k.find("*0*");
			string::size_type pos2 = k.find("-0*");
			string::size_type pos3 = k.find("+0*");
			if (pos1 != string::npos || pos2 != string::npos || pos3 != string::npos ) {
				bypass = true;
			}

			// do we bypass this term?
			if (bypass) {
				int index = NULL_POS;
				string coeff = " ";
				coeArray.push_back(coeff);
				indexArray.push_back(index);
				continue;
			}
		}

		// finally, add + to these positive terms
		if (item > 0 && k[0] != '-') {
			k = "+" + k;
		}

		// take care of the operator information
		int newO = oper;             
		if (!isDerivWork() && !isHRRWork() && intOperChanged(oper)) {
			newO = intTypes[item];
		}

		// take care of the M value
		// M value is only in VRR
		int newM = I.getMValue(); 
		if (isVRRWork() && hasMValueDefined(oper)) {
			newM = newM + mVals[item];
		}

		// take care of the division
		// we will carry division information from LHS
		// to RHS for all of modules
		long long newDivision = I.getDivision();

		// now let's take care of the variables
		// we note, that for OS type of algorithm the maximum
		// change number is two and we use this information
		// here
		int vp1  = lenVars[item];
		int v1   = vars[vp1];     // the first changing variable
		int ang1 = angs[vp1];     // its changing angular momentum 
		int varNum = lenVars[item+1] - lenVars[item];          

		// the next changing variable is set on top of the previous one
		int v2   = NULL_POS;      // the second changing variable
		int ang2 = NULL_POS;      // its changing angular momentum
		if (varNum>1) {
			int vp2  = vp1+1;
			v2   = vars[vp2];
			ang2 = angs[vp2];
		}

		// now let's analyze the corresponding value
		Basis b1(I.getBasis(BRA1));
		Basis b2(I.getBasis(BRA2));
		Basis b3(I.getBasis(KET1));
		Basis b4(I.getBasis(KET2));
		bool hasNewI = true;
		for(int n=0; n<varNum; n++) {

			// get the variable and its angular momentum infor
			int newVar = v1;
			int angInc = ang1;
			if (n == 1) {
				newVar = v2;
				angInc = ang2;
			}

			// get the basis going to be changed
			Basis oldBas;
			if (newVar == BRA1) {
				oldBas = b1;
			}else if(newVar == BRA2){
				oldBas = b2;
			}else if(newVar == KET1){
				oldBas = b3;
			}else {
				oldBas = b4;
			}

			// identify the changing xyz for th basis set
			string changeXYZ = "NONE";
			if (isRRWork()) {
				changeXYZ = xyz;
			} else if (isDerivWork()) {
				if (derivOrder == 1) {
					changeXYZ = xyz1;
				}else if (derivOrder == 2) {
					if (derivOnSamePosAndDir()) {
						changeXYZ = xyz1;
					}else{
						// now it depends on which term
						// v1 is always associated with firstDerivDir
						// v2 is always associated with secondDerivDir
						changeXYZ = xyz1;
						if (n == 1) {
							changeXYZ = xyz2;
						}
					}
				}
			}

			// now get the basis set change
			Basis newBas = oldBas;
			if (angInc != 0) {
				newBas = getChangeBasisSet(oldBas,changeXYZ,angInc);
			}

			// replace the original basis set with new one
			// we note that if the newBas is none, then this integral term 
			// does not exist therefore we proceed to the next expression
			if (newBas.isnull()) {
				hasNewI = false;
				break;
			}else{
				if (newVar == BRA1) {
					b1 = newBas;
				}else if(newVar == BRA2){
					b2 = newBas;
				}else if(newVar == KET1){
					b3 = newBas;
				}else {
					b4 = newBas;
				}
			}
		}

		// now if all of newV exisits, then integral exists; this integral exists
		if (hasNewI) {
			Integral newI(b1,b2,b3,b4,newO,newM,newDivision);

			// we need to take care of exponential factor multipler
			if (!isRRWork()) {

				// for the deriv work, for each item we have 
				// number of deriv order possible multiplers
				// exps in dimension of derivOrder*nitems
				if (isDerivWork()) {
					for(int i=0; i<derivOrder; i++) {
						int pos = exps[derivOrder*item+i];
						if (pos == BRA1 || pos == BRA2 || pos == KET1 || pos == KET2) {
							newI.addExpFac(pos); 
						}
					}
				}
			}

			// for the HRR case, we need to pass the same exponential
			// information from the LHS sq to the RHS
			if (I.withExpFac() && rrType == HRR) {
				newI.copyExpFac(I);
			}

			// now push in the information
			int index = newI.getIndex();
			coeArray.push_back(k);
			indexArray.push_back(index);
		}else{
			int index = NULL_POS;
			string coeff = " ";
			coeArray.push_back(coeff);
			indexArray.push_back(index);
		}
	}
}

void RRBuild::build3BodyKIInt(const Integral& I, const int& direction, 
		vector<string>& coeArray, vector<int>& indexArray) const
{
	//
	// this function is modified from buildRRInt 
	// for the RR and deriv case you can not use it
	// but we do not check inside
	//
	// this function only applies to 3 body KI right now
	//
	
	// set the xyz according to the direction
	string xyz = "NONE";
	if (direction == DERIV_X) {
		xyz = "X";
	}else if (direction == DERIV_Y) {
		xyz = "Y";
	}else if (direction == DERIV_Z) {
		xyz = "Z";
	}else{
		crash(true,"the input direction in build3BodyKIInt is not correct");
	}

	// we determine the uncertainNiA etc.
	string NiA = "0";
	string NiB = "0";
	string NiC = "0";
	string NiD = "0";
	determineNi(I,xyz,NiA,NiB,NiC,NiD);

	// loop over the coefficients and expressions in the RHS
	int nItems = getNItems();
	for(int item=0; item<nItems; item++) {

		// first anaylze the coefficient of this item
		string k = coes[item];
		string::size_type pos;
		if ((pos = k.find(unangA)) != string::npos) {
			k.replace(pos,unangA.size(),NiA);
		} 
		if ((pos = k.find(unangB)) != string::npos) {
			k.replace(pos,unangB.size(),NiB);
		} 
		if ((pos = k.find(unangC)) != string::npos) {
			k.replace(pos,unangC.size(),NiC);
		} 
		if ((pos = k.find(unangD)) != string::npos) {
			k.replace(pos,unangD.size(),NiD);
		}

		//
		// check wether this item has zero coefficient
		// since the coefficients is merely a series of 
		// multiplication, therefore the 0 could appear
		// in the head, in the middle or in the last position.
		// based on this point, we have rules following:
		// 1  for the coefficient, the first element is 0
		//    and second element is *
		// 2  if it's not 1, then we check the following 
		//    situations:
		//    2.1   *0* 
		//    2.2   -0*
		//    2.3   +0*
		//    2.4   *0, 0 is in the last position
		//
		if (k[0] == '0' && k[1] == '*') {
			int index = NULL_POS;
			string coeff = " ";
			coeArray.push_back(coeff);
			indexArray.push_back(index);
			continue;
		}else {

			// consider the other cases
			bool bypass = false;

			// whether 0 is in the last position?
			int len = k.length();
			if (len >= 2 && k[len-1] == '0' && k[len-2] == '*') bypass = true;

			// now consider that 0 is in the middle
			string::size_type pos1 = k.find("*0*");
			string::size_type pos2 = k.find("-0*");
			string::size_type pos3 = k.find("+0*");
			if (pos1 != string::npos || pos2 != string::npos || pos3 != string::npos ) {
				bypass = true;
			}

			// do we bypass this term?
			if (bypass) {
				int index = NULL_POS;
				string coeff = " ";
				coeArray.push_back(coeff);
				indexArray.push_back(index);
				continue;
			}
		}

		// finally, add + to these positive terms
		if (item > 0 && k[0] != '-') {
			k = "+" + k;
		}

		// take care of the operator information
		int newO = oper;             
		if (!isDerivWork() && !isHRRWork() && intOperChanged(oper)) {
			newO = intTypes[item];
		}

		// take care of the M value
		// M value is only in VRR
		int newM = I.getMValue(); 
		if (isVRRWork() && hasMValueDefined(oper)) {
			newM = newM + mVals[item];
		}

		// take care of the division
		// division information will be carried from LHS
		// to RHS for all of modules
		long long newDivision = I.getDivision();

		// now let's take care of the variables
		// we note, that for OS type of algorithm the maximum
		// change number is two and we use this information
		// here
		int vp1  = lenVars[item];
		int v1   = vars[vp1];     // the first changing variable
		int ang1 = angs[vp1];     // its changing angular momentum 
		int varNum = lenVars[item+1] - lenVars[item];          

		// the next changing variable is set on top of the previous one
		int v2   = NULL_POS;      // the second changing variable
		int ang2 = NULL_POS;      // its changing angular momentum
		if (varNum>1) {
			int vp2  = vp1+1;
			v2   = vars[vp2];
			ang2 = angs[vp2];
		}

		// now let's analyze the corresponding value
		Basis b1(I.getBasis(BRA1));
		Basis b2(I.getBasis(BRA2));
		Basis b3(I.getBasis(KET1));
		Basis b4(I.getBasis(KET2));
		bool hasNewI = true;
		for(int n=0; n<varNum; n++) {

			// get the variable and its angular momentum infor
			int newVar = v1;
			int angInc = ang1;
			if (n == 1) {
				newVar = v2;
				angInc = ang2;
			}

			// now get the basis set change
			Basis newBas;
			const Basis& tmpBas = I.getBasis(newVar); 
			if (angInc == 0) {
				newBas = tmpBas;
			}else {
				newBas = getChangeBasisSet(tmpBas,xyz,angInc);
			}

			// replace the original basis set with new one
			// we note that if the newBas is none, then this integral term 
			// does not exist therefore we proceed to the next expression
			if (newBas.isnull()) {
				hasNewI = false;
				break;
			}else{
				if (newVar == BRA1) {
					b1 = newBas;
				}else if(newVar == BRA2){
					b2 = newBas;
				}else if(newVar == KET1){
					b3 = newBas;
				}else {
					b4 = newBas;
				}
			}
		}

		// now if all of newV exisits, then integral exists; this integral exists
		if (hasNewI) {
			Integral newI(b1,b2,b3,b4,newO,newM,newDivision);

			// finally, we need to consider the possible exp factor adding in
			if (isNonRRNonDerivWork()) {
				if (oper == THREEBODYKI) {
					int order = getNExpFacFor3BodyKi();
					for(int i=0; i<order; i++) {
						int pos = exps[order*item+i];
						if (pos == BRA1 || pos == BRA2 || pos == KET1) {
							newI.addExpFac(pos); 
						}
					}
				}
			}

			// now push in the information
			int index = newI.getIndex();
			coeArray.push_back(k);
			indexArray.push_back(index);
		}else{
			int index = NULL_POS;
			string coeff = " ";
			coeArray.push_back(coeff);
			indexArray.push_back(index);
		}
	}
}

void RRBuild::buildRRSQ(const ShellQuartet& sq, vector<ShellQuartet>& rrSQList) const
{
	//
	// originally buildRRSQ only applies to the RR situation
	// then this function is extended to all of situations,
	// including derivatives and non-RR, non-deriv 
	// 
	
	// for producing the RHS shell quartet we do not need 
	// the coefficients information
	int nItems = getNItems();
	for(int item=0; item<nItems; item++) {

		// take care of the operator information
		// in default, operator is not change, the only
		// exception now is the kinetic energy calculation
		int newO = oper;             
		if (!isDerivWork() && !isHRRWork() && intOperChanged(oper)) {
			newO = intTypes[item];
		}

		// take care of the M value
		// M value is only in VRR
		int newM = sq.getM(); 
		if (isVRRWork() && hasMValueDefined(oper)) {
			newM = newM + mVals[item];
		}

		// take care of the division
		// division information will be carried from LHS
		// to RHS for all of modules
		long long newDivision = sq.getDivision();

		// now let's take care of the variables
		// we note, that for OS type of algorithm the maximum
		// change number is two. therefore, we just use two 
		// variables
		int vp1  = lenVars[item];
		int v1   = vars[vp1];     // the first changing variable
		int ang1 = angs[vp1];     // its changing angular momentum 
		int varNum = lenVars[item+1] - lenVars[item];          

		// the next changing variable is set on top of the previous one
		int v2   = NULL_POS;      // the second changing variable
		int ang2 = NULL_POS;      // its changing angular momentum
		if (varNum>1) {
			int vp2  = vp1+1;
			v2   = vars[vp2];
			ang2 = angs[vp2];
		}

		// now this is an exception that needs particular check
		// that for the deriv order = 2, and same derivatives
		// positions, the second term 2*alpha*l*I(l+0), if the original
		// shell is s; this term does not appear (see the definition
		// of formular above)
		// for integral derivation it's very clear. Because the l=0
		// in the coefficients. However, for rrsq we need to identify
		// this term
		if (derivOrder == 2 && derivOnSamePosAndDir()) {
			const Shell& s = sq.getShell(v1);  // this case v1 is same as v2
			if (item == 1 && s.getL() == 0) {
				ShellQuartet newSQ; // default empty shell quartet
				rrSQList.push_back(newSQ);
				continue;
			}
		}

		// now let's analyze the corresponding value
		Shell s1(sq.getShell(BRA1));
		Shell s2(sq.getShell(BRA2));
		Shell s3(sq.getShell(KET1));
		Shell s4(sq.getShell(KET2));
		bool hasNewI = true;
		for(int n=0; n<varNum; n++) {

			// get the variable and its angular momentum infor
			int newVar = v1;
			int angInc = ang1;
			if (n == 1) {
				newVar = v2;
				angInc = ang2;
			}

			// get the shell going to be change
			Shell oldShell;
			if (newVar == BRA1) {
				oldShell = s1;
			}else if(newVar == BRA2){
				oldShell = s2;
			}else if(newVar == KET1){
				oldShell = s3;
			}else {
				oldShell = s4;
			}
			int newL = oldShell.getL();
			newL += angInc;

			// replace the original shell with new one
			// we note that if the new shell is none, then this term 
			// does not exist 
			if (newL < 0) {
				hasNewI = false;
				break;
			}else{
				Shell newShell(newL);
				if (newVar == BRA1) {
					s1 = newShell;
				}else if(newVar == BRA2){
					s2 = newShell;
				}else if(newVar == KET1){
					s3 = newShell;
				}else {
					s4 = newShell;
				}
			}
		}

		// now if all of newV exisits, then integral exists; this integral exists
		if (hasNewI) {

			// here we have new shell quartet forming
			ShellQuartet newSQ(s1,s2,s3,s4,newO,newM,newDivision);

			// finally, for deriv work or the non-RR-deriv work 
			// we need to take care of exponential factor multipler
			if (isDerivWork()) {

				// for the deriv work, for each item we have 
				// number of deriv order possible multiplers
				// exps in dimension of derivOrder*nitems
				for(int i=0; i<derivOrder; i++) {
					int pos = exps[derivOrder*item+i];
					if (pos == BRA1 || pos == BRA2 || pos == KET1 || pos == KET2) {
						newSQ.addExpFac(pos); 
					}
				}
			}

			// now consider the non-deriv non-RR case
			if (isNonRRNonDerivWork()) {
				if (oper == THREEBODYKI) {
					int order = getNExpFacFor3BodyKi();
					for(int i=0; i<order; i++) {
						int pos = exps[order*item+i];
						if (pos == BRA1 || pos == BRA2 || pos == KET1) {
							newSQ.addExpFac(pos); 
						}
					}
				}
			}

			// for the HRR case, we need to pass the same exponential
			// information from the LHS sq to the RHS
			if (sq.withExpFac() && rrType == HRR) {
				newSQ.copyExpFac(sq);
			}

			// push the result newSQ in 
			rrSQList.push_back(newSQ);
		}else{
			ShellQuartet newSQ; // default empty shell quartet
			rrSQList.push_back(newSQ);
		}
	}
}
	
void RRBuild::print() const
{
	string operName = getOperStringName(oper);
	cout << "Operator is " << operName << endl; 
	if (isRRWork()) {
		cout << "RR expansion position is " << position << endl;
	}else if (isDerivWork()) {
		cout << "deriv order is " << derivOrder << endl;
		cout << "first  order derivative position is " << firstDerivPos << endl;
		cout << "second order derivative position is " << secondDerivPos << endl;
		cout << "first  order derivative direction is " << firstDerivDir << endl;
		cout << "second order derivative direction is " << secondDerivDir << endl;
	}

	int nitems = getNItems();
	for(int i=0; i<nitems; i++) {
		cout << "For the item          " << i+1 << endl;
		cout << "coefficients          " << coes[i] << endl;

		// does the integral type changes?
		if (intOperChanged(oper)){
			string oName = getOperStringName(intTypes[i]); 
			cout << "Operator              " << oName << endl;
		}
	
		// does M value are defined?
		// only for VRR process
		if (isVRRWork() && hasMValueDefined(oper)) {
			cout << "M value's Change      " << mVals[i] << endl;
		}

		// does the exp factor defined?
		if (! isRRWork()) {

			// number of terms
			int n = -1;
			if (isDerivWork()) n = derivOrder;
			if (isNonRRNonDerivWork()) {
				if (oper == THREEBODYKI) n = getNExpFacFor3BodyKi();
			}

			// now print them out
			for(int j=0; j<n; j++) {
				int pos = exps[n*i+j];
				if (pos == BRA1) {
					cout << "add in alpha exp factor" << endl;
				}else if (pos == BRA2) {
					cout << "add in beta  exp factor" << endl;
				}else if (pos == KET1) {
					cout << "add in gamma exp factor" << endl;
				}else if (pos == KET2) {
					cout << "add in delta exp factor" << endl;
				}
			}
		}

		// now change of variables
		int vp1  = lenVars[i];
		int v1   = vars[vp1];     // the first changing variable
		int ang1 = angs[vp1];     // its changing angular momentum 
		int varNum = lenVars[i+1] - lenVars[i];          
		if (varNum>1) {
			int v2   = NULL_POS;      // the second changing variable
			int ang2 = NULL_POS;      // its changing angular momentum
			int vp2  = vp1+1;
			v2   = vars[vp2];
			ang2 = angs[vp2];
			cout << "var1 " << v1 << " ang change " << ang1 << 
				" var2 " << v2 << " ang change " << ang2 << endl;
		}else{
			cout << "var1 " << v1 << " ang change " << ang1 << endl;
		}
		cout << endl;
	}
}

