//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2014 Fenglai Liu
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
#include "boost/lexical_cast.hpp"
#include "inttype.h"
using namespace inttype;

//
// we note, that the two body eri and three body
// eri appears only as "names", however in real
// code generation they will be treated as real
// ERI codes
// ESP integrals in fact is same to the NAI
// therefore, in the operator part we only 
// provide the name access function for them
//

int inttype::getOperOrder(const int& name) {
	int order = -1;
	if (name == TWOBODYOVERLAP){
		order = 2;
	}else if (name == THREEBODYOVERLAP) {
		order = 3;
	}else if (name == THREEBODYKI) {
		order = 3;
	}else if (name == KINETIC) {
		order = 2;
	}else if (name == NAI) {
		order = 2;
	}else if (name == ESP) {
		order = 2;
	}else if (name == MOM) {
		order = 2;
	}else if (name == ERI) {
		order = 4;
	}else{
		crash(true, "Incorrect integral type name given in getOrder");
	}
	return order;
}

int inttype::getRROrder(const int& name) {
	int order = -1;
	if (name == TWOBODYOVERLAP){
		order = 2;
	}else if (name == THREEBODYOVERLAP) {
		order = 3;
	}else if (name == THREEBODYKI) {
		order = 3;
	}else if (name == KINETIC) {
		order = 2;
	}else if (name == NAI) {
		order = 2;
	}else if (name == ESP) {
		order = 2;
	}else if (name == MOM) {
		order = 3;
	}else if (name == ERI) {
		order = 4;
	}else{
		crash(true, "Incorrect integral type name given in getRROrder");
	}
	return order;
}

bool inttype::hasKetSide(const int& name) {
	int order = getOperOrder(name);
	if (order > 2) return true;
	return false;
}


string inttype::getOperStringName(const int& name) {
	if (name == TWOBODYOVERLAP){
		return "TWOBODYOVERLAP";
	}else if (name == THREEBODYOVERLAP) {
		return "THREEBODYOVERLAP";
	}else if (name == THREEBODYKI) {
		return "THREEBODYKI";
	}else if (name == KINETIC) {
		return "KINETIC";
	}else if (name == NAI) {
		return "NAI";
	}else if (name == ESP) {
		return "ESP";
	}else if (name == MOM) {
		return "MOM";
	}else if (name == ERI) {
		return "ERI";
	}else if (name == TWOBODYERI) {
		return "TWOBODYERI";
	}else if (name == THREEBODYERI) {
		return "THREEBODYERI";
	}else{
		crash(true, "Incorrect integral type name given in getStringName");
		return "NONE";
	}
}

int inttype::getOperIntName(const string& name) {
	if (name == "TWOBODYOVERLAP"){
		return TWOBODYOVERLAP;
	}else if (name == "THREEBODYOVERLAP") {
		return THREEBODYOVERLAP;
	}else if (name == "THREEBODYKI") {
		return THREEBODYKI;
	}else if (name == "KINETIC") {
		return KINETIC;
	}else if (name == "MOM") {
		return MOM;
	}else if (name == "NAI") {
		return NAI;
	}else if (name == "ESP") {
		return ESP;
	}else if (name == "ERI") {
		return ERI;
	}else if (name == "TWOBODYERI") {
		return TWOBODYERI;
	}else if (name == "THREEBODYERI") {
		return THREEBODYERI;
	}else{
		crash(true, "Incorrect integral type name given in getOperIntName");
		return -1;
	}
}

int inttype::compareOper(const int& O1, const int& O2) {

	// we will compare the given two operators
	// currently we just support the situation below
	// that one operator has gradient operator
	// the other one belongs to overlap type
	
	// let's check that O1 is gradient operator
	if (O1 == KINETIC || O1 == THREEBODYKI) {

		// whether O2 is overlap type operator?
		if (O2 == TWOBODYOVERLAP || O2 == THREEBODYOVERLAP) {
			return LT;
		}else{
			return FAIL_COMPARE;
		}

	// now the situation is that O1 is overlap type of operator
	}else if (O1 == TWOBODYOVERLAP || O1 == THREEBODYOVERLAP){

		// whether O2 is gradient type operator?
		if (O2 == KINETIC || O2 == THREEBODYKI) {
			return GT;
		}else{
			return FAIL_COMPARE;
		}

	}else{
		return FAIL_COMPARE;
	}
}

void inttype::selectOperListInRR(const int& O, vector<int>& oList) {

	if(O == KINETIC) {
		oList.reserve(2);
		oList.push_back(KINETIC);
		oList.push_back(TWOBODYOVERLAP);
	}else{
		cout << "Operator is " << getOperStringName(O) << endl;
		crash(true, "Invalid operator passed in selectOperListInRR");
	}
}

bool inttype::hasMValueDefined(const int& O) 
{
	if(O == NAI || O == ERI || O == ESP) return true;
	return false;
}

bool inttype::intOperChanged(const int& O) 
{
	if(O == THREEBODYKI || O == KINETIC) return true;
	return false;
}

bool inttype::useFmt(const int& O) 
{
	if(O == NAI || O == ERI || O == ESP) return true;
	return false;
}

bool inttype::sigCheck(const int& O) 
{
	if(O == ERI) return true;
	return false;
}

bool inttype::canDOHRR(const int& O) 
{
	if(O == KINETIC || O == THREEBODYKI || O == ESP) return false;
	return true;
}

int inttype::getNSpaceByOper(const int& oper)
{
	// order of the operator determines that 
	// how many loops we may have
	// we note, here we should have the operator order
	// rather than the rr order
	// for example, the MOM, only have one for loop over the bra side
	// even though it's RR order is 3
	int nSpace= 2;
	int nBody = getOperOrder(oper);
	int additionalSpace = 2;
	if(nBody>2) additionalSpace = 4; 
	if (oper == NAI || oper == ESP) {
		additionalSpace += 2; // this is because we have loop over natoms/ngrids
	}
	return nSpace + additionalSpace;
}

bool inttype::isNullOper(const int& oper) 
{
	if (oper<0) return true;
	return false;
}

int inttype::rrProp(const int& oper) 
{
	if (oper == ERI || oper == NAI || oper == ESP) {
		return SQ_ON_L_M;
	}else if (oper == KINETIC) {
		return SQ_ON_L_OPER;
	}else if (oper == TWOBODYOVERLAP || oper == THREEBODYOVERLAP || oper == MOM) {
		return SQ_ON_L;
	}else{
		crash(true,"Incorrect operator passed in rrProp function in inttype");
	}
	return -1;
}

bool inttype::complicatedBottomIntegrals(const int& oper)
{
	if (oper == MOM) return true;
	return false;
}

bool inttype::needIntDeriv(const int& oper)
{
	if (oper == THREEBODYKI) return true;
	return false;
}

int inttype::getDerivNumberFromOper(const int& oper, const int& jobOrder)
{
	if (jobOrder == 0) {
		if (oper == THREEBODYKI) {
			return 3;  // for x y and z
		}else{
			crash(true, "oper is not supported in getDerivNumberFromOper for job order = 0");
		}
	}else{
		crash(true, "job order> 0 is not supported in getDerivNumberFromOper");
	}
	return -1;
}

int inttype::getDerivOrderFromOper(const int& oper, const int& jobOrder)
{
	if (jobOrder == 0) {
		if (oper == THREEBODYKI) {
			return 1;  
		}else{
			crash(true, "oper is not supported in getDerivOrderFromOper for job order = 0");
		}
	}else{
		crash(true, "job order> 0 is not supported in getDerivOrderFromOper");
	}
	return -1;
}

bool inttype::resultIntegralHasAdditionalOffset(const int& oper)
{
	if (oper == ESP) return true;
	return false;
}

string inttype::determineAdditionalOffset(const int& oper, const int& nInts)
{
	if (oper == ESP) {
		if (nInts<=0) {
			crash(true, "Invalid nInts in determineAdditionalOffset");
		}
		string offset = "iGrid*";
		offset = offset + boost::lexical_cast<string>(nInts);
		return offset;
	}else{
		crash(true, "Invalid operation here in determineAdditionalOffset");
	}
	return "NONE";
}
