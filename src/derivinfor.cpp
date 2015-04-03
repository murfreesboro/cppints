//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2015 Fenglai Liu
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
#include "derivinfor.h"
using namespace derivinfor;

int derivinfor::getJobOrder(const int& deriv) {
	int order = -1;
	if (deriv == NO_DERIV) {
		order = 0;
	}else if (deriv == DERIV_X || deriv == DERIV_Y || deriv == DERIV_Z) {
		order = 1;
	}else if (deriv == DERIV_XX || deriv == DERIV_XY || deriv == DERIV_XZ || 
			deriv == DERIV_YY || deriv == DERIV_YZ || deriv == DERIV_ZZ ) {
		order = 2;
	}else{
		crash(true, "Currently we do not have job order > 2 for getJobOrder in derivinfor");
	}
	return order;
}

string derivinfor::getDerivDir(const int& deriv) {
	if (deriv == DERIV_X ) {
		return "X";
	}else if (deriv == DERIV_Y) {
		return "Y";
	}else if (deriv == DERIV_Z) {
		return "Z";
	}else{
		cout << "It could be only DERIV_X(Y or Z)" << endl;
		crash(true, "Illegal deriv provided in getDerivDir");
	}
	return "NONE";
}

bool derivinfor::isDerivInfor(const int& pos) {
	if (pos == DERIV_X || pos == DERIV_Y || pos == DERIV_Z) return true;
	return false;
}

bool derivinfor::derivCompare(const int& pos1, const int& pos2) {
	// here we simply compare the macro of the deriv 
	return (pos1 < pos2);
}

int derivinfor::getNDeriv(const int& order) {
	if (order == 1) return 3;
	if (order == 2) return 6;
	if (order == 3) return 10;
	if (order == 4) return 15;
	if (order>4 || order<1) {
		cout << "order is " << order << endl;
		crash(true, "input order in getNDeriv is not valid");
	}
	return -1;
}

void derivinfor::formDerivOrderArray(const int& order, vector<int>& derivOrder) {
	if (order == 1) {
		derivOrder[0] = DERIV_X;
		derivOrder[1] = DERIV_Y;
		derivOrder[2] = DERIV_Z;
	}else{
		crash(true,"Currently we do not have order > 1 in formDerivOrderArray");
	}
}

string derivinfor::symTransform(const int& var) {
	if (var == DERIV_X) {
		return "DX";
	} else if (var == DERIV_Y) {
		return "DY";
	} else if (var == DERIV_Z) {
		return "DZ";
	}else{
		crash(true,"Currently we do not support this var in symtransform");
	}
	return "NONE";
}

int derivinfor::getDerivVar(const int& index, const int& order) {
	if (order == 1) {
		if (index == 0) return DERIV_X;
		if (index == 1) return DERIV_Y;
		if (index == 2) return DERIV_Z;
	}else{
		crash(true,"Currently we do not support this order in getDerivVar");
	}
	return NULL_POS;
}
