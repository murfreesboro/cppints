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
#include "shellsymbol.h"
#include "basisutil.h"
#include "basis.h"
#include "shell.h"
using namespace basis;
using namespace basisutil;
using namespace shell;

void Shell::getBasisSetFromIndex(const int& index, int& l, int& m, int& n) const {
	BasisUtil bu;
	bu.getLMNVal(L,index,l,m,n);
}

bool Shell::hasThisBasisSet(const Basis& b) const {
	int totalL = b.getL();
	if (totalL >= 0 && totalL == L) return true;
	return false;
}

void Shell::getBasis(vector<Basis>& basis) const {
	int num = getBasisSetNumber();
	basis.reserve(num);
	BasisUtil bu;
	int l = -1;
	int m = -1;
	int n = -1;
	for(int i=0; i<num; i++) {
		bu.getLMNVal(L,i,l,m,n);
		Basis b(l,m,n);
		basis.push_back(b);
	}
}

string Shell::getName() const {
	string name = getShellSymbol(L);
	return name;
}
