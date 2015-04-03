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
#include "basis.h"
#include "shellsymbol.h"
#include "boost/lexical_cast.hpp"
#include "basisutil.h"
using boost::lexical_cast;
using namespace basisutil;
using namespace basis;

string Basis::getName() const {
	int L = l + m + n;
	string name = getShellSymbol(L);
	if (l>0) {
		if (l == 1) {
			name += "x";
		}else{
			string nx = lexical_cast<string>(l);
			name += nx + "x";
		}
	}
	if (m>0) {
		if (m == 1) {
			name += "y";
		}else{
			string ny = lexical_cast<string>(m);
			name += ny + "y";
		}
	}
	if (n>0) {
		if (n == 1) {
			name += "z";
		}else{
			string nz = lexical_cast<string>(n);
			name += nz + "z";
		}
	}
	return name;
}

int Basis::getLocalIndex() const
{
	BasisUtil bu;
	int index = bu.getLocalBasisSetIndex(l,m,n);
	return index;
}
