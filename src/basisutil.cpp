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
#include "basisutil.h"
#include "shellsymbol.h"
#include <boost/algorithm/string.hpp> 
using namespace basisutil;
using namespace boost;

void BasisUtil::getLMNVal(const int& L, const int& index,
		int& l, int& m, int& n) const 
{
	if (L<=MAX_L) {
		int ii = L*(L+1)*(L+2)/6 + index;
		l = LIBINT_BASIS_SET_ORDER[3*ii  ];
		m = LIBINT_BASIS_SET_ORDER[3*ii+1];
		n = LIBINT_BASIS_SET_ORDER[3*ii+2];
	}else{
		// now we have to use the libint order to generate th l m n
		int nx = -1;
		int ny = -1;
		int nz = -1;
		int count = 0;
		for(int x=0; x<=L; x++) {
			nx = L - x;
			for(int y=0; y<=x; y++) {
				ny = x-y; 
				nz = y;   
				if (count == index) {
					l = nx;
					m = ny;
					n = nz;
					return;
				}
				count++;
			}
		}
	}
}

int BasisUtil::getLocalBasisSetIndex(const int& l, const int& m,
		const int& n) const 
{
	// actually we note that since here we calculate the local index
	// within the shell, therefore we only need the m and n value actually
	// however, we still pass in the l,m,n in full since we do not want 
	// to mistake about whether it's l,m passed or m,n passed
	int index = -1;
	int x =  m + n;
	index = (x*(x+1))/2 + n;
	return index;
}

