/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2015 Fenglai Liu
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
 * \file    compositeshell.h
 * \brief   parsing the input composite shells information
 * \author  Fenglai Liu
 */
#ifndef COMPOSITESHELL_H
#define COMPOSITESHELL_H
#include "general.h"

namespace compositeshell {

	/**
	 * possible input shells defined in the program
	 * maximum angular momentum is I
	 * 10 is the SP shell
	 * the sequence means the shell order
	 * for example, SP<P<D etc.
	 */
	const int SHELL_ANG_MOM_CODE[] = {
		0, 10, 1, 2, 3, 4, 5, 6
	};

	/**
	 * number of input shell types
	 */
	const int MAX_INPUT_SHELL_TYPES = 8;

	///
	/// get the minimum L for this composite shell
	///
	inline int getLMin(const int& code) {
		if (code < 10) {
			return code;
		}else{
			int lmax = code/10;
			return code - lmax*10;
		}
	};

	///
	/// get the maximum L for this composite shell
	///
	inline int getLMax(const int& code) {
		if (code < 10) {
			return code;
		}else{
			return code/10;
		}
	};

	///
	/// for the given shell code, we count number of single shells inside
	/// just list all of situations here
	///
	inline int getNShells(const int& code) {
		if (code == 10) return 2;
		return 1;
	};

	///
	/// for the given group of shell code (in vector form), we calculate
	/// the number of single shell quartet here
	///
	inline int getNSQ(const vector<int>& codeList) {
		int n = 1;
		for(int iCode=0; iCode<(int)codeList.size(); iCode++) {
			n *= getNShells(codeList[iCode]);
		}
		return n;
	};


}


#endif

