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
 * \file   printing.h
 * \brief  inline function for printing
 * \author Fenglai Liu 
 */
#ifndef PRINTING_H
#define PRINTING_H
#include "general.h"

namespace printing {

	/**
	 * print the given line with the file stream
	 * print it with the given number of space
	 */
	inline void printLine(const int& nSpace, const string& line, ofstream& file) {
		switch(nSpace) {
			case 0: 
				file << line << endl;
				break;
			case 2: 
				file << "  " << line << endl;
				break;
			case 4: 
				file << "    " << line << endl;
				break;
			case 6: 
				file << "      " << line << endl;
				break;
			case 8: 
				file << "        " << line << endl;
				break;
			case 10: 
				file << "          " << line << endl;
				break;
			case 12: 
				file << "            " << line << endl;
				break;
			case 14: 
				file << "              " << line << endl;
				break;
			default: 
				crash(true, "invalid nspace given in printLine");
				break;
		}
	};

}


#endif

