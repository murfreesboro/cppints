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
 * \file    derivinfor.h
 * \brief   parsing the derivatives information 
 * \author  Fenglai Liu
 */
#ifndef DERIVINFOR_H
#define DERIVINFOR_H
#include "general.h"

namespace derivinfor {

	///
	/// according to the input of deriv infor (defined in general.h)
	/// let's say what is the derivatives order for it
	///
	int getJobOrder(const int& deriv);

	///
	/// get the derivative direction
	/// could be only X, Y or Z
	///
	string getDerivDir(const int& deriv);

	///
	/// whether this is derivatives symbol defined in general.h? 
	///
	bool isDerivInfor(const int& pos);

	///
	/// compare two derivative symbols
	/// we do not do checking inside!!!
	///
	bool derivCompare(const int& pos1, const int& pos2); 

	///
	/// get the number of derivatives for the given order
	///
	int getNDeriv(const int& order);

	///
	/// form the derivatives order array for the given order
	///
	void formDerivOrderArray(const int& order, vector<int>& derivOrder);

	///
	/// transform the integer derivatives into it's symbol form
	///
	string symTransform(const int& var); 

	///
	/// add in a small function to get var from given index
	/// \param order  the derivatives order
	/// \param index  the derivative variable index
	///
	int getDerivVar(const int& index, const int& order);
}

#endif
