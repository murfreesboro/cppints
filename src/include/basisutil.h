/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2012-2014 Fenglai Liu
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
 * \file    basisutil.h
 * \brief   this is used to generate the basis set orders
 * \author  Fenglai Liu
 * \date    Nov, 2012
 * \note
 *
 * In this file, we define the ordering of basis sets used in the integral
 * generation. We note that for basis set order, you need to give the arrange
 * orders in explicit way. Now in the program we only have the libint order.
 *
 * The basis set order is used to map the connection between basis set local
 * order and the l,m,n value for each basis set(see basis.h). By specifically
 * listed all of basis set orders(up to 20), we avoid to get the l,m,n through
 * for loop so save some time.
 */
#ifndef BASISUTIL_H
#define BASISUTIL_H
#include "general.h"

namespace basisutil {

	//
	// this is generated according to the MAX_L defined in shellsymbol.h
	//
	const int LIBINT_BASIS_SET_ORDER[ ] = {
		0,  0,  0,  /*L = 0*/
		1,  0,  0,  /*L = 1*/
		0,  1,  0,
		0,  0,  1,
		2,  0,  0,  /*L = 2*/
		1,  1,  0,
		1,  0,  1,
		0,  2,  0,
		0,  1,  1,
		0,  0,  2,
		3,  0,  0,  /*L = 3*/
		2,  1,  0,
		2,  0,  1,
		1,  2,  0,
		1,  1,  1,
		1,  0,  2,
		0,  3,  0,
		0,  2,  1,
		0,  1,  2,
		0,  0,  3,
		4,  0,  0,  /*L = 4*/
		3,  1,  0,
		3,  0,  1,
		2,  2,  0,
		2,  1,  1,
		2,  0,  2,
		1,  3,  0,
		1,  2,  1,
		1,  1,  2,
		1,  0,  3,
		0,  4,  0,
		0,  3,  1,
		0,  2,  2,
		0,  1,  3,
		0,  0,  4,
		5,  0,  0,  /*L = 5*/
		4,  1,  0,
		4,  0,  1,
		3,  2,  0,
		3,  1,  1,
		3,  0,  2,
		2,  3,  0,
		2,  2,  1,
		2,  1,  2,
		2,  0,  3,
		1,  4,  0,
		1,  3,  1,
		1,  2,  2,
		1,  1,  3,
		1,  0,  4,
		0,  5,  0,
		0,  4,  1,
		0,  3,  2,
		0,  2,  3,
		0,  1,  4,
		0,  0,  5,
		6,  0,  0,  /*L = 6*/
		5,  1,  0,
		5,  0,  1,
		4,  2,  0,
		4,  1,  1,
		4,  0,  2,
		3,  3,  0,
		3,  2,  1,
		3,  1,  2,
		3,  0,  3,
		2,  4,  0,
		2,  3,  1,
		2,  2,  2,
		2,  1,  3,
		2,  0,  4,
		1,  5,  0,
		1,  4,  1,
		1,  3,  2,
		1,  2,  3,
		1,  1,  4,
		1,  0,  5,
		0,  6,  0,
		0,  5,  1,
		0,  4,  2,
		0,  3,  3,
		0,  2,  4,
		0,  1,  5,
		0,  0,  6,
		7,  0,  0,  /*L = 7*/
		6,  1,  0,
		6,  0,  1,
		5,  2,  0,
		5,  1,  1,
		5,  0,  2,
		4,  3,  0,
		4,  2,  1,
		4,  1,  2,
		4,  0,  3,
		3,  4,  0,
		3,  3,  1,
		3,  2,  2,
		3,  1,  3,
		3,  0,  4,
		2,  5,  0,
		2,  4,  1,
		2,  3,  2,
		2,  2,  3,
		2,  1,  4,
		2,  0,  5,
		1,  6,  0,
		1,  5,  1,
		1,  4,  2,
		1,  3,  3,
		1,  2,  4,
		1,  1,  5,
		1,  0,  6,
		0,  7,  0,
		0,  6,  1,
		0,  5,  2,
		0,  4,  3,
		0,  3,  4,
		0,  2,  5,
		0,  1,  6,
		0,  0,  7,
		8,  0,  0,  /*L = 8*/
		7,  1,  0,
		7,  0,  1,
		6,  2,  0,
		6,  1,  1,
		6,  0,  2,
		5,  3,  0,
		5,  2,  1,
		5,  1,  2,
		5,  0,  3,
		4,  4,  0,
		4,  3,  1,
		4,  2,  2,
		4,  1,  3,
		4,  0,  4,
		3,  5,  0,
		3,  4,  1,
		3,  3,  2,
		3,  2,  3,
		3,  1,  4,
		3,  0,  5,
		2,  6,  0,
		2,  5,  1,
		2,  4,  2,
		2,  3,  3,
		2,  2,  4,
		2,  1,  5,
		2,  0,  6,
		1,  7,  0,
		1,  6,  1,
		1,  5,  2,
		1,  4,  3,
		1,  3,  4,
		1,  2,  5,
		1,  1,  6,
		1,  0,  7,
		0,  8,  0,
		0,  7,  1,
		0,  6,  2,
		0,  5,  3,
		0,  4,  4,
		0,  3,  5,
		0,  2,  6,
		0,  1,  7,
		0,  0,  8,
		9,  0,  0,  /*L = 9*/
		8,  1,  0,
		8,  0,  1,
		7,  2,  0,
		7,  1,  1,
		7,  0,  2,
		6,  3,  0,
		6,  2,  1,
		6,  1,  2,
		6,  0,  3,
		5,  4,  0,
		5,  3,  1,
		5,  2,  2,
		5,  1,  3,
		5,  0,  4,
		4,  5,  0,
		4,  4,  1,
		4,  3,  2,
		4,  2,  3,
		4,  1,  4,
		4,  0,  5,
		3,  6,  0,
		3,  5,  1,
		3,  4,  2,
		3,  3,  3,
		3,  2,  4,
		3,  1,  5,
		3,  0,  6,
		2,  7,  0,
		2,  6,  1,
		2,  5,  2,
		2,  4,  3,
		2,  3,  4,
		2,  2,  5,
		2,  1,  6,
		2,  0,  7,
		1,  8,  0,
		1,  7,  1,
		1,  6,  2,
		1,  5,  3,
		1,  4,  4,
		1,  3,  5,
		1,  2,  6,
		1,  1,  7,
		1,  0,  8,
		0,  9,  0,
		0,  8,  1,
		0,  7,  2,
		0,  6,  3,
		0,  5,  4,
		0,  4,  5,
		0,  3,  6,
		0,  2,  7,
		0,  1,  8,
		0,  0,  9,
		10,  0,  0,  /*L = 10*/
		9,  1,  0,
		9,  0,  1,
		8,  2,  0,
		8,  1,  1,
		8,  0,  2,
		7,  3,  0,
		7,  2,  1,
		7,  1,  2,
		7,  0,  3,
		6,  4,  0,
		6,  3,  1,
		6,  2,  2,
		6,  1,  3,
		6,  0,  4,
		5,  5,  0,
		5,  4,  1,
		5,  3,  2,
		5,  2,  3,
		5,  1,  4,
		5,  0,  5,
		4,  6,  0,
		4,  5,  1,
		4,  4,  2,
		4,  3,  3,
		4,  2,  4,
		4,  1,  5,
		4,  0,  6,
		3,  7,  0,
		3,  6,  1,
		3,  5,  2,
		3,  4,  3,
		3,  3,  4,
		3,  2,  5,
		3,  1,  6,
		3,  0,  7,
		2,  8,  0,
		2,  7,  1,
		2,  6,  2,
		2,  5,  3,
		2,  4,  4,
		2,  3,  5,
		2,  2,  6,
		2,  1,  7,
		2,  0,  8,
		1,  9,  0,
		1,  8,  1,
		1,  7,  2,
		1,  6,  3,
		1,  5,  4,
		1,  4,  5,
		1,  3,  6,
		1,  2,  7,
		1,  1,  8,
		1,  0,  9,
		0,  10,  0,
		0,  9,  1,
		0,  8,  2,
		0,  7,  3,
		0,  6,  4,
		0,  5,  5,
		0,  4,  6,
		0,  3,  7,
		0,  2,  8,
		0,  1,  9,
		0,  0,  10,
		11,  0,  0,  /*L = 11*/
		10,  1,  0,
		10,  0,  1,
		9,  2,  0,
		9,  1,  1,
		9,  0,  2,
		8,  3,  0,
		8,  2,  1,
		8,  1,  2,
		8,  0,  3,
		7,  4,  0,
		7,  3,  1,
		7,  2,  2,
		7,  1,  3,
		7,  0,  4,
		6,  5,  0,
		6,  4,  1,
		6,  3,  2,
		6,  2,  3,
		6,  1,  4,
		6,  0,  5,
		5,  6,  0,
		5,  5,  1,
		5,  4,  2,
		5,  3,  3,
		5,  2,  4,
		5,  1,  5,
		5,  0,  6,
		4,  7,  0,
		4,  6,  1,
		4,  5,  2,
		4,  4,  3,
		4,  3,  4,
		4,  2,  5,
		4,  1,  6,
		4,  0,  7,
		3,  8,  0,
		3,  7,  1,
		3,  6,  2,
		3,  5,  3,
		3,  4,  4,
		3,  3,  5,
		3,  2,  6,
		3,  1,  7,
		3,  0,  8,
		2,  9,  0,
		2,  8,  1,
		2,  7,  2,
		2,  6,  3,
		2,  5,  4,
		2,  4,  5,
		2,  3,  6,
		2,  2,  7,
		2,  1,  8,
		2,  0,  9,
		1,  10,  0,
		1,  9,  1,
		1,  8,  2,
		1,  7,  3,
		1,  6,  4,
		1,  5,  5,
		1,  4,  6,
		1,  3,  7,
		1,  2,  8,
		1,  1,  9,
		1,  0,  10,
		0,  11,  0,
		0,  10,  1,
		0,  9,  2,
		0,  8,  3,
		0,  7,  4,
		0,  6,  5,
		0,  5,  6,
		0,  4,  7,
		0,  3,  8,
		0,  2,  9,
		0,  1,  10,
		0,  0,  11,
		12,  0,  0,  /*L = 12*/
		11,  1,  0,
		11,  0,  1,
		10,  2,  0,
		10,  1,  1,
		10,  0,  2,
		9,  3,  0,
		9,  2,  1,
		9,  1,  2,
		9,  0,  3,
		8,  4,  0,
		8,  3,  1,
		8,  2,  2,
		8,  1,  3,
		8,  0,  4,
		7,  5,  0,
		7,  4,  1,
		7,  3,  2,
		7,  2,  3,
		7,  1,  4,
		7,  0,  5,
		6,  6,  0,
		6,  5,  1,
		6,  4,  2,
		6,  3,  3,
		6,  2,  4,
		6,  1,  5,
		6,  0,  6,
		5,  7,  0,
		5,  6,  1,
		5,  5,  2,
		5,  4,  3,
		5,  3,  4,
		5,  2,  5,
		5,  1,  6,
		5,  0,  7,
		4,  8,  0,
		4,  7,  1,
		4,  6,  2,
		4,  5,  3,
		4,  4,  4,
		4,  3,  5,
		4,  2,  6,
		4,  1,  7,
		4,  0,  8,
		3,  9,  0,
		3,  8,  1,
		3,  7,  2,
		3,  6,  3,
		3,  5,  4,
		3,  4,  5,
		3,  3,  6,
		3,  2,  7,
		3,  1,  8,
		3,  0,  9,
		2,  10,  0,
		2,  9,  1,
		2,  8,  2,
		2,  7,  3,
		2,  6,  4,
		2,  5,  5,
		2,  4,  6,
		2,  3,  7,
		2,  2,  8,
		2,  1,  9,
		2,  0,  10,
		1,  11,  0,
		1,  10,  1,
		1,  9,  2,
		1,  8,  3,
		1,  7,  4,
		1,  6,  5,
		1,  5,  6,
		1,  4,  7,
		1,  3,  8,
		1,  2,  9,
		1,  1,  10,
		1,  0,  11,
		0,  12,  0,
		0,  11,  1,
		0,  10,  2,
		0,  9,  3,
		0,  8,  4,
		0,  7,  5,
		0,  6,  6,
		0,  5,  7,
		0,  4,  8,
		0,  3,  9,
		0,  2,  10,
		0,  1,  11,
		0,  0,  12,
		13,  0,  0,  /*L = 13*/
		12,  1,  0,
		12,  0,  1,
		11,  2,  0,
		11,  1,  1,
		11,  0,  2,
		10,  3,  0,
		10,  2,  1,
		10,  1,  2,
		10,  0,  3,
		9,  4,  0,
		9,  3,  1,
		9,  2,  2,
		9,  1,  3,
		9,  0,  4,
		8,  5,  0,
		8,  4,  1,
		8,  3,  2,
		8,  2,  3,
		8,  1,  4,
		8,  0,  5,
		7,  6,  0,
		7,  5,  1,
		7,  4,  2,
		7,  3,  3,
		7,  2,  4,
		7,  1,  5,
		7,  0,  6,
		6,  7,  0,
		6,  6,  1,
		6,  5,  2,
		6,  4,  3,
		6,  3,  4,
		6,  2,  5,
		6,  1,  6,
		6,  0,  7,
		5,  8,  0,
		5,  7,  1,
		5,  6,  2,
		5,  5,  3,
		5,  4,  4,
		5,  3,  5,
		5,  2,  6,
		5,  1,  7,
		5,  0,  8,
		4,  9,  0,
		4,  8,  1,
		4,  7,  2,
		4,  6,  3,
		4,  5,  4,
		4,  4,  5,
		4,  3,  6,
		4,  2,  7,
		4,  1,  8,
		4,  0,  9,
		3,  10,  0,
		3,  9,  1,
		3,  8,  2,
		3,  7,  3,
		3,  6,  4,
		3,  5,  5,
		3,  4,  6,
		3,  3,  7,
		3,  2,  8,
		3,  1,  9,
		3,  0,  10,
		2,  11,  0,
		2,  10,  1,
		2,  9,  2,
		2,  8,  3,
		2,  7,  4,
		2,  6,  5,
		2,  5,  6,
		2,  4,  7,
		2,  3,  8,
		2,  2,  9,
		2,  1,  10,
		2,  0,  11,
		1,  12,  0,
		1,  11,  1,
		1,  10,  2,
		1,  9,  3,
		1,  8,  4,
		1,  7,  5,
		1,  6,  6,
		1,  5,  7,
		1,  4,  8,
		1,  3,  9,
		1,  2,  10,
		1,  1,  11,
		1,  0,  12,
		0,  13,  0,
		0,  12,  1,
		0,  11,  2,
		0,  10,  3,
		0,  9,  4,
		0,  8,  5,
		0,  7,  6,
		0,  6,  7,
		0,  5,  8,
		0,  4,  9,
		0,  3,  10,
		0,  2,  11,
		0,  1,  12,
		0,  0,  13,
		14,  0,  0,  /*L = 14*/
		13,  1,  0,
		13,  0,  1,
		12,  2,  0,
		12,  1,  1,
		12,  0,  2,
		11,  3,  0,
		11,  2,  1,
		11,  1,  2,
		11,  0,  3,
		10,  4,  0,
		10,  3,  1,
		10,  2,  2,
		10,  1,  3,
		10,  0,  4,
		9,  5,  0,
		9,  4,  1,
		9,  3,  2,
		9,  2,  3,
		9,  1,  4,
		9,  0,  5,
		8,  6,  0,
		8,  5,  1,
		8,  4,  2,
		8,  3,  3,
		8,  2,  4,
		8,  1,  5,
		8,  0,  6,
		7,  7,  0,
		7,  6,  1,
		7,  5,  2,
		7,  4,  3,
		7,  3,  4,
		7,  2,  5,
		7,  1,  6,
		7,  0,  7,
		6,  8,  0,
		6,  7,  1,
		6,  6,  2,
		6,  5,  3,
		6,  4,  4,
		6,  3,  5,
		6,  2,  6,
		6,  1,  7,
		6,  0,  8,
		5,  9,  0,
		5,  8,  1,
		5,  7,  2,
		5,  6,  3,
		5,  5,  4,
		5,  4,  5,
		5,  3,  6,
		5,  2,  7,
		5,  1,  8,
		5,  0,  9,
		4,  10,  0,
		4,  9,  1,
		4,  8,  2,
		4,  7,  3,
		4,  6,  4,
		4,  5,  5,
		4,  4,  6,
		4,  3,  7,
		4,  2,  8,
		4,  1,  9,
		4,  0,  10,
		3,  11,  0,
		3,  10,  1,
		3,  9,  2,
		3,  8,  3,
		3,  7,  4,
		3,  6,  5,
		3,  5,  6,
		3,  4,  7,
		3,  3,  8,
		3,  2,  9,
		3,  1,  10,
		3,  0,  11,
		2,  12,  0,
		2,  11,  1,
		2,  10,  2,
		2,  9,  3,
		2,  8,  4,
		2,  7,  5,
		2,  6,  6,
		2,  5,  7,
		2,  4,  8,
		2,  3,  9,
		2,  2,  10,
		2,  1,  11,
		2,  0,  12,
		1,  13,  0,
		1,  12,  1,
		1,  11,  2,
		1,  10,  3,
		1,  9,  4,
		1,  8,  5,
		1,  7,  6,
		1,  6,  7,
		1,  5,  8,
		1,  4,  9,
		1,  3,  10,
		1,  2,  11,
		1,  1,  12,
		1,  0,  13,
		0,  14,  0,
		0,  13,  1,
		0,  12,  2,
		0,  11,  3,
		0,  10,  4,
		0,  9,  5,
		0,  8,  6,
		0,  7,  7,
		0,  6,  8,
		0,  5,  9,
		0,  4,  10,
		0,  3,  11,
		0,  2,  12,
		0,  1,  13,
		0,  0,  14,
		15,  0,  0,  /*L = 15*/
		14,  1,  0,
		14,  0,  1,
		13,  2,  0,
		13,  1,  1,
		13,  0,  2,
		12,  3,  0,
		12,  2,  1,
		12,  1,  2,
		12,  0,  3,
		11,  4,  0,
		11,  3,  1,
		11,  2,  2,
		11,  1,  3,
		11,  0,  4,
		10,  5,  0,
		10,  4,  1,
		10,  3,  2,
		10,  2,  3,
		10,  1,  4,
		10,  0,  5,
		9,  6,  0,
		9,  5,  1,
		9,  4,  2,
		9,  3,  3,
		9,  2,  4,
		9,  1,  5,
		9,  0,  6,
		8,  7,  0,
		8,  6,  1,
		8,  5,  2,
		8,  4,  3,
		8,  3,  4,
		8,  2,  5,
		8,  1,  6,
		8,  0,  7,
		7,  8,  0,
		7,  7,  1,
		7,  6,  2,
		7,  5,  3,
		7,  4,  4,
		7,  3,  5,
		7,  2,  6,
		7,  1,  7,
		7,  0,  8,
		6,  9,  0,
		6,  8,  1,
		6,  7,  2,
		6,  6,  3,
		6,  5,  4,
		6,  4,  5,
		6,  3,  6,
		6,  2,  7,
		6,  1,  8,
		6,  0,  9,
		5,  10,  0,
		5,  9,  1,
		5,  8,  2,
		5,  7,  3,
		5,  6,  4,
		5,  5,  5,
		5,  4,  6,
		5,  3,  7,
		5,  2,  8,
		5,  1,  9,
		5,  0,  10,
		4,  11,  0,
		4,  10,  1,
		4,  9,  2,
		4,  8,  3,
		4,  7,  4,
		4,  6,  5,
		4,  5,  6,
		4,  4,  7,
		4,  3,  8,
		4,  2,  9,
		4,  1,  10,
		4,  0,  11,
		3,  12,  0,
		3,  11,  1,
		3,  10,  2,
		3,  9,  3,
		3,  8,  4,
		3,  7,  5,
		3,  6,  6,
		3,  5,  7,
		3,  4,  8,
		3,  3,  9,
		3,  2,  10,
		3,  1,  11,
		3,  0,  12,
		2,  13,  0,
		2,  12,  1,
		2,  11,  2,
		2,  10,  3,
		2,  9,  4,
		2,  8,  5,
		2,  7,  6,
		2,  6,  7,
		2,  5,  8,
		2,  4,  9,
		2,  3,  10,
		2,  2,  11,
		2,  1,  12,
		2,  0,  13,
		1,  14,  0,
		1,  13,  1,
		1,  12,  2,
		1,  11,  3,
		1,  10,  4,
		1,  9,  5,
		1,  8,  6,
		1,  7,  7,
		1,  6,  8,
		1,  5,  9,
		1,  4,  10,
		1,  3,  11,
		1,  2,  12,
		1,  1,  13,
		1,  0,  14,
		0,  15,  0,
		0,  14,  1,
		0,  13,  2,
		0,  12,  3,
		0,  11,  4,
		0,  10,  5,
		0,  9,  6,
		0,  8,  7,
		0,  7,  8,
		0,  6,  9,
		0,  5,  10,
		0,  4,  11,
		0,  3,  12,
		0,  2,  13,
		0,  1,  14,
		0,  0,  15,
		16,  0,  0,  /*L = 16*/
		15,  1,  0,
		15,  0,  1,
		14,  2,  0,
		14,  1,  1,
		14,  0,  2,
		13,  3,  0,
		13,  2,  1,
		13,  1,  2,
		13,  0,  3,
		12,  4,  0,
		12,  3,  1,
		12,  2,  2,
		12,  1,  3,
		12,  0,  4,
		11,  5,  0,
		11,  4,  1,
		11,  3,  2,
		11,  2,  3,
		11,  1,  4,
		11,  0,  5,
		10,  6,  0,
		10,  5,  1,
		10,  4,  2,
		10,  3,  3,
		10,  2,  4,
		10,  1,  5,
		10,  0,  6,
		9,  7,  0,
		9,  6,  1,
		9,  5,  2,
		9,  4,  3,
		9,  3,  4,
		9,  2,  5,
		9,  1,  6,
		9,  0,  7,
		8,  8,  0,
		8,  7,  1,
		8,  6,  2,
		8,  5,  3,
		8,  4,  4,
		8,  3,  5,
		8,  2,  6,
		8,  1,  7,
		8,  0,  8,
		7,  9,  0,
		7,  8,  1,
		7,  7,  2,
		7,  6,  3,
		7,  5,  4,
		7,  4,  5,
		7,  3,  6,
		7,  2,  7,
		7,  1,  8,
		7,  0,  9,
		6,  10,  0,
		6,  9,  1,
		6,  8,  2,
		6,  7,  3,
		6,  6,  4,
		6,  5,  5,
		6,  4,  6,
		6,  3,  7,
		6,  2,  8,
		6,  1,  9,
		6,  0,  10,
		5,  11,  0,
		5,  10,  1,
		5,  9,  2,
		5,  8,  3,
		5,  7,  4,
		5,  6,  5,
		5,  5,  6,
		5,  4,  7,
		5,  3,  8,
		5,  2,  9,
		5,  1,  10,
		5,  0,  11,
		4,  12,  0,
		4,  11,  1,
		4,  10,  2,
		4,  9,  3,
		4,  8,  4,
		4,  7,  5,
		4,  6,  6,
		4,  5,  7,
		4,  4,  8,
		4,  3,  9,
		4,  2,  10,
		4,  1,  11,
		4,  0,  12,
		3,  13,  0,
		3,  12,  1,
		3,  11,  2,
		3,  10,  3,
		3,  9,  4,
		3,  8,  5,
		3,  7,  6,
		3,  6,  7,
		3,  5,  8,
		3,  4,  9,
		3,  3,  10,
		3,  2,  11,
		3,  1,  12,
		3,  0,  13,
		2,  14,  0,
		2,  13,  1,
		2,  12,  2,
		2,  11,  3,
		2,  10,  4,
		2,  9,  5,
		2,  8,  6,
		2,  7,  7,
		2,  6,  8,
		2,  5,  9,
		2,  4,  10,
		2,  3,  11,
		2,  2,  12,
		2,  1,  13,
		2,  0,  14,
		1,  15,  0,
		1,  14,  1,
		1,  13,  2,
		1,  12,  3,
		1,  11,  4,
		1,  10,  5,
		1,  9,  6,
		1,  8,  7,
		1,  7,  8,
		1,  6,  9,
		1,  5,  10,
		1,  4,  11,
		1,  3,  12,
		1,  2,  13,
		1,  1,  14,
		1,  0,  15,
		0,  16,  0,
		0,  15,  1,
		0,  14,  2,
		0,  13,  3,
		0,  12,  4,
		0,  11,  5,
		0,  10,  6,
		0,  9,  7,
		0,  8,  8,
		0,  7,  9,
		0,  6,  10,
		0,  5,  11,
		0,  4,  12,
		0,  3,  13,
		0,  2,  14,
		0,  1,  15,
		0,  0,  16,
		17,  0,  0,  /*L = 17*/
		16,  1,  0,
		16,  0,  1,
		15,  2,  0,
		15,  1,  1,
		15,  0,  2,
		14,  3,  0,
		14,  2,  1,
		14,  1,  2,
		14,  0,  3,
		13,  4,  0,
		13,  3,  1,
		13,  2,  2,
		13,  1,  3,
		13,  0,  4,
		12,  5,  0,
		12,  4,  1,
		12,  3,  2,
		12,  2,  3,
		12,  1,  4,
		12,  0,  5,
		11,  6,  0,
		11,  5,  1,
		11,  4,  2,
		11,  3,  3,
		11,  2,  4,
		11,  1,  5,
		11,  0,  6,
		10,  7,  0,
		10,  6,  1,
		10,  5,  2,
		10,  4,  3,
		10,  3,  4,
		10,  2,  5,
		10,  1,  6,
		10,  0,  7,
		9,  8,  0,
		9,  7,  1,
		9,  6,  2,
		9,  5,  3,
		9,  4,  4,
		9,  3,  5,
		9,  2,  6,
		9,  1,  7,
		9,  0,  8,
		8,  9,  0,
		8,  8,  1,
		8,  7,  2,
		8,  6,  3,
		8,  5,  4,
		8,  4,  5,
		8,  3,  6,
		8,  2,  7,
		8,  1,  8,
		8,  0,  9,
		7,  10,  0,
		7,  9,  1,
		7,  8,  2,
		7,  7,  3,
		7,  6,  4,
		7,  5,  5,
		7,  4,  6,
		7,  3,  7,
		7,  2,  8,
		7,  1,  9,
		7,  0,  10,
		6,  11,  0,
		6,  10,  1,
		6,  9,  2,
		6,  8,  3,
		6,  7,  4,
		6,  6,  5,
		6,  5,  6,
		6,  4,  7,
		6,  3,  8,
		6,  2,  9,
		6,  1,  10,
		6,  0,  11,
		5,  12,  0,
		5,  11,  1,
		5,  10,  2,
		5,  9,  3,
		5,  8,  4,
		5,  7,  5,
		5,  6,  6,
		5,  5,  7,
		5,  4,  8,
		5,  3,  9,
		5,  2,  10,
		5,  1,  11,
		5,  0,  12,
		4,  13,  0,
		4,  12,  1,
		4,  11,  2,
		4,  10,  3,
		4,  9,  4,
		4,  8,  5,
		4,  7,  6,
		4,  6,  7,
		4,  5,  8,
		4,  4,  9,
		4,  3,  10,
		4,  2,  11,
		4,  1,  12,
		4,  0,  13,
		3,  14,  0,
		3,  13,  1,
		3,  12,  2,
		3,  11,  3,
		3,  10,  4,
		3,  9,  5,
		3,  8,  6,
		3,  7,  7,
		3,  6,  8,
		3,  5,  9,
		3,  4,  10,
		3,  3,  11,
		3,  2,  12,
		3,  1,  13,
		3,  0,  14,
		2,  15,  0,
		2,  14,  1,
		2,  13,  2,
		2,  12,  3,
		2,  11,  4,
		2,  10,  5,
		2,  9,  6,
		2,  8,  7,
		2,  7,  8,
		2,  6,  9,
		2,  5,  10,
		2,  4,  11,
		2,  3,  12,
		2,  2,  13,
		2,  1,  14,
		2,  0,  15,
		1,  16,  0,
		1,  15,  1,
		1,  14,  2,
		1,  13,  3,
		1,  12,  4,
		1,  11,  5,
		1,  10,  6,
		1,  9,  7,
		1,  8,  8,
		1,  7,  9,
		1,  6,  10,
		1,  5,  11,
		1,  4,  12,
		1,  3,  13,
		1,  2,  14,
		1,  1,  15,
		1,  0,  16,
		0,  17,  0,
		0,  16,  1,
		0,  15,  2,
		0,  14,  3,
		0,  13,  4,
		0,  12,  5,
		0,  11,  6,
		0,  10,  7,
		0,  9,  8,
		0,  8,  9,
		0,  7,  10,
		0,  6,  11,
		0,  5,  12,
		0,  4,  13,
		0,  3,  14,
		0,  2,  15,
		0,  1,  16,
		0,  0,  17,
		18,  0,  0,  /*L = 18*/
		17,  1,  0,
		17,  0,  1,
		16,  2,  0,
		16,  1,  1,
		16,  0,  2,
		15,  3,  0,
		15,  2,  1,
		15,  1,  2,
		15,  0,  3,
		14,  4,  0,
		14,  3,  1,
		14,  2,  2,
		14,  1,  3,
		14,  0,  4,
		13,  5,  0,
		13,  4,  1,
		13,  3,  2,
		13,  2,  3,
		13,  1,  4,
		13,  0,  5,
		12,  6,  0,
		12,  5,  1,
		12,  4,  2,
		12,  3,  3,
		12,  2,  4,
		12,  1,  5,
		12,  0,  6,
		11,  7,  0,
		11,  6,  1,
		11,  5,  2,
		11,  4,  3,
		11,  3,  4,
		11,  2,  5,
		11,  1,  6,
		11,  0,  7,
		10,  8,  0,
		10,  7,  1,
		10,  6,  2,
		10,  5,  3,
		10,  4,  4,
		10,  3,  5,
		10,  2,  6,
		10,  1,  7,
		10,  0,  8,
		9,  9,  0,
		9,  8,  1,
		9,  7,  2,
		9,  6,  3,
		9,  5,  4,
		9,  4,  5,
		9,  3,  6,
		9,  2,  7,
		9,  1,  8,
		9,  0,  9,
		8,  10,  0,
		8,  9,  1,
		8,  8,  2,
		8,  7,  3,
		8,  6,  4,
		8,  5,  5,
		8,  4,  6,
		8,  3,  7,
		8,  2,  8,
		8,  1,  9,
		8,  0,  10,
		7,  11,  0,
		7,  10,  1,
		7,  9,  2,
		7,  8,  3,
		7,  7,  4,
		7,  6,  5,
		7,  5,  6,
		7,  4,  7,
		7,  3,  8,
		7,  2,  9,
		7,  1,  10,
		7,  0,  11,
		6,  12,  0,
		6,  11,  1,
		6,  10,  2,
		6,  9,  3,
		6,  8,  4,
		6,  7,  5,
		6,  6,  6,
		6,  5,  7,
		6,  4,  8,
		6,  3,  9,
		6,  2,  10,
		6,  1,  11,
		6,  0,  12,
		5,  13,  0,
		5,  12,  1,
		5,  11,  2,
		5,  10,  3,
		5,  9,  4,
		5,  8,  5,
		5,  7,  6,
		5,  6,  7,
		5,  5,  8,
		5,  4,  9,
		5,  3,  10,
		5,  2,  11,
		5,  1,  12,
		5,  0,  13,
		4,  14,  0,
		4,  13,  1,
		4,  12,  2,
		4,  11,  3,
		4,  10,  4,
		4,  9,  5,
		4,  8,  6,
		4,  7,  7,
		4,  6,  8,
		4,  5,  9,
		4,  4,  10,
		4,  3,  11,
		4,  2,  12,
		4,  1,  13,
		4,  0,  14,
		3,  15,  0,
		3,  14,  1,
		3,  13,  2,
		3,  12,  3,
		3,  11,  4,
		3,  10,  5,
		3,  9,  6,
		3,  8,  7,
		3,  7,  8,
		3,  6,  9,
		3,  5,  10,
		3,  4,  11,
		3,  3,  12,
		3,  2,  13,
		3,  1,  14,
		3,  0,  15,
		2,  16,  0,
		2,  15,  1,
		2,  14,  2,
		2,  13,  3,
		2,  12,  4,
		2,  11,  5,
		2,  10,  6,
		2,  9,  7,
		2,  8,  8,
		2,  7,  9,
		2,  6,  10,
		2,  5,  11,
		2,  4,  12,
		2,  3,  13,
		2,  2,  14,
		2,  1,  15,
		2,  0,  16,
		1,  17,  0,
		1,  16,  1,
		1,  15,  2,
		1,  14,  3,
		1,  13,  4,
		1,  12,  5,
		1,  11,  6,
		1,  10,  7,
		1,  9,  8,
		1,  8,  9,
		1,  7,  10,
		1,  6,  11,
		1,  5,  12,
		1,  4,  13,
		1,  3,  14,
		1,  2,  15,
		1,  1,  16,
		1,  0,  17,
		0,  18,  0,
		0,  17,  1,
		0,  16,  2,
		0,  15,  3,
		0,  14,  4,
		0,  13,  5,
		0,  12,  6,
		0,  11,  7,
		0,  10,  8,
		0,  9,  9,
		0,  8,  10,
		0,  7,  11,
		0,  6,  12,
		0,  5,  13,
		0,  4,  14,
		0,  3,  15,
		0,  2,  16,
		0,  1,  17,
		0,  0,  18,
		19,  0,  0,  /*L = 19*/
		18,  1,  0,
		18,  0,  1,
		17,  2,  0,
		17,  1,  1,
		17,  0,  2,
		16,  3,  0,
		16,  2,  1,
		16,  1,  2,
		16,  0,  3,
		15,  4,  0,
		15,  3,  1,
		15,  2,  2,
		15,  1,  3,
		15,  0,  4,
		14,  5,  0,
		14,  4,  1,
		14,  3,  2,
		14,  2,  3,
		14,  1,  4,
		14,  0,  5,
		13,  6,  0,
		13,  5,  1,
		13,  4,  2,
		13,  3,  3,
		13,  2,  4,
		13,  1,  5,
		13,  0,  6,
		12,  7,  0,
		12,  6,  1,
		12,  5,  2,
		12,  4,  3,
		12,  3,  4,
		12,  2,  5,
		12,  1,  6,
		12,  0,  7,
		11,  8,  0,
		11,  7,  1,
		11,  6,  2,
		11,  5,  3,
		11,  4,  4,
		11,  3,  5,
		11,  2,  6,
		11,  1,  7,
		11,  0,  8,
		10,  9,  0,
		10,  8,  1,
		10,  7,  2,
		10,  6,  3,
		10,  5,  4,
		10,  4,  5,
		10,  3,  6,
		10,  2,  7,
		10,  1,  8,
		10,  0,  9,
		9,  10,  0,
		9,  9,  1,
		9,  8,  2,
		9,  7,  3,
		9,  6,  4,
		9,  5,  5,
		9,  4,  6,
		9,  3,  7,
		9,  2,  8,
		9,  1,  9,
		9,  0,  10,
		8,  11,  0,
		8,  10,  1,
		8,  9,  2,
		8,  8,  3,
		8,  7,  4,
		8,  6,  5,
		8,  5,  6,
		8,  4,  7,
		8,  3,  8,
		8,  2,  9,
		8,  1,  10,
		8,  0,  11,
		7,  12,  0,
		7,  11,  1,
		7,  10,  2,
		7,  9,  3,
		7,  8,  4,
		7,  7,  5,
		7,  6,  6,
		7,  5,  7,
		7,  4,  8,
		7,  3,  9,
		7,  2,  10,
		7,  1,  11,
		7,  0,  12,
		6,  13,  0,
		6,  12,  1,
		6,  11,  2,
		6,  10,  3,
		6,  9,  4,
		6,  8,  5,
		6,  7,  6,
		6,  6,  7,
		6,  5,  8,
		6,  4,  9,
		6,  3,  10,
		6,  2,  11,
		6,  1,  12,
		6,  0,  13,
		5,  14,  0,
		5,  13,  1,
		5,  12,  2,
		5,  11,  3,
		5,  10,  4,
		5,  9,  5,
		5,  8,  6,
		5,  7,  7,
		5,  6,  8,
		5,  5,  9,
		5,  4,  10,
		5,  3,  11,
		5,  2,  12,
		5,  1,  13,
		5,  0,  14,
		4,  15,  0,
		4,  14,  1,
		4,  13,  2,
		4,  12,  3,
		4,  11,  4,
		4,  10,  5,
		4,  9,  6,
		4,  8,  7,
		4,  7,  8,
		4,  6,  9,
		4,  5,  10,
		4,  4,  11,
		4,  3,  12,
		4,  2,  13,
		4,  1,  14,
		4,  0,  15,
		3,  16,  0,
		3,  15,  1,
		3,  14,  2,
		3,  13,  3,
		3,  12,  4,
		3,  11,  5,
		3,  10,  6,
		3,  9,  7,
		3,  8,  8,
		3,  7,  9,
		3,  6,  10,
		3,  5,  11,
		3,  4,  12,
		3,  3,  13,
		3,  2,  14,
		3,  1,  15,
		3,  0,  16,
		2,  17,  0,
		2,  16,  1,
		2,  15,  2,
		2,  14,  3,
		2,  13,  4,
		2,  12,  5,
		2,  11,  6,
		2,  10,  7,
		2,  9,  8,
		2,  8,  9,
		2,  7,  10,
		2,  6,  11,
		2,  5,  12,
		2,  4,  13,
		2,  3,  14,
		2,  2,  15,
		2,  1,  16,
		2,  0,  17,
		1,  18,  0,
		1,  17,  1,
		1,  16,  2,
		1,  15,  3,
		1,  14,  4,
		1,  13,  5,
		1,  12,  6,
		1,  11,  7,
		1,  10,  8,
		1,  9,  9,
		1,  8,  10,
		1,  7,  11,
		1,  6,  12,
		1,  5,  13,
		1,  4,  14,
		1,  3,  15,
		1,  2,  16,
		1,  1,  17,
		1,  0,  18,
		0,  19,  0,
		0,  18,  1,
		0,  17,  2,
		0,  16,  3,
		0,  15,  4,
		0,  14,  5,
		0,  13,  6,
		0,  12,  7,
		0,  11,  8,
		0,  10,  9,
		0,  9,  10,
		0,  8,  11,
		0,  7,  12,
		0,  6,  13,
		0,  5,  14,
		0,  4,  15,
		0,  3,  16,
		0,  2,  17,
		0,  1,  18,
		0,  0,  19,
		20,  0,  0,  /*L = 20*/
		19,  1,  0,
		19,  0,  1,
		18,  2,  0,
		18,  1,  1,
		18,  0,  2,
		17,  3,  0,
		17,  2,  1,
		17,  1,  2,
		17,  0,  3,
		16,  4,  0,
		16,  3,  1,
		16,  2,  2,
		16,  1,  3,
		16,  0,  4,
		15,  5,  0,
		15,  4,  1,
		15,  3,  2,
		15,  2,  3,
		15,  1,  4,
		15,  0,  5,
		14,  6,  0,
		14,  5,  1,
		14,  4,  2,
		14,  3,  3,
		14,  2,  4,
		14,  1,  5,
		14,  0,  6,
		13,  7,  0,
		13,  6,  1,
		13,  5,  2,
		13,  4,  3,
		13,  3,  4,
		13,  2,  5,
		13,  1,  6,
		13,  0,  7,
		12,  8,  0,
		12,  7,  1,
		12,  6,  2,
		12,  5,  3,
		12,  4,  4,
		12,  3,  5,
		12,  2,  6,
		12,  1,  7,
		12,  0,  8,
		11,  9,  0,
		11,  8,  1,
		11,  7,  2,
		11,  6,  3,
		11,  5,  4,
		11,  4,  5,
		11,  3,  6,
		11,  2,  7,
		11,  1,  8,
		11,  0,  9,
		10,  10,  0,
		10,  9,  1,
		10,  8,  2,
		10,  7,  3,
		10,  6,  4,
		10,  5,  5,
		10,  4,  6,
		10,  3,  7,
		10,  2,  8,
		10,  1,  9,
		10,  0,  10,
		9,  11,  0,
		9,  10,  1,
		9,  9,  2,
		9,  8,  3,
		9,  7,  4,
		9,  6,  5,
		9,  5,  6,
		9,  4,  7,
		9,  3,  8,
		9,  2,  9,
		9,  1,  10,
		9,  0,  11,
		8,  12,  0,
		8,  11,  1,
		8,  10,  2,
		8,  9,  3,
		8,  8,  4,
		8,  7,  5,
		8,  6,  6,
		8,  5,  7,
		8,  4,  8,
		8,  3,  9,
		8,  2,  10,
		8,  1,  11,
		8,  0,  12,
		7,  13,  0,
		7,  12,  1,
		7,  11,  2,
		7,  10,  3,
		7,  9,  4,
		7,  8,  5,
		7,  7,  6,
		7,  6,  7,
		7,  5,  8,
		7,  4,  9,
		7,  3,  10,
		7,  2,  11,
		7,  1,  12,
		7,  0,  13,
		6,  14,  0,
		6,  13,  1,
		6,  12,  2,
		6,  11,  3,
		6,  10,  4,
		6,  9,  5,
		6,  8,  6,
		6,  7,  7,
		6,  6,  8,
		6,  5,  9,
		6,  4,  10,
		6,  3,  11,
		6,  2,  12,
		6,  1,  13,
		6,  0,  14,
		5,  15,  0,
		5,  14,  1,
		5,  13,  2,
		5,  12,  3,
		5,  11,  4,
		5,  10,  5,
		5,  9,  6,
		5,  8,  7,
		5,  7,  8,
		5,  6,  9,
		5,  5,  10,
		5,  4,  11,
		5,  3,  12,
		5,  2,  13,
		5,  1,  14,
		5,  0,  15,
		4,  16,  0,
		4,  15,  1,
		4,  14,  2,
		4,  13,  3,
		4,  12,  4,
		4,  11,  5,
		4,  10,  6,
		4,  9,  7,
		4,  8,  8,
		4,  7,  9,
		4,  6,  10,
		4,  5,  11,
		4,  4,  12,
		4,  3,  13,
		4,  2,  14,
		4,  1,  15,
		4,  0,  16,
		3,  17,  0,
		3,  16,  1,
		3,  15,  2,
		3,  14,  3,
		3,  13,  4,
		3,  12,  5,
		3,  11,  6,
		3,  10,  7,
		3,  9,  8,
		3,  8,  9,
		3,  7,  10,
		3,  6,  11,
		3,  5,  12,
		3,  4,  13,
		3,  3,  14,
		3,  2,  15,
		3,  1,  16,
		3,  0,  17,
		2,  18,  0,
		2,  17,  1,
		2,  16,  2,
		2,  15,  3,
		2,  14,  4,
		2,  13,  5,
		2,  12,  6,
		2,  11,  7,
		2,  10,  8,
		2,  9,  9,
		2,  8,  10,
		2,  7,  11,
		2,  6,  12,
		2,  5,  13,
		2,  4,  14,
		2,  3,  15,
		2,  2,  16,
		2,  1,  17,
		2,  0,  18,
		1,  19,  0,
		1,  18,  1,
		1,  17,  2,
		1,  16,  3,
		1,  15,  4,
		1,  14,  5,
		1,  13,  6,
		1,  12,  7,
		1,  11,  8,
		1,  10,  9,
		1,  9,  10,
		1,  8,  11,
		1,  7,  12,
		1,  6,  13,
		1,  5,  14,
		1,  4,  15,
		1,  3,  16,
		1,  2,  17,
		1,  1,  18,
		1,  0,  19,
		0,  20,  0,
		0,  19,  1,
		0,  18,  2,
		0,  17,  3,
		0,  16,  4,
		0,  15,  5,
		0,  14,  6,
		0,  13,  7,
		0,  12,  8,
		0,  11,  9,
		0,  10,  10,
		0,  9,  11,
		0,  8,  12,
		0,  7,  13,
		0,  6,  14,
		0,  5,  15,
		0,  4,  16,
		0,  3,  17,
		0,  2,  18,
		0,  1,  19,
		0,  0,  20
	};

	/**
	 * this is a utility class used to build the mapping between index
	 * and basis set for a given basis set order type
	 */
	class BasisUtil {

		public:

			BasisUtil() { };
			~BasisUtil() { };

			// get the l,m and n value for the given local index for basis set
			// in the shell L
			void getLMNVal(const int& L, const int& index, int& l, int& m, int& n) const;

			// give the local basis set index by given the l,m and n value
			int getLocalBasisSetIndex(const int& l, const int& m, const int& n) const;

	};

}

#endif

