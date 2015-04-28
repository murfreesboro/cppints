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
 * \file    inttype.h
 * \brief   describing the integral type in quantum chemistry
 * \author  Fenglai Liu
 * \note
 *
 * This module is used define the integral operator, and through
 * a bunch of functions we are able to handle the operator information
 */
#ifndef INTTYPE_H
#define INTTYPE_H
#include "general.h"

namespace inttype {

	///
	/// according to the name of integral, we give the order of 
	/// integrals
	/// order = 1 -> only bra1 presents
	/// order = 2 -> only bra1, bra2 present
	/// order = 3 -> bra1, bra2 and ket1 present
	/// order = 4 -> bra1, bra2, ket1 and ket2 present
	///
	int getOperOrder(const int& name);

	///
	/// get the order used in RR expansion
	/// usually this order is equal to the above getOperOrder
	/// however, in some cases it's not.
	/// For example, in the MOM integral, even though it's 
	/// two body integral, but the ket1 part (operator infor)
	/// is also appearing in the RR expansion. Therefore,
	/// for MOM it's RR order is 3
	///
	int getRROrder(const int& name);

	/**
	 * whether the integral has ket side?
	 */
	bool hasKetSide(const int& name);

	/**
	 * get the operator string ID
	 */
	string getOperStringName(const int& name);

	/**
	 * get the operator int ID
	 */
	int getOperIntName(const string& name);

	/**
	 * compare the two given operators
	 */
	int compareOper(const int& O1, const int& O2);

	/**
	 * for the given operator O, we try to get a list
	 * of operators that appears in the RR expression 
	 * characterized by O.
	 */
	void selectOperListInRR(const int& O, vector<int>& oList);

	/**
	 * for the given integral type, whether we have M value really
	 * defined in RR
	 */
	bool hasMValueDefined(const int& O); 

	/**
	 * whether the operator will generate the fmt function for 
	 * bottom integrals
	 */
	bool useFmt(const int& O);

	/**
	 * do we use significance check in the code?
	 */
	bool sigCheck(const int& O);

	/**
	 * this function return the type of RR property for the 
	 * given operator
	 * the RR property determines the property in RR process
	 * for example, NAI and ERI use both L and M value in
	 * forming the RR expression, therefore we return
	 * SQ_ON_L_M
	 * KINETIC use L and operator (two body overlap integrals
	 * appearing in RHS of its RR), therefore it's 
	 * SQ_ON_L_OPER
	 *
	 * the return value is defined in general.h
	 */
	int rrProp(const int& oper);

	/**
	 * for some operator, we can not do HRR process
	 */
	bool canDOHRR(const int& O); 

	/**
	 * get the number of space by operator
	 */
	int getNSpaceByOper(const int& oper);

	/**
	 * whether the operator itself is null
	 */
	bool isNullOper(const int& oper);

	/**
	 * for the RR/derivatives expression, whether the RHS has different
	 * operator other than the input one?
	 */
	bool intOperChanged(const int& O);

	/**
	 * for the given operator, whethe it's RR demands some 
	 * complicated bottom integrals?
	 */
	bool complicatedBottomIntegrals(const int& oper);

	/**
	 * does the oper involves direct integral derivatives 
	 * calculation?
	 */
	bool needIntDeriv(const int& oper);

	/**
	 * get the deriv number from the given operator as well as job order
	 */
	int getDerivNumberFromOper(const int& oper, const int& jobOrder);

	/**
	 * get the deriv order from the given operator as well as job order
	 */
	int getDerivOrderFromOper(const int& oper, const int& jobOrder);

	/**
	 * does the result integrals has additional offset?
	 * actually it means that besides the number of basis set dimension,
	 * whether the result integrals has other dimensions? for example,
	 * ESP has an additional dimension for number of grid points
	 */
	bool resultIntegralHasAdditionalOffset(const int& oper);

	/**
	 * if the result integrals has complicated offset, then here
	 * we will determine it and return it in string format
	 * 
	 * We note, that the additional dimension should be after the 
	 * integrals. For example, ESP result is in dimension of 
	 * (nBas,nBas,NGrids). Therefore the additional dimension of 
	 * nGrids is last dimension. that's the reason why we pass
	 * the nInts.
	 *
	 * \param  nInts number of total integrals
	 */
	string determineAdditionalOffset(const int& oper, const int& nInts);
}

#endif
