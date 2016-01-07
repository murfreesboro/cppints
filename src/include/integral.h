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
 * \file    integral.h
 * \brief   describing the integral
 * \author  Fenglai Liu
 * \note
 *
 * This module is used describe the class related to the "integral".
 * The integral is just the explicit electron integrals based on the 
 * basis set functions. For example; the ERI(two electron 
 * repulsion integral), overlap integral etc. For each integral,
 * we have bra1, bra2, ket1 and ket2 as its input basis set functions.
 *
 * There's an additional comment about m value used in the integral,
 * this is related to the auxiliary integral used in OS framework. 
 * We note that M = 0 is the true integral.
 * 
 * See the equation 38 for description of m value in the reference below:
 * Obara, S. and Saika, A.
 * Efficient recursive computation of molecular integrals over Cartesian 
 * Gaussian functions
 * The Journal of chemical physics
 * 1986,84,3963
 */
#ifndef INTEGRAL_H
#define INTEGRAL_H
#include "general.h"
#include "basis.h"
using namespace basis;

namespace shellquartet {
	class ShellQuartet;
}

namespace integral {

	using namespace shellquartet;

	class Integral {

		private:

			Basis bra1;   ///< bra1 basis
			Basis bra2;   ///< bra2 basis
			Basis ket1;   ///< ket1 basis
			Basis ket2;   ///< ket2 basis
			int      O;   ///< integral operator
			int mvalue;   ///< M value needed in the ERI integral etc.
			long long division; ///< this is associated with shell quartets
			int expFacListLen;  ///< the real length of expFacList
			int expFacList[MAX_EXP_FAC_LIST];  ///< possible exponential factor multiply with integral

			// now it's the derivative information
			int firstDerivPos;    ///< the first derivative position 
			int secondDerivPos;   ///< the second derivative position
			int firstDerivDir;    ///< the first derivative direction
			int secondDerivDir;   ///< the second derivative direction

		public:

			/**
			 * this constructor builds the integral directly with 
			 * the given basis sets
			 *
			 * we note, that in default we do not have any derivatives
			 * information, these will be added through function
			 * addDerivInfor
			 */
			Integral(const Basis& oribra1, const Basis& oribra2,
					const Basis& oriket1, const Basis& oriket2, int oriO,
					int m = 0, long long div = NULL_POS):bra1(oribra1),bra2(oribra2),
			ket1(oriket1),ket2(oriket2),O(oriO),mvalue(m),division(div),expFacListLen(0) {
				crash(mvalue<0, "In Integral constructor m value is < 0!");
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = NULL_POS;
				firstDerivPos = NULL_POS;  
				secondDerivPos = NULL_POS; 
				firstDerivDir  = NO_DERIV;   
				secondDerivDir = NO_DERIV;  
			};

			/**
			 * this is the default construcor, used to build the 
			 * NULL integral
			 */
			Integral( ):O(NULL_POS),mvalue(NULL_POS),division(NULL_POS),
			expFacListLen(0) {
				Basis null(NULL_POS,NULL_POS,NULL_POS);
				bra1 = null;
				bra2 = null;
				ket1 = null;
				ket2 = null;
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = NULL_POS;
				firstDerivPos = NULL_POS;  
				secondDerivPos = NULL_POS; 
				firstDerivDir  = NO_DERIV;   
				secondDerivDir = NO_DERIV;  
			};

			/**
			 * constructor builds integral through shell quartet and 
			 * index - therefore, this integral is in the given position
			 * of the integral list for the given shell quartet
			 */
			Integral(const ShellQuartet& sq, const int& index);

			/**
			 * destructor
			 */
			~Integral() { };

			/*
			 *operator 
			 */
			int getOper() const { return O; };

			/**
			 * get M value
			 */
			int getMValue() const {return mvalue;};

			/**
			 * get division
			 */
			long long getDivision() const {return division;};

			/**
			 * remove the division
			 */
			void destroyDivision()  {division = NULL_POS;};

			/**
			 * remove the exponetial factor modifier
			 */
			void destroyExpFac()  {
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = NULL_POS;
				expFacListLen = 0;
			};

			/**
			 * destroy all of the multipilers information
			 * such as exponential factors etc.
			 */
			void destroyMultipliers() {
				destroyDivision();
				destroyExpFac(); 
			};

			/**
			 * destroy the derivative information
			 */
			void destroyDerivInfor() {
				firstDerivPos = NULL_POS;  
				secondDerivPos = NULL_POS; 
				firstDerivDir  = NO_DERIV;   
				secondDerivDir = NO_DERIV;  
			};

			/**
			 * get the 1st order deriv pos
			 */
			int get1stDerivPos() const { return firstDerivPos; };  

			/**
			 * get the 2ed order deriv pos
			 */
			int get2edDerivPos() const { return secondDerivPos; }; 

			/**
			 * get the 1st order deriv direction, x, y or z?
			 */
			int get1stDerivDir() const { return firstDerivDir; };   

			/**
			 * get the 2ed order deriv direction, x, y or z?
			 */
			int get2edDerivDir() const { return secondDerivDir; };  

			/**
			 * add deriv information for the given position and direction
			 */
			void addDerivInfor(int derivPos, int derivDir) {
				firstDerivPos = derivPos;
				firstDerivDir = derivDir;
			};

			/**
			 * add deriv information for the given positions and directions
			 */
			void addDerivInfor(int pos1, int pos2, int dir1, int dir2) {
				firstDerivPos  = pos1;
				firstDerivDir  = dir1;
				secondDerivPos = pos2;
				secondDerivDir = dir2;
			};

			/**
			 * get the derivatives job order
			 */
			int getDerivJobOrder() const {
				int jobOrder = 0;
				if (firstDerivPos != NULL_POS) jobOrder++;
				if (secondDerivPos != NULL_POS) jobOrder++;
				return jobOrder;
			};

			/**
			 * whether it's null integral?
			 * we note, that bra1 will never be null if 
			 * the integral is meaningful
			 */
			bool isnull() const {
				if (bra1.isnull()) return true;
				return false;
			};

			/**
			 * get the basis set
			 */
			const Basis& getBasis(const int& pos) const;

			// eq 
			bool operator==(const Integral& I) const;

			/**
			 * get the name for the shell quartet
			 */
			string getName() const;

			/**
			 * based on the integral name (see above),
			 * here we can further form it's variable name for RR in the final codes
			 */
			string formVarName(const int& rrType) const;

			/**
			 * through the basis set information, we are able to 
			 * locate the index for this integral in the corresponding
			 * shell quartet - this location we note is unique
			 */
			int getIndex() const;

			/**
			 * adding the expotential factor to the integral
			 * here for the expotential factors like alpha, beta 
			 * etc. multiplied with integrals, the most important
			 * thing is that they are exchangable.
			 *
			 * For example, two expFacList array:
			 * expFacList1 = [ BRA1(alpha) KET1(gamma) -1 -1]
			 * expFacList2 = [ KET1(gamma) BRA1(alpha) -1 -1]
			 * they are actually same. Therefore, in adding the 
			 * exponential factor, we will sort the list and 
			 * make it arrange in the order that:
			 *
			 * bra1 -> bra2 -> ket1 -> ket2
			 */
			void addExpFac(int pos); 

			/**
			 * copy the exponential factor information from the 
			 * input shell quartet to this one
			 */
			void copyExpFac(const Integral& I) {
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = I.expFacList[i];
				expFacListLen = I.expFacListLen;
			};

			/**
			 * get the given order of exp factor value
			 */
			int getExpFacVal(int order) const {
				if (order<0 || order>=expFacListLen) {
					crash(true,"Invalid order pass into the getExpVal function in Integral class");
				}
				return expFacList[order];
			};

			/**
			 * return the array of exp factor list
			 */
			const int* getExpFacList() const { return expFacList; };

			/**
			 * get the len of exp factor list
			 */
			int getExpFacListLen() const { 
				return expFacListLen;
			};

			/**
			 * whether the integral with exponential factors?
			 */
			bool withExpFac() const {
				if (expFacListLen == 0) return false;
				return true;
			};

			/**
			 * whether the integral with division information defined?
			 */
			bool withDivision() const {
				if(division == NULL_POS) return false;
				return true;
			};

			/**
			 * whether the integral is with modifier information?
			 */
			bool withModifier() const {
				if (withDivision()) return true;
				if (withExpFac()) return true;
				return false;
			};

			/**
			 * whether the given integral is S type integral(bottom integral)?
			 * it seems this function is not used, so comment it out 
			 */
			//bool isSTypeIntegral() const;
	};

}

#endif

