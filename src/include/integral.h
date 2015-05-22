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

		public:

			/**
			 * this constructor builds the integral directly with 
			 * the given basis sets
			 */
			Integral(const Basis& oribra1, const Basis& oribra2,
					const Basis& oriket1, const Basis& oriket2, int oriO,
					int m = 0, long long div = NULL_POS):bra1(oribra1),bra2(oribra2),
			ket1(oriket1),ket2(oriket2),O(oriO),mvalue(m),division(div) {
				crash(mvalue<0, "In Integral constructor m value is < 0!");
			};

			/**
			 * this is the default construcor, used to build the 
			 * NULL integral
			 */
			Integral( ):O(NULL_POS),mvalue(NULL_POS),division(NULL_POS){
				Basis null(NULL_POS,NULL_POS,NULL_POS);
				bra1 = null;
				bra2 = null;
				ket1 = null;
				ket2 = null;
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
			 * here we can further form it's variable name in the final codes
			 * such forming will depending on the following information:
			 * - rrType      : HRR or VRR? (VRR matters)
			 * - status      : whether it's tmp result, module result, final result?
			 * - pos         : the RR expanding position
			 * - withModifier: this only works for VRR module result. Since the VRR
			 *                 module result is also the input for HRR, therefore;
			 *                 whether the corresponding VRR result sq is with modifiier?
			 *                 for example, division (composite shell quartet), and/or
			 *                 the exponents?
			 */
			string formVarName(const int& rrType, const int& status, 
					const int& pos, bool withModifier = false) const;

			/**
			 * through the basis set information, we are able to 
			 * locate the index for this integral in the corresponding
			 * shell quartet - this location we note is unique
			 */
			int getIndex() const;
	};

}

#endif

