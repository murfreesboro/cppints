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
 * \file    shell.h
 * \brief   describing the shell quartet
 * \author  Fenglai Liu
 * \note
 *
 * This module is used describe the class related to the "Shell quartet". 
 * Shell quartet is a group of integrals based on shells.
 * For example; the ERI(two electron repulsion integral) is: (l1l2|1/r12|l3l4)
 * each l representing a "Shell". 
 *
 * There's an additional comment about m value used in the integral evaluation,
 * this is related to the auxiliary integral used in OS framework. 
 * We note that M = 0 is the true integral.
 *
 * See the equation 38 for description of m value in the reference below:
 * Obara, S. and Saika, A.
 * Efficient recursive computation of molecular integrals over Cartesian 
 * Gaussian functions
 * The Journal of chemical physics
 * 1986,84,3963
 *
 * important note about shell quartet:
 *
 * 1  in the code we have introduce an important assumption, that if the operator
 * is introduced as a "shell component" (for example, in the momentum integrals,
 * the operation is introduced in the ket1 position as a "shell". See the code in
 * the codegen.cpp of MOM part), then this special shell component is always "APPENDING"
 * to the end of shells. For example, im MOM integrals bra1 and bra2 is defined;
 * and ket1 is the operator. After ket1, ket2 is null shell. It seems that in 
 * out code we stick to this rule. So be careful in the future!
 *
 * 2 in the hrr less than operation, the L(bra1) >= L(bra2) and L(ket1) >= L(ket2)
 * rules has been used. Therefore in the HRR process, you need to stick to this
 * to generate correct ordering of shell quartets. See the comment inside
 *
 */
#ifndef SHELLQUARTET_H
#define SHELLQUARTET_H
#include "shell.h"
#include "general.h"
using namespace shell;

namespace integral {
	class Integral;
}
namespace shellquartet {

	using namespace integral;

	class ShellQuartet {

		private:

			Shell bra1;   ///< bra1 shell
			Shell bra2;   ///< bra2 shell
			Shell ket1;   ///< ket1 shell
			Shell ket2;   ///< ket2 shell
			int      O;   ///< integral operator
			int mvalue;   ///< M value needed in the ERI integral etc.
			long long division; ///< a identifier for this single sq with respect 
			                    ///< to its oringinal composite sq    

		public:

			///
			/// constructors and destructors
			/// here the passed in m value could have several meanings:
			/// first, if the integral do not use m value, we just use the default 0;
			/// second, for the integrals that has M defined but m value is not set;
			/// this indicates the integral has M value of 0
			/// third, if m value passed in is some integer and larger than zero
			/// (could also =0), then the value passed in is the M value generated
			/// from recursive process
			///
			ShellQuartet(const Shell& oribra1, const Shell& oribra2,
					const Shell& oriket1, const Shell& oriket2, int oriO, 
					int m = 0, long long div = NULL_POS):bra1(oribra1),bra2(oribra2),
			ket1(oriket1),ket2(oriket2),O(oriO),mvalue(m),division(div) {
#ifdef DEBUG
				crash(mvalue<0, "In ShellQuartet constructor m value is < 0!");
#endif
			};

			///
			/// this is the default constructor to build the null shell quartet
			///
			ShellQuartet():O(NULL_POS),mvalue(NULL_POS),division(-1) {
				Shell nullshell;
				bra1 = nullshell;
				bra2 = nullshell;
				ket1 = nullshell;
				ket2 = nullshell;
			};

			///
			/// destructor
			///
			~ShellQuartet() { };

			/**
			 * whether this shell quartet is null?
			 */
			bool isnull() const;

			/**
			 * copy constructor
			 */
			ShellQuartet(const ShellQuartet& sq):bra1(sq.bra1),bra2(sq.bra2), 
			ket1(sq.ket1),ket2(sq.ket2),O(sq.O),mvalue(sq.mvalue),division(sq.division){ };

			/**
			 * get the division
			 */
			long long getDivision() const { return division; };

			/**
			 * get the maximum L from all of shell components
			 */
			//int getMaxL() const;

			/**
			 * get the number of full integrals generated from this shell quartet
			 */
			int getNInts() const;

			/**
			 * get the number of integrals for the given side
			 */
			int getNInts(const int& side) const;

			/**
			 * whether the given shell quartet matches the given L and M combination?
			 */
			bool matchLM(const int& l, const int& m) const;

			/**
			 * whether the given shell quartet matches the given L (consider LSum)?
			 */
			bool matchL(const int& l) const; 

			/**
			 * whether the given shell quartet matches the given L and oper combination?
			 */
			bool matchLOper(const int& l, const int& oper) const;

			/**
			 * get the L sum form all shells
			 */
			int getLSum() const;

			/**
			 * get the L sum for the given side, bra or ket
			 */
			int getLSum(const int& side) const;

			/**
			 * get the number of shells in the given shell quartet
			 * these shells are non-null shells
			 * not used
			 */
			//int getNShell() const;

			/**
			 * get the number of non-S type of shells in this shell quartet
			 */
			int getNNonSShell() const;

			/**
			 * get the first non-S shell position
			 * return BRA1, BRA2 etc.
			 */
			int getNonSShellPos() const;

			/**
			 * whether the given position corresponds to non-S type of shell
			 */
			bool isNonSShell(const int& pos) const;

			/**
			 *less than relation used in the HRR step
			 */
			bool hrrLessThan(const ShellQuartet& sq, const int& side) const;

			/**
			 * core algorithm to define the less than relation 
			 * for HRR 
			 */
			int hrrCompare(const ShellQuartet& sq, bool isBraSide) const;

			/**
			 *less than relation used to identify the shell quartet 
			 *position in the RR process
			 */
			bool operator<(const ShellQuartet& sq) const;

			/*
			 * operator information 
			 */
			int getOper() const { return O;};

			/**
			 * get M value
			 */
			int getM() const {return mvalue;};

			/**
			 * get the shell 
			 */
			const Shell& getShell(const int& pos) const;

			// eq 
			bool operator==(const ShellQuartet& sq) const;

			/**
			 * get the name for the shell quartet
			 */
			string getName() const;

			/**
			 * based on the shell quartet name (see above),
			 * here we can further form it's array name in the final codes
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
			string formArrayName(const int& rrType, const int& status, 
					const int& pos, bool withModifier = false) const;

			/**
			 * get the full integral list for this shell quartet
			 */
			void getIntegralList(set<int>& list) const;

			/**
			 * test that whether this is pure S type of shell quartet
			 */
			bool isSTypeSQ() const;

			/**
			 * whether this shell quartet has the given integral
			 */
			bool hasThisIntegral(const Integral& I) const;

			/**
			 * determine the bra/ket side in HRR
			 * for the HRR work, it's either working on the bra side
			 * or on the ket side. Sometimes, for the given shell quartet
			 * we even do not have any side(for example, (DS|DS)). therefore;
			 * we need to determine it before all of real work
			 */
			bool canDoHRR(const int& side) const;

			/**
			 * division information is only necessary for HRR process
			 * for VRR, it's better to destroy it before in use
			 * therefore we provide an function here to destroy the 
			 * shell quartet's division information
			 */
			void destroyDivision() { division = NULL_POS; };

			/**
			 * sometimes we can determine the RR expandable position from
			 * the shell quartet itself. In this function we will see 
			 * whether it's doable
			 */
			bool canDODirectRRPosSearch() const; 

			/**
			 * if the RR position direct search is doable, then in this
			 * function we will return this position
			 * only bra1, bra2, ket1 and ket2 value are returned
			 */
			int getRRPos() const; 
	};

}

#endif

