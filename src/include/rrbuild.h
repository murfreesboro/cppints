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
 * \file    rrbuild.h
 * \brief   build the basic recursive relation 
 * \author  Fenglai Liu
 * \note
 *
 * This module is used to build the general recursive relation,
 * which is called as "RR". The RR could be "VRR"(vertical type)
 * or the HRR(horizontal type).
 *
 * There are a lot of symbols in this module, they are divided 
 * into two kinds:
 *
 * the first group of symbol is defined in the inttype.h. They 
 * are related to the integral types. On the other hand, we use 
 * bra1, bra2, ket1, ket2 name to define the position in the 
 * integral/shell quartet. They are everywhere in this program, too.
 *
 * the second group of symbol is the names of constants or variables 
 * in the real calculation. For the meanings of the symbols, please 
 * refer to the comments on the top of source code.
 *
 * For each general recursive relation(it's named as general expression
 * here), it's composed by many items. For each item, it's composed by
 * "coeffieicnt" and the potential "integral/shell quartets" which 
 * could be referred from the paper below. 
 *
 * For each potential "integral/shell quartet", here we only record
 * its changing part(in comparison with the integral/shell quartet
 * on the LHS). For example, bra1 may arise 1 on X, and ket2 may 
 * lower down 1 on Z. We call the changing part as "variable" in
 * the codes. The variable is referring to whether the changes is 
 * on bra1, bra2 etc.
 *
 * See the reference below to get the original idea about VRR(OS type):
 * Obara, S. and Saika, A.
 * Efficient recursive computation of molecular integrals over Cartesian 
 * Gaussian functions
 * The Journal of chemical physics
 * 1986,84,3963
 *
 * By the way, besides the OS type of VRR, there could be other type of
 * VRR used; see the review for more information:
 * Gill, P.M.W.
 * Molecular integrals over Gaussian basis functions
 * Advances in quantum chemistry
 * 1994,25,141
 *
 * For the original idea about HRR, see the reference below:
 * Head-Gordon, M. and Pople, J.A.
 * A method for two-electron Gaussian integral and integral derivative
 * evaluation using recurrence relations
 * The Journal of chemical physics
 * 1988,89,5777
 *
 * important note for RRBuild:
 *
 * 1 in the RR expression, now we assume that each term(shell quartet 
 * or integral) in the RHS has maximum two changing positions in terms
 * of bra1, bra2, ket1 and ket2. This is used in buildRRSQ and buildRRInt.
 * See the code part generating v1/ang1 and v2/ang2. You need to remember this.
 *
 */
#ifndef RRBUILD_H
#define RRBUILD_H
#include "general.h"

namespace integral {
	class Integral;
};

namespace basis {
	class Basis;
};

namespace shellquartet {
	class ShellQuartet;
};

namespace rrbuild {

	using namespace basis;
	using namespace integral;
	using namespace shellquartet;


	//////////////////////////////////////////////////////////////////////////////
	//                                                                          //
	//  first section: constant data and utility functions shared by classes    //
	//                                                                          //
	//////////////////////////////////////////////////////////////////////////////

	// uni is used to characterize the non-determined components on the x, y or z
	// axis, for example; PAi could be PAx, PAy or PAz 
	const string uni = "uncertainI";

	// unang is used to characterize the non-determined angular momentum numbers(Ni)
	// defined in equation 6 of OS paper
	const string unangA = "uncertainNiA";
	const string unangB = "uncertainNiB";
	const string unangC = "uncertainNiC";
	const string unangD = "uncertainNiD";

	// unM is used to label the un-determined M value in the codes
	const string unM = "uncertainM";

	//////////////////////////////////////////////////////////////////////////////
	//                                                                          //
	//                     second section: RRBuild class                        //
	//                                                                          //
	//////////////////////////////////////////////////////////////////////////////

	/**
	 * \class RRBuild
	 *
	 * This class has two folds of purpose. The first, that it is used to build 
	 * the general RR expressions, for the given RR position (bra1, ket1 etc.)
	 * it's able to return the the RHS shell quartet list for an input LHS
	 * shell quartet or for a given integral it returns the RHS coefficients 
	 * as well as RHS integrals. All of these functions are used in VRR/HRR
	 * process.
	 *
	 * On the other hand, RRbuild also performs complicated contraction work
	 * between the input shell quartets and RR shell quartets. In the integral
	 * forming process, usually the input shell quartet list (which is also the
	 * final results) could be obtained directly from its RR work. For example,
	 * in ERI/NAI etc. integrals (job order == 0), we do RR; then simply adding
	 * the RR sq into the result sq who share the same name; like this:
	 * SQ_ERI_P_S_P_S[0] += SQ_ERI_P_S_P_S_vrr[0]; 
	 * SQ_ERI_P_S_P_S[1] += SQ_ERI_P_S_P_S_vrr[1]; 
	 * SQ_ERI_P_S_P_S[2] += SQ_ERI_P_S_P_S_vrr[2]; etc.
	 * this is called ``simple contraction''.
	 *
	 * However, the input sq list may not be the direct input for RR work. For example,
	 * the three body kinetic integrals, or the situations that jobOrder>0; we need
	 * to lower/raise the angular momentum in the input shell quartet so that to 
	 * perform the derivatives calculation. Therefore it's necessary to have a class
	 * to build the rr sq list from the input sq, this is what this class really wants
	 * to do.
	 *
	 * additionally, if we have the rrsqlist != inputsqlist, then it implies that 
	 * we will has a linear combination between the rr sq so that to make the final
	 * input sq. Usually we need another thing, that is the coefficients. Therefore
	 * we also have a coefficient array for the contraction work.
	 *
	 * Compared with the simple contraction, this is the ``complicated contraction''
	 * step. We use RRBuild to perform this function as well. The difference between
	 * RR building work is that the position is different. In contraction work,
	 * position will be x, y, z etc. defined in general.h
	 *
	 */
	class RRBuild {

		private:

			// input data
			int rrType;            ///< RR type
			int oper;              ///< integral type
			int jobOrder;          ///< job order: 0 energy, 1 2 etc. drivatives calculation
			int position;          ///< in which position we build the general RR
			                       ///< or in which derivatives position we form contraction

			// data of the general RR/contraction expression
			vector<string> coes;   ///< coefficients for each item in the general expression
			vector<int> lenVars;   ///< used to locate the length of variable in each item 
			vector<int> intTypes;  ///< for each item, it's integral type (only for changed)
			vector<int> angs;      ///< angular momentum changes on variable for each item
			vector<int> vars;      ///< variable changes on variable for each item(bra1 etc.)
			vector<int> mVals;     ///< M value for the integral/shell quartets

			/**
			 * build the general vertical RR expression
			 */
			void buildGeneralVRR(); 

			/**
			 * build the general HRR relation
			 */
			void buildGeneralHRR();

			/**
			 * build the general derivatives expression
			 */
			void buildDerivExpression(); 

			/**
			 * checking compatibility for the given input data
			 */
			bool compatibilityCheck() const;

			/**
			 * by giving a basis set, we determine the 
			 * direction to lower down it's angular momentum
			 * on x, y or z?
			 */
			string determineDirection(const Basis& bas) const;

			/**
			 * derive the basis set by raising up/lowering down the angular momentum
			 * i is the direction to lower down/raising up angular momentum
			 * inc is the incremental value for angular momentum reduction
			 */
			Basis getChangeBasisSet(const Basis& bas, const string& i, const int& inc) const; 

			/**
			 * get the Ni value in the recursive relation
			 */
			string Ni(const Basis& bas,const string& state) const;

			/**
			 * according to the direction and position given, 
			 * for a given integral we could evaluate
			 * the Ni for A,B,C and D
			 * i is the direction to lower down/raising up 
			 * angular momentum
			 */
			void determineNi(const Integral& I, const string& i,
					string& NiA, string& NiB, string& NiC, string& NiD) const;

			/**
			 * whether this is RR work?
			 */
			bool isRRWork() const {
				if (position == BRA1 || position == BRA2 ||
						position == KET1 || position == KET2 ) return true;
				return false;
			};

			/**
			 * whether this is VRR work?
			 */
			bool isVRRWork() const {
				if (isRRWork() && rrType != HRR ) return true;
				return false;
			};

			/**
			 * whether this is HRR work?
			 */
			bool isHRRWork() const {
				if (isRRWork() && rrType == HRR ) return true;
				return false;
			};

			/**
			 * return the length of items in the VRR general expression
			 * We note that they are fixed number. See the paper above
			 * we note, that this is for normal VRR forming
			 */
			int getVRRLenItems() const;

			/**
			 * return the length of items in the HRR general expression
			 * We note that this number is fixed for all situations
			 * we note, that this is for normal HRR forming
			 */
			int getHRRLenItems() const {return 2; };

			/**
			 * if the input position is derivatives, then we will 
			 * search its length of items defined
			 */
			int getDerivLenItems() const;

		public:

			/**
			 * constructor 
			 * \param rrType   RR working algorithm
			 * \param oper     operator
			 * \param pos      position information
			 * \param jobOrder job order, 0 is energy; 1 is gradient, 2 is 2ed derivatives etc.
			 *
			 * position could be two kinds of information. One is when you do RR, then
			 * position is BRA1, BRA2, KET1 and KET2. On the other hand, if you need to
			 * do derivatives, then position could be DERIV_X etc.
			 *
			 */
			RRBuild(const int& rrType0, const int& oper0, const int& pos0, int jobOrder = 0);

			/**
			 * destructor
			 */
			~RRBuild() { };

			/**
			 * used for debug print
			 */
			void print() const;

			/**
			 * return the length of items 
			 */
			int getNItems() const; 

			/**
			 * build the concrete RRInt for a given integral of I
			 */
			void buildRRInt(const Integral& I, vector<string>& coeArray,
					vector<int>& indexArray) const;

			/**
			 * build the RHS SQ list for a given shell quartet
			 */
			void buildRRSQ(const ShellQuartet& sq, vector<ShellQuartet>& rrSQList) const;
	};
}

#endif
