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

	// this is used to express the undetermined Ni for derivative expression etc.
	// 1 means the first deriv posistion
	// 2 means the second deriv posistion
	const string unang1   = "uncertainL1";
	const string unangm1  = "uncertainLOneMinus";
	const string unangp1  = "uncertainLOnePlus";
	const string unang2   = "uncertainL2";

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
	 * RRBuild is used to build a general expression for future use. For using 
	 * with RR, in terms of a given RR position (bra1, ket1 etc.)
	 * it's able to return the the RHS shell quartet list for an input LHS
	 * shell quartet or for a given integral it returns the RHS coefficients 
	 * as well as RHS integrals. All of these functions are used in VRR/HRR
	 * process.
	 *
	 * RRbuild not only build the VRR and HRR expressions. Additionally, it's
	 * also used to build the derivative expression(for example, the ERI
	 * 1st or 2ed derivatives), and the other type of expression(for example,
	 * to compute three body kinetic integral from the three body overlap
	 * integral for the given contracted shell quartets). In summary, RRBuild
	 * provides the following three type of functions:
	 * - build general RR expression (VRR and HRR);
	 * - build general derivative expression;
	 * - build a specific type of expression which is non-RR and non-derivative
	 *
	 * the RR expression build will be specified by rrType and the variable of 
	 * position. If position is a NULL value, then this is not RR expression.
	 *
	 * the derivative expression build is specified by derivative positions. 
	 * If all of derivative positions are NULL, then this is not the derivative
	 * expression build.
	 *
	 * The non-RR and non-derivative expression is specified by the operator.
	 * For example, the three body KI calculation is just a non-RR and non-derivative
	 * expression.
	 *
	 * Please remember, for the general derivative expression it's not working 
	 * on the electron coordinates, it's working on the nuclear coordinate.
	 *
	 */
	class RRBuild {

		private:

			// input data
			int rrType;            ///< RR type
			int oper;              ///< integral type
			int position;          ///< in which position we build the general RR
			int firstDerivPos;     ///< the first derivatives position
			int secondDerivPos;    ///< the second derivatives position
			int firstDerivDir;     ///< the first derivatives direction on x, y or z
			int secondDerivDir;    ///< the second derivatives direction on x, y or z
			int derivOrder;        ///< the practical derivatives order

			// data of the general RR/contraction expression
			vector<string> coes;   ///< coefficients for each item in the general expression
			vector<int> lenVars;   ///< used to locate the length of variable in each item 
			vector<int> intTypes;  ///< for each item, it's integral type (only for changed)
			vector<int> angs;      ///< angular momentum changes on variable for each item
			vector<int> vars;      ///< variable changes on variable for each item(bra1 etc.)
			vector<int> mVals;     ///< M value for the integral/shell quartets
			vector<int> exps;      ///< exponetial factors(in number of derivOrder) adding into the terms

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
			 * build the expression for three body KI
			 */
			void buildThreeBodyKIExpression(); 

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
			 * consider the possible number of exponetial factors
			 * adding in for three body KI 
			 * this is a constant
			 */
			int getNExpFacFor3BodyKi() const { return 2; };

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
			 * constructor to build the RR expression
			 * \param rrType   RR working algorithm
			 * \param oper     operator
			 * \param pos      position information, only BRA1, BRA2, KET1, KET2 are allowed
			 */
			RRBuild(const int& rrType0, const int& oper0, const int& pos0);

			/**
			 * constructor to build the derivative expression
			 * \param oper     operator
			 * \param pos1     first derivative position
			 * \param pos2     second derivative position
			 * \param dir1     first derivative direction
			 * \param dir2     second derivative direction
			 */
			RRBuild(const int& oper0, const int& pos1, const int& pos2,
					const int& dir1, const int& dir2);

			/**
			 * constructor to build the non-derivative non-RR expression
			 * for this type of expression only operator is involved
			 * \param oper     operator
			 */
			RRBuild(const int& oper0);

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
			 * build the integral expression expression in terms of 
			 * three body integral for direction xyz(should be x, y or z)
			 */
			void build3BodyKIInt(const Integral& I, const int& direction, 
					vector<string>& coeArray, vector<int>& indexArray) const;

			/**
			 * build the RHS SQ list for a given shell quartet
			 */
			void buildRRSQ(const ShellQuartet& sq, vector<ShellQuartet>& rrSQList) const;

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
			 * whether this is deriv work?
			 */
			bool isDerivWork() const {
				if (derivOrder>0) return true;
				return false;
			};

			/**
			 * whether the derivative work is on same nuclear and same direction? 
			 * for the situation that direction is null, it's assumed to be false
			 */
			bool derivOnSamePosAndDir() const {
				if (firstDerivPos == secondDerivPos && 
						firstDerivDir == secondDerivDir && 
						firstDerivDir != NO_DERIV) return true;
				return false;
			}

			/**
			 * whether this is non-deriv work and non-RR work
			 */
			bool isNonRRNonDerivWork() const {
				if (isDerivWork()) return false;
				if (isRRWork()) return false;
				return true;
			};

	};
}

#endif
