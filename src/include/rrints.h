/**
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
 * \file    rrints.h
 * \brief   describing the RR results
 * \author  Fenglai Liu
 * \note
 *
 * Here we need to explain more about two concepts. One is the "Integral index", the other
 * is the "Array index". Both of them are used in the following classes.
 *
 * Integral index is actually the position for the given integral in the "full" integral
 * list generated by the corresponding shell quartet. For example, the integral
 * <P_y|O|D_zz> has the index of "11" (basis set is in libint order, see basisutil.h). 
 * the full integrals for <P|O|D> shell quartet are listed below:
 * <P_x|O|D_xx> <P_x|O|D_xy> <P_x|O|D_xz> <P_x|O|D_yy> <P_x|O|D_yz> <P_x|O|D_zz> 
 * <P_y|O|D_xx> <P_y|O|D_xy> <P_y|O|D_xz> <P_y|O|D_yy> <P_y|O|D_yz> <P_y|O|D_zz> 
 * <P_z|O|D_xx> <P_z|O|D_xy> <P_z|O|D_xz> <P_z|O|D_yy> <P_z|O|D_yz> <P_z|O|D_zz> 
 * It's clear that <P_y|O|D_zz> has position of 11(starting from 0).
 *
 * The integral index is used to set up a mapping between an integral and its shell
 * quartet. Therefore, from this index we can retrieve the basis set information
 * so that to rebuild the integral. Therefore, we do not need to store the integral
 * in the rrints.
 *
 * On the other hand, array index indicates the practical position for the given 
 * integral in the final RR expression. Still take the <P_y|O|D_zz> as an example,
 * even though it's ranking 11 in the full integral list, however; in the practical
 * RR generation process some of the integrals may not appear therefore for the <P|O|D>
 * shell quartet, it may only have 15 integrals(rather than 18 in the full list). 
 * therefore, the position for <P_y|O|D_zz> maybe different from 11. This is the 
 * meaning of array index.
 *
 * We may use the array index in the final code printing.  
 *
 * Finally, it's worthy to note that if the input shell quartet(oriSQ) is either 
 * the bottom shell quartet for RR process or the result shell quartet for RR
 * process, then it's clear the array index has no difference with the integral
 * index. since they are both full integral list.
 *
 * Here we mainly use set as container to hold the data. The reasons
 * is listed below:
 * - for set, data is automatically sorted
 * - inserting element is automatically processed
 *
 * The RRSQ class is originally only used for RR process. Now we have extended it
 * to the non-RR process like derivatives expression.
 *
 */
#ifndef RRINTS_H
#define RRINTS_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet;

namespace rrbuild {
	class RRBuild;
}

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace rrints {

	using namespace rrbuild;
	using namespace sqintsinfor;

	/**
	 * \class RRSQ
	 *
	 *	This class is used to set up the practical RR search process. The data
	 *	below is in two dimension:
	 *
	 *	first dimension is nInts, that is number of integral in RR for each sq;
	 *	second dimension is nSQ, that is number of sq in RHS of RR
	 *
	 *	For first dimension, we use list since frequently insert function may be
	 *	needed. For second dimension, we use vector since we want directly access
	 *	to any given RHS sq.
	 *
	 *	Additionally, the input/output unsolvedSQList we use set. The benefit for
	 *	set is that we do not need to care about the repeat elements in the RHS
	 *	integrals.
	 *
	 */
	class RRSQ {

		private:

			// basic information
			int rrType;                    ///< rr type
			int oper;                      ///< integral type
			int position;                  ///< where we are going to expand the RR?
			int direction;                 ///< current x, y or z direction (use for 3 body KI etc.)

			// RR integrals expression data
			ShellQuartet oriSQ;            ///< original shell quartet, appearing on LHS
			vector<int> sqPosList;         ///< record each non-null sq position
			vector<ShellQuartet> sqlist;   ///< the shell quartets on the RHS, length is nitems
			list<int>          LHS;        ///< LHS integral index, length is nInts
			vector<list<int> > RHS;        ///< RHS integrals (nInts, nitems)
			vector<list<string> > coe;     ///< coefficients for each terms (nInts, nitems)

			///
			/// for the original shell quartet and unsolved integral list, 
			/// set up the base expression
			///
			void setupExpression(const set<int>& unsolvedIntegralList,const RRBuild& generalRR);

			///
			/// this is the working function to update the LHS expression
			/// for the given unsolvedIntList
			///
			/// if there is LHS integral in unsolvedIntList did not defined;
			/// we will update the expression here by adding this LHS integral
			///
			void updateLHSExpression(const set<int>& unsolvedIntList, const RRBuild& generalRR); 

			///
			/// whether this is derivative work?
			/// when we do derivatives, the first order derivative position 
			/// is at least in some meaningful value
			///
			bool isDerivWork() const {
				int derivOrder = oriSQ.getDerivJobOrder();
				if (derivOrder>0) return true;
				return false;
			};

			/**
			 * whether this is RR work?
			 */
			bool isRRWork() const {
				if (position == BRA1 || position == BRA2 ||
						position == KET1 || position == KET2 ) return true;
				return false;
			};

			/**
			 * whether this is non-deriv work and non-RR work
			 */
			bool isNonRRNonDerivWork() const {
				if (isDerivWork()) return false;
				if (isRRWork()) return false;
				return true;
			};

		public:

			/**
			 * build the RRSQ in the practical RR/derivative process
			 * \param rrtype  what kind of rr, HRR or OS etc.? could be NULL for deriv job
			 * \param pos     position to expand the RR, be NULL for deriv job
			 * \param sq      the input shell quartet
			 * \param unsolvedIntegralList integral index list. The number of integrals 
			 *        contained in the list may less than the full integral number in sq
			 * \param dir  the direction information used by 3 body KI expression etc.
			 *        in default it's NO_DERIV
			 */
			RRSQ(const int& rrType0, const int& pos, const ShellQuartet& sq,
					const set<int>& unsolvedIntegralList, int dir = NO_DERIV);

			///
			/// destructor
			///
			~RRSQ() { };

			///
			/// operator ==
			/// basically, we need to check every information that set up the 
			/// expression
			///
			bool operator==(const RRSQ& rrsq) const {
				return (oriSQ == rrsq.oriSQ && position == rrsq.position && 
						rrType == rrsq.rrType && direction == rrsq.direction);
			};

			///
			/// operator <
			///
			bool operator<(const RRSQ& rrsq) const; 

			/**
			 * get the number of items by counting the sq in sqlist
			 */
			int getNItems() const { return sqlist.size(); };

			/**
			 * return the position
			 */
			int getPosition() const { return position; };

			/**
			 * get the LHS sq
			 */
			const ShellQuartet& getLHSSQ() const { return oriSQ; };

			/**
			 * get the LHS index array
			 */
			const list<int>& getLHSIndexArray() const { return LHS; };

			/**
			 * return the LHS index array - in set form
			 */
			void formLHSIndexSet(set<int>& intList) const;

			/** 
			 * return the unsolved Integral List for the given RHS sq index
			 */
			void getUnsolvedIntList(const int& sqindex, set<int>& unsolvedList) const;

			/**
			 * return the unsolved integral list for the given RHS sq
			 * however, for this function we only search the given set of lhs integral
			 * index which is shown in lhsIndex array
			 */
			void getUnsolvedIntList(const int& sqindex, const set<int>& lhsIndex, 
					set<int>& unsolvedList) const;

			/**
			 * return the RHS shell quartet for the given index
			 */
			const ShellQuartet& getRHSSQ(const int& sqindex) const {
				return sqlist.at(sqindex);
			};

			/**
			 * get the RHS index array for the given item
			 */
			const list<int>& getRHSIndexArray(const int& item) const { return RHS[item]; };

			/**
			 * update the LHS of this rrsq with the given 
			 * unsolvedIntList. We do not check the LHS 
			 * SQ here since it should be realdy checked
			 */
			void updateLHS(const set<int>& unsolvedIntList); 

			/**
			 * check whether the rhs of this rrsq is complete
			 * if it's complete, we return true
			 * else the missing integrals will be pushed into
			 * the missingLHS set for future use
			 */
			bool checkCompleteness(const ShellQuartet& sq, 
					const set<int>& sqLHS, set<int>& missingLHS) const;

			/**
			 * transform the integral index into array index for the rhs
			 * here the rhs term corresponding to a local rrsq in the module
			 */
			void rhsArrayIndexTransform(const RRSQ& rrsq);

			/**
			 * transform the integral index into array index for the rhs
			 * the rhs corresponding to the bottom integrals which is 
			 * defined in the previous module
			 */
			void rhsArrayIndexTransform(const vector<ShellQuartet>& bottomSQList,
					const vector<set<int> >& unsolvedIntList);

			/**
			 * transform the integral index into array index for the lhs
			 * this should be done after call to rhsArrayIndexTransform!
			 */
			void lhsArrayIndexTransform();

			///
			/// this function is used to count every RHS integral (repeat one
			/// also included inside) for the given recurrence relation
			/// on every LHS
			///
			/// it's only that when RHS is NULL (value<0) that we do not count
			///
			int countRHSIntegralNum() const;

			/**
			 * printing the RRSQ(HRR section) to the given code file
			 * \param module the name of the module, must be HRR1 or HRR2
			 * \param infor  provide general sqints information for printing
			 * \param status the result status(module result, or final result etc.)
			 * \param nspace the number of space for printing
			 * \param file   the output file stream
			 */
			void hrrPrint(const int& module, const int& nSpace, const int& status, 
					const vector<ShellQuartet>& moduleResultList, const SQIntsInfor& infor, ofstream& file) const;

			/**
			 * printing the RRSQ(VRR section) to the given code file
			 * \param nspace the number of space for printing
			 * \param file   the output file stream
			 */
			void vrrPrint(const int& nSpace, ofstream& file) const;

			/**
			 * printing the RRSQ(non-RR section) to the given code file
			 * \param module the name of the module, must be DERIV or NON_RR
			 * \param infor  provide general information for printing
			 * \param nspace the number of space for printing
			 * \param file   the output file stream
			 */
			void nonRRPrint(const int& module, const int& nSpace, 
					const SQIntsInfor& infor, ofstream& file) const; 

			/**
			 * printing for debug purpose
			 * we note that this function can not be used after the 
			 * integral index transformed into array index
			 */
			void debug_print() const; 
	};

}


#endif

