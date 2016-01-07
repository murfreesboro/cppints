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
 * \file    nonrr.h
 * \brief   this is to do the non-RR part of integral formula expansion, generate the RRSQ list
 * \author  Fenglai Liu
 *
 * the NON-RR work is to do the integral generation so that to form the RRSQ list, for these non-RR
 * work. For example, the derivatives generation, three body kinetic integrals generation. For 
 * these kinds of work, it's different from RR work by following features:
 * -  it does not have recursive loop, the generation only take place once;
 * -  it does not have completeness check and updateRRSQlist
 *
 * Therefore, the non-RR work could be viewed as a simple version of RR work. We gather all of these
 * kinds of work together to form this nonrr.cpp.
 */
#ifndef NONRR_H
#define NONRR_H
#include "general.h"
#include "rrints.h"
#include "shellquartet.h"
using namespace rrints;
using namespace shellquartet;

namespace nonrr {

	/**
	 * this class is used to calculate the derivatives etc. non-RR work for the given shell quartets
	 */
	class NONRR {

		private:

			int codeSec;                       ///< code section name, DERIV or NON_RR?
			int derivOrder;                    ///< deriv job starts with 1
			int oper;                          ///< the operator of input SQ list
			vector<int> directions;            ///< used by non-RR non-Deriv job
			vector<ShellQuartet> inputSQList;  ///< the input shell quartets with derivatives infor inside
			vector<set<int> > inputIntList;    ///< the integral list corresponding to input SQ list
			vector<ShellQuartet> resultSQList; ///< the output bottom shell quartet list for next stage
			vector<set<int> > unsolvedIntList; ///< the integral list corresponding to result SQ list
			list<RRSQ> rrsqList;               ///< collecting the rrsq result

		public:

			///
			/// contruct the non-RR integral generation 
			/// with input input shell quartet list as well as it's associated
			/// integral list
			///
			NONRR(const vector<ShellQuartet>& sqlist, const vector<set<int> >& intList);

			///
			/// destructor
			///
			~NONRR() { };

			///
			/// this is the driver function to build the whole RRSQList
			///
			void buildRRSQList();

			///
			/// return the result sq list
			///
			const vector<ShellQuartet>& getResultSQList() const { return  resultSQList; };

			///
			/// return the unsolved integrals list for result sq
			///
			const vector<set<int> >& getResultIntSQList() const { return  unsolvedIntList; };

			///
			/// print out the rrsqlist for non rr
			///
			void nonRRPrint(const SQIntsInfor& infor) const;

			///
			/// print out the rrsqlist for deriv
			///
			void derivPrint(const SQIntsInfor& infor) const;

			///
			/// perform possible array index transformation to the rrsq list
			/// 
			/// for the RHS of the non-RR and deriv, it will use the array/variable
			/// form defined for the module. 
			/// 
			/// for the LHS of non-RR and deriv, they are either the final results;
			/// or the input of next code section; so we will refer to the next
			/// code section for array/var information for LHS
			///
			void arrayIndexTransformation(const SQIntsInfor& infor);

			///
			/// checking the LHS integral number for the given batch of input 
			/// shell quartets, this is used to make sure that under the current
			/// deriv position arrangement, how many integrals (LHS) will be 
			/// generated
			///
			size_t evalDerivIntProcess(const SQIntsInfor& infor);

			///
			/// this is used to count the LHS integral numbers for the reslt RRSQ list
			/// 
			int countLHSIntNumbers() const;
	};
}


#endif

