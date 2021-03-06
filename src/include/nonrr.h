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

namespace nonrrinfor {
	class NONRRInfor;
}

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace nonrr {

	using namespace nonrrinfor;
	using namespace sqintsinfor;

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
			/// return the rrsq list
			///
			const list<RRSQ>& getRRSQList() const { return  rrsqList; };

			///
			/// return the result sq list, which is the bottom shell quartet list
			///
			const vector<ShellQuartet>& getBottomSQList() const { return  resultSQList; };

			///
			/// return the unsolved integrals list for bottom sq
			///
			const vector<set<int> >& getBottomIntList() const { return  unsolvedIntList; };

			///
			/// return the input result shell quartet list
			///
			const vector<ShellQuartet>& getNonRRResultSQList() const { return inputSQList; }; 

			///
			/// return the input unsolved integral list
			///
			const vector<set<int> >& getNonRRResultIntList() const { return inputIntList; }; 

			///
			/// print out the rrsqlist for non rr
			///
			void print(const SQIntsInfor& infor, const NONRRInfor& nonrrinfor);

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

			///
			/// this function count all of RHS integrals (including the 
			/// repeat one) for the non-rr code section
			///
			int countRHSIntNumbers() const;

			///
			/// return the section name
			///
			int getSection() const { return codeSec; };

	};
}


#endif

