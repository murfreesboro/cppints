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
 * \file    rr.h
 * \brief   describing the practical RR process
 * \author  Fenglai Liu
 */
#ifndef RR_H
#define RR_H
#include "general.h"
#include "rrints.h"
#include "shellquartet.h"
using namespace rrints;
using namespace shellquartet;

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace rr {

	using namespace sqintsinfor;

	/**
	 * RR is responsible to form vertical RR for the given group 
	 * of shell quartets
	 */
	class RR {

		private:

			int rrType;                        ///< rr type
			int jobOrder;                      ///< 0:energy, 1,2: derivatives calculation etc.
			int side;                          ///< for HRR, we do RR on bra/ket
			vector<ShellQuartet> inputSQList;  ///< input sq list for RR work
			vector<int> posList;               ///< for each sq, where the RR expansion did
			vector<ShellQuartet> optRRList;    ///< the opt shell quartet list from rrsqsearch
			list<RRSQ> rrsqList;               ///< collecting the rrsq result

			////////////////////////////////////
			// functions shared by HRR and VRR//
			////////////////////////////////////

			///
			/// search the optimum position for the given shell quartet
			///
			int searchPos(const ShellQuartet& sq) const;

			///
			/// check that whether for the given sq, the 
			/// set of intList will appear in the LHS of 
			/// corresponding rrsq. If it's old (already 
			/// there), we delete it from the set.
			///
			void lhsIntegralCheck(const ShellQuartet& sq, 
					set<int>& intList) const;

			///
			/// build/update the rrsq list for the given batch 
			/// of shell quartets and its corresponding unsolved
			/// integral list
			///
			bool buildRRSQList(vector<ShellQuartet>& sqlist, 
					vector<set<int> >& unsolvedIntList);

			///
			/// form the entire rrsqlist
			/// do building/updating work, sorting, checking completeness
			///
			void formRRSQList();

			///
			/// after completeness check, we may still have some missing 
			/// integrals that they only appear on the RHS of RR but 
			/// without define them on the LHS of RR. Therefore, we 
			/// will use this function to recursively generate all of 
			/// missing integrals for the already built RR tree
			///
			void rrUpdating(vector<ShellQuartet>& sqlist, 
					vector<set<int> >& unsolvedIntList);

			///
			/// this is the general function to perform completeness check
			///
			void completenessCheck(); 

			///
			/// this function is used to print head part for each rrsq section
			///
			//void printHead(const int& nSpace, const int& status, const SQIntsInfor& infor, 
			//		const RRSQ& rrsq, ofstream& file) const; 

			///
			/// forming the status of the given LHS sq for a rrsq section
			/// for the status meaning, please see the rrints.h's head
			///
			int sqStatusCheck(const SQIntsInfor& infor, const ShellQuartet& sq) const; 

			///
			/// perform array index transformation
			/// convert the original integral index to array index as required
			///
			void arrayIndexTransformation();

		public:

			/**
			 * initilize the VRR process by calling the RRSQSearch 
			 * \param rrtype      what kind of vrr?
			 * \param inputSQList the input raw shell quartets
			 * \param side        bra/ket for the HRR, for VRR is set to NULL
			 * \param infor       provides additional array forming information
			 */
			RR(const int& rrType0, const vector<ShellQuartet>& inputSQList0,
				const SQIntsInfor& infor,int side0 = NULL_POS);

			///
			/// destructor
			///
			~RR() { };

			///
			/// get bottom sq list for this rr process
			///
			void getBottomSQList(vector<ShellQuartet>& sqlist) const; 

			///
			/// printing the whole rrsq list for the RR body to the given file
			/// we note that withArrayIndex should be consistent with the 
			/// doArrayIndexTransform option passed into the constructor
			///
			/// in theory this function should be constant member function.
			/// however, since we use reverse iterator inside for rrsqlist,
			/// it can not be constant anymore. However, still safe enough
			/// for use.
			///
			void print(const SQIntsInfor& infor, const string& filename); 
	};
}


#endif

