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

namespace vrrinfor {
	class VRRInfor;
}

namespace rr {

	using namespace sqintsinfor;
	using namespace vrrinfor;

	/**
	 * RR is responsible to form vertical/horizontal RR for the given group 
	 * of shell quartets
	 *
	 * here in the RR we have two shell quartet list, one is the input SQ list
	 * and the other one is work SQ list. Basically, the two SQ lists are same
	 * for HRR work. However, for VRR it's different. The work SQ list do not
	 * contains the composite shell modifier(division data), and the possible
	 * expoential factors adding etc. are all removed in the workSQList. 
	 *
	 * the workSQList is a necessary copy. For example, in checking the status
	 * of the shell quartet we need the workSQList to see whether this is the 
	 * module result(for VRR).
	 *
	 * the input shell quartet is also necessary. In the contraction step we 
	 * will use the inputSQList and the initial unsolved integral list to 
	 * form the output of VRR. so that's the reason why we also store a copy 
	 * of input unsolved integral list (named as initUnsolvedIntList). 
	 * Together with the inputSQList, such data is used to initialize 
	 * the the RR section's result, also do possible contraction work.
	 *
	 * The bottom sq list and corresponding integral list, is used to preserve
	 * the output shell quartets as well as the integral list for the next 
	 * step code generation. Such data will be generated on the fly during 
	 * the RRSQ list generation.
	 *
	 */
	class RR {

		private:

			int codeSec;                            ///< RR, HRR1 or HRR2
			int rrType;                             ///< RR type of work
			int side;                               ///< for HRR it should be BRA or KET, VRR is not used
			vector<ShellQuartet> inputSQList;       ///< input sq list for RR work
			vector<ShellQuartet> workSQList;        ///< working sq list for RR work
			vector<set<int> > initUnsolvedIntList;  ///< input unsolved integral list
			vector<int> posList;                    ///< for each sq, where the RR expansion did
			vector<ShellQuartet> optRRList;         ///< the opt shell quartet list from rrsqsearch
			vector<ShellQuartet> bottomSQList;      ///< the bottom sq list for HRR
			vector<set<int> > bottomIntList;        ///< integral list correspondiing to bottom sq list for HRR
			list<RRSQ> rrsqList;                    ///< collecting the rrsq result

			////////////////////////////////////
			// functions shared by HRR and VRR//
			////////////////////////////////////

			///
			/// this function is used to update the bottom shell quartet list
			/// as well as the corresponding integral list for HRR module
			///
			void updateBottomSQIntsList(const ShellQuartet& sq, const set<int>& intList); 

			///
			/// search the optimum position for the given shell quartet
			///
			int searchPos(const ShellQuartet& sq) const;

			///
			/// check that whether for the given sq, the 
			/// set of intList already appeared in the LHS of 
			/// corresponding rrsq. If it's old (already 
			/// there), we delete it from the set.
			///
			/// this is used in rr sq lists updation process
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
			/// form the entire rrsqlist from the input unsolved integral list
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
			/// forming the status of the given LHS sq for a rrsq section
			/// for the status meaning, please see the rrints.h's head
			///
			int sqStatusCheck(const SQIntsInfor& infor, const ShellQuartet& sq) const; 

		public:

			/**
			 * constructor
			 * \param codeSec     VRR, HRR1 or HRR2
			 * \param rrtype      what kind of rr, VRR or HRR?
			 * \param inputSQList the result shell quartets on the LHS
			 * \param inputUnsolvedIntList for the result shell quartets, it's corresponding unsolved list
			 */
			RR(const int& codeSec0, const int& rrType0, const vector<ShellQuartet>& inputSQList0,
			const vector<set<int> >& inputUnsolvedIntList);

			///
			/// destructor
			///
			~RR() { };

			///
			/// this is the driver function to generate the RRSQ list
			/// \param side0  bra/ket for the HRR, for VRR is set to NULL, not used
			///
			void generateRRSQList(const int& side0);

			///
			/// for VRR work, we may need to remove all of modifiers (division
			/// information etc.). Because we remove it, it's possible that
			/// the work sq list may become duplicate so as the initial
			/// integral list. This function will perform merging work
			/// for the input sqlist as well as the input unsolvedIntList
			///
			/// we note, that this function may be used in the other places
			/// like VRR Infor class
			///
			void updateSQIntListForVRR(vector<ShellQuartet>& sqlist, 
					vector<set<int> >& unsolvedIntList) const;

			///
			/// this is a simple function to count the integral numbers
			/// for the LHS
			///
			int countLHSIntNumbers() const;

			///
			/// perform possible array index transformation if infor requires
			/// convert the original integral index to array index as required
			///
			/// for the VRR module, because it never uses the array form, therefore
			/// it does not apply to VRR case (will check inside)
			///
			/// for the HRR module, it's a bit of different. RHS is always local
			/// to the module, however; the LHS may be the module result; which
			/// is also the input for the next code section. For these LHS integrals,
			/// we will see whether it uses array index in terms of next code section.
			///
			/// we note that before this function is called, all of index in the 
			/// result RRSQ are integral index.
			///
			void arrayIndexTransformation(const SQIntsInfor& infor);

			///
			/// output bottom sq list for the HRR part
			/// VRR part it's empty
			///
			const vector<ShellQuartet>& getHRRBottomSQList() const { return bottomSQList; }; 

			///
			/// output bottom integrals list for the HRR part
			/// VRR part it's empty
			///
			const vector<set<int> >& getHRRBottomIntList() const { return bottomIntList; }; 

			///
			/// return the input result shell quartet list
			///
			const vector<ShellQuartet>& getRRResultSQList() const { return inputSQList; }; 

			///
			/// return the input unsolved integral list
			///
			const vector<set<int> >& getRRUnsolvedIntList() const { return initUnsolvedIntList; }; 

			///
			/// for VRR module, print out the wholerrsq list to the file
			///
			void vrrPrint(const SQIntsInfor& infor, const VRRInfor& vrrinfor) const;

			///
			/// print HRR codes
			///
			void hrrPrint(const SQIntsInfor& infor) const;

			///
			/// for HRR, determines first side and second side for HRR process
			///
			void sideDeterminationInHRR(int& firstSide, int& secondSide) const;
	};
}


#endif

