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
 * \file    sqintsprint.h
 * \brief   control the printing for code generation
 * \author  Fenglai Liu
 */
#ifndef SQINTSPRINT_H
#define SQINTSPRINT_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet;

namespace sqintsprint {

	/**
	 * \class SubFileRecord
	 *
	 * This is the class to collect shell quartet status in terms of file split
	 * for each sub file used in all code sections.
	 */
	class SubFileRecord {

		private:

			vector<int> LHSSQStatus;           ///< LHS shell quartet status in this sub file
			vector<int> RHSSQStatus;           ///< RHS shell quartet status in this sub file
			vector<ShellQuartet> LHSSQList;    ///< LHS shell quartet list in this sub file
			vector<ShellQuartet> RHSSQList;    ///< RHS shell quartet list in this sub file
			vector<ShellQuartet> inputSQList;  ///< the input shell quartet list for this sub file
			vector<ShellQuartet> outputSQList; ///< the output shell quartet list for this sub file

		public:

			/**
			 * contructor - in default all of elements are empty
			 */
			SubFileRecord() { };

			/**
			 * destructor
			 */
			~SubFileRecord() { };

			/**
			 * initialize the content, reserve space for later push back
			 */
			void init() {
				int n = 500;
				LHSSQStatus.reserve(n);  
				RHSSQStatus.reserve(n); 
				LHSSQList.reserve(n);   
				RHSSQList.reserve(n);   
				inputSQList.reserve(n); 
				outputSQList.reserve(n);
			};

			/**
			 * clear it's content
			 */
			void clear() {
				LHSSQStatus.clear();  
				RHSSQStatus.clear(); 
				LHSSQList.clear();   
				RHSSQList.clear();   
				inputSQList.clear(); 
				outputSQList.clear();
			};

			/**
			 * update the LHS and RHS from RRSQ
			 *
			 * here status will be non-determined, leave for later
			 * determination
			 */
			void updateFromRRSQ(const RRSQ& rrsq);

			/**
			 * update shell quartet status from section input/output  
			 */
			void updateFromSection(const vector<ShellQuartet>& sectionInput, 
					const vector<ShellQuartet>& sectionOutput);

			/**
			 * update the sub output form the following sub files
			 */
			void updateSubOutput(const SubFileRecord& followingRecord);

			/**
			 * update the sub input form the previous sub files
			 */
			void updateSubInput(const SubFileRecord& previousRecord);

	};

	/**
	 * \class VRRSectionRecord
	 *
	 * This is to form the printing information in terms of VRR section
	 *
	 * several things to note:
	 * - VRR does not have input shell quartet list. All of input for VRR
	 *   should be directly calculated (VRR bottom integrals)
	 */
	class VRRSectionRecord {

		private:

			vector<ShellQuartet> outputSQList;       ///< the output shell quartet list for VRR
			vector<set<int> > outputIntegralList;    ///< the output integral index for each output shell quartet
			vector<SubFileRecord> subFunctionsList;  ///< divide the whole VRR into different sub functions

		public:

			/**
			 * contructor 
			 */
			VRRSectionRecord(const vector<ShellQuartet>& outputList,
					const vector<set<int> >& integralList):outputSQList(outputList),
			outputIntegralList(integralList) { };

			/**
			 * destructor
			 */
			~VRRSectionRecord() { };

			/**
			 * form the sub files from input RRSQ list
			 */
			void updateFromRRSQ(const RRSQ& rrsq);

			/**
			 * update shell quartet status from section input/output  
			 */
			void updateFromSection(const vector<ShellQuartet>& sectionInput, 
					const vector<ShellQuartet>& sectionOutput);

			/**
			 * update the sub output form the following sub files
			 */
			void updateSubOutput(const SubFileRecord& followingRecord);

			/**
			 * update the sub input form the previous sub files
			 */
			void updateSubInput(const SubFileRecord& previousRecord);

	};

}


#endif

