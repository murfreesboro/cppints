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
 * \file    hrrinfor.h
 * \brief   utility class to help non-RR modules to print out the codes
 * \author  Fenglai Liu
 */
#ifndef NONRRINFOR_H
#define NONRRINFOR_H
#include "general.h"
#include "infor.h"
#include "shellquartet.h"
#include "subfilerecord.h"
using namespace infor;
using namespace shellquartet;
using namespace subfilerecord;

namespace nonrr {
	class NONRR;
}

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace nonrrinfor {

	using namespace nonrr;
	using namespace sqintsinfor;

	/**
	 * \class NONRRInfor
	 */
	class NONRRInfor : public Infor {

		private:

			bool nonrrFileSplit;               ///< does the nonrr part of code in sub files?
			int section;                       ///< is this the non-RR or DERIV?
			int nextSection;                   ///< which is the next section of VRR
			vector<ShellQuartet> inputSQList;  ///< input shell quartet list 
			vector<ShellQuartet> outputSQList; ///< output shell quartet list 
			vector<int> inputSQStatus;         ///< the input shell quartet for HRR, in array or variable form?
			vector<int> outputSQStatus;        ///< the output shell quartet for HRR, in array or variable form?
			vector<int> outputSQIntNumList;    ///< the integral number list for the ouput sq list

			//
			// sub files record
			//
			vector<SubFileRecord> subFilesList;  ///< divide the whole NONRR into different sub functions

		public:

			///
			/// constructor to form nonrr information, initialize the 
			/// content
			///
			NONRRInfor(const SQIntsInfor& infor, const NONRR& nonrr);

			///
			/// destructor
			///
			~NONRRInfor() { };

			///
			/// whether this is the last section?
			///
         bool isLastSection() const {
				if (nextSection == NULL_POS) return true;
				return false;
			};

			///
			/// declare the array for the given NON-RR section
			///
			void declareArray(const SQIntsInfor& infor) const;

			///
			/// this is the working function to form sub file records 
			///
			void formSubFiles(const SQIntsInfor& infor, const NONRR& nonrr);

			///
			/// update the NONRR output shell quartet list against
			/// other modules input. This does not apply to DERIV module
			///
			void updateOutputSQInArray(const vector<ShellQuartet>& moduleInput);

			///
			/// update the NONRR input shell quartet list against
			/// other modules output, the input sq in the vector 
			/// are all in array form 
			///
			void updateInputSQInArray(const vector<ShellQuartet>& moduleOutput);

			///
			/// whether NONRR is in file split mode?
			///
			bool fileSplit() const { return nonrrFileSplit; };

			/// 
			/// return the output sq list
			///
			const vector<ShellQuartet>& getOutputSQList() const { return outputSQList; };

			/// 
			/// return the output sq status list
			///
			const vector<int>& getOutputSQStatus() const { return outputSQStatus; };

			///
			/// return the input sq status list
			///
			const vector<int>& getInputSQStatus() const { return inputSQStatus; };

			///
			/// get the number of sub file records
			///
			int getNSubFiles() const { return subFilesList.size(); };

			///
			/// return the given sub file record
			///
			const SubFileRecord& getSubFileRecord(int i) const { return subFilesList[i]; };
	};

}


#endif

