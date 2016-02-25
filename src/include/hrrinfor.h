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
 * \brief   utility class to help hrr to print out the codes
 * \author  Fenglai Liu
 */
#ifndef HRRINFOR_H
#define HRRINFOR_H
#include "general.h"
#include "infor.h"
#include "shellquartet.h"
#include "subfilerecord.h"
using namespace infor;
using namespace shellquartet;
using namespace subfilerecord;

namespace rr {
	class RR;
}

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace hrrinfor {

	using namespace rr;
	using namespace sqintsinfor;

	/**
	 * \class HRRInfor
	 */
	class HRRInfor : public Infor {

		private:

			bool hrrFileSplit;                 ///< does the hrr part of code in sub files?
			int section;                       ///< is this the HRR1 or HRR2?
			int nextSection;                   ///< which is the next section of VRR
			int oper;                          ///< we still need operator information
			int side;                          ///< HRR working side, bra or ket?
			vector<ShellQuartet> inputSQList;  ///< input shell quartet list 
			vector<ShellQuartet> outputSQList; ///< output shell quartet list 
			vector<int> inputSQStatus;         ///< the input shell quartet for HRR, in array or variable form?
			vector<int> outputSQStatus;        ///< the output shell quartet for HRR, in array or variable form?
			vector<int> outputSQIntNumList;    ///< the integral number list for the ouput sq list

			//
			// sub files record
			//
			vector<SubFileRecord> subFilesList;  ///< divide the whole HRR into different sub functions

		public:

			///
			/// constructor to form vrr information, initialize the 
			/// content
			///
			HRRInfor(const SQIntsInfor& infor, const RR& hrr);

			///
			/// destructor
			///
			~HRRInfor() { };

			///
			/// whether this is the last section?
			///
         bool isLastSection() const {
				if (nextSection == NULL_POS) return true;
				return false;
			};

			///
			/// return the current section
			///
			int getSection() const { return section; };

			///
			/// declare the array for the given HRR section
			///
			void declareArray(const SQIntsInfor& infor) const;

			///
			/// this is the working function to form sub file records 
			///
			void formSubFiles(const SQIntsInfor& infor, const RR& hrr);

			///
			/// reset the file split status to be true
			///
			void updateFileSplit();

			///
			/// update the HRR output shell quartet list against
			/// other modules input, the input sq in the vector 
			/// are all in array form 
			///
			void updateOutputSQInArray(const vector<ShellQuartet>& moduleInput);

			///
			/// update the HRR input shell quartet list against
			/// other modules output, the input sq in the vector 
			/// are all in array form 
			///
			void updateInputSQInArray(const vector<ShellQuartet>& moduleOutput);

			///
			/// whether HRR is in file split mode?
			///
			bool fileSplit() const { return hrrFileSplit; };

			/// 
			/// return the output sq list
			///
			const vector<ShellQuartet>& getOutputSQList() const { return outputSQList; };

			/// 
			/// return the in sq list
			///
			const vector<ShellQuartet>& getInputSQList() const { return inputSQList; };

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

