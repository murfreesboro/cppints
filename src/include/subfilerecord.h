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
 * \file    subfilerecord.h
 * \brief   record for sub file forming
 * \author  Fenglai Liu
 */
#ifndef SUBFILERECORD_H
#define SUBFILERECORD_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet;

namespace subfilerecord {

	/**
	 * \class SubFileRecord
	 *
	 * When the given code section is split, we will have one or multiple 
	 * sub files. Each sub file has a working function defined to print
	 * a chunk of rrsqlist of the whole code section.
	 *
	 * sub file record is used to record the information about the sub file,
	 * including in/out, how many RRSQ terms inside; and most importantly,
	 * the array/variable form the LHS/RHS shell quartets.
	 */
	class SubFileRecord {

		private:

			bool hasABCD;                      ///< whether this sub file output the final results?
			int moduleName;                    ///< the module name
			vector<int> LHSSQStatus;           ///< LHS shell quartet array/var status in this sub file
			vector<int> LHSSQIntNum;           ///< the number of integrals generated for this LHS shell quartet 
			vector<int> RHSSQStatus;           ///< RHS shell quartet array/var status in this sub file
			vector<ShellQuartet> LHSSQList;    ///< LHS shell quartet list in this sub file
			vector<ShellQuartet> RHSSQList;    ///< RHS shell quartet list in this sub file
			vector<ShellQuartet> inputSQList;  ///< the input shell quartet list for this sub file
			vector<ShellQuartet> outputSQList; ///< the output shell quartet list for this sub file

		public:

			/**
			 * contructor - in default all of elements are empty
			 */
			SubFileRecord(int name):hasABCD(false),moduleName(name) { };

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
				LHSSQIntNum.reserve(n); 
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
				LHSSQIntNum.calr(); 
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
			 * return the sub file output shell quartets
			 */
			const vector<ShellQuartet>& getFunctionOutput() const {
				return outputSQList; 
			};

			/**
			 * return the sub file input shell quartets
			 */
			const vector<ShellQuartet>& getFunctionInput() const {
				return inputSQList; 
			};

			/**
			 * for the given shell quartet, return it's integral number 
			 * used for array form declare
			 *
			 * the input sq must be the module output 
			 */
			int getLHSIntNum(const ShellQuartet& sq) const;

			/**
			 * the input is the module output shell quartet list
			 * we will compare this list against the LHS to see 
			 * which one is the module output
			 *
			 * For module output, it must be in the function output list
			 * so it must choose either array form or it's final result
			 *
			 * we use the input infor class to judge whether this is 
			 * global result
			 */
			void updateModuleOutput(const SQIntsInfor& infor, const vector<ShellQuartet>& sqlist);

			/**
			 * for VRR, the update of output list and lhs sq status
			 * comparing with the input module output (sqlist) is 
			 * a bit of different from other modules
			 */
			void updateVRROutput(bool destroyMultiplerInfor,
					const SQIntsInfor& infor, const vector<ShellQuartet>& sqlist);
	};

}


#endif

