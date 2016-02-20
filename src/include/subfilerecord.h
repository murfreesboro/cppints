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

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace subfilerecord {

	using namespace sqintsinfor;

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

			int moduleName;                    ///< the module name
			vector<int> LHSSQStatus;           ///< LHS shell quartet status in this sub file
			vector<int> LHSSQIntNum;           ///< the number of integrals generated for this LHS shell quartet 
			vector<int> RHSSQStatus;           ///< RHS shell quartet status in this sub file
			vector<ShellQuartet> LHSSQList;    ///< LHS shell quartet list in this sub file
			vector<ShellQuartet> RHSSQList;    ///< RHS shell quartet list in this sub file

		public:

			/**
			 * contructor - in default all of elements are empty
			 */
			SubFileRecord(int name):moduleName(name) { };

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
			};

			/**
			 * update the LHS and RHS from RRSQ
			 *
			 * here status will be non-determined, leave for later
			 * determination
			 */
			void updateFromRRSQ(const RRSQ& rrsq);

			/**
			 * return the sub file lhs shell quartets
			 */
			const vector<ShellQuartet>& getLHSSQList() const {
				return LHSSQList; 
			};

			/**
			 * return the lhs shell quartets status
			 */
			const vector<int>& getLHSSQStatus() const {
				return LHSSQStatus; 
			};

			/**
			 * return the sub file RHS shell quartets
			 */
			const vector<ShellQuartet>& getRHSSQList() const {
				return RHSSQList; 
			};

			/**
			 * return the rhs shell quartets status
			 */
			const vector<int>& getRHSSQStatus() const {
				return RHSSQStatus; 
			};

			/**
			 * for the given shell quartet, return it's integral number 
			 * used for array form declare
			 */
			int getLHSIntNum(const ShellQuartet& sq) const;

			/**
			 * to get the M value limit for the LHS shell quartets
			 */
			void getMValueLimit(int& lowerM, int& upperM) const;

			/**
			 * the input is the module output shell quartet list
			 * we will compare this list against the LHS to see 
			 * which one is the module output
			 *
			 * we use the input infor class to judge whether this is 
			 * global result
			 */
			void updateModuleOutput(const SQIntsInfor& infor, const vector<ShellQuartet>& sqlist);

			/**
			 * the input is the module input shell quartet list
			 * we will compare this list against the RHS to see 
			 * which one is the module input
			 */
			void updateModuleInput(const vector<ShellQuartet>& sqlist);

			/**
			 * according to the previous sub file record, update the input 
			 */
			void updateInput(const SubFileRecord& record);

			/**
			 * according to the next sub file record, update the output sq
			 */
			void updateOutput(const SubFileRecord& record);

			///
			/// update all of LHS with the given status
			///
			void updateLHSSQStatus(int status) {
				LHSSQStatus.assign(LHSSQStatus.size(),status);
			};
	};

}


#endif

