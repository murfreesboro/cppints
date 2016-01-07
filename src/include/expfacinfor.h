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
 * \file    expfacinforinfor.h
 * \brief   parsing the exponential factors information
 * \author  Fenglai Liu
 */
#ifndef EXPFACINFOR_H
#define EXPFACINFOR_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet;

namespace expfacinfor {

	///
	/// this small class is used to collect the basic exponential
	/// factor information for a given group of shell quartets
	///
	class ExpFacInfor {

		private:

			int expFacListLen;                 ///< the real length of expFacList
			int expFacList[MAX_EXP_FAC_LIST];  ///< possible exponential factors 

		public:

			///
			/// constructor
			///
			ExpFacInfor(const int* expFacList0, const int& expFacListLen0):expFacListLen(expFacListLen0) {
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = expFacList0[i];
			};

			///
			/// destructor
			///
			~ExpFacInfor() { };

			///
			/// whether two given exp factor information class are same?
			///
			bool operator==(const ExpFacInfor& infor0) const {
				if (expFacListLen != infor0.expFacListLen) return false;
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) {
					if (expFacList[i] != infor0.expFacList[i]) return false;
				}
				return true;
			};

			///
			/// whether two given exp factor information class are same?
			///
			bool operator!=(const ExpFacInfor& infor0) const {
				if (expFacListLen != infor0.expFacListLen) return true;
				for(int i=0; i<MAX_EXP_FAC_LIST; i++) {
					if (expFacList[i] != infor0.expFacList[i]) return true;
				}
				return false;
			};

			///
			/// sort the input shell quartet according to the exp factor information
			/// we will pick up all of shell quartets which contains the same exp factor
			/// infor and put it into the output list
			///
			void sortInputSQList(const vector<ShellQuartet>& inputList, 
					vector<ShellQuartet>& outputList) const {
				for(int iSQ=0; iSQ<(int)inputList.size(); iSQ++) {
					ExpFacInfor newInfor(inputList[iSQ].getExpFacList(),inputList[iSQ].getExpFacListLen());
					if (*this == newInfor) outputList.push_back(inputList[iSQ]);
				}
			};
	};

}

#endif
