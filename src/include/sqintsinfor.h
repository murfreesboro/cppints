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
 * \file    sqintsinfor.h
 * \brief   processing the input information for the whole program
 * \author  Fenglai Liu
 */
#ifndef SQINTSINFOR_H
#define SQINTSINFOR_H
#include "general.h"
#include "infor.h"
#include "shellquartet.h"
using namespace infor;
using namespace shellquartet;

namespace sqintsinfor {


	/**
	 * \class SQIntsInfor
	 *
	 * This class is the messager for SQInts. In general, SQInts is the commander
	 * and all of data in the SQInts are actually folded here so that SQInts
	 * could send out this infor class to direct the movement of working 
	 * modules.
	 *
	 * For example, SQInts use this infor class to direct the printing in
	 * the SQIntsPrint class, and use this infor class to direct the RR
	 * printing in RR class.
	 */
	class SQIntsInfor : public Infor {

		private:

			bool vrrWithArrayIndex;            ///< do we use array form in vrr part?
			bool hrrWithArrayIndex;            ///< do we use array form in hrr part?
			bool splitFile;                    ///< do we split file into VRR and HRR parts?
			int firstSide;                     ///< HRR work on the first side (BRA/KET)
			int secondSide;                    ///< HRR work on the second side(BRA/KET)
			int jobOrder;                      ///< 0:integral calculation, 1,2: derivatives order
			int oper;                          ///< operator for the input shell quartet
			vector<int> inputShellCodes;       ///< keep a copy of input shell codes
			vector<int> sqOffsetList;          ///< for each single sq, it's offset
			vector<int> braCoeOffset;          ///< for BRA side composite shell quartet 
			                                   ///< record coe offset
			vector<int> ketCoeOffset;          ///< for KET side composite shell quartet 
			                                   ///< record coe offset
			vector<ShellQuartet> inputSQList;  ///< decompose the input sq into single ones 

			///
			/// for the input shell codes, we form the rest of other data here
			/// 
			void formSQInfor();

			///
			/// a small utility function to determine the HRR expanding position
			/// for the given two LCodes
			/// \param LCode1  LCode for BRA1 or KET1
			/// \param LCode2  LCode for BRA2 or KET2
			/// \param side    BRA or KET
			///
			int determineHRRPos(const int& side, const int& LCode1, const int LCode2) const;

			///
			/// for the input shell code, check that for the given
			/// side whether it has S type of shell
			///
			bool hasSTypeSQ(int side) const;

			///
			/// determine the first side and second side through the given shell codes
			/// by using the above two utility functions
			///
			void sideDeterminationInHRR(); 

		public:

			///
			/// constructor for infor class
			///
			SQIntsInfor(const int& job, const int& oper0, 
					const Infor& infor, const int& codeBra1, 
					const int& codeBra2, const int& codeKet1,
					const int& codeKet2);

			///
			/// default destructor
			///
			~SQIntsInfor() { };

			///
			/// whether the given sq is in the result list?
			///
			bool isResult(const ShellQuartet& sq) const;

			///
			/// whether is the composite sq?
			///
			bool isComSQ() const {
				if (inputSQList.size() > 1) return true;
				return false;
			};

			///
			/// whether the input SQ are all bottom ones?
			///
			bool areAllBottomSQ() const;

			///
			/// return the total number of integrals given by the input
			/// shell quartets
			///
			int nInts() const;

			///
			/// for the result shell quartet and it's index, we get its offset
			///
			int getOffset(const ShellQuartet& sq, const int& index) const;

			///
			/// for a given shell quartet, get it's coefficient 
			/// array offset
			/// we note that the shell quartet may not need to 
			/// be the input shell quartet. What we do here is 
			/// to use division information to get the offset 
			/// for coefficients
			///
			void getCoeOffset(const ShellQuartet& sq, 
					int& ic2Offset, int& jc2Offset) const;

			///
			/// do we do split file?
			/// we note, that the maxL should be inputed here
			/// for judgement
			///
			bool splitCPPFile() const { return splitFile; }; 

			///
			/// return the array index status
			///
			bool withArrayIndex(const int& rrType) const { 
				if (rrType == HRR) return hrrWithArrayIndex; 
				return vrrWithArrayIndex;
			};

			///
			/// get the file name for sqints
			/// in default with temp work dir should be set true
			///
			string getFileName(bool withTmpWorkDir = true) const;

			///
			/// get the function name for sqints
			///
			string getFuncName() const;

			///
			/// get the operator
			///
			int getOper() const { return oper; };

			///
			/// get the job order
			///
			int getJobOrder() const { return jobOrder; };

			///
			/// get the side information
			///
			int get1stSide() const { return firstSide; };

			///
			/// get the side information
			///
			int get2edSide() const { return secondSide; };

			///
			/// can we really do HRR?
			///
			bool hasHrr() const {
				bool userChooseHRR = hasHRR();
				if (userChooseHRR) {
					if (firstSide<0 && secondSide<0) return false;
					return true;
				}else{
					return false;
				}
			};

			/**
			 * return the input sq list
			 */
			const vector<ShellQuartet>& getInputSQList() const {
				return inputSQList;
			};

			/**
			 * determine the total length of coefficient array for bra/ket side
			 */
			int getCoeArrayLength(const int& side) const;

			/**
			 * whether the code will be with scr form of vector?
			 */
			bool withSCRVec() const {
				if (hrrWithArrayIndex || vrrWithArrayIndex) {
					if (useSCRVec()) return true;
				}
				return false;
			};

			/**
			 * whether the code will be with TBB or STD form of vector?
			 */
			bool withDoubleVec() const {
				if (hrrWithArrayIndex || vrrWithArrayIndex) {
					if (! useSCRVec()) return true;
				}
				return false;
			};

			/**
			 * whether the code will use the boost library for imcomplete gamma function?
			 */
			bool useBoostGamma() const;
	};

}


#endif

